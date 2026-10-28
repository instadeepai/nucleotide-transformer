# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests that the window plan, including the tail shift-back, reaches the output."""

from pathlib import Path

from nucleotide_transformer_v3.gene_annotation.annotate_genome import run
from nucleotide_transformer_v3.gene_annotation.annotate_utils import select_windows

_SHIPPED_FILTERS: dict[str, object] = {
    "decoder": "hmm_constrained",
    "apply_biological_filters": True,
    "min_cds_length": 201,
}


def _gene_rows(text: str) -> list[list[str]]:
    """Every ``gene`` row of a GFF3, split into columns."""
    return [
        cols
        for cols in (
            line.split("\t")
            for line in text.splitlines()
            if line and not line.startswith("#")
        )
        if cols[2] == "gene"
    ]


def _strands_over(
    genes: list[list[str]], seqid: str, span: tuple[int, int]
) -> list[str]:
    """Strands of the genes on ``seqid`` meeting the 1-based closed ``span``."""
    start, end = span
    return [
        cols[6]
        for cols in genes
        if cols[0] == seqid and int(cols[3]) <= end and int(cols[4]) >= start
    ]


def test_select_windows_pulls_the_last_region_onto_the_contig_end(
    production_config,
) -> None:
    """The last kept region is shifted back to end on the contig's final base."""
    config = production_config()
    contig_length = 4_096
    windows, crop = select_windows(
        sequence_length=config.sequence_length,
        keep_target_center_fraction=config.keep_target_center_fraction,
        chrom_lengths={"chrShift": contig_length},
        chrom="chrShift",
    )
    stride = config.sequence_length - 2 * crop

    assert (crop, stride) == (512, 3_072), "the golden geometry moved"
    kept_starts = [window.start + crop for window in windows]
    assert kept_starts == [0, contig_length - stride]


def test_the_window_plan_decides_what_the_tail_of_a_contig_emits(
    probe_genome: list[str],
    probe_spans: dict[str, tuple[int, int]],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """The shift-back puts ``chrShift``'s probe gene on the minus strand.

    ``chrProbe`` is checked first, so a probe that stopped emitting fails there.
    """
    output = tmp_path / "probe.gff3"
    run(
        production_config(
            fasta_paths=probe_genome,
            output_path=str(output),
            **_SHIPPED_FILTERS,
        )
    )
    genes = _gene_rows(output.read_text())

    control = _strands_over(genes, "chrProbe", probe_spans["chrProbe"])
    assert control == ["+"], f"the probe itself is broken: {control}"

    tail = _strands_over(genes, "chrShift", probe_spans["chrShift"])
    assert tail == ["-"], f"the tail window was not shifted back: {tail}"
