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

"""Tests for the output-format dispatch and GTF ``exon_number`` in ``formats.py``.

``exon_number`` runs in transcript order (Ensembl/GENCODE), so on the minus
strand it counts down as genomic coordinates rise.
"""

import pytest

from nucleotide_transformer_v3.gene_annotation.formats import (
    FORMAT_BY_SUFFIX,
    FORMAT_GFF3,
    FORMAT_GTF,
    gff3_lines_to_gtf,
    resolve_output_format,
)

# Every extension the module advertises, with the format it must select spelled
# out rather than read back from FORMAT_BY_SUFFIX: a test that consults the same
# table the resolver does would pass on any table.
_EXPECTED_FORMAT_BY_SUFFIX: dict[str, str] = {
    ".gff3": FORMAT_GFF3,
    ".gff": FORMAT_GFF3,
    ".gtf": FORMAT_GTF,
}

# Extensions that must be refused rather than guessed at. ``.fa`` and ``.txt``
# are plausible typos, ``.gff2`` is a near miss on an accepted one, ``.gz``
# alone leaves nothing behind once the suffix is stripped, and a bare name has
# no extension at all.
_REJECTED_PATHS = (
    "out.fa",
    "out.txt",
    "out.gff2",
    "out.gz",
    "out.bed.gz",
    "out",
    "",
)


def test_the_expected_table_covers_every_advertised_extension() -> None:
    """Every extension ``FORMAT_BY_SUFFIX`` advertises has a case here."""
    assert set(_EXPECTED_FORMAT_BY_SUFFIX) == set(FORMAT_BY_SUFFIX)


@pytest.mark.parametrize("suffix", sorted(_EXPECTED_FORMAT_BY_SUFFIX))
@pytest.mark.parametrize("gzipped", [False, True])
@pytest.mark.parametrize("upper", [False, True])
def test_resolve_output_format_selects_the_format_the_extension_names(
    suffix: str, gzipped: bool, upper: bool
) -> None:
    """Each accepted extension selects its own format, in any case, gzipped or not."""
    path = f"/tmp/annotation{suffix.upper() if upper else suffix}"
    if gzipped:
        path += ".GZ" if upper else ".gz"

    assert resolve_output_format(path) == _EXPECTED_FORMAT_BY_SUFFIX[suffix]


@pytest.mark.parametrize("path", _REJECTED_PATHS)
def test_resolve_output_format_rejects_an_unsupported_extension(path: str) -> None:
    """An unsupported extension is refused rather than guessed at."""
    with pytest.raises(ValueError, match="Cannot infer an output format"):
        resolve_output_format(path)


def _gff3_gene(gene_id: str, strand: str, blocks: list[tuple[int, int]]) -> list[str]:
    """One GFF3 gene row plus a CDS row per block, as ``find_genes`` writes them.

    Args:
        gene_id: The gene's ``ID``, which its CDS rows name as ``Parent``.
        strand: ``"+"`` or ``"-"``.
        blocks: CDS blocks as 1-based inclusive ``(start, end)``, ascending.

    Returns:
        The GFF3 feature lines, gene row first.
    """
    start, end = blocks[0][0], blocks[-1][1]
    lines = [f"chrTest\tNTv3_HMM\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id}"]
    lines += [
        f"chrTest\tNTv3_HMM\tCDS\t{lo}\t{hi}\t.\t{strand}\t0\tParent={gene_id}"
        for lo, hi in blocks
    ]
    return lines


def _exon_numbers_by_start(lines: list[str], gene_id: str) -> dict[int, int]:
    """``exon`` rows of one gene as genomic start -> ``exon_number``."""
    numbers: dict[int, int] = {}
    for line in lines:
        cols = line.split("\t")
        if cols[2] != "exon" or f'gene_id "{gene_id}";' not in cols[8]:
            continue
        number = cols[8].split('exon_number "')[1].split('"')[0]
        numbers[int(cols[3])] = int(number)
    return numbers


def test_minus_strand_exon_number_counts_down_in_genomic_order() -> None:
    """Minus-strand exons number from the highest coordinate, plus from the lowest."""
    blocks = [(100, 199), (300, 399), (500, 599)]
    lines = gff3_lines_to_gtf(
        _gff3_gene("gene_m", "-", blocks) + _gff3_gene("gene_p", "+", blocks)
    )

    assert _exon_numbers_by_start(lines, "gene_m") == {100: 3, 300: 2, 500: 1}
    assert _exon_numbers_by_start(lines, "gene_p") == {100: 1, 300: 2, 500: 3}


def test_a_cds_row_carries_the_exon_number_of_the_exon_it_pairs_with() -> None:
    """Each block's ``exon`` and ``CDS`` rows agree on the number."""
    lines = gff3_lines_to_gtf(_gff3_gene("gene_m", "-", [(100, 199), (300, 399)]))
    numbers = {
        (cols[2], cols[3]): cols[8].split('exon_number "')[1].split('"')[0]
        for cols in (line.split("\t") for line in lines)
        if cols[2] in {"exon", "CDS"}
    }

    assert numbers[("exon", "100")] == numbers[("CDS", "100")] == "2"
    assert numbers[("exon", "300")] == numbers[("CDS", "300")] == "1"
