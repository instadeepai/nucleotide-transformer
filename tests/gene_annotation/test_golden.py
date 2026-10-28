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

"""Tier 1: production's GFF3 feature rows are pinned line for line.

Every non-comment line is compared exactly; the ``#`` lines are stripped, so
the ``##ntv3-*`` provenance header is outside the pin rather than part of it.

This is the regression contract every later efficiency change is measured
against. Multi-device sharding, a larger batch size and a different Viterbi
implementation must all leave these files identical, and if one of them does
not, that is the finding. It drives the real ``run()`` end to end (window
planning, batching, the real decoders, stitching, the filters and the streaming
writer) against a deterministic stub head, so it needs no accelerator, no
network and no trained weights, and it is exact rather than tolerance-bounded.

**Why four goldens.** ``hmm_constrained`` is the shipped default, and it vetoes
states whose sequence context lacks the required start, stop and splice motifs,
so pinning it only pins something worth pinning if the synthetic genome carries
those motifs. It does: the stub reads soft-mask case and the decoder reads the
base (see ``conftest``), so the same bases spell a real ``ATG``, real
``gt``..``ag`` and ``at``..``ac`` intron boundaries and a real ``TAA`` while
still telling the stub which positions are CDS, which are intron and which are
intergenic. The constrained pair therefore carries multi-exon genes on both
strands, and both post-processing filters remove genes from it. The ``hmm``
pair is kept beside it because the unconstrained decoder is a real code path
with a different answer, and because turning the filters on breaks
comparability with ``annotate/``'s raw output.

``argmax`` is deliberately not pinned. It reads each position on its own, so on
this input it returns the stub's token groups back: one ``CDS_p0`` run per exon
block, phase 0 throughout, no start, stop or splice state anywhere. That would
pin the stub rather than the pipeline, which matches ``argmax``'s documented
status as a diagnostic rather than an output.

The goldens are not regenerated from this package. They came from the research
pipeline, which is what makes agreeing with them evidence of anything.
A failure here means this pipeline's output moved away from that reference.
"""

from pathlib import Path

import pytest

from nucleotide_transformer_v3.gene_annotation.annotate_genome import run

GOLDEN_DIR = Path(__file__).parent / "data"

_MIN_CDS_LENGTH = 201
_FILTERS_ON = {"apply_biological_filters": True, "min_cds_length": _MIN_CDS_LENGTH}
_FILTERS_OFF = {"apply_biological_filters": False, "min_cds_length": 0}

# The pinned configurations, one filtered / unfiltered pair per decoder.
_GOLDENS: dict[str, dict[str, object]] = {
    "hmm_unfiltered": {"decoder": "hmm", **_FILTERS_OFF},
    "hmm_filtered": {"decoder": "hmm", **_FILTERS_ON},
    "hmm_constrained_unfiltered": {"decoder": "hmm_constrained", **_FILTERS_OFF},
    "hmm_constrained_filtered": {"decoder": "hmm_constrained", **_FILTERS_ON},
}

_DECODERS = ("hmm", "hmm_constrained")

# The one configuration also pinned in GTF: the shipped decoder and filters.
_GTF_GOLDEN = "hmm_constrained_filtered"

_CONTIGS = {"chrShort", "chrSeam", "chrTail", "chrHead"}

# Splice-site pairs the constrained decoder is allowed to open an intron on,
# read in transcribed orientation: the U2 donors with their AG acceptor, and the
# U12 lane's AT..AC. Anything else means the motif mask let through an intron the
# grammar has no state for.
_LEGAL_SPLICE_PAIRS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
_U12_SPLICE_PAIR = ("AT", "AC")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _feature_lines(text: str) -> list[str]:
    """Non-comment, non-blank lines, i.e. everything provenance-independent."""
    return [ln for ln in text.splitlines() if ln and not ln.startswith("#")]


def _golden(name: str) -> list[str]:
    """Read one committed golden's feature lines."""
    return _feature_lines((GOLDEN_DIR / f"golden_{name}.gff3").read_text())


def _genes(lines: list[str]) -> list[str]:
    """Only the ``gene`` rows."""
    return [ln for ln in lines if ln.split("\t")[2] == "gene"]


def _gene_spans(lines: list[str]) -> set[tuple[str, int, int, str]]:
    """Every gene as ``(seqid, start, end, strand)``, ignoring its id.

    Gene ids renumber whenever a filter drops an earlier gene, so they cannot be
    compared across configurations; the span and strand can.
    """
    return {
        (cols[0], int(cols[3]), int(cols[4]), cols[6])
        for cols in (line.split("\t") for line in lines)
        if cols[2] == "gene"
    }


def _cds_length_by_gene(lines: list[str]) -> dict[str, int]:
    """Total CDS bp per gene id, which is what ``min_cds_length`` thresholds."""
    totals: dict[str, int] = {}
    for line in lines:
        cols = line.split("\t")
        if cols[2] != "CDS":
            continue
        parent = cols[8].split("Parent=")[1]
        totals[parent] = totals.get(parent, 0) + int(cols[4]) - int(cols[3]) + 1
    return totals


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_production_output_matches_the_golden(
    name: str,
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """The full pipeline's GFF3 feature lines match the committed golden."""
    output = tmp_path / f"{name}.gff3"
    run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            **_GOLDENS[name],
        )
    )

    # Deliberately no regeneration switch. The one property that makes these
    # goldens worth having is that they came from a different implementation;
    # a flag that rewrites them from this pipeline's own output turns the
    # suite into a self-comparison, and nothing in CI would notice. They are
    # regenerated upstream or not at all.
    assert _feature_lines(output.read_text()) == _golden(name)


def _golden_gtf() -> list[str]:
    """Read the committed GTF golden's feature lines."""
    return _feature_lines((GOLDEN_DIR / f"golden_{_GTF_GOLDEN}.gtf").read_text())


def test_production_gtf_matches_the_golden(
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """Asking for ``.gtf`` yields GTF, identical to the committed golden."""
    output = tmp_path / f"{_GTF_GOLDEN}.gtf"
    run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            **_GOLDENS[_GTF_GOLDEN],
        )
    )

    assert _feature_lines(output.read_text()) == _golden_gtf()


def test_the_gtf_golden_carries_a_multi_exon_minus_strand_gene() -> None:
    """The GTF golden exercises minus-strand ``exon_number``, counting down."""
    numbers: dict[str, list[tuple[int, int]]] = {}
    for line in _golden_gtf():
        cols = line.split("\t")
        if cols[2] != "exon" or cols[6] != "-":
            continue
        gene_id = cols[8].split('gene_id "')[1].split('"')[0]
        number = int(cols[8].split('exon_number "')[1].split('"')[0])
        numbers.setdefault(gene_id, []).append((int(cols[3]), number))

    descending = [
        gene_id
        for gene_id, blocks in numbers.items()
        if len(blocks) > 1
        and [number for _, number in sorted(blocks)]
        == sorted((number for _, number in blocks), reverse=True)
    ]
    assert descending, (
        "no minus-strand gene in the GTF golden has more than one exon numbered "
        "in descending genomic order, so the reversal is untested"
    )


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_the_golden_is_not_empty(name: str) -> None:
    """Guard the guard: an empty golden would pin nothing at all."""
    assert len(_genes(_golden(name))) >= 5


@pytest.mark.parametrize("decoder", _DECODERS)
def test_the_filters_actually_remove_genes(
    decoder: str, are_goldens_present: None
) -> None:
    """Guard the guard: identical goldens leave the filtered case untested.

    If the synthetic genome produced no gene the filters remove, both files
    would be the same and the filtered configuration would be pinned only by
    accident. Asserted for the shipped decoder as well as for ``hmm``, because
    a filter is only pinned on the decoder whose output it was run on.
    """
    filtered = _golden(f"{decoder}_filtered")
    unfiltered = _golden(f"{decoder}_unfiltered")
    assert filtered != unfiltered
    assert len(_genes(filtered)) < len(_genes(unfiltered))


def test_each_filter_removes_genes_under_the_shipped_decoder(
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
    are_goldens_present: None,
) -> None:
    """Both filters bite on ``hmm_constrained``, not just one of them.

    The two pinned configurations only bracket the filters together, so one of
    them removing everything would look the same as both working. Running the
    in-frame-stop filter on its own splits the difference: the gene set has to
    shrink once when it is enabled and again when ``min_cds_length`` joins it.
    """
    output = tmp_path / "stop_filter_only.gff3"
    run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            decoder="hmm_constrained",
            apply_biological_filters=True,
            min_cds_length=0,
        )
    )

    unfiltered_lines = _golden("hmm_constrained_unfiltered")
    unfiltered = _gene_spans(unfiltered_lines)
    stop_filtered = _gene_spans(_feature_lines(output.read_text()))
    both_filtered = _gene_spans(_golden("hmm_constrained_filtered"))

    assert stop_filtered < unfiltered, "the in-frame-stop filter removed nothing"
    assert both_filtered < stop_filtered, "min_cds_length removed nothing"

    # The in-frame-stop filter is the only one that can reach a gene this long,
    # so the gene it removed is not one min_cds_length would have removed anyway.
    cds_length = _cds_length_by_gene(unfiltered_lines)
    id_by_span = {
        (cols[0], int(cols[3]), int(cols[4]), cols[6]): cols[8].removeprefix("ID=")
        for cols in (line.split("\t") for line in _genes(unfiltered_lines))
    }
    removed = {id_by_span[span] for span in unfiltered - stop_filtered}
    assert any(cds_length[gene_id] >= _MIN_CDS_LENGTH for gene_id in removed)


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_the_golden_covers_intron_subfeatures(
    name: str, are_goldens_present: None
) -> None:
    """``find_genes``' sub-feature partition must be on the pinned path.

    A golden of single-CDS genes would not notice a change in how runs of equal
    GFF3 type are split.
    """
    types = {line.split("\t")[2] for line in _golden(name)}
    assert types == {"gene", "CDS", "intron"}


@pytest.mark.parametrize(
    "name", ["hmm_constrained_unfiltered", "hmm_constrained_filtered"]
)
def test_constrained_introns_sit_on_real_splice_motifs(
    name: str,
    synthetic_records: dict[str, str],
    are_goldens_present: None,
) -> None:
    """Every pinned intron has legal boundaries, and one of them is U12.

    This is what the constrained decoder is for, so it is checked against the
    sequence rather than assumed. The U12 (``AT``..``AC``) case matters twice
    over: it has no path through the 15-state grammar at all, so it exercises
    the paired extended lane, which nothing else in the suite reaches.
    """
    pairs = [
        _splice_pair(synthetic_records[cols[0]], int(cols[3]), int(cols[4]), cols[6])
        for cols in (line.split("\t") for line in _golden(name))
        if cols[2] == "intron"
    ]
    assert pairs, "no intron to check"
    illegal = [pair for pair in pairs if pair not in _LEGAL_SPLICE_PAIRS]
    assert not illegal, f"introns on impossible splice sites: {illegal}"
    assert _U12_SPLICE_PAIR in pairs, "the U12 decode lane is not covered"


def _splice_pair(sequence: str, start: int, end: int, strand: str) -> tuple[str, str]:
    """The (donor, acceptor) dinucleotides of one intron, in transcribed order.

    Args:
        sequence: The contig the intron sits on.
        start: GFF3 1-based inclusive start.
        end: GFF3 1-based inclusive end.
        strand: ``"+"`` or ``"-"``; a minus-strand intron is read 3' to 5' in
            genomic coordinates, so its donor is at the genomic end.

    Returns:
        The two dinucleotides, uppercased (soft-mask case is not a motif).
    """
    left, right = sequence[start - 1 : start + 1], sequence[end - 2 : end]
    if strand == "+":
        return left.upper(), right.upper()
    return (
        right.translate(_COMPLEMENT)[::-1].upper(),
        left.translate(_COMPLEMENT)[::-1].upper(),
    )


def test_every_contig_appears_in_the_golden(are_goldens_present: None) -> None:
    """All four geometry branches must reach the output, not just the easy ones.

    A contig silently dropped by the window planner would otherwise shrink the
    golden without failing anything.
    """
    for name in _GOLDENS:
        seqids = {line.split("\t")[0] for line in _golden(name)}
        assert seqids == _CONTIGS, f"golden_{name}.gff3 covers only {seqids}"


@pytest.mark.parametrize("name", sorted(_GOLDENS))
def test_golden_coordinates_and_phases_are_well_formed(name: str) -> None:
    """Structural invariants of the pinned files, checked rather than assumed.

    A golden generated from a bug is a bug that is now pinned, so the properties
    that would reveal one are asserted here instead of trusted.

    Nesting is checked against each sub-feature's ``Parent``, not against the
    preceding gene row. ``write_gff3`` sorts globally by ``(seqid, start)``, so a
    gene's sub-features are interleaved with any other gene overlapping the same
    span, and on the opposite strand that is common.
    """
    lines = _golden(name)
    gene_spans: dict[str, tuple[int, int]] = {}
    for line in lines:
        cols = line.split("\t")
        if cols[2] == "gene":
            gene_id = cols[8].removeprefix("ID=")
            assert gene_id not in gene_spans, f"duplicate gene id: {line}"
            gene_spans[gene_id] = (int(cols[3]), int(cols[4]))

    for line in lines:
        cols = line.split("\t")
        start, end = int(cols[3]), int(cols[4])
        assert start >= 1, f"GFF3 is 1-based, got start {start}: {line}"
        assert end >= start, f"end before start: {line}"
        assert cols[6] in {"+", "-"}, f"bad strand: {line}"

        if cols[2] != "gene":
            parent = cols[8].split("Parent=")[1]
            assert parent in gene_spans, f"sub-feature with no gene row: {line}"
            lo, hi = gene_spans[parent]
            assert lo <= start and end <= hi, (
                f"sub-feature outside its parent {parent} ({lo}, {hi}): {line}"
            )

        if cols[2] == "CDS":
            assert cols[7] in {"0", "1", "2"}, f"bad CDS phase: {line}"
        else:
            assert cols[7] == ".", f"phase set on a non-CDS feature: {line}"


@pytest.fixture
def are_goldens_present() -> None:
    """Skip the cross-golden guards while the goldens are being regenerated."""
    missing = [
        name for name in _GOLDENS if not (GOLDEN_DIR / f"golden_{name}.gff3").exists()
    ]
    if missing:
        pytest.skip(f"goldens not generated yet: {missing}")


@pytest.mark.parametrize("batch_size", [1, 3, 5])
def test_output_is_invariant_to_batch_size(
    batch_size: int,
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """Batching must not change the annotation, only how fast it is produced.

    This is the property PR 3's efficiency work depends on, and it is not
    implied by the goldens: those all run at one batch size, so a batching bug
    in the retained-slice bookkeeping would pin itself happily. Asserting
    against the committed golden rather than against another run also means a
    batch size that changes the result cannot agree with itself and pass.
    """
    output = tmp_path / f"batch_{batch_size}.gff3"
    run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            batch_size=batch_size,
            **_GOLDENS["hmm_constrained_filtered"],
        )
    )
    assert _feature_lines(output.read_text()) == _golden("hmm_constrained_filtered")


@pytest.mark.parametrize("chunk_len", [64, 128, 512, 4096, 100_000])
def test_output_is_invariant_to_the_viterbi_block_length(
    chunk_len: int,
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """The blocked Viterbi must be exact, not an approximation of block length.

    ``decode_chunk_len`` splits the state path into blocks decoded in parallel
    and reconciled across their boundaries. If that reconciliation were
    approximate, the decoded path would depend on where the blocks fall, and
    the tuning knob would quietly be an accuracy knob. Asserting against the
    committed golden rather than against another chunk length means a chunking
    that changed the answer cannot agree with itself and pass.

    The PR claimed this invariance before anything checked it; it is checked
    now.
    """
    output = tmp_path / f"chunk_{chunk_len}.gff3"
    run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            decode_chunk_len=chunk_len,
            **_GOLDENS["hmm_constrained_filtered"],
        )
    )
    assert _feature_lines(output.read_text()) == _golden("hmm_constrained_filtered")


def test_the_viterbi_block_length_reaches_the_decoder(
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """``decode_chunk_len`` must be the block length the decoder actually uses.

    The invariance test above cannot establish this, and that is not a quibble
    about rigour: pinning the value at both call sites in ``annotate_genome``
    leaves all five of its parametrizations green, because five runs that
    secretly share one chunk length do all reproduce the golden. So the
    invariance assertion proves the decode is exact and says nothing about
    whether the knob is wired up.

    The general form, worth stating once: **an invariance test over a parameter
    cannot also prove the parameter is read.** Every such parameter needs a
    second, separate test that observes the value arriving somewhere. This is
    the same gap that was found and fixed for ``batch_size``, then
    reintroduced one commit later in the test written to close it.

    So this spies on the block length handed to ``viterbi_jax_parallel``, which
    is where the setting is finally consumed.
    """
    from nucleotide_transformer_v3.gene_annotation import annotate_utils

    real = annotate_utils.viterbi_jax_parallel
    for chunk_len in (64, 4096):
        seen: list[int] = []

        def _spy(emissions, block_len, *rest, _sink=seen, _real=real, **kwargs):
            _sink.append(int(block_len))
            return _real(emissions, block_len, *rest, **kwargs)

        annotate_utils.viterbi_jax_parallel = _spy  # type: ignore[assignment]
        try:
            run(
                production_config(
                    fasta_paths=synthetic_genome,
                    output_path=str(tmp_path / f"plumbed_{chunk_len}.gff3"),
                    decode_chunk_len=chunk_len,
                    **_GOLDENS["hmm_constrained_filtered"],
                )
            )
        finally:
            annotate_utils.viterbi_jax_parallel = real  # type: ignore[assignment]

        assert seen, (
            f"decode_chunk_len={chunk_len}: viterbi_jax_parallel was never "
            "called, so this test cannot observe the block length at all."
        )
        assert set(seen) == {chunk_len}, (
            f"decode_chunk_len={chunk_len} but the decoder was handed "
            f"{sorted(set(seen))}. The setting is not reaching the decoder, so "
            "the invariance test above is comparing runs that share one block "
            "length."
        )


def test_a_gene_truncated_at_a_run_start_gets_the_stop_codon_phase_right() -> None:
    """A leading ``stop_codon`` is the third base of its codon, so phase is 1.

    The state is reachable at index 0 of a run, where a gene is truncated by the
    contig or window boundary. Phase 0 there frameshifts the partial CDS.
    """
    import numpy as np

    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        CLASS_NAME_TO_INDEX,
    )
    from nucleotide_transformer_v3.gene_annotation.gff_utils import find_genes

    states = np.array(
        [CLASS_NAME_TO_INDEX["stop_codon"]] + [CLASS_NAME_TO_INDEX["intergenic"]] * 5
    )
    lines, _ = find_genes(states, "plus", 0, chrom="chr1", genomic_start=0)

    cds = [line for line in lines if line.split("\t")[2] == "CDS"]
    assert len(cds) == 1
    assert cds[0].split("\t")[7] == "1"
