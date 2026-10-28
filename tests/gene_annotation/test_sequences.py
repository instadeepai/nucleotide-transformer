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

"""Tests for the protein and CDS FASTA path, from one codon to a merged file."""

import gzip
from pathlib import Path

import pytest

from nucleotide_transformer_v3.gene_annotation.sequences import (
    STOP_SYMBOL,
    UNKNOWN_SYMBOL,
    GeneSequence,
    gene_sequences,
    merge_fasta_parts,
    render_fasta_records,
    spliced_coding_sequence,
    translate,
    translate_codon,
    write_fasta_part,
)

# The standard genetic code (NCBI translation table 1), written as one codon
# list per amino acid rather than as the module's positional string, so the two
# encodings share no data and a corrupted entry cannot agree with itself.
_STANDARD_CODE: dict[str, tuple[str, ...]] = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "N": ("AAT", "AAC"),
    "D": ("GAT", "GAC"),
    "C": ("TGT", "TGC"),
    "Q": ("CAA", "CAG"),
    "E": ("GAA", "GAG"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "K": ("AAA", "AAG"),
    "M": ("ATG",),
    "F": ("TTT", "TTC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "*": ("TAA", "TAG", "TGA"),
}

_EXPECTED_AMINO_ACID: dict[str, str] = {
    codon: amino_acid
    for amino_acid, codons in _STANDARD_CODE.items()
    for codon in codons
}

_STOP_CODONS = ("TAA", "TAG", "TGA")

# Complement table of this test's own, so the minus-strand expectations are not
# derived from the helper they are checking.
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# One synthetic two-exon gene on a 42 bp contig, laid out so every coordinate
# below is readable: 1-based inclusive, as GFF3 counts.
#
#   1-4    TTTT              intergenic
#   5-13   ATGAAACCC         CDS block 1
#   14-29  GTAAGTTTTTTTTTAG  intron, real gt..ag boundaries
#   30-38  GGGTGGTAA         CDS block 2
#   39-42  CCCC              intergenic
_CONTIG = "TTTT" "ATGAAACCC" "GTAAGTTTTTTTTTAG" "GGGTGGTAA" "CCCC"
_EXON_ONE = (5, 13, 0)
_EXON_TWO = (30, 38, 0)
_GENE_START, _GENE_END = 5, 38
_SPLICED = "ATGAAACCCGGGTGGTAA"
_PROTEIN = "MKPGW*"


def _revcomp(sequence: str) -> str:
    """Reverse-complement a plain ACGT sequence, without the module's helper."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _mirror(start: int, end: int, phase: int = 0) -> tuple[int, int, int]:
    """The coordinates of one block once ``_CONTIG`` is reverse-complemented.

    Args:
        start: 1-based inclusive start on the plus-strand contig.
        end: 1-based inclusive end on the plus-strand contig.
        phase: GFF3 phase to carry over unchanged.

    Returns:
        The block as ``(start, end, phase)`` on the reverse-complemented contig.
    """
    length = len(_CONTIG)
    return (length - end + 1, length - start + 1, phase)


def _feature(
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    attributes: str,
    phase: str = ".",
    seqid: str = "chrTest",
) -> str:
    """One GFF3 feature line, shaped exactly as the pipeline writes them."""
    return "\t".join(
        [
            seqid,
            "NTv3_HMM",
            feature_type,
            str(start),
            str(end),
            ".",
            strand,
            phase,
            attributes,
        ]
    )


def _gene_lines(strand: str, blocks: list[tuple[int, int, int]]) -> list[str]:
    """A gene row plus its CDS rows, sorted by start as the writer sorts them."""
    gene = _feature("gene", _GENE_START, _GENE_END, strand, "ID=gene_1")
    cds = [
        _feature(
            "CDS",
            start,
            end,
            strand,
            f"ID=gene_1_CDS_{start};Parent=gene_1",
            phase=str(phase),
        )
        for start, end, phase in blocks
    ]
    return sorted([gene, *cds], key=lambda line: int(line.split("\t")[3]))


# ---------------------------------------------------------------------------
# The genetic code
# ---------------------------------------------------------------------------


def test_the_reference_genetic_code_is_complete_and_unambiguous() -> None:
    """The reference lists all 64 codons exactly once, since it drives the cases."""
    listed = [codon for codons in _STANDARD_CODE.values() for codon in codons]
    all_codons = {
        f"{first}{second}{third}"
        for first in "ACGT"
        for second in "ACGT"
        for third in "ACGT"
    }
    assert len(listed) == 64, f"reference lists {len(listed)} codons, not 64"
    assert set(listed) == all_codons
    assert len(set(listed)) == 64, "a codon is listed under two amino acids"


def test_the_reference_genetic_code_agrees_with_biopython() -> None:
    """The hand-written reference matches Biopython's table 1 (a test-only dep)."""
    bio_data = pytest.importorskip("Bio.Data.CodonTable")
    table = bio_data.unambiguous_dna_by_id[1]
    oracle = dict(table.forward_table)
    for codon in table.stop_codons:
        oracle[codon] = STOP_SYMBOL
    assert oracle == _EXPECTED_AMINO_ACID


@pytest.mark.parametrize("codon", sorted(_EXPECTED_AMINO_ACID))
def test_translate_codon_matches_the_standard_genetic_code(codon: str) -> None:
    """Every one of the 64 codons translates to its published amino acid."""
    assert translate_codon(codon) == _EXPECTED_AMINO_ACID[codon]


@pytest.mark.parametrize("codon", _STOP_CODONS)
def test_a_stop_codon_translates_to_the_stop_symbol(codon: str) -> None:
    """The three stops translate to the stop symbol."""
    assert translate_codon(codon) == STOP_SYMBOL


@pytest.mark.parametrize("codon", sorted(_EXPECTED_AMINO_ACID))
def test_translate_codon_reads_a_soft_masked_codon(codon: str) -> None:
    """Lower case is soft-masking, not an unknown codon."""
    assert translate_codon(codon.lower()) == _EXPECTED_AMINO_ACID[codon]


@pytest.mark.parametrize(
    ("codon", "expected"),
    [
        ("GGN", "G"),  # every resolution is glycine
        ("TRA", STOP_SYMBOL),  # TAA and TGA, both stops
        ("CAR", "Q"),
        ("RAT", "B"),  # Asn or Asp, which IUPAC calls B
        ("SAA", "Z"),  # Gln or Glu, which IUPAC calls Z
        ("MTA", "J"),  # Ile or Leu, which IUPAC calls J
        ("ATN", UNKNOWN_SYMBOL),  # Ile or Met, no letter for the pair
        ("TAM", UNKNOWN_SYMBOL),  # a possible stop is not a stop
        ("NNN", UNKNOWN_SYMBOL),
        ("AT-", UNKNOWN_SYMBOL),  # outside the IUPAC alphabet
        ("AT", UNKNOWN_SYMBOL),  # not a codon at all
    ],
)
def test_translate_codon_resolves_an_ambiguous_codon(codon: str, expected: str) -> None:
    """An ambiguity code is expanded; a possible stop is not a certain one."""
    assert translate_codon(codon) == expected


# ---------------------------------------------------------------------------
# translate
# ---------------------------------------------------------------------------


def test_translate_reads_a_spliced_coding_sequence_in_frame() -> None:
    """A whole CDS translates codon by codon, terminal stop included."""
    assert translate(_SPLICED) == _PROTEIN


def test_translate_reads_an_internal_stop_rather_than_halting() -> None:
    """An internal stop is kept as ``*``; the biological filters act on it."""
    assert translate("ATGTGAAAATAA") == "M*K*"


@pytest.mark.parametrize("trailing", ["A", "AA"])
def test_translate_drops_a_trailing_partial_codon(trailing: str) -> None:
    """One or two leftover bases are dropped rather than guessed at."""
    assert translate(_SPLICED + trailing) == _PROTEIN


def test_translate_of_an_empty_sequence_is_empty() -> None:
    """No coding sequence means no protein, not a one-residue one."""
    assert translate("") == ""


def test_translate_reads_a_soft_masked_coding_sequence() -> None:
    """Soft-mask case must not change the protein."""
    assert translate(_SPLICED.lower()) == _PROTEIN


# ---------------------------------------------------------------------------
# spliced_coding_sequence
# ---------------------------------------------------------------------------


def test_spliced_coding_sequence_removes_the_intron() -> None:
    """The two CDS blocks are concatenated and the intron between them is gone."""
    spliced = spliced_coding_sequence([_EXON_ONE, _EXON_TWO], "+", _CONTIG)
    assert spliced == _SPLICED
    assert "GTAAGT" not in spliced, "the intron survived the splice"


def test_spliced_coding_sequence_ignores_the_order_blocks_arrive_in() -> None:
    """Blocks are sorted here, so a caller's ordering cannot frameshift a gene."""
    assert spliced_coding_sequence([_EXON_TWO, _EXON_ONE], "+", _CONTIG) == _SPLICED


@pytest.mark.parametrize("phase", [0, 1, 2])
def test_spliced_coding_sequence_trims_the_phase_on_the_plus_strand(
    phase: int,
) -> None:
    """The first block's phase is removed from the front of the spliced CDS."""
    start, end, _ = _EXON_ONE
    blocks = [(start, end, phase), _EXON_TWO]
    assert spliced_coding_sequence(blocks, "+", _CONTIG) == _SPLICED[phase:]


def test_spliced_coding_sequence_reads_a_minus_strand_gene_in_transcribed_order() -> (
    None
):
    """A minus-strand gene splices 3' to 5' genomically and reverse-complements."""
    contig = _revcomp(_CONTIG)
    blocks = [_mirror(*_EXON_ONE[:2]), _mirror(*_EXON_TWO[:2])]
    spliced = spliced_coding_sequence(blocks, "-", contig)
    assert spliced == _SPLICED
    assert translate(spliced) == _PROTEIN


@pytest.mark.parametrize("phase", [0, 1, 2])
def test_spliced_coding_sequence_trims_the_phase_of_the_first_transcribed_block(
    phase: int,
) -> None:
    """On the minus strand that is the *last* block in genomic order.

    The blocks carry different phases, so reading the wrong one fails.
    """
    contig = _revcomp(_CONTIG)
    first_transcribed = _mirror(*_EXON_ONE[:2], phase)
    second_transcribed = _mirror(*_EXON_TWO[:2], (phase + 1) % 3)
    assert first_transcribed[0] > second_transcribed[0], "mirroring went the wrong way"
    blocks = [second_transcribed, first_transcribed]
    assert spliced_coding_sequence(blocks, "-", contig) == _SPLICED[phase:]


def test_spliced_coding_sequence_of_no_blocks_is_empty() -> None:
    """A gene with no CDS has no coding sequence, rather than the whole contig."""
    assert spliced_coding_sequence([], "+", _CONTIG) == ""
    assert spliced_coding_sequence([], "-", _CONTIG) == ""


@pytest.mark.parametrize("strand", [".", "?", "", "1", "plus"])
def test_spliced_coding_sequence_rejects_an_unknown_strand(strand: str) -> None:
    """An unrecognised strand must raise rather than default to the plus strand."""
    with pytest.raises(ValueError, match="strand must be"):
        spliced_coding_sequence([_EXON_ONE], strand, _CONTIG)


def test_spliced_coding_sequence_upper_cases_a_soft_masked_contig() -> None:
    """Soft-masked input yields the same upper-case CDS."""
    assert spliced_coding_sequence([_EXON_ONE, _EXON_TWO], "+", _CONTIG.lower()) == (
        _SPLICED
    )


# ---------------------------------------------------------------------------
# gene_sequences
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("strand", ["+", "-"])
def test_gene_sequences_recovers_the_cds_and_protein_of_a_synthetic_gene(
    strand: str,
) -> None:
    """End to end: GFF3 lines plus a contig in, CDS and protein out."""
    if strand == "+":
        contig, blocks = _CONTIG, [_EXON_ONE, _EXON_TWO]
    else:
        contig = _revcomp(_CONTIG)
        blocks = [_mirror(*_EXON_ONE[:2]), _mirror(*_EXON_TWO[:2])]

    genes = gene_sequences(_gene_lines(strand, blocks), contig)

    assert genes == [
        GeneSequence(
            gene_id="gene_1",
            seqid="chrTest",
            start=_GENE_START,
            end=_GENE_END,
            strand=strand,
            coding=_SPLICED,
        )
    ]
    assert translate(genes[0].coding) == _PROTEIN
    assert genes[0].description == f"chrTest:{_GENE_START}-{_GENE_END}({strand})"


def test_gene_sequences_renders_both_fasta_outputs() -> None:
    """The rendered records carry the gene's id, locus and payload."""
    genes = gene_sequences(_gene_lines("+", [_EXON_ONE, _EXON_TWO]), _CONTIG)
    coding = "".join(render_fasta_records(genes, translate_to_protein=False))
    protein = "".join(render_fasta_records(genes, translate_to_protein=True))
    assert coding == f">gene_1 chrTest:5-38(+)\n{_SPLICED}\n"
    assert protein == f">gene_1 chrTest:5-38(+)\n{_PROTEIN}\n"


def test_gene_sequences_matches_blocks_through_parent_not_adjacency() -> None:
    """Interleaved rows of two overlapping genes are matched by ``Parent``."""
    plus = _gene_lines("+", [_EXON_ONE, _EXON_TWO])
    minus = [
        _feature("gene", 6, 37, "-", "ID=gene_2"),
        _feature("CDS", 6, 11, "-", "ID=gene_2_CDS_6;Parent=gene_2", phase="0"),
    ]
    interleaved = sorted([*plus, *minus], key=lambda line: int(line.split("\t")[3]))

    genes = {gene.gene_id: gene.coding for gene in gene_sequences(interleaved, _CONTIG)}
    assert genes["gene_1"] == _SPLICED
    assert genes["gene_2"] == _revcomp(_CONTIG[5:11])


def test_gene_sequences_skips_a_gene_with_no_cds() -> None:
    """An all-intron gene contributes nothing, not an empty FASTA record."""
    lines = [
        _feature("gene", _GENE_START, _GENE_END, "+", "ID=gene_1"),
        _feature("intron", 14, 29, "+", "ID=gene_1_intron_14;Parent=gene_1"),
    ]
    assert gene_sequences(lines, _CONTIG) == []


def test_gene_sequences_preserves_the_order_the_gene_rows_appear_in() -> None:
    """Output order follows the gene rows, so the FASTA tracks the GFF3."""
    lines = [
        _feature("gene", 5, 13, "+", "ID=gene_b"),
        _feature("CDS", 5, 13, "+", "ID=gene_b_CDS_5;Parent=gene_b", phase="0"),
        _feature("gene", 30, 38, "+", "ID=gene_a"),
        _feature("CDS", 30, 38, "+", "ID=gene_a_CDS_30;Parent=gene_a", phase="0"),
    ]
    assert [gene.gene_id for gene in gene_sequences(lines, _CONTIG)] == [
        "gene_b",
        "gene_a",
    ]


def test_gene_sequences_rejects_a_cds_whose_parent_has_no_gene_row() -> None:
    """An orphan CDS must raise; its sequence would otherwise vanish silently."""
    lines = [
        _feature("gene", _GENE_START, _GENE_END, "+", "ID=gene_1"),
        _feature("CDS", 5, 13, "+", "ID=x;Parent=gene_missing", phase="0"),
    ]
    with pytest.raises(ValueError, match="gene_missing"):
        gene_sequences(lines, _CONTIG)


# ---------------------------------------------------------------------------
# merge_fasta_parts
# ---------------------------------------------------------------------------


def _record(record_id: str) -> str:
    """One trivial FASTA record, enough to be counted and ordered."""
    return f">{record_id}_gene chr:1-9(+)\nATGAAATAA\n"


def test_merge_fasta_parts_concatenates_the_parts_in_record_id_order(
    tmp_path: Path, parts_dir: Path
) -> None:
    """Record order follows the record ids, as in the GFF3 merge."""
    for record_id in ("chrC", "chrA", "chrB"):
        write_fasta_part(parts_dir, record_id, [_record(record_id)])

    output = tmp_path / "merged.fa"
    written = merge_fasta_parts(output, parts_dir, ["chrC", "chrA", "chrB"])

    assert written == 3
    headers = [
        line for line in output.read_text().splitlines() if line.startswith(">")
    ]
    assert headers == [f">{name}_gene chr:1-9(+)" for name in ("chrA", "chrB", "chrC")]


def test_merge_fasta_parts_keeps_every_record_of_a_multi_gene_part(
    tmp_path: Path, parts_dir: Path
) -> None:
    """A part holding several genes contributes all of them, counted correctly."""
    write_fasta_part(parts_dir, "chrA", [_record("a1"), _record("a2")])
    write_fasta_part(parts_dir, "chrB", [_record("b1")])

    output = tmp_path / "merged.fa"
    assert merge_fasta_parts(output, parts_dir, ["chrA", "chrB"]) == 3
    assert output.read_text().count(">") == 3


def test_merge_fasta_parts_accepts_a_record_that_produced_nothing(
    tmp_path: Path, parts_dir: Path
) -> None:
    """A gene-free contig's empty part counts as finished and adds no record."""
    write_fasta_part(parts_dir, "chrA", [_record("a1")])
    empty = write_fasta_part(parts_dir, "chrEmpty", [])

    assert empty.exists() and empty.read_text() == ""
    output = tmp_path / "merged.fa"
    assert merge_fasta_parts(output, parts_dir, ["chrA", "chrEmpty"]) == 1
    assert output.read_text() == _record("a1")


def test_merge_fasta_parts_raises_when_a_part_is_missing(
    tmp_path: Path, parts_dir: Path
) -> None:
    """A record with no part means the FASTA would be silently incomplete."""
    write_fasta_part(parts_dir, "chrA", [_record("a1")])
    output = tmp_path / "merged.fa"

    with pytest.raises(RuntimeError, match="chrGone"):
        merge_fasta_parts(output, parts_dir, ["chrA", "chrGone"])
    assert not output.exists(), "a partial FASTA was written before the failure"


def test_merge_fasta_parts_ignores_a_part_left_by_an_earlier_job(
    tmp_path: Path, parts_dir: Path
) -> None:
    """Only this job's records are merged, so stale work cannot leak in."""
    write_fasta_part(parts_dir, "chrA", [_record("a1")])
    write_fasta_part(parts_dir, "chrStale", [_record("stale")])

    output = tmp_path / "merged.fa"
    assert merge_fasta_parts(output, parts_dir, ["chrA"]) == 1
    assert "stale" not in output.read_text()


@pytest.mark.parametrize("name", ["merged.fa.gz", "merged.FA.GZ"])
def test_merge_fasta_parts_gzips_a_gz_output(
    name: str, tmp_path: Path, parts_dir: Path
) -> None:
    """The ``.gz`` branch is chosen case-insensitively, as for the GFF3 writer."""
    write_fasta_part(parts_dir, "chrA", [_record("a1")])
    output = tmp_path / name
    assert merge_fasta_parts(output, parts_dir, ["chrA"]) == 1
    assert gzip.decompress(output.read_bytes()).decode() == _record("a1")


def test_merge_fasta_parts_handles_a_record_id_needing_escaping(
    tmp_path: Path, parts_dir: Path
) -> None:
    """A seqid holding ``/`` or ``|`` round-trips through the part filename."""
    record_ids = ["scaffold/1", "gi|12345|ref"]
    for record_id in record_ids:
        write_fasta_part(parts_dir, record_id, [_record("g")])

    output = tmp_path / "merged.fa"
    assert merge_fasta_parts(output, parts_dir, record_ids) == 2
