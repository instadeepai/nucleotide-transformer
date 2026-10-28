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

"""Spliced coding sequences and their translations, written as FASTA.

The two extra artifacts a gene caller is expected to produce alongside its
annotation, matching Tiberius' ``--codingseq`` and ``--protseq``: the spliced CDS
of every predicted gene, and the protein it encodes.

Both are derived from the GFF3 lines the pipeline already produced plus the
contig they were called on, so they describe exactly the genes in the annotation
and cannot drift from it. Three things make that derivation less obvious than it
looks, and each is tested:

**Transcribed order.** A gene's CDS blocks are concatenated 5' to 3' *of the
transcript*, which for a minus-strand gene means descending genomic coordinate
and a reverse complement. Concatenating in genomic order instead would splice
the exons backwards and translate nonsense.

**Phase.** GFF3 phase is the number of bases to remove from the start of a
feature to reach its first complete codon, so the first block's phase (first in
transcribed order, i.e. the *last* genomic block on the minus strand) is trimmed
off the front of the concatenated CDS. Getting this wrong frameshifts the whole
protein rather than failing visibly.

**Ambiguity.** Genome assemblies carry ``N`` and the other IUPAC codes, and a
codon holding one is not automatically untranslatable: ``GGN`` is glycine
whichever base it hides, while ``ATN`` could be three different residues. The
translator expands the ambiguity, returns the residue when the expansions agree,
returns ``B`` / ``Z`` / ``J`` when they are one of the pairs IUPAC has a letter
for, and only then falls back to ``X``. That is what Biopython and the NCBI
ambiguous codon tables do, and the tests pin it against Biopython directly
rather than restating this module's own table.

Production may not import Biopython (see the import-boundary test), so the
genetic code lives here as the NCBI translation-table-1 string it is published
as, and the tests check all 64 codons against the oracle.
"""

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from urllib.parse import quote

from nucleotide_transformer_v3.gene_annotation.formats import parse_gff3_attributes
from nucleotide_transformer_v3.logging import get_logger

logger = get_logger()

# NCBI genetic code table 1 (the standard code), in the form NCBI publishes it:
# the amino acid of each codon, with bases enumerated TCAG for the first,
# second and third position in that nesting order.
_TRANSLATION_TABLE_1 = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
_TABLE_BASES = "TCAG"

STOP_SYMBOL = "*"
UNKNOWN_SYMBOL = "X"

# The three IUPAC ambiguous amino acids. An ambiguous codon whose resolutions
# disagree is not automatically unknown: if they are exactly one of these pairs,
# the pair has its own letter, which is what NCBI's ambiguous codon tables and
# Biopython emit. Any other disagreement, a possible-but-not-certain stop
# included, is X.
_AMBIGUOUS_AMINO_ACIDS: dict[frozenset[str], str] = {
    frozenset("DN"): "B",  # Asx
    frozenset("EQ"): "Z",  # Glx
    frozenset("IL"): "J",  # Xle
}

CODON_TABLE: dict[str, str] = {
    f"{first}{second}{third}": _TRANSLATION_TABLE_1[index]
    for index, (first, second, third) in enumerate(
        (a, b, c) for a in _TABLE_BASES for b in _TABLE_BASES for c in _TABLE_BASES
    )
}

# The DNA bases each IUPAC code stands for. A codon holding one of these is
# translated by resolving every combination it can spell.
IUPAC_BASES: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

# Complement of every IUPAC code, so a minus-strand gene keeps its ambiguity
# codes meaningful instead of having them silently pass through unflipped.
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# FASTA sequence lines are wrapped at this width, the usual convention.
FASTA_LINE_WIDTH = 60

# Suffix of a finished per-record FASTA part, mirroring the GFF3 parts: written
# by an atomic rename, so a killed job leaves nothing a resumed one mistakes for
# finished work.
FASTA_PART_SUFFIX = ".fapart"


def reverse_complement(sequence: str) -> str:
    """Reverse-complement a DNA sequence, resolving IUPAC ambiguity codes.

    Case-insensitive, and the result is upper case. ``_COMPLEMENT`` covers only
    the upper-case alphabet, so lower-case input used to pass through the
    translation untouched and come back *reversed but not complemented*:
    ``reverse_complement("acgt")`` returned ``"tgca"`` rather than ``"acgt"``.
    Inside this module that was latent, because ``spliced_coding_sequence``
    upper-cases each piece before reverse-complementing it. It is not latent for
    an external caller, which is what this module is shaped for, and soft-masked
    genomes are lower case by convention - so the obvious helper returned
    silently wrong sequence rather than an error.

    Args:
        sequence: DNA in either case, possibly with IUPAC ambiguity codes.

    Returns:
        The upper-case reverse complement. Characters outside the IUPAC alphabet
        are reversed but left as they are.
    """
    return sequence.upper().translate(_COMPLEMENT)[::-1]


@cache
def translate_codon(codon: str) -> str:
    """Translate one codon under the standard genetic code, in either case.

    An unambiguous codon is a table lookup. An ambiguous one is expanded into
    every codon it can spell: if they all encode the same residue (or all encode
    a stop) that residue is returned, and if they encode exactly one of the IUPAC
    ambiguous pairs the pair's letter is returned. Anything else is unknown. So
    ``GGN`` is ``G``, ``TRA`` is ``*``, ``RAT`` is ``B``, and ``ATN`` is ``X``,
    as is ``TAM``, which may or may not be a stop.

    Args:
        codon: Three upper-case characters.

    Returns:
        The one-letter amino acid, ``*`` for a stop, ``B`` / ``Z`` / ``J`` for
        an ambiguity NCBI has a letter for, or ``X`` when the codon is ambiguous
        enough to encode more than one thing, or holds a character outside the
        IUPAC alphabet.
    """
    # Upper-cased first, for the same reason `reverse_complement` is: both
    # tables are upper-case only, so `translate_codon("atg")` returned "X"
    # instead of "M". Latent inside this module, where the caller upper-cases,
    # and not latent for the external caller this module exists to serve, since
    # a soft-masked genome is lower case by convention.
    codon = codon.upper()
    amino_acid = CODON_TABLE.get(codon)
    if amino_acid is not None:
        return amino_acid
    options = [IUPAC_BASES.get(base, "") for base in codon]
    if len(codon) != 3 or not all(options):
        return UNKNOWN_SYMBOL
    resolved = {
        CODON_TABLE[f"{first}{second}{third}"]
        for first in options[0]
        for second in options[1]
        for third in options[2]
    }
    if len(resolved) == 1:
        return resolved.pop()
    return _AMBIGUOUS_AMINO_ACIDS.get(frozenset(resolved), UNKNOWN_SYMBOL)


def translate(coding_sequence: str) -> str:
    """Translate a spliced coding sequence under the standard genetic code.

    The sequence is read from its first base, so any phase has to have been
    trimmed already (:func:`spliced_coding_sequence` does that). A trailing
    partial codon is dropped rather than guessed at, and a stop codon is
    translated to ``*`` wherever it falls, including the terminal one that the
    CDS features carry.

    Args:
        coding_sequence: The spliced CDS, in transcribed orientation.

    Returns:
        The protein sequence, one letter per complete codon.
    """
    sequence = coding_sequence.upper()
    complete = len(sequence) - len(sequence) % 3
    return "".join(
        translate_codon(sequence[start : start + 3]) for start in range(0, complete, 3)
    )


@dataclass(frozen=True)
class GeneSequence:
    """One predicted gene's spliced CDS, with where it came from.

    Attributes:
        gene_id: The gene's GFF3 ``ID``, so the FASTA, the GFF3 and the GTF all
            name the same object.
        seqid: The contig the gene was called on.
        start: The gene's 1-based inclusive genomic start.
        end: The gene's 1-based inclusive genomic end.
        strand: ``"+"`` or ``"-"``.
        coding: The spliced CDS in transcribed orientation, upper-cased and with
            the first block's phase already trimmed.
    """

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    coding: str

    @property
    def description(self) -> str:
        """The FASTA description: where on the assembly this gene sits."""
        return f"{self.seqid}:{self.start}-{self.end}({self.strand})"


def spliced_coding_sequence(
    blocks: list[tuple[int, int, int]], strand: str, sequence: str
) -> str:
    """Splice one gene's CDS blocks into a coding sequence.

    Args:
        blocks: One ``(start, end, phase)`` per CDS feature, with GFF3's 1-based
            inclusive coordinates and GFF3's phase. Order is irrelevant; they
            are sorted here.
        strand: ``"+"`` or ``"-"``.
        sequence: The contig the blocks index into.

    Returns:
        The spliced CDS in transcribed orientation, upper-cased, with the first
        block's phase trimmed from the front. Empty when there are no blocks.

    Raises:
        ValueError: If ``strand`` is neither ``"+"`` nor ``"-"``.
    """
    if strand not in {"+", "-"}:
        raise ValueError(f"strand must be '+' or '-', got {strand!r}.")
    if not blocks:
        return ""
    # Transcribed order: ascending genomic on the plus strand, descending on the
    # minus strand, where each block is also read as its own reverse complement.
    ordered = sorted(blocks, reverse=strand == "-")
    pieces = [sequence[start - 1 : end].upper() for start, end, _ in ordered]
    if strand == "-":
        pieces = [reverse_complement(piece) for piece in pieces]
    # The first block in transcribed order carries the phase of the whole CDS.
    phase = ordered[0][2]
    return "".join(pieces)[phase:]


def gene_sequences(feature_lines: list[str], sequence: str) -> list[GeneSequence]:
    """Extract every gene's spliced CDS from one record's GFF3 lines.

    Genes are matched to their CDS blocks through ``Parent``, never by
    adjacency: the lines are sorted by start, so a gene's rows are routinely
    interleaved with those of a gene overlapping it on the other strand.

    A gene with no CDS block contributes nothing. ``find_genes`` can emit one
    when an unconstrained decode calls a run that is entirely intron, and an
    empty FASTA record would be worse than no record.

    Args:
        feature_lines: GFF3 feature lines for one record. Comment lines must
            already be excluded.
        sequence: That record's sequence, exactly as it was annotated.

    Returns:
        One :class:`GeneSequence` per gene that has a CDS, in the order the gene
        rows appear.

    Raises:
        ValueError: If a CDS feature names a ``Parent`` with no gene row.
    """
    gene_order: list[str] = []
    gene_row: dict[str, list[str]] = {}
    blocks: dict[str, list[tuple[int, int, int]]] = {}

    for line in feature_lines:
        cols = line.split("\t")
        attributes = parse_gff3_attributes(cols[8])
        if cols[2] == "gene":
            gene_id = attributes["ID"]
            gene_order.append(gene_id)
            gene_row[gene_id] = cols
            blocks.setdefault(gene_id, [])
        elif cols[2] == "CDS":
            blocks.setdefault(attributes["Parent"], []).append(
                (int(cols[3]), int(cols[4]), int(cols[7]))
            )

    orphans = sorted(set(blocks) - set(gene_row))
    if orphans:
        raise ValueError(
            f"{len(orphans)} CDS feature(s) name a Parent with no gene row "
            f"(first 5: {orphans[:5]}); their coding sequence would belong to "
            "no gene."
        )

    extracted: list[GeneSequence] = []
    for gene_id in gene_order:
        cols = gene_row[gene_id]
        coding = spliced_coding_sequence(blocks[gene_id], cols[6], sequence)
        if not coding:
            continue
        extracted.append(
            GeneSequence(
                gene_id=gene_id,
                seqid=cols[0],
                start=int(cols[3]),
                end=int(cols[4]),
                strand=cols[6],
                coding=coding,
            )
        )
    return extracted


def format_fasta_record(record_id: str, sequence: str, description: str = "") -> str:
    """Render one FASTA record, wrapped at :data:`FASTA_LINE_WIDTH`.

    Args:
        record_id: The identifier, which becomes the header's first field.
        sequence: The residues.
        description: Free text after the identifier, or empty for none.

    Returns:
        The record, header included, ending in a newline.
    """
    header = f">{record_id} {description}" if description else f">{record_id}"
    wrapped = [
        sequence[start : start + FASTA_LINE_WIDTH]
        for start in range(0, len(sequence), FASTA_LINE_WIDTH)
    ] or [""]
    return "".join(f"{line}\n" for line in [header, *wrapped])


def render_fasta_records(
    genes: list[GeneSequence], translate_to_protein: bool
) -> Iterator[str]:
    """Render one FASTA record per gene, coding or translated.

    Args:
        genes: The genes to render.
        translate_to_protein: Emit the protein rather than the coding sequence.

    Yields:
        Rendered FASTA records.
    """
    for gene in genes:
        payload = translate(gene.coding) if translate_to_protein else gene.coding
        yield format_fasta_record(gene.gene_id, payload, gene.description)


def fasta_part_path_for(parts_dir: Path, record_id: str) -> Path:
    """Path of one record's FASTA part file.

    The record id is percent-escaped for the same reason the GFF3 parts escape
    it: a FASTA seqid may legally contain ``/`` or ``|``.

    Args:
        parts_dir: Directory the FASTA parts live in.
        record_id: The record's seqid.

    Returns:
        The part file's path, which exists only once the record is finished.
    """
    return parts_dir / f"{quote(record_id, safe='')}{FASTA_PART_SUFFIX}"


def write_fasta_part(
    parts_dir: Path, record_id: str, records: Iterable[str]
) -> Path:
    """Write one record's FASTA entries to its part file.

    Mirrors ``write_record_part``: created by an atomic rename from a ``.tmp``
    sibling, and written even when the record has no gene, so a resumed run does
    not redo it forever.

    Args:
        parts_dir: Directory the FASTA parts live in.
        record_id: The record's seqid.
        records: Rendered FASTA records, each ending in a newline.

    Returns:
        The path written.
    """
    parts_dir.mkdir(parents=True, exist_ok=True)
    final = fasta_part_path_for(parts_dir, record_id)
    staging = final.parent / f"{final.name}.tmp"
    staging.write_text("".join(records))
    staging.replace(final)
    return final


def merge_fasta_parts(
    output_path: str | Path, parts_dir: Path, record_ids: list[str]
) -> int:
    """Concatenate this job's FASTA parts into one file, in record-id order.

    Ordered and filtered exactly as ``merge_parts`` orders and filters the GFF3
    parts, so the FASTA's records follow the same contig order as the
    annotation's, and a leftover part from an earlier job cannot leak in.

    Args:
        output_path: Destination FASTA; gzipped when it ends in ``.gz``.
        parts_dir: Directory holding the FASTA parts.
        record_ids: Every record this job covered, annotated or resumed.

    Returns:
        The number of FASTA records written.

    Raises:
        RuntimeError: If a named record has no part, which means the output
            would be silently incomplete.
    """
    missing = [
        record_id
        for record_id in record_ids
        if not fasta_part_path_for(parts_dir, record_id).exists()
    ]
    if missing:
        raise RuntimeError(
            f"{len(missing)} record(s) have no FASTA part file in {parts_dir} and "
            f"would be missing from {output_path}: {sorted(missing)[:5]}."
        )
    parts = sorted(
        (record_id, fasta_part_path_for(parts_dir, record_id))
        for record_id in record_ids
    )
    body = [path.read_text() for _, path in parts]
    text = "".join(body)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Case-insensitive, as for the GFF3/GTF writer.
    if str(out).lower().endswith(".gz"):
        out.write_bytes(gzip.compress(text.encode("utf-8")))
    else:
        out.write_text(text)

    num_records = sum(
        1 for chunk in body for line in chunk.splitlines() if line.startswith(">")
    )
    logger.info(f"Wrote {num_records} sequence(s) from {len(parts)} record(s) to {out}")
    return num_records
