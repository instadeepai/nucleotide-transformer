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

"""Output-format selection, and the GTF re-expression of production's GFF3.

Two formats are offered because the gene-annotation community reads two, and a
tool that emits only one gets wrapped in somebody's conversion script. Which one
is written is inferred from ``output_path``'s extension, as Tiberius does:
``.gff3`` / ``.gff`` give GFF3 and ``.gtf`` gives GTF. Anything else is rejected
rather than guessed at, so a typo cannot silently produce a file whose contents
do not match its name.

GFF3 stays the canonical intermediate. ``find_genes`` emits GFF3 lines, the
per-record part files hold GFF3 lines, and GTF is produced from them at merge
time. So the GTF is by construction a re-expression of exactly the same
intervals, and the GFF3 path is untouched by this module existing.

**Which GTF dialect.** The one Ensembl and GENCODE write, sometimes called
GTF2.5: the nine tab-separated columns of GFF with a GTF attribute column, plus
the ``gene`` and ``transcript`` feature rows that GTF2.2 itself does not define
but that every modern consumer (gffread, StringTie, featureCounts, IGV) expects
to group a transcript's exons. Specifically:

* Column 8 is ``frame``, and GTF's frame is GFF3's phase: "0" means the first
  base of the feature is the first base of a codon. The number is carried over
  from the GFF3 line unchanged, not recomputed.
* Column 9 follows the GTF2.2 attribute rule -- ``gene_id`` first,
  ``transcript_id`` second, each attribute ``key "value"`` ending in a semicolon
  and separated from the next by exactly one space, values double-quoted.
* One transcript per gene, ``<gene_id>.t1``. The decoder emits a single state
  path per strand, so there are no isoforms to name; inventing an isoform
  numbering would claim more than the model predicts.
* ``exon_number`` counts from the 5' end in transcribed order, so a
  minus-strand gene numbers from its highest coordinate.

**What is deliberately not emitted.**

* ``start_codon`` and ``stop_codon`` rows. GTF2.2 requires them of a complete
  gene, but the decoder's 15-state path has already been merged into maximal
  runs of one GFF3 type by ``find_genes``, so the start and stop positions are
  no longer separable from the CDS runs that contain them. Recovering them here
  would mean re-deriving from sequence what the GFF3 does not record, i.e.
  asserting a boundary the annotation itself does not carry.
* ``intron`` rows, which the GFF3 does carry. GTF2.2 has no ``intron`` feature
  type (only ``intron_CNS``, which means something else), and in GTF an intron
  is the gap between consecutive exons of a transcript. No information is lost.
* UTR rows (``5UTR`` / ``3UTR`` / ``UTR``), ``Selenocysteine``, ``inter`` and
  ``inter_CNS``. The model predicts none of them.
* A score. The GFF3 has none either; ``.`` is written through.

**One dialect difference worth stating.** The CDS rows include the terminal stop
codon, because the GFF3 does: the decoder's ``stop_codon`` state maps to GFF3
type ``CDS``, and the two outputs must describe the same intervals. GTF2.2 and
Ensembl exclude the stop codon from the CDS; AUGUSTUS and Tiberius include it by
default. A consumer that assumes the Ensembl convention will read the last three
CDS bases as a trailing codon rather than as the stop.
"""

FORMAT_GFF3 = "gff3"
FORMAT_GTF = "gtf"

# Output extension -> format, matching Tiberius' own extension dispatch. A
# ``.gz`` suffix is stripped before the lookup, so ``genome.gtf.gz`` is a
# gzipped GTF rather than an unknown format.
FORMAT_BY_SUFFIX: dict[str, str] = {
    ".gff3": FORMAT_GFF3,
    ".gff": FORMAT_GFF3,
    ".gtf": FORMAT_GTF,
}

# GFF3 sub-feature types that survive into GTF, and what they become. ``intron``
# is absent on purpose: see the module docstring.
_GTF_TYPE_BY_GFF3_TYPE: dict[str, str] = {"CDS": "CDS"}

# Suffix appended to a gene id to name its single transcript.
TRANSCRIPT_SUFFIX = ".t1"


def resolve_output_format(output_path: str) -> str:
    """Infer the annotation output format from a path's extension.

    Args:
        output_path: The configured output path, optionally gzipped.

    Returns:
        ``"gff3"`` or ``"gtf"``.

    Raises:
        ValueError: If the extension is not one of the accepted ones. Guessing
            would produce a file whose contents contradict its name.
    """
    # Case-insensitive, to match is_fasta_name and _open_text on the input side.
    # They accepted chr4.fna.GZ while this rejected out.gff3.GZ, so the two ends
    # of the same pipeline disagreed about what a gzip suffix looks like.
    stem = (
        output_path[: -len(".gz")]
        if output_path.lower().endswith(".gz")
        else output_path
    )
    _, dot, suffix = stem.rpartition(".")
    resolved = FORMAT_BY_SUFFIX.get(f".{suffix.lower()}") if dot else None
    if resolved is None:
        accepted = ", ".join(sorted(FORMAT_BY_SUFFIX))
        raise ValueError(
            f"Cannot infer an output format from output_path={output_path!r}. "
            f"The extension must be one of {accepted} (optionally followed by "
            "'.gz'), which is what selects GFF3 or GTF."
        )
    return resolved


def _attributes(*pairs: tuple[str, str]) -> str:
    """Render a GTF attribute column from ordered key/value pairs.

    Follows the GTF2.2 rule: every attribute ends in a semicolon and is
    separated from the next by exactly one space, and values are double-quoted.

    Args:
        pairs: Attribute keys and values, in the order they should appear.
            ``gene_id`` must come first and ``transcript_id`` second.

    Returns:
        The rendered column, ending in a semicolon.
    """
    return " ".join(f'{key} "{value}";' for key, value in pairs)


def parse_gff3_attributes(column: str) -> dict[str, str]:
    """Parse a GFF3 attribute column into its tag/value pairs.

    Args:
        column: The ninth column of a GFF3 feature line.

    Returns:
        Tag to value. Malformed entries (no ``=``) are dropped rather than
        raising, since they cannot name a gene either way.
    """
    parsed: dict[str, str] = {}
    for item in column.split(";"):
        tag, sep, value = item.partition("=")
        if sep:
            parsed[tag.strip()] = value
    return parsed


def gff3_lines_to_gtf(lines: list[str]) -> list[str]:
    """Re-express one record's GFF3 feature lines as GTF.

    Genes are recognised by their ``gene`` row and their sub-features by
    ``Parent``, never by adjacency: ``write_record_part`` sorts globally by
    start, so a gene's rows are routinely interleaved with those of a gene
    overlapping it on the other strand.

    Each gene yields a ``gene`` row, one ``transcript`` row, and an ``exon`` /
    ``CDS`` pair per CDS block. The output is sorted by start with a stable sort,
    so a gene row still precedes its transcript row and both precede the
    sub-features that share their start, exactly as the GFF3 orders them.

    Args:
        lines: GFF3 feature lines for one record, as ``find_genes`` emits them.
            Comment lines must already be excluded.

    Returns:
        The GTF lines, sorted by start.

    Raises:
        ValueError: If a sub-feature names a ``Parent`` with no ``gene`` row,
            which would produce a transcript belonging to nothing.
    """
    gene_rows: dict[str, list[str]] = {}
    gene_order: list[str] = []
    cds_by_gene: dict[str, list[list[str]]] = {}

    for line in lines:
        cols = line.split("\t")
        attributes = parse_gff3_attributes(cols[8])
        if cols[2] == "gene":
            gene_id = attributes["ID"]
            gene_rows[gene_id] = cols
            gene_order.append(gene_id)
            cds_by_gene.setdefault(gene_id, [])
        elif cols[2] in _GTF_TYPE_BY_GFF3_TYPE:
            cds_by_gene.setdefault(attributes["Parent"], []).append(cols)

    orphans = sorted(set(cds_by_gene) - set(gene_rows))
    if orphans:
        raise ValueError(
            f"{len(orphans)} CDS feature(s) name a Parent with no gene row "
            f"(first 5: {orphans[:5]}). The GTF would carry a transcript "
            "belonging to no gene."
        )

    out: list[str] = []
    for gene_id in gene_order:
        # A gene with no CDS block is skipped, exactly as `sequences.py` skips it
        # when building the CDS and protein FASTAs (`if not coding: continue`).
        # Emitting it produced a `transcript` row with no exon and no CDS under
        # it, which is malformed GTF that gffread, StringTie and agat drop or
        # error on, and it made the GTF and the FASTAs disagree on gene count.
        #
        # Not reachable at the shipped `min_cds_length: 201`, which filters such
        # a gene out first. It is reachable at `min_cds_length=0`, which is
        # `annotate/`'s default and what `run_equivalence_experiment.sh` runs
        # under, so it is reachable in the arm the equivalence numbers come from.
        if not cds_by_gene[gene_id]:
            continue
        gene_cols = gene_rows[gene_id]
        transcript_id = f"{gene_id}{TRANSCRIPT_SUFFIX}"
        seqid, source, _, start, end, score, strand, _, _ = gene_cols
        prefix = f"{seqid}\t{source}"
        span = f"{start}\t{end}\t{score}\t{strand}"
        out.append(f"{prefix}\tgene\t{span}\t.\t{_attributes(('gene_id', gene_id))}")
        out.append(
            f"{prefix}\ttranscript\t{span}\t.\t"
            f"{_attributes(('gene_id', gene_id), ('transcript_id', transcript_id))}"
        )
        # Transcribed order, so exon 1 is the 5'-most block on either strand.
        blocks = sorted(cds_by_gene[gene_id], key=lambda cols: int(cols[3]))
        if strand == "-":
            blocks.reverse()
        for exon_number, cols in enumerate(blocks, start=1):
            attributes = _attributes(
                ("gene_id", gene_id),
                ("transcript_id", transcript_id),
                ("exon_number", str(exon_number)),
            )
            block = f"{cols[3]}\t{cols[4]}\t{cols[5]}\t{cols[6]}"
            out.append(f"{prefix}\texon\t{block}\t.\t{attributes}")
            gtf_type = _GTF_TYPE_BY_GFF3_TYPE[cols[2]]
            # GTF's frame column and GFF3's phase column mean the same thing,
            # so it is carried over rather than recomputed.
            out.append(f"{prefix}\t{gtf_type}\t{block}\t{cols[7]}\t{attributes}")

    return sorted(out, key=lambda line: int(line.split("\t")[3]))
