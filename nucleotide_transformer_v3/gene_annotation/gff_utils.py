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

"""GFF3 writer for 15-class genome-annotation state paths."""

import numpy as np

from nucleotide_transformer_v3.gene_annotation.constants import (
    CLASS_NAME_TO_INDEX,
    GENOME_ANNOTATION_CLASSES,
    NUC_A,
    NUC_G,
    NUC_N,
    NUC_T,
    STRAND_NAMES,
)
from nucleotide_transformer_v3.logging import get_logger

logger = get_logger()

_STRAND_TO_GFF3: dict[str, str] = dict(zip(STRAND_NAMES, ("+", "-"), strict=True))

_INTERGENIC_IDX = CLASS_NAME_TO_INDEX["intergenic"]


# GFF3 feature type per class index. Order matches GENOME_ANNOTATION_CLASSES;
# None for intergenic so every state index is covered.
_GFF_TYPE_BY_INDEX: list[str | None] = [
    None,  # 0  intergenic
    "intron",  # 1  intron_p0
    "intron",  # 2  intron_p1
    "intron",  # 3  intron_p2
    "CDS",  # 4  CDS_p0
    "CDS",  # 5  CDS_p1
    "CDS",  # 6  CDS_p2
    "CDS",  # 7  splice_acceptor_p0
    "CDS",  # 8  splice_acceptor_p1
    "CDS",  # 9  splice_acceptor_p2
    "CDS",  # 10 splice_donor_p0
    "CDS",  # 11 splice_donor_p1
    "CDS",  # 12 splice_donor_p2
    "CDS",  # 13 start_codon
    "CDS",  # 14 stop_codon
]
assert len(_GFF_TYPE_BY_INDEX) == len(GENOME_ANNOTATION_CLASSES)

_GFF_TYPE_NAMES: list[str] = sorted({t for t in _GFF_TYPE_BY_INDEX if t})

# The same mapping as integer codes, so a whole state path can be turned into
# feature types with one vectorised LUT index. Only the equality pattern matters
# (0 = intergenic); find_genes reads the type name back out of the list above.
_GFF_TYPE_CODE_BY_INDEX: np.ndarray = np.array(
    [0 if t is None else _GFF_TYPE_NAMES.index(t) + 1 for t in _GFF_TYPE_BY_INDEX],
    dtype=np.uint8,
)

# Codon position of every class that can begin OR end a CDS run.
# Plus strand reads the 5' end (run start); minus strand reads the 3' end
# (run end after reversal). CDS_p{N} carries codon position N directly,
# splice_donor_p{N} mirrors splice_acceptor_p{N}, start_codon is the first base
# of its codon and stop_codon the third.
_CODON_POS_BY_CLASS: dict[str, int] = {
    "splice_acceptor_p0": 0,
    "splice_acceptor_p1": 1,
    "splice_acceptor_p2": 2,
    "splice_donor_p0": 0,
    "splice_donor_p1": 1,
    "splice_donor_p2": 2,
    "CDS_p0": 0,
    "CDS_p1": 1,
    "CDS_p2": 2,
    "start_codon": 0,
    "stop_codon": 2,  # The only predecessor of stop_codon is CDS_p1
}
_CODON_POS_BY_INDEX: dict[int, int] = {
    CLASS_NAME_TO_INDEX[name]: pos for name, pos in _CODON_POS_BY_CLASS.items()
}

# Class-index ranges used by the in-frame-stop filter. Introns are the three
# intron_p{0,1,2} classes; the spliced CDS is everything non-intergenic minus
# the introns.
_INTRON_IDX_MIN = CLASS_NAME_TO_INDEX["intron_p0"]
_INTRON_IDX_MAX = CLASS_NAME_TO_INDEX["intron_p2"]

# Stop codons as base-code triples — TAA, TAG, TGA.
_STOP_TRIPLES: frozenset[tuple[int, int, int]] = frozenset(
    {(NUC_T, NUC_A, NUC_A), (NUC_T, NUC_A, NUC_G), (NUC_T, NUC_G, NUC_A)}
)


def _rc_base_codes_np(arr: np.ndarray) -> np.ndarray:
    """Reverse-complement an ACGTN-encoded numpy array (A=0,C=1,G=2,T=3,N=4)."""
    return np.where(arr == NUC_N, NUC_N, NUC_T - arr)[::-1]


def _has_internal_inframe_stop(
    is_intron: np.ndarray, gene_base_codes: np.ndarray, strand: str
) -> bool:
    """Return True if the gene's spliced CDS contains an internal in-frame stop.

    The spliced CDS is collected from all non-intron positions of the gene
    (``is_intron`` is the per-position intron mask the caller already computed).
    For minus-strand genes it is reverse-complemented to transcribed direction
    before codon scanning. The terminal stop codon is excluded from the check.
    """
    cds_codes = gene_base_codes[~is_intron]
    if strand == "minus":
        cds_codes = _rc_base_codes_np(cds_codes)
    n_full_codons = len(cds_codes) // 3
    if n_full_codons < 2:
        return False
    for i in range(0, (n_full_codons - 1) * 3, 3):
        triple = (int(cds_codes[i]), int(cds_codes[i + 1]), int(cds_codes[i + 2]))
        if triple in _STOP_TRIPLES:
            return True
    return False


def find_genes(
    states: np.ndarray,
    strand: str,
    gene_id_offset: int,
    *,
    chrom: str,
    genomic_start: int,
    base_codes: np.ndarray | None = None,
    min_cds_length: int = 0,
) -> tuple[list[str], int]:
    """State path -> GFF3 gene + sub-feature lines for one strand.

    A gene is any contiguous run of non-intergenic positions; sub-features
    merge by GFF3 type. CDS phase is derived from the codon position of the
    starting class (or the ending class on minus strand).

    Two post-decode filters are applied when their inputs are provided: genes
    with total CDS < ``min_cds_length`` are dropped, and (when ``base_codes``
    is provided) genes whose spliced CDS contains an internal in-frame stop
    codon are also dropped.

    Args:
        states: 1-D integer array of class indices (matching
            GENOME_ANNOTATION_CLASSES) for every base-pair position on one
            strand.
        strand: "plus" or "minus" (STRAND_NAMES convention). Controls the
            GFF3 strand column and which end of a CDS run determines the
            phase.
        gene_id_offset: Counter seed for gene numbering. The first gene
            discovered is numbered gene_id_offset + 1, allowing sequential
            calls (e.g. plus then minus) to produce unique IDs.
        chrom: Chromosome name emitted in the GFF3 seqid column.
        genomic_start: 0-based genomic offset of states[0]; coordinates are
            converted to 1-based closed intervals for GFF3.
        base_codes: Optional uint8 array aligned to ``states`` carrying per-base
            codes (A=0, C=1, G=2, T=3, N=4). When provided, the in-frame-stop
            filter drops any gene whose spliced CDS has an internal in-frame
            stop codon.
        min_cds_length: Drop any gene whose total CDS length is strictly less
            than this. Default 0 disables the filter.

    Returns:
        Tuple of (lines, last_gene_id) where lines is a list of tab-separated
        GFF3 feature strings (gene rows followed by their sub-features) and
        last_gene_id is the highest gene counter assigned, suitable as
        gene_id_offset for a subsequent call.
    """
    lines: list[str] = []
    gene_id = gene_id_offset
    strand_col = _STRAND_TO_GFF3[strand]
    strand_tag = strand[0]  # "plus" -> "p", "minus" -> "m"
    states = np.asarray(states)
    if states.size == 0:
        return lines, gene_id

    # Use vectorised diffs to find the gene partition: maximal runs of
    # non-intergenic positions. The zero pad at each end lets a run touching a
    # boundary produce its +1/-1 edge like any other.
    non_intergenic = (states != _INTERGENIC_IDX).astype(np.int8)
    edges = np.diff(non_intergenic, prepend=np.int8(0), append=np.int8(0))
    gene_starts = np.flatnonzero(edges == 1)
    if gene_starts.size == 0:
        return lines, gene_id
    gene_ends = np.flatnonzero(edges == -1)

    # Same for the sub-feature partition: maximal runs of equal GFF3 *type*, so the
    # boundaries come from the type codes rather than the states (several states
    # share a type). Computed once for the whole path; each gene slices its own.
    # Comparing the two slices keeps the mask boolean, which lets flatnonzero take
    # its zero-skipping path — an arithmetic diff would hand it a uint8 array and
    # cost ~10x more on this step.
    type_codes = _GFF_TYPE_CODE_BY_INDEX[states]
    type_changes = np.flatnonzero(type_codes[1:] != type_codes[:-1]) + 1

    for gene_start, gene_end in zip(
        gene_starts.tolist(), gene_ends.tolist(), strict=True
    ):
        if min_cds_length > 0 or base_codes is not None:
            gene_states = states[gene_start:gene_end]
            is_intron = (gene_states >= _INTRON_IDX_MIN) & (
                gene_states <= _INTRON_IDX_MAX
            )
            cds_length = int((~is_intron).sum())
            if min_cds_length > 0 and cds_length < min_cds_length:
                continue
            if base_codes is not None and _has_internal_inframe_stop(
                is_intron, base_codes[gene_start:gene_end], strand
            ):
                continue

        gene_id += 1
        gid = f"gene_{gene_id}_{strand_tag}"
        gene_start_coord = genomic_start + gene_start + 1  # GFF3 is 1-based
        gene_end_coord = genomic_start + gene_end
        lines.append(
            f"{chrom}\tNTv3_HMM\tgene\t{gene_start_coord}\t{gene_end_coord}"
            f"\t.\t{strand_col}\t.\tID={gid}"
        )

        # Type changes strictly inside (gene_start, gene_end) — both gene ends are
        # themselves changes, to/from intergenic, hence the exclusive sides.
        lo = np.searchsorted(type_changes, gene_start, "right")
        hi = np.searchsorted(type_changes, gene_end, "left")
        bounds = [gene_start, *type_changes[lo:hi].tolist(), gene_end]
        for i in range(len(bounds) - 1):
            feat_pos, feat_end = bounds[i], bounds[i + 1]
            feat_state = int(states[feat_pos])
            # Every position in a gene run is non-intergenic, so never None.
            out_type = _GFF_TYPE_BY_INDEX[feat_state]
            assert out_type is not None
            feat_start_coord = genomic_start + feat_pos + 1
            feat_end_coord = genomic_start + feat_end
            feat_id = f"{gid}_{out_type}_{feat_start_coord}"

            # GFF3 phase column: defined for CDS, '.' otherwise.
            # GFF phase = (3 - codon_position) % 3. Every CDS-class state is
            # in _CODON_POS_BY_INDEX, so the direct lookup is safe and fails
            # loudly if a new CDS class is added without updating the table.
            phase_col = "."
            if out_type == "CDS":
                # On the minus strand the 5' end of the CDS is the feature end.
                frame_state = (
                    int(states[feat_end - 1]) if strand == "minus" else feat_state
                )
                codon_pos = _CODON_POS_BY_INDEX[frame_state]
                phase_col = str((3 - codon_pos) % 3)
            lines.append(
                f"{chrom}\tNTv3_HMM\t{out_type}\t{feat_start_coord}\t{feat_end_coord}"
                f"\t.\t{strand_col}\t{phase_col}\tID={feat_id};Parent={gid}"
            )

    return lines, gene_id
