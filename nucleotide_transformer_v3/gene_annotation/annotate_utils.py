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

"""Window geometry, the 15-state HMM and its decoders.
"""

import dataclasses
from collections.abc import Callable, Mapping
from typing import Protocol, cast
from dataclasses import dataclass, field

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx
from typing_extensions import Self

from nucleotide_transformer_v3.logging import get_logger
from nucleotide_transformer_v3.gene_annotation.constants import (
    CLASS_NAME_TO_INDEX,
    GENOME_ANNOTATION_CLASSES,
    NUC_A,
    NUC_C,
    NUC_G,
    NUC_N,
    NUC_T,
)
from nucleotide_transformer_v3.tokenizers import (
    BaseTokenizer,
    NucleotideTokenizer,
    get_ntv3_tokenizer,
)

logger = get_logger()

# ``select_windows`` returns these.
@dataclass
class GenomicInterval:
    """A genomic interval defined by a chromosome name, start and end position. The
    file assembly (optional) can be used to specify the interval's genome assembly.

    This is a base class for working with genomic intervals and can be inherited to
    create more specific interval classes.

    Args:
        chr_name: The name of the chromosome.
        start: The start position of the interval.
        end: The end position of the interval.
        is_reverse_complement: If True, sources that support reverse-complement read
            this flag and flip their output (sequence reversed/complemented, tracks
            read with reversed position axis and swapped strand columns).
        file_assembly: The file assembly of the interval (optional).

    """

    chr_name: str
    start: int
    end: int
    is_reverse_complement: bool = field(default=False, kw_only=True)

    def __str__(self) -> str:
        return f"{self.chr_name}:{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return int(self.end - self.start)

    def shrink(self, keep_center_fraction: float) -> Self:
        """Shrink a genomic interval by keeping the center fraction of the interval."""
        offset = int(self.length * (1 - keep_center_fraction) / 2)
        return dataclasses.replace(
            self,
            start=self.start + offset,
            end=self.end - offset,
        )

    def extend(self, new_length: int) -> Self:
        """Extend a genomic interval by the given length.
        The interval is extended to the new length by adding padding on both sides.
        """
        if self.length > new_length:
            raise ValueError(f"Cannot extend: {self.length=} > {new_length=}")
        if new_length % 2 != self.length % 2:
            raise ValueError(
                f"Parity mismatch: {self.length=} and {new_length=} must both be odd or both even"
            )
        offset = (new_length - self.length) // 2
        if self.start - offset < 0:
            raise ValueError(f"Cannot extend: {self.start=} - {offset=} < 0")
        return dataclasses.replace(
            self,
            start=self.start - offset,
            end=self.end + offset,
        )

    def __repr__(self) -> str:
        return f"GenomicInterval(chr_name={self.chr_name}, start={self.start}, end={self.end})"


NUM_STATES = len(GENOME_ANNOTATION_CLASSES)


# Tiberius mean lengths: intergenic ~10 kb, intron ~4.5 kb, exon ~200 bp.
# TODO: add additional transition probabilities for use in plants
_P_IR = 1 / 10_000


_P_I = 1 / 4_500


_P_E = 1 / 200


# Allowed transitions per source state; self-loops fill in to sum to 1, omitted
# destinations are zero. start_codon / stop_codon / splice_* are 1-bp states by
# construction (self-loop = 0): multi-bp ground-truth runs collapse under Viterbi
# (e.g. a 3-bp start_codon decodes as start_codon -> CDS_p1 -> CDS_p2), but CDS
# boundaries are preserved.
_TRANSITIONS: dict[str, dict[str, float]] = {
    "intergenic": {"start_codon": _P_IR},
    "intron_p0": {"splice_acceptor_p1": _P_I},
    "intron_p1": {"splice_acceptor_p2": _P_I},
    "intron_p2": {"splice_acceptor_p0": _P_I},
    "CDS_p0": {"CDS_p1": 1 - _P_E, "splice_donor_p1": _P_E},
    "CDS_p1": {"CDS_p2": 1 - 2 * _P_E, "splice_donor_p2": _P_E, "stop_codon": _P_E},
    "CDS_p2": {"CDS_p0": 1 - _P_E, "splice_donor_p0": _P_E},
    "splice_acceptor_p0": {"CDS_p1": 1.0},
    "splice_acceptor_p1": {"CDS_p2": 1.0},
    "splice_acceptor_p2": {"CDS_p0": 1.0},
    "splice_donor_p0": {"intron_p0": 1.0},
    "splice_donor_p1": {"intron_p1": 1.0},
    "splice_donor_p2": {"intron_p2": 1.0},
    "start_codon": {"CDS_p1": 1.0},
    "stop_codon": {"intergenic": 1.0},
}


def _build_transition_matrix(
    transitions: dict[str, dict[str, float]] | None = None,
    index_of: dict[str, int] | None = None,
) -> list[list[float]]:
    """Build the HMM transition matrix from a {from: {to: prob}} grammar.

    Each row sums to 1.0: explicit destination probabilities are taken from
    ``transitions``, and the remainder (1.0 - outflow) is assigned to the
    self-loop on the diagonal. States absent from a row's destination dict get
    probability 0.0. ``index_of`` maps state names to their row/column and its
    size sets the matrix size. Both default to the canonical 15-state grammar;
    the U12 decode passes the extended pair.
    """
    transitions = _TRANSITIONS if transitions is None else transitions
    index_of = CLASS_NAME_TO_INDEX if index_of is None else index_of
    num_states = len(index_of)
    matrix = [[0.0] * num_states for _ in range(num_states)]
    for from_name, to_dict in transitions.items():
        i = index_of[from_name]
        # clamped so a rate change pushing outflow past 1.0 forbids the
        # self-loop (log 0 = -inf) rather than yielding NaN
        matrix[i][i] = max(0.0, 1.0 - sum(to_dict.values()))
        for to_name, prob in to_dict.items():
            matrix[i][index_of[to_name]] = prob
    return matrix


# Disallowed transitions are 0.0 → log = -inf (the intended sentinel); silence
# the expected divide-by-zero warning. NumPy (not jnp) so importing this module
# does not initialise the JAX backend before the platform default is set.
with np.errstate(divide="ignore"):
    _LOG_TRANS_PROBS = np.log(np.array(_build_transition_matrix(), dtype=np.float32))


# Named state-index aliases used by the constraint mask. Single source of truth
# is still CLASS_NAME_TO_INDEX — the symbolic names are convenience aliases.
_ASS0 = CLASS_NAME_TO_INDEX["splice_acceptor_p0"]


_ASS1 = CLASS_NAME_TO_INDEX["splice_acceptor_p1"]


_ASS2 = CLASS_NAME_TO_INDEX["splice_acceptor_p2"]


_DSS0 = CLASS_NAME_TO_INDEX["splice_donor_p0"]


_DSS1 = CLASS_NAME_TO_INDEX["splice_donor_p1"]


_DSS2 = CLASS_NAME_TO_INDEX["splice_donor_p2"]


_START = CLASS_NAME_TO_INDEX["start_codon"]


_STOP = CLASS_NAME_TO_INDEX["stop_codon"]


# Token-id -> base code (A=0, C=1, G=2, T=3, N=4), case-insensitive. Built from
# the soft-masked tokenizer so both upper- and lower-case nucleotide tokens
# (RepeatMasker soft-mask) decode to the same base code; the non-soft-masked
# vocab is a prefix of it, so this single LUT serves both model types. Soft-mask
# case marks repeats but is orthogonal to the start/splice/stop motifs the
# constrained Viterbi checks, so repeat bases are decoded by identity rather than
# waived. Any non-ACGTN token id (pad / unk / ...) maps to N.
def _build_token_to_base_code() -> np.ndarray:
    tokenizer = get_ntv3_tokenizer(include_softmasked_tokens=True)
    lut = np.full(tokenizer.vocabulary_size, NUC_N, dtype=np.uint8)
    for base, code in (
        ("A", NUC_A),
        ("C", NUC_C),
        ("G", NUC_G),
        ("T", NUC_T),
        ("N", NUC_N),
    ):
        lut[tokenizer.token_to_id(base)] = code
        lut[tokenizer.token_to_id(base.lower())] = code
    return lut


TOKEN_TO_BASE_CODE: np.ndarray = _build_token_to_base_code()


def build_rc_token_lut(tokenizer: NucleotideTokenizer) -> np.ndarray:
    """Token-id -> reverse-complement token-id permutation.

    Each single-nucleotide token maps to its complement with case preserved, so
    soft-masked lowercase tokens stay soft-masked; every other token (specials,
    padding) maps to itself. Indexing tokens through this LUT and reversing along
    the sequence axis reverse-complements a tokenized window without
    re-tokenizing — used by single-strand models to recover the minus strand via
    a second forward pass.
    """
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",
        "n": "n",
    }
    lut = np.arange(tokenizer.vocabulary_size, dtype=np.int32)
    for base, comp in complement.items():
        if base in tokenizer.standard_tokens and comp in tokenizer.standard_tokens:
            lut[tokenizer.token_to_id(base)] = tokenizer.token_to_id(comp)
    return lut


def reverse_complement_base_codes(base_codes: jax.Array) -> jax.Array:
    """Reverse-complement an A=0/C=1/G=2/T=3/N=4 base-code array along the last axis."""
    is_n = base_codes == NUC_N
    complement = jnp.where(is_n, NUC_N, NUC_T - base_codes)
    return complement[..., ::-1]


_NEG_INF = np.float32(-np.inf)


@dataclass(frozen=True)
class _MotifRule:
    """One accepted motif pattern and the log-score it awards when matched.

    offsets are positions relative to t (e.g. (0, 1, 2) for an ATG starting at
    t); bases are the expected base codes at those offsets; log_score is the
    additive log-probability applied to the constrained state on a match.
    """

    offsets: tuple[int, ...]
    bases: tuple[int, ...]
    log_score: float


# Motif log-scores (float32 to match the model's log-emission dtype). Scores
# follow Tiberius §3.9.2 — a small TAG stop-codon preference and the
# canonical/rare GT/GC splice-donor split.
_LOG_034 = float(np.log(np.float32(0.34)))


_LOG_033 = float(np.log(np.float32(0.33)))


_LOG_099 = float(np.log(np.float32(0.99)))


_LOG_001 = float(np.log(np.float32(0.01)))


# Biological constraints, one entry per constrained motif. Each motif lists its
# mutually-exclusive accepted patterns; a position matching none is vetoed
# (-inf).
_MOTIF_RULES: dict[str, list[_MotifRule]] = {
    "start": [_MotifRule((0, 1, 2), (NUC_A, NUC_T, NUC_G), 0.0)],
    "stop": [
        _MotifRule((-2, -1, 0), (NUC_T, NUC_A, NUC_G), _LOG_034),
        _MotifRule((-2, -1, 0), (NUC_T, NUC_A, NUC_A), _LOG_033),
        _MotifRule((-2, -1, 0), (NUC_T, NUC_G, NUC_A), _LOG_033),
    ],
    "donor": [
        _MotifRule((1, 2), (NUC_G, NUC_T), _LOG_099),
        _MotifRule((1, 2), (NUC_G, NUC_C), _LOG_001),
    ],
    "acceptor": [_MotifRule((-2, -1), (NUC_A, NUC_G), 0.0)],
}


# HMM states each motif constrains.
_MOTIF_TO_STATES: dict[str, list[int]] = {
    "start": [_START],
    "stop": [_STOP],
    "donor": [_DSS0, _DSS1, _DSS2],
    "acceptor": [_ASS0, _ASS1, _ASS2],
}


def _compute_motif_mask(base_codes: jax.Array, rules: list[_MotifRule]) -> jax.Array:
    """Score one motif's accepted patterns at every position.

    Returns a (B, T) additive log-mask: log_score of the matched pattern, 0 if
    any base in the motif's context is N (constraint waived for assembly gaps),
    and -inf where no pattern matches or the required context falls outside the
    window.
    """
    _, seq_len = base_codes.shape
    pos = jnp.arange(seq_len)
    padded = jnp.pad(base_codes, ((0, 0), (2, 2)), constant_values=NUC_N)

    all_offsets = {off for rule in rules for off in rule.offsets}
    min_off, max_off = min(all_offsets), max(all_offsets)
    in_bounds = (pos >= -min_off) & (pos < seq_len - max_off)

    def window(off: int) -> jax.Array:
        return padded[:, 2 + off : 2 + off + seq_len]

    has_n = jnp.zeros((base_codes.shape[0], seq_len), dtype=bool)
    for off in all_offsets:
        has_n = has_n | (window(off) == NUC_N)

    score = jnp.full((base_codes.shape[0], seq_len), _NEG_INF)
    for rule in rules:
        match = jnp.ones((base_codes.shape[0], seq_len), dtype=bool)
        for off, base in zip(rule.offsets, rule.bases, strict=True):
            match = match & (window(off) == base)
        score = jnp.where(match, jnp.float32(rule.log_score), score)

    return jnp.where(~in_bounds, _NEG_INF, jnp.where(has_n, jnp.float32(0.0), score))


def build_state_constraint_log_mask(base_codes: jax.Array) -> jax.Array:
    """Additive (B, T, 15) log-emission mask vetoing start/splice/stop states.

    A state is forbidden wherever the local sequence can't form the motif it
    requires. States are constrained in the sequence-direction frame (the
    caller passes reverse-complemented codes for the minus strand). Mask values
    are added to the model's log-emissions before Viterbi:

    - start codon: hard ATG at seq[t:t+3].
    - stop codon: seq[t-2:t+1] ∈ {TAG, TAA, TGA} with log-probs log(0.34/0.33/0.33).
    - splice donor: seq[t+1:t+3] ∈ {GT, GC} with log-probs log(0.99/0.01).
    - splice acceptor: hard AG at seq[t-2:t].

    Any N in the relevant context bypasses the constraint (0 = allowed) so
    assembly gaps don't suppress valid states. The codes are already cropped to
    the kept region, so its outermost ≤ 2 positions lack the context a motif
    needs and forbid these states (-inf) whatever the sequence holds.
    """
    batch_size, seq_len = base_codes.shape
    zeros = jnp.zeros((batch_size, seq_len), dtype=jnp.float32)

    column_for: dict[int, jax.Array] = {}
    for motif_name, rules in _MOTIF_RULES.items():
        mask = _compute_motif_mask(base_codes, rules)
        for state_idx in _MOTIF_TO_STATES[motif_name]:
            column_for[state_idx] = mask

    columns = [column_for.get(i, zeros) for i in range(NUM_STATES)]
    return jnp.stack(columns, axis=-1)


# The U2 splice chain the lane mirrors, one state per phase.
_U12_MIRRORED: tuple[str, ...] = tuple(
    f"{kind}_p{phase}"
    for kind in ("splice_donor", "intron", "splice_acceptor")
    for phase in range(3)
)


def _u12(state_name: str) -> str:
    """A mirrored U2 state name -> its lane twin (splice_donor_p1 -> splice_donor_u12_p1)."""
    kind, phase = state_name.rsplit("_p", 1)
    return f"{kind}_u12_p{phase}"


EXT_STATE_NAMES: tuple[str, ...] = tuple(GENOME_ANNOTATION_CLASSES) + tuple(
    _u12(name) for name in _U12_MIRRORED
)


NUM_EXT_STATES = len(EXT_STATE_NAMES)


_EXT_INDEX: dict[str, int] = {name: i for i, name in enumerate(EXT_STATE_NAMES)}


# Share of introns taking the minor spliceosome, from the census above. Priced
# on the CDS -> donor transition rather than in the motif score, so the U12
# motifs stay hard and the frequency is not counted twice. Mammalian, like the
# _P_IR / _P_I / _P_E mean lengths and their plant TODO above.
U12_INTRON_FRACTION = 0.00105


# Extended state -> the emission column it reads, and equally the canonical
# state it collapses to after decode. The canonical states lead EXT_STATE_NAMES,
# so they keep their index and this is the identity on them.
_EXT_TO_CANONICAL = np.arange(NUM_EXT_STATES, dtype=np.int32)


for _u2_name in _U12_MIRRORED:
    _EXT_TO_CANONICAL[_EXT_INDEX[_u12(_u2_name)]] = CLASS_NAME_TO_INDEX[_u2_name]


def _build_ext_transition_matrix() -> list[list[float]]:
    """The 15-state grammar plus a paired U12 intron lane, NUM_EXT x NUM_EXT.

    The lane is the U2 splice chain with every within-lane destination renamed,
    so it inherits its counterpart's phase bookkeeping and rejoins the canonical
    CDS lane exactly where the U2 acceptor would. Each CDS state's exon-exit
    mass is then split between the two donors by U12_INTRON_FRACTION.
    """
    transitions = {name: dict(dests) for name, dests in _TRANSITIONS.items()}
    for u2_name in _U12_MIRRORED:
        transitions[_u12(u2_name)] = {
            (_u12(dest) if dest in _U12_MIRRORED else dest): prob
            for dest, prob in _TRANSITIONS[u2_name].items()
        }
    for phase in range(3):
        cds = f"CDS_p{phase}"
        donor = next(d for d in transitions[cds] if d.startswith("splice_donor_"))
        exit_mass = transitions[cds].pop(donor)
        transitions[cds][donor] = exit_mass * (1.0 - U12_INTRON_FRACTION)
        transitions[cds][_u12(donor)] = exit_mass * U12_INTRON_FRACTION
    return _build_transition_matrix(transitions, _EXT_INDEX)


with np.errstate(divide="ignore"):
    _LOG_TRANS_PROBS_EXT = np.log(
        np.array(_build_ext_transition_matrix(), dtype=np.float32)
    )


# The U12 ends are hard single motifs; their rarity lives in the transition.
_U12_MOTIF_RULES: dict[str, list[_MotifRule]] = {
    "splice_donor": [_MotifRule((1, 2), (NUC_A, NUC_T), 0.0)],
    "intron": [],
    "splice_acceptor": [_MotifRule((-2, -1), (NUC_A, NUC_C), 0.0)],
}


def build_ext_state_constraint_log_mask(base_codes: jax.Array) -> jax.Array:
    """The 15-state mask with the nine U12 columns appended.

    The canonical columns are ``build_state_constraint_log_mask`` unchanged. The
    U12 donor columns require AT at ``seq[t+1:t+3]``, the acceptor columns AC at
    ``seq[t-2:t]``, and the intron columns are unconstrained; they are ordered by
    _U12_MIRRORED to match EXT_STATE_NAMES.
    """
    unconstrained = jnp.zeros(base_codes.shape, dtype=jnp.float32)
    per_kind = {
        kind: _compute_motif_mask(base_codes, rules) if rules else unconstrained
        for kind, rules in _U12_MOTIF_RULES.items()
    }
    lane = [per_kind[name.rsplit("_p", 1)[0]] for name in _U12_MIRRORED]
    return jnp.concatenate(
        [build_state_constraint_log_mask(base_codes), jnp.stack(lane, axis=-1)],
        axis=-1,
    )


def compute_window_geometry(
    sequence_length: int,
    keep_target_center_fraction: float,
) -> int:
    """Return the per-side model-context margin: the bp the head does not predict.

    It is the whole geometry of a window — the kept region is what is left between
    the two margins — and it is why covering a contig end to end needs a window
    outside it at each end (see `select_windows`).
    """
    prediction_length = int(sequence_length * keep_target_center_fraction)
    return (sequence_length - prediction_length) // 2


def select_windows(
    sequence_length: int,
    keep_target_center_fraction: float,
    chrom_lengths: dict[str, int],
    chrom: str,
) -> tuple[list[GenomicInterval], int]:
    """Tile a chromosome so the kept regions cover `[0, chrom_length)` exactly.

    A window's kept region is the central `stride` bp of its row — what the head
    predicts — so the window covering a region starts `center_offset` bp ahead of
    it. The regions are a stride grid anchored at bp 0, which puts a window outside
    the contig at each end:

    - The first region starts at bp 0, so its window starts at `-center_offset`,
      returned as such: `build_dataloader` reads it as `left_pad` N bases followed by
      the contig. Keeping more of the window instead would score positions the head
      was never trained on.
    - The last region is shifted back to end on the contig's final base rather than
      tiled past the end and N-padded, which would leave it part synthetic. Its
      window still overhangs, but only by `center_offset`, into the right context
      margin.

    A contig shorter than one kept region tiles to a single mostly-N window rather than
    being dropped.

    That shifted last region is the only place two kept regions overlap, and
    `plan_contig_windows` gives it to the earlier of the two.

    Args:
        sequence_length: Window length fed to the model (bp).
        keep_target_center_fraction: Fraction of `sequence_length` kept as
            the prediction region; the remainder is split as model-context
            margin on both sides (see `compute_window_geometry`). Must leave a
            positive prediction region.
        chrom_lengths: Map chromosome name -> chromosome length (bp).
        chrom: Chromosome to tile; must be a key in `chrom_lengths`.

    Returns:
        windows: `GenomicInterval`s of length `sequence_length`, ascending by
            start, at least one however short the contig. The leading ones may start
            before 0 and the last may end past `chrom_lengths[chrom]`, by at most
            `center_offset` either way.
        center_offset: bp cropped at each side as model-context margin.

    """
    center_offset = compute_window_geometry(
        sequence_length, keep_target_center_fraction
    )
    stride = sequence_length - 2 * center_offset
    if stride <= 0:
        raise ValueError(
            f"keep_target_center_fraction={keep_target_center_fraction} leaves no "
            f"kept region for sequence_length={sequence_length}."
        )
    chrom_length = chrom_lengths[chrom]
    # Kept regions are a stride grid anchored at bp 0, with the last one pulled back
    # onto the contig end: clamping every start to it both shifts that region and
    # absorbs the grid points it now covers.
    last_start = max(chrom_length - stride, 0)
    kept_starts = [min(start, last_start) for start in range(0, chrom_length, stride)]
    windows = [
        GenomicInterval(
            chrom, start - center_offset, start - center_offset + sequence_length
        )
        for start in kept_starts
    ]
    return windows, center_offset


def _advance_stable(
    scores: jax.Array,
    log_emissions_t: jax.Array,
    log_trans: jax.Array,
    *,
    matrix_carry: bool,
) -> jax.Array:
    """One tropical (max-plus) HMM step, rescaled by the carry's own global max.

    Two different carries run through this step, and they differ only in how many
    trailing state axes they have — hence ``matrix_carry``: pass 1 advances a full
    ``(..., K, K)`` transfer matrix (best score from every chunk-entry state to
    every current state), pass 3 a plain ``(..., K)`` score vector.

    Subtracting one scalar (the max over those state axes) from every entry leaves
    all later cross-state argmax comparisons invariant — the shift cancels — so
    the decoded path is exactly unchanged, while running magnitudes stay bounded
    near 0 instead of growing ~linearly with sequence length. That unbounded
    growth is what makes a plain float32 Viterbi accumulation lose the
    true-optimal path at ~1 Mb+ scale.
    """
    stepped = jnp.max(scores[..., :, None] + log_trans, axis=-2) + log_emissions_t

    # Reduce over the state axes only, so each chunk and sequence keeps its own
    # offset (a per-entry-state offset would not cancel in the later argmaxes).
    state_axes: tuple[int, ...]
    if matrix_carry:
        state_axes = (-2, -1)  # (..., K, K): entry state and current state.
    else:
        state_axes = (-1,)  # (..., K): current state.

    offset = jnp.maximum(
        jnp.max(stepped, axis=state_axes, keepdims=True), jnp.float32(-1e30)
    )
    return stepped - offset


def viterbi_jax_parallel(
    log_emissions: jax.Array,
    chunk_len: int = 512,
    log_trans_probs: np.ndarray | None = None,
) -> jax.Array:
    """Parallel (blocked) Viterbi decode, returning the exact optimal path.

    ``(batch, seq_len, num_states)`` log-emissions -> ``(batch, seq_len)`` int8
    states. The caller is responsible for strand reversal: minus-strand log
    emissions must be reversed along the time axis before calling and the
    returned states reversed back into genomic order.

    Where a textbook Viterbi is one length-``T`` sequential ``jax.lax.scan``
    (critical path ``O(T)`` steps), this decoder splits the ``T-1`` transitions
    into local blocks of ``chunk_len`` positions and does three shallow passes:

      1. ``extend_transfer_matrices``: per-chunk transfer matrices (K x K), all
         chunks advanced in parallel (a scan of depth ``chunk_len``, not ``T``);
      2. ``advance_across_chunk``: an exclusive scan over the chunk matrices
         (depth ``n_chunks = T/chunk_len``) for each chunk's entry scores;
      3. ``replay_with_backpointers``: a per-chunk backpointer replay (depth
         ``chunk_len``), then an ``O(log T)`` associative-scan compose of
         backpointers to read the path.

    Sequential depth drops from ``T`` to ``~max(chunk_len, T/chunk_len)`` (minimized
    near ``sqrt(T)``), which is much faster on GPU at genome scale. Each accumulation
    step is rescaled by its global max (``_advance_stable``); this is a uniform shift
    that cancels in every later argmax, so the path is exactly preserved while running
    magnitudes stay bounded (the sequential float32 scan instead drifts off the
    optimum past ~1 Mb as its log-scores grow unboundedly).

    Args:
        log_emissions: ``(batch, seq_len, num_states)`` log-emission scores.
        chunk_len: Local-Viterbi block length; the parallelism granularity.
        log_trans_probs: ``(num_states, num_states)`` log transition matrix;
            defaults to the canonical 15-state grammar. The U12-paired decode
            passes the extended matrix (see _build_ext_transition_matrix).

    Note: this materializes ``n_chunks`` K x K matrices and has heavy parallel
    fan-out, so it is meant for GPU/TPU. On CPU the fan-out is pure overhead and
    only ``chunk_len >= seq_len`` (one chunk) stays close to a sequential scan.
    The caller decodes bounded per-window / few-Mb dispatches rather than a whole
    chromosome in one call.
    """
    log_trans = jnp.asarray(
        _LOG_TRANS_PROBS if log_trans_probs is None else log_trans_probs
    )
    batch, seq_len, num_states = log_emissions.shape

    init_scores = log_emissions[:, 0, :]
    n_steps = seq_len - 1
    n_chunks = -(-n_steps // chunk_len)
    n_pad = n_chunks * chunk_len - n_steps

    # (chunk_len, n_chunks, batch, num_states): time-major inside the per-chunk scans.
    padded = jnp.pad(log_emissions[:, 1:, :], ((0, 0), (0, n_pad), (0, 0)))
    chunked_emissions = padded.reshape(
        batch, n_chunks, chunk_len, num_states
    ).transpose(2, 1, 0, 3)
    valid = (jnp.arange(n_chunks * chunk_len) < n_steps).reshape(n_chunks, chunk_len).T

    # Pass 1: per-chunk transfer matrices, all chunks in parallel. Carry (C, B, K, K).
    def extend_transfer_matrices(
        carry: jax.Array, xs: tuple[jax.Array, jax.Array]
    ) -> tuple[jax.Array, None]:
        """Extend every chunk's entry-state -> current-state matrix by one position."""
        log_emissions_t, valid_t = xs
        stepped = _advance_stable(
            carry, log_emissions_t[:, :, None, :], log_trans, matrix_carry=True
        )
        return jnp.where(valid_t[:, None, None, None], stepped, carry), None

    tropical_eye = jnp.where(jnp.eye(num_states, dtype=bool), 0.0, -jnp.inf)
    eye = jnp.broadcast_to(tropical_eye, (n_chunks, batch, num_states, num_states))
    chunk_mats, _ = jax.lax.scan(
        extend_transfer_matrices, eye, (chunked_emissions, valid)
    )

    # Pass 2: exclusive scan over chunk matrices -> per-chunk entry scores.
    def advance_across_chunk(
        scores: jax.Array, mat: jax.Array
    ) -> tuple[jax.Array, jax.Array]:
        """Jump the running scores over one whole chunk, emitting its entry scores."""
        stepped = jnp.max(scores[:, :, None] + mat, axis=-2)
        offset = jnp.maximum(
            jnp.max(stepped, axis=-1, keepdims=True), jnp.float32(-1e30)
        )
        return stepped - offset, scores

    _, entry_scores = jax.lax.scan(advance_across_chunk, init_scores, chunk_mats)

    # Pass 3: replay each chunk from its entry scores, emitting backpointers.
    def replay_with_backpointers(
        scores: jax.Array, xs: tuple[jax.Array, jax.Array]
    ) -> tuple[jax.Array, jax.Array]:
        """Re-run one position of every chunk, recording each state's predecessor."""
        log_emissions_t, valid_t = xs
        backpointers_t = jnp.argmax(scores[..., :, None] + log_trans, axis=-2).astype(
            jnp.int8
        )
        stepped = _advance_stable(
            scores, log_emissions_t, log_trans, matrix_carry=False
        )
        return jnp.where(valid_t[:, None, None], stepped, scores), backpointers_t

    final_scores, backpointers = jax.lax.scan(
        replay_with_backpointers, entry_scores, (chunked_emissions, valid)
    )

    backpointers = backpointers.transpose(1, 0, 2, 3).reshape(-1, batch, num_states)[
        :n_steps
    ]
    final_state = jnp.argmax(final_scores[-1], axis=-1).astype(jnp.int8)

    # Backpointer maps compose like functions: (a o b)[..., x] = b[..., a[..., x]].
    composed = jax.lax.associative_scan(
        lambda a, b: jnp.take_along_axis(b, a, axis=-1), backpointers[::-1], axis=0
    )
    head = jnp.take_along_axis(composed, final_state[None, :, None], axis=-1)[..., 0]
    path = jnp.concatenate([head[::-1], final_state[None]], axis=0)
    return path.transpose(1, 0).astype(jnp.int8)


def decode_log_emissions_argmax(
    log_emissions_plus: jax.Array, log_emissions_minus: jax.Array
) -> jax.Array:
    """Per-position argmax over log emissions. Returns (B, T, 2) int8; axis -1 is strand (0=+, 1=-)."""
    arg_p = jnp.argmax(log_emissions_plus, axis=-1).astype(jnp.int8)
    arg_m = jnp.argmax(log_emissions_minus, axis=-1).astype(jnp.int8)
    return jnp.stack([arg_p, arg_m], axis=-1)


def decode_log_emissions_hmm(
    log_emissions_plus: jax.Array,
    log_emissions_minus: jax.Array,
    decode_chunk_len: int,
) -> jax.Array:
    """Per-strand Viterbi decode (both strands batched into one scan). Returns (B, T, 2) int8.

    ``decode_chunk_len`` is the decoder's block length (see viterbi_jax_parallel).
    """
    log_emissions_minus_rev = log_emissions_minus[:, ::-1, :]
    stacked = jnp.concatenate([log_emissions_plus, log_emissions_minus_rev], axis=0)
    states = viterbi_jax_parallel(stacked, decode_chunk_len)
    batch = log_emissions_plus.shape[0]
    hmm_p = states[:batch]
    hmm_m = states[batch:][:, ::-1]
    return jnp.stack([hmm_p, hmm_m], axis=-1)


def decode_log_emissions_hmm_constrained(
    log_emissions_plus: jax.Array,
    log_emissions_minus: jax.Array,
    base_codes: jax.Array,
    decode_chunk_len: int,
) -> jax.Array:
    """Constraint-masked per-strand Viterbi decode. Returns (B, T, 2) int8.

    Adds ``build_ext_state_constraint_log_mask`` to each strand's log-emissions
    before Viterbi. The minus mask is built from the reverse-complement of
    ``base_codes`` so constraints apply in the transcribed direction.

    Decoding runs on the extended state space carrying the paired U12 (AT..AC)
    intron lane; its states read their U2 counterparts' emission column and
    collapse back onto them before returning, so the returned states are the 15
    canonical ones.

    ``decode_chunk_len`` is the decoder's block length (see viterbi_jax_parallel).
    """
    gather = jnp.asarray(_EXT_TO_CANONICAL)
    log_emissions_minus_rev = log_emissions_minus[:, ::-1, :]
    mask_plus = build_ext_state_constraint_log_mask(base_codes)
    mask_minus_rev = build_ext_state_constraint_log_mask(
        reverse_complement_base_codes(base_codes)
    )
    stacked = jnp.concatenate(
        [
            log_emissions_plus[..., gather] + mask_plus,
            log_emissions_minus_rev[..., gather] + mask_minus_rev,
        ],
        axis=0,
    )
    states = viterbi_jax_parallel(stacked, decode_chunk_len, _LOG_TRANS_PROBS_EXT)
    # U12 states are decode-only; nothing downstream knows about them
    states = jnp.asarray(_EXT_TO_CANONICAL, dtype=states.dtype)[states]
    batch = log_emissions_plus.shape[0]
    hc_p = states[:batch]
    hc_m = states[batch:][:, ::-1]
    return jnp.stack([hc_p, hc_m], axis=-1)


def resolve_condition_inputs(
    conditions_names: list[str],
    name_to_tokenizer: Mapping[str, BaseTokenizer],
    common_name: str,
    *,
    mask_species_condition: bool,
    data_sharding: jax.sharding.NamedSharding,
) -> tuple[list[jax.Array], list[jax.Array] | None]:
    """Resolve device-side species condition tokens and optional masks.

    Encodes common_name with each condition's tokenizer as a (1,)-shaped int
    array (broadcast to the batch inside the jitted inference fn). When
    mask_species_condition is set the species token is zeroed and a boolean mask
    list (False for species, True elsewhere) is returned so the model drops the
    species conditioning; otherwise the mask is None.

    These (1,)-shaped arrays are placed with a replicated sharding rather than
    ``data_sharding``: they are broadcast to the batch inside the jit, so every
    device needs the full value. Sharding a size-1 leading dimension along the
    data axis crashes on more than one device (``data_sharding`` requires the
    dimension be divisible by the device count), so only the mesh is reused.

    When conditioning is off (``conditions_names == []``) there is nothing to
    place, so this returns empty lists without touching ``data_sharding`` (the
    unconditioned path never sets up a mesh).
    """
    if not conditions_names:
        return [], None
    replicated_sharding = jax.sharding.NamedSharding(
        data_sharding.mesh, jax.sharding.PartitionSpec()
    )
    conditions_tokens_1d: list[jax.Array] = []
    for name in conditions_names:
        if mask_species_condition and name == "species":
            token = np.zeros((1,), dtype=np.int32)
        else:
            # [:1] keeps shape (1,): np_tokenize falls back to char-level for
            # OOV names, which would otherwise yield (N>1,) and break the
            # broadcast at inference.
            token = name_to_tokenizer[name].np_tokenize(common_name)[:1]  # type: ignore[attr-defined]
        conditions_tokens_1d.append(
            jax.device_put(jnp.asarray(token), replicated_sharding)
        )

    conditions_masks_1d: list[jax.Array] | None = None
    if mask_species_condition and "species" in conditions_names:
        # True keeps the condition active, False drops it. Only species is
        # droppable — taxonomy conditions are always active.
        conditions_masks_1d = [
            jax.device_put(
                jnp.array([name != "species"], dtype=jnp.bool_), replicated_sharding
            )
            for name in conditions_names
        ]
    return conditions_tokens_1d, conditions_masks_1d


# ---------------------------------------------------------------------------
# Model inference
#
# From scripts/ntv3_post_training/annotate/utils.py. The old-arch branch is
# gone: only the NTv3 BED-head family loads here, so there is no skip_lm_head
# to pass conditionally and no LongRangeNT isinstance check.
# ---------------------------------------------------------------------------


def run_annotation_inference(
    model: nnx.Module,
    tokens: jax.Array,
    conditions_tokens: list[jax.Array],
    conditions_masks: list[jax.Array] | None = None,
) -> jax.Array:
    """Forward the model and return per-position log-emissions over the 15 classes.

    Annotation logits live under ``bed_tracks_logits``, shape (B, S, 1, 15).
    ``file_assembly_indices=None`` is passed for signature compatibility; the
    released annotation models build no bigwig head (``number_bigwig_tracks``
    is 0), so nothing reads it.
    """
    outputs = model(
        tokens=tokens,
        conditions_tokens=conditions_tokens,
        conditions_masks=conditions_masks,
        file_assembly_indices=None,
    )
    return jax.nn.log_softmax(outputs["bed_tracks_logits"], axis=-1)


def commit_model_state(model: nnx.Module) -> None:
    """Place every uncommitted state leaf where the model's parameters live.

    ``model.eval()`` writes each ``RMSBatchNorm.use_running_average`` flag as an
    uncommitted host array. ``nnx.jit`` hands those flags back committed after the
    first call, so the second call sees new input shardings and recompiles the whole
    forward pass (~24 s on the 700M model).
    """
    state = nnx.state(model)
    leaves = [x for x in jax.tree.leaves(state) if isinstance(x, jax.Array)]
    committed = [x.sharding for x in leaves if x.committed]
    if len(committed) == len(leaves):
        return
    # A loaded model's parameters are committed, so the flags join them, replicated
    # on their mesh. When nothing is committed (a stub model nothing has placed),
    # each leaf is committed where it already lives, so no array moves.
    reference = committed[0] if committed else None
    if isinstance(reference, jax.sharding.NamedSharding):
        target = jax.sharding.NamedSharding(
            reference.mesh, jax.sharding.PartitionSpec()
        )
    else:
        target = reference

    def commit(x: jax.Array) -> jax.Array:
        return jax.device_put(x, target or x.sharding)

    nnx.update(
        model,
        jax.tree.map(
            lambda x: commit(x) if isinstance(x, jax.Array) and not x.committed else x,
            state,
        ),
    )


def build_annotation_inference_fn(
    model: nnx.Module,
    *,
    decoders: list[str],
    center_offset: int,
    rc_token_lut: jax.Array,
    decode_chunk_len: int | None = None,
) -> Callable[..., dict[str, jax.Array]]:
    """Build the fused forward -> log-softmax -> decode closure.

    The plus strand comes from the forward pass; the minus strand from a second
    pass on the reverse-complemented tokens, flipped back to genomic order.

    Returns a cached, jitted callable taking
    ``(tokens, conditions_tokens_1d, conditions_masks_1d)`` and returning a dict
    keyed by enabled decoder name. The ``(1,)``-shaped condition arrays are
    broadcast to the batch inside the jit.

    Args:
        model: The loaded annotation model.
        decoders: Enabled decoder names (argmax / hmm / hmm_constrained).
        center_offset: Context margin from the window geometry; locates the kept
            region when deriving hmm_constrained base codes.
        rc_token_lut: Token-id -> complement-token-id LUT for the minus pass.
        decode_chunk_len: Viterbi block length, required when hmm or
            hmm_constrained is enabled and ignored otherwise.
    """
    run_argmax = "argmax" in decoders
    run_hmm = "hmm" in decoders
    run_hmm_constrained = "hmm_constrained" in decoders
    if (run_hmm or run_hmm_constrained) and decode_chunk_len is None:
        raise ValueError("hmm / hmm_constrained decoding requires decode_chunk_len.")
    chunk_len = cast(int, decode_chunk_len)
    base_code_lut = jnp.asarray(TOKEN_TO_BASE_CODE)

    def inference_and_decode(
        model: nnx.Module,
        tokens: jax.Array,
        conditions_tokens_1d: list[jax.Array],
        conditions_masks_1d: list[jax.Array] | None,
    ) -> dict[str, jax.Array]:
        batch_size = tokens.shape[0]
        conditions_tokens = [
            jnp.broadcast_to(t, (batch_size,)) for t in conditions_tokens_1d
        ]
        conditions_masks = (
            [jnp.broadcast_to(m, (batch_size,)) for m in conditions_masks_1d]
            if conditions_masks_1d is not None
            else None
        )
        # Single-strand head: forward gives plus-strand emissions; a second pass on
        # the reverse complement gives minus-strand emissions in RC'd order, flipped
        # back to genomic order. A boundary row carries N padding at one end — ahead of
        # bp 0, or past the contig end — which the RC relocates to the other. Alignment
        # holds either way: the flip spans the whole row, so position i maps to L-1-i.
        # A front pad sits in the cropped margin; a tail pad is dropped by accumulate.
        rc_tokens = rc_token_lut[tokens][:, ::-1]
        # Both strands through one lax.map body, so XLA compiles the network once;
        # two inlined calls doubled the compile.
        graphdef, state = nnx.split(model)

        def strand_pass(strand_tokens: jax.Array) -> jax.Array:
            return run_annotation_inference(
                nnx.merge(graphdef, state),
                strand_tokens,
                conditions_tokens,
                conditions_masks,
            )

        both = jax.lax.map(strand_pass, jnp.stack([tokens, rc_tokens]))
        log_probs, log_probs_rc = both[0], both[1]
        log_emissions_plus = log_probs[:, :, 0, :]
        log_emissions_minus = log_probs_rc[:, ::-1, 0, :]
        out: dict[str, jax.Array] = {}
        if run_argmax:
            out["argmax"] = decode_log_emissions_argmax(
                log_emissions_plus, log_emissions_minus
            )
        if run_hmm:
            out["hmm"] = decode_log_emissions_hmm(
                log_emissions_plus, log_emissions_minus, chunk_len
            )
        if run_hmm_constrained:
            base_codes = base_code_lut[tokens][
                :, center_offset : tokens.shape[1] - center_offset
            ]
            out["hmm_constrained"] = decode_log_emissions_hmm_constrained(
                log_emissions_plus,
                log_emissions_minus,
                base_codes,
                chunk_len,
            )
        return out

    commit_model_state(model)
    return cast(
        Callable[..., dict[str, jax.Array]],
        nnx.cached_partial(nnx.jit(inference_and_decode), model),
    )

# ---------------------------------------------------------------------------
# Crop resolution, from scripts/ntv3_genome_annotation/annotate/config.py
# ---------------------------------------------------------------------------


class _HasCropFraction(Protocol):
    """Any model config exposing the head's centre-crop fraction."""

    keep_target_center_fraction: float


DEFAULT_KEEP_TARGET_CENTER_FRACTION = 0.75


def pin_keep_target_center_fraction(
    model_config: _HasCropFraction,
    requested: float | None,
    *,
    trained: float | None = None,
) -> float:
    """Resolve the crop to annotate at and pin it onto the model config.

    The head crops the embedding with ``model_config.keep_target_center_fraction`` at
    forward time, so writing the resolved value back is what keeps the head and the
    window geometry in step. Call before building the inference fn, which bakes the
    value in at trace time.

    Args:
        model_config: Loaded model config, mutated in place when the crop narrows.
        requested: Crop to annotate at, or None for
            :data:`DEFAULT_KEEP_TARGET_CENTER_FRACTION` capped by the trained crop.
            Narrowing only increases window overlap; widening asks the head for
            window-edge positions it was never supervised on.
        trained: The crop the head was **trained** at, which is not always the
            one its config carries. A released export substitutes the *served*
            crop into ``keep_target_center_fraction`` so that a bare forward
            pass produces the served geometry, and records the trained value
            under ``provenance.trained_keep_target_center_fraction``. Reading
            the config alone then binds ``trained`` to 0.75 for a head trained
            at 1.0, which lowers the ceiling, states a falsehood in the error,
            and makes the narrowing warning below unreachable because
            ``min(0.75, 0.75)`` always equals it. None falls back to the config
            for an export that records nothing.

    Returns:
        The crop the window geometry must be built from.

    Raises:
        ValueError: If ``requested`` is outside ``(0, trained]``.
    """
    trained = (
        model_config.keep_target_center_fraction if trained is None else trained
    )
    if requested is None:
        # Capped by the trained crop rather than applied flat. A head trained
        # below the default was never supervised on the window-edge positions a
        # wider crop asks for, so the default narrows to what that head supports
        # instead of raising on a model it was not chosen for.
        requested = min(DEFAULT_KEEP_TARGET_CENTER_FRACTION, trained)
    if not 0.0 < requested <= trained:
        raise ValueError(
            f"keep_target_center_fraction={requested} must be in (0, {trained}] — "
            f"the crop this head was trained with. Above it, the head is asked for "
            "window-edge positions it was never supervised on."
        )
    if requested != trained:
        logger.warning(
            f"Annotating at keep_target_center_fraction={requested} instead of "
            f"the trained {trained}: more windows per contig, each kept position "
            "further from a window edge."
        )
    model_config.keep_target_center_fraction = requested
    return requested
