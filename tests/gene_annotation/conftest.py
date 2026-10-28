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

"""Fixtures for the production annotation tests.

The golden test needs production's output to be a deterministic function of its
input, with no trained weights and no accelerator. That is what the stub head
here provides: its emissions are a pure lookup on the token id, so the whole
pipeline downstream of the model (window planning, batching, retained-slice
bookkeeping, the real decoder, ``find_genes``, the filters and the writer) is
exercised for real while the output stays exactly reproducible.

**How the synthetic genome signals its states.** Not through base identity,
which would spend the bases the biological motifs need, but through soft-mask
case, which ``hmm_constrained`` cannot see:

* uppercase ``ACGT`` emits the eleven exonic (GFF3 ``CDS``) states,
* lowercase ``acgt`` emits the three ``intron_p*`` states,
* ``N``, ``n`` and every other token (padding, folded IUPAC codes) emit
  intergenic.

The checkpoint's soft-masked vocabulary gives the two cases distinct token ids,
while ``TOKEN_TO_BASE_CODE`` folds them onto the same base code. So the stub
reads the case and the constrained decoder reads the base, and the bases stay
free to spell a real ``ATG`` start, real ``gt``..``ag`` and ``at``..``ac``
(U12) intron boundaries and a real ``TAA`` stop. Each group's emission is
uniform over its states, so which member a position takes is settled by the
grammar and by those motifs rather than by the stub: the constrained decoder
does the structural work and the golden pins its answer.

All three token sets are closed under reverse complement, so a contig reads the
same way on both strands. Genes are placed a spacer apart, and a spacer is
longer than the exon beside it, so no path can pay its way from one gene to the
next: every gene decodes on its own.

**Why each gene appears twice.** ``_compute_motif_mask`` waives a motif whose
context holds an ``N``, so that an assembly gap cannot suppress a valid state.
An exon block flanked by the intergenic ``N`` spacers is therefore openable
from either strand, whatever its bases say, and each planted gene comes back
with a mirror image of itself on the opposite strand, a base or two wider.
That is the waiver made visible and is worth pinning: nothing else in the suite
reaches it. The mirror is not free-form, though. It has to find real splice
sites of its own to place its introns on, which is exactly what
``test_constrained_introns_sit_on_real_splice_motifs`` checks.

Inside a block the sequence still decides. A gene is written without ``C``, so
no antisense ``ATG`` sits inside it and its mirror cannot open anywhere but at
the waived flank; ``chrTail``'s minus-strand gene is the reverse complement of
such a gene, so its plus strand carries no ``G`` and the same holds there.

**The probe base** ``n`` is the stub's one context dependence: it emits exonic in
the leading half of a window's kept region and intergenic in the trailing half.
The minus pass runs over the reversed row, so the strand of the gene a probe run
yields says which half it sat in. No golden contig contains ``n``.
"""

from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from flax import nnx

from nucleotide_transformer_v3.tokenizers import get_ntv3_tokenizer

from nucleotide_transformer_v3.gene_annotation.constants import CLASS_NAME_TO_INDEX
from nucleotide_transformer_v3.gene_annotation.config import (
    ProductionAnnotationConfig,
)
from nucleotide_transformer_v3.gene_annotation.model_registry import (
    FIXED_KEEP_TARGET_CENTER_FRACTION,
)

NUM_ANNOTATION_CLASSES = 15

# Window geometry for the golden. Small enough to stay a CPU-second, cropped
# (keep fraction below 1.0) so the left-padded head window and the shift-back
# tail are both on the tested path.
#
# The crop TRACKS THE PRODUCTION PIN rather than naming a value. It is passed to
# `ProductionAnnotationConfig`, and production refuses any crop but its own, so a
# literal here means every config-building test in this directory breaks the day
# the pin moves - which is how 44 of them failed at once when it went to 0.75.
# 0.75 still satisfies the "below 1.0" requirement above.
GOLDEN_SEQUENCE_LENGTH = 4_096
GOLDEN_KEEP_FRACTION = FIXED_KEEP_TARGET_CENTER_FRACTION
GOLDEN_DECODE_CHUNK_LEN = 512

_INTERGENIC = CLASS_NAME_TO_INDEX["intergenic"]
_INTRON_STATES: tuple[int, ...] = tuple(
    CLASS_NAME_TO_INDEX[f"intron_p{phase}"] for phase in range(3)
)
# Every state ``find_genes`` writes out as a CDS sub-feature. The stub cannot
# tell them apart, which is the point: the grammar and the motif mask do.
_EXON_STATES: tuple[int, ...] = tuple(
    CLASS_NAME_TO_INDEX[name]
    for name in (
        "CDS_p0",
        "CDS_p1",
        "CDS_p2",
        "splice_donor_p0",
        "splice_donor_p1",
        "splice_donor_p2",
        "splice_acceptor_p0",
        "splice_acceptor_p1",
        "splice_acceptor_p2",
        "start_codon",
        "stop_codon",
    )
)

# Logit awarded to every state of a token's group. After the log-softmax the
# pipeline applies, a state outside the group costs this many nats, which
# outweighs the grammar's dearest transition (intergenic -> start_codon, about
# 9.2 nats) and so keeps the group assignment, not the transition prices, in
# charge of where a gene begins and ends.
_GROUP_LOGIT = 10.0


def _group_row(states: tuple[int, ...]) -> np.ndarray:
    """Per-class logits awarding ``_GROUP_LOGIT`` to every state of one group."""
    row = np.zeros(NUM_ANNOTATION_CLASSES, dtype=np.float32)
    row[list(states)] = _GROUP_LOGIT
    return row


def _build_token_logits() -> np.ndarray:
    """Token id -> per-class logits, uniform over the token's state group.

    Built against the soft-masked vocabulary, which is the one the stub declares,
    so the two cases of a base have distinct ids and can carry distinct groups.
    Anything not an ``ACGT`` nucleotide, ``N`` and the folded IUPAC codes
    included, is intergenic, so a stray ambiguity code cannot invent a gene.
    """
    tokenizer = get_ntv3_tokenizer(include_softmasked_tokens=True)
    logits = np.zeros((len(tokenizer.vocabulary), NUM_ANNOTATION_CLASSES), np.float32)
    logits[:, _INTERGENIC] = _GROUP_LOGIT
    for base in "ACGT":
        for token, states in ((base, _EXON_STATES), (base.lower(), _INTRON_STATES)):
            logits[tokenizer.token_to_id(token)] = _group_row(states)
    return logits


# The one token whose emission depends on its position in the window. Lowercase
# ``n`` folds onto ``N``, so the constrained decoder waives every motif across it.
_PROBE_BASE = "n"


def _golden_center_offset() -> int:
    """Context margin the golden geometry crops from each side of a window."""
    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        compute_window_geometry,
    )

    return compute_window_geometry(GOLDEN_SEQUENCE_LENGTH, GOLDEN_KEEP_FRACTION)


# Module-level rather than instance attributes: NNX rejects a raw array leaf on
# a Module, which is why the proven stub in test_annotate_ga_equivalence.py
# carries no attributes either. Closed-over constants also keep the stub
# parameter-free, so there is nothing to initialise and nothing to seed.
_TOKEN_LOGITS = jnp.asarray(_build_token_logits())
_CENTER_OFFSET = _golden_center_offset()
_PROBE_TOKEN_ID = get_ntv3_tokenizer(include_softmasked_tokens=True).token_to_id(
    _PROBE_BASE
)
_PROBE_LEADING_LOGITS = jnp.asarray(_group_row(_EXON_STATES))


class DeterministicAnnotationStub(nnx.Module):
    """Annotation head whose emissions are a token lookup plus one probe base.

    Emits ``_TOKEN_LOGITS[token]`` under both the ``logits`` and
    ``bed_tracks_logits`` keys, so either inference builder accepts it, with
    shape ``(B, S, 1, 15)`` as a single-strand head does.

    Emissions are cropped to the *kept* region, because that is the contract a
    real post-training annotation head has: it predicts only the window centre,
    and ``build_annotation_inference_fn`` consumes ``bed_tracks_logits`` at that
    resolution while cropping the base codes it derives from the full-length
    tokens by ``center_offset``.

    The probe base (``_PROBE_BASE``) is the one exception: exonic in the leading
    half of the kept region, intergenic in the trailing half.
    """

    def __call__(
        self,
        tokens: jax.Array,
        conditions_tokens: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        file_assembly_indices: jax.Array | None = None,
        **kwargs: object,
    ) -> dict[str, jax.Array]:
        kept = tokens[:, _CENTER_OFFSET : tokens.shape[1] - _CENTER_OFFSET]
        in_leading_half = jnp.arange(kept.shape[1]) < kept.shape[1] // 2
        is_leading_probe = (kept == _PROBE_TOKEN_ID) & in_leading_half
        per_position = jnp.where(
            is_leading_probe[:, :, None], _PROBE_LEADING_LOGITS, _TOKEN_LOGITS[kept]
        )
        logits = per_position[:, :, None, :]
        return {"logits": logits, "bed_tracks_logits": logits}


# ---------------------------------------------------------------------------
# The synthetic genome
# ---------------------------------------------------------------------------

# Codons the exon bodies cycle through. No ``C``, so no antisense ``ATG``; no
# ``ATG`` across their joins either, so a gene has exactly one place to open;
# and none of them is a stop, so the only in-frame stop is the one planted.
_BODY_CODONS: tuple[str, ...] = ("GGA", "TTG", "AAG", "GAT", "TGG", "ATT")
_START_CODON = "ATG"
_STOP_CODON = "TAA"
_PLANTED_STOP_CODON = "TGA"

# Intron boundaries, lowercase so the stub reads them as intron whatever their
# bases spell. The filler carries neither ``a`` nor ``c``, so the only ``ag`` /
# ``ac`` an intron holds is its own acceptor and the decoder cannot end it early.
_INTRON_ENDS: dict[str, tuple[str, str]] = {"u2": ("gt", "ag"), "u12": ("at", "ac")}
_INTRON_FILLER = "ttgt"

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

_LONG_CODONS = 84  # 252 bp of CDS, comfortably over the 201 bp filter threshold
_SHORT_CODONS = 20  # 60 bp of CDS, which min_cds_length removes
_INTRON_BP = 48
# Every gap is this long, which is more than the longest exon beside it: a
# decode that strayed into a spacer would pay ten nats a base for it, so no
# gene can reach across one.
_SPACER_BP = 384


def _reverse_complement(sequence: str) -> str:
    """Reverse-complement a sequence, preserving soft-mask case."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _intron(kind: str, length: int) -> str:
    """One lowercase intron of ``length`` bp with real ``kind`` boundaries.

    Args:
        kind: ``"u2"`` for ``gt``..``ag``, ``"u12"`` for the minor spliceosome's
            ``at``..``ac``, which only the extended decode lane can accept.
        length: Total intron length in bp, boundaries included.

    Returns:
        The intron sequence.
    """
    donor, acceptor = _INTRON_ENDS[kind]
    body_length = length - len(donor) - len(acceptor)
    filler = (_INTRON_FILLER * length)[:body_length]
    return f"{donor}{filler}{acceptor}"


def _gene(
    *,
    codons: int,
    intron: tuple[int, str, int] | None = None,
    planted_stop_codon: int | None = None,
) -> str:
    """One plus-strand gene: ``ATG``, a body, a stop, optionally split by an intron.

    Args:
        codons: CDS length in codons, start and stop codons included.
        intron: ``(cds_offset, kind, length)`` to splice into the CDS at
            ``cds_offset`` bp, or None for a single-exon gene. The offset need
            not be codon-aligned; the phase bookkeeping is the decoder's job.
        planted_stop_codon: Index into the body of a codon to replace with a
            stop, which makes the gene one ``apply_biological_filters`` removes.

    Returns:
        The sequence, uppercase apart from the intron.
    """
    body = [_BODY_CODONS[i % len(_BODY_CODONS)] for i in range(codons - 2)]
    if planted_stop_codon is not None:
        body[planted_stop_codon] = _PLANTED_STOP_CODON
    cds = f"{_START_CODON}{''.join(body)}{_STOP_CODON}"
    if intron is None:
        return cds
    offset, kind, length = intron
    return f"{cds[:offset]}{_intron(kind, length)}{cds[offset:]}"


def _contig(length: int, blocks: dict[int, str]) -> str:
    """Place ``blocks`` at their genomic starts in an intergenic ``N`` background.

    Args:
        length: Contig length in bp.
        blocks: Genomic start -> sequence, none of them overlapping.

    Returns:
        The contig sequence.

    Raises:
        ValueError: If a block overruns the contig or another block.
    """
    pieces: list[str] = []
    cursor = 0
    for start, block in sorted(blocks.items()):
        if start < cursor or start + len(block) > length:
            raise ValueError(f"block at {start} overruns the contig or another block")
        pieces.append("N" * (start - cursor))
        pieces.append(block)
        cursor = start + len(block)
    pieces.append("N" * (length - cursor))
    return "".join(pieces)


def _synthetic_records(stride: int) -> dict[str, str]:
    """Contigs covering every geometry branch the pipeline has.

    Every gene below is designed, and the constrained decoder is expected to
    recover exactly it. Coordinates are 0-based starts.

    Args:
        stride: Decoded bp per window, i.e. the kept-region length.

    Returns:
        Record id to sequence.
    """
    long_u2 = _gene(codons=_LONG_CODONS, intron=(120, "u2", _INTRON_BP))
    long_u12 = _gene(codons=_LONG_CODONS, intron=(122, "u12", _INTRON_BP))
    head_gene = _gene(codons=_LONG_CODONS, intron=(121, "u2", _INTRON_BP))
    seam_gene = _gene(codons=_LONG_CODONS, intron=(180, "u2", _INTRON_BP))
    short_gene = _gene(codons=_SHORT_CODONS)
    # Survives min_cds_length at 252 bp of CDS, so only the in-frame-stop filter
    # can remove it.
    stop_gene = _gene(codons=_LONG_CODONS, planted_stop_codon=40)
    minus_gene = _reverse_complement(
        _gene(codons=_LONG_CODONS, intron=(150, "u2", _INTRON_BP))
    )
    long_bp = len(long_u2)

    return {
        # Shorter than one kept region: a single N-padded window.
        "chrShort": _contig(
            3 * _SPACER_BP + long_bp + len(short_gene),
            {_SPACER_BP: long_u2, 2 * _SPACER_BP + long_bp: short_gene},
        ),
        # A gene straddling the first seam between two kept regions, which is
        # where an off-by-crop error in the retained-slice bookkeeping would
        # surface. The second window decodes it whole, having started far enough
        # ahead of it, and the first window has to hand over at the seam.
        "chrSeam": _contig(
            stride + 2 * _SPACER_BP + long_bp + len(short_gene) - 96,
            {
                _SPACER_BP: long_u12,
                stride - 96: seam_gene,
                stride + _SPACER_BP + long_bp - 96: short_gene,
            },
        ),
        # A single window at the golden stride. Carries the only minus-strand
        # gene and the only gene with a planted in-frame stop.
        "chrTail": _contig(
            4 * _SPACER_BP + len(stop_gene) + len(minus_gene) + len(short_gene),
            {
                _SPACER_BP: stop_gene,
                2 * _SPACER_BP + len(stop_gene): minus_gene,
                3 * _SPACER_BP + len(stop_gene) + len(minus_gene): short_gene,
            },
        ),
        # A gene at bp 0, only reachable through the left-padded head window
        # under a cropping geometry, plus a run of folded IUPAC codes, which
        # must read as intergenic rather than invent anything.
        "chrHead": _contig(
            2 * _SPACER_BP + 2 * long_bp,
            {
                0: head_gene,
                _SPACER_BP + long_bp: long_u12,
                _SPACER_BP + 2 * long_bp + 200: "RYKMSWBDHV",
            },
        ),
    }


# The probe run and how far past the grid point it starts. 512 bp clears
# ``min_cds_length`` (201 bp), so the gene the probe raises survives the shipped
# filters and reaches the output; the 128 bp lead keeps the run clear of the
# kept-region edge, whose outermost bases the motif mask vetoes outright for
# want of context.
_PROBE_BP = 512
_PROBE_LEAD_BP = 128


def _probe_records(stride: int) -> dict[str, str]:
    """Two contigs that tell the correct window plan from an unclamped one.

    ``chrProbe`` is a single window with its probe run in the leading half, so
    it yields a plus-strand gene: the liveness control. ``chrShift`` is
    ``stride + stride // 3`` bp, so its last window is shifted back, and its
    probe run sits just past the first grid point: in the trailing half of the
    shifted window (minus-strand gene), or the leading half without the shift
    (plus-strand gene).

    Args:
        stride: Decoded bp per window, i.e. the kept-region length.

    Returns:
        Record id to sequence.
    """
    probe = _PROBE_BASE * _PROBE_BP
    tail_bp = stride // 3
    return {
        "chrProbe": _contig(stride // 2, {_PROBE_LEAD_BP: probe}),
        "chrShift": _contig(stride + tail_bp, {stride + _PROBE_LEAD_BP: probe}),
    }


@pytest.fixture
def golden_geometry() -> tuple[int, int]:
    """The golden's ``(crop, stride)``, derived from the real geometry helper."""
    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        compute_window_geometry,
    )

    crop = compute_window_geometry(GOLDEN_SEQUENCE_LENGTH, GOLDEN_KEEP_FRACTION)
    return crop, GOLDEN_SEQUENCE_LENGTH - 2 * crop


@pytest.fixture
def synthetic_records(golden_geometry: tuple[int, int]) -> dict[str, str]:
    """The synthetic genome's sequences, keyed by record id."""
    _, stride = golden_geometry
    records = _synthetic_records(stride)
    # The goldens stay unchanged only while no golden contig holds the probe.
    assert not any(_PROBE_BASE in sequence for sequence in records.values())
    return records


@pytest.fixture
def synthetic_genome(tmp_path: Path, synthetic_records: dict[str, str]) -> list[str]:
    """Write the synthetic genome and return the configured input paths.

    Two contigs go in one multi-record file and two in a directory of
    per-contig FASTAs, so ``resolve_fasta_paths``' file and directory branches
    are both on the golden's path.
    """
    multi = tmp_path / "pair.fa"
    multi.write_text(
        "".join(
            f">{name} some description\n{synthetic_records[name]}\n"
            for name in ("chrShort", "chrSeam")
        )
    )
    per_contig = tmp_path / "contigs"
    per_contig.mkdir()
    for name in ("chrTail", "chrHead"):
        (per_contig / f"{name}.fasta").write_text(
            f">{name}\n{synthetic_records[name]}\n"
        )

    return [str(multi), str(per_contig)]


@pytest.fixture
def probe_records(golden_geometry: tuple[int, int]) -> dict[str, str]:
    """The window-plan probe contigs, keyed by record id."""
    _, stride = golden_geometry
    return _probe_records(stride)


@pytest.fixture
def probe_spans(golden_geometry: tuple[int, int]) -> dict[str, tuple[int, int]]:
    """Each probe run's 1-based inclusive GFF3 span, keyed by record id."""
    _, stride = golden_geometry
    return {
        "chrProbe": (_PROBE_LEAD_BP + 1, _PROBE_LEAD_BP + _PROBE_BP),
        "chrShift": (
            stride + _PROBE_LEAD_BP + 1,
            stride + _PROBE_LEAD_BP + _PROBE_BP,
        ),
    }


@pytest.fixture
def probe_genome(tmp_path: Path, probe_records: dict[str, str]) -> list[str]:
    """Write the probe contigs and return the configured input paths."""
    fasta = tmp_path / "probe.fa"
    fasta.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in probe_records.items())
    )
    return [str(fasta)]


@pytest.fixture
def stub_model_load(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace the checkpoint load with the deterministic stub head.

    Only the load is replaced. Everything the entry point does afterwards runs
    unmodified: the inference builder, the real Viterbi, stitching, filters and
    the writer. That is what makes the golden a regression test on the pipeline
    rather than on a mock. ``load_annotation_model`` itself gets its real
    coverage in the reproducibility tests, which load a genuine checkpoint.
    """
    from nucleotide_transformer_v3.gene_annotation import annotate_genome

    tokenizer = get_ntv3_tokenizer(include_softmasked_tokens=True)
    model = DeterministicAnnotationStub()
    def _fake_load(config, sharding):
        keep_fraction = config.keep_target_center_fraction or GOLDEN_KEEP_FRACTION
        return model, tokenizer, [], keep_fraction

    monkeypatch.setattr(annotate_genome, "load_annotation_model", _fake_load)


@pytest.fixture
def production_config(tmp_path: Path):
    """Build a golden-geometry config with the given overrides."""

    def _build(**overrides: object) -> ProductionAnnotationConfig:
        base: dict[str, object] = {
            "fasta_paths": [],
            # Local-looking: the loader is monkeypatched, but the config still goes
        # through the real validation, which rejects a URI here.
        "weights_dir": "/stub/checkpoint",
            "output_path": str(tmp_path / "out.gff3"),
            "decoder": "hmm_constrained",
            "sequence_length": GOLDEN_SEQUENCE_LENGTH,
            "keep_target_center_fraction": GOLDEN_KEEP_FRACTION,
            "decode_chunk_len": GOLDEN_DECODE_CHUNK_LEN,
            "batch_size": 2,
            "use_bf16_compute": False,
        }
        base.update(overrides)
        return ProductionAnnotationConfig(**base)  # type: ignore[arg-type]

    return _build


@pytest.fixture
def parts_dir(tmp_path: Path) -> Path:
    """A clean part-file staging directory."""
    return Path(tmp_path / "parts")
