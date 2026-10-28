"""Pure policy and artifact I/O for single-SNV track disruption."""

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ntv3_tools.errors import ArtifactError, RequestError

_DNA_BASES = "ACGT"
_ARTIFACT_NAMES = (
    "variant.json",
    "ref_tracks.npy",
    "alt_tracks.npy",
    "delta_tracks.npy",
    "disrupt_top.tsv",
    "effect_percentile.json",
)


@dataclass(frozen=True)
class Snv:
    """A biallelic single-nucleotide variant.

    Attributes:
        chrom: Chromosome name.
        position: One-based genomic position.
        ref: Reference allele.
        alt: Alternate allele.
    """

    chrom: str
    position: int
    ref: str
    alt: str


@dataclass(frozen=True)
class NullSamplingRequest:
    """Configuration for deterministic null SNV sampling.

    Attributes:
        sequence: UCSC reference sequence, shape (window_bp,).
        chrom: Approved window chromosome.
        window_start: Approved zero-based window start.
        query: Query SNV to exclude.
        null_n: Number of null SNVs.
        seed: Explicit NumPy random seed.
    """

    sequence: str
    chrom: str
    window_start: int
    query: Snv
    null_n: int
    seed: int


@dataclass(frozen=True)
class DisruptArtifacts:
    """Values written under a run's disrupt directory.

    Attributes:
        variant: Variant and inference provenance.
        ref_tracks: Query reference tracks, shape (n_tracks, n_bins).
        alt_tracks: Query alternate tracks, shape (n_tracks, n_bins).
        delta_tracks: Alternate-minus-reference tracks, shape (n_tracks, n_bins).
        top_rows: Track effect rows sorted by absolute effect.
        effect_percentile: Empirical null percentile metadata.
        log2fc_png: Paper-style per-bin log-fold-change PNG bytes.
        log2fc_pdf: Paper-style per-bin log-fold-change PDF bytes.
    """

    variant: dict[str, Any]
    ref_tracks: np.ndarray
    alt_tracks: np.ndarray
    delta_tracks: np.ndarray
    top_rows: list[dict[str, Any]]
    effect_percentile: dict[str, Any]
    log2fc_png: bytes
    log2fc_pdf: bytes


def validate_snv(
    variant: Snv,
    window_chrom: str,
    window_start: int,
    sequence: str,
) -> int:
    """Validate an SNV against an approved window and reference sequence.

    Args:
        variant: One-based query SNV.
        window_chrom: Approved window chromosome.
        window_start: Approved zero-based window start.
        sequence: UCSC reference sequence, shape (window_bp,).

    Returns:
        Zero-based variant offset within the window.

    Raises:
        RequestError: The variant is invalid or does not match the reference.
    """
    ref = variant.ref.upper()
    alt = variant.alt.upper()
    if variant.chrom != window_chrom:
        raise RequestError(
            f"Variant chromosome {variant.chrom} does not match the approved "
            f"window chromosome {window_chrom}."
        )
    if len(ref) != 1 or len(alt) != 1:
        raise RequestError("Disrupt requires a single-base ref and alt allele.")
    if ref not in _DNA_BASES or alt not in _DNA_BASES:
        raise RequestError("Disrupt alleles must each be one of A/C/G/T.")
    if ref == alt:
        raise RequestError("Reference and alternate alleles must be different.")
    genomic_position = variant.position - 1
    window_end = window_start + len(sequence)
    if genomic_position < window_start or genomic_position >= window_end:
        raise RequestError(
            f"Variant position {variant.position} is outside the approved "
            f"half-open window [{window_start}, {window_end}) after conversion "
            "to zero-based coordinates."
        )
    offset = genomic_position - window_start
    observed = sequence[offset].upper()
    if observed != ref:
        raise RequestError(
            f"Variant reference {ref} does not match UCSC reference {observed} "
            f"at {variant.chrom}:{variant.position}."
        )
    return offset


def apply_snv(sequence: str, offset: int, alt: str) -> str:
    """Return a sequence with one alternate allele, in the reference base's case.

    Args:
        sequence: Reference DNA sequence, soft mask kept, shape (sequence_length,).
        offset: Zero-based mutation offset.
        alt: Alternate A/C/G/T base.

    Returns:
        New DNA sequence, shape (sequence_length,).
    """
    if offset < 0 or offset >= len(sequence):
        raise RequestError(f"Mutation offset {offset} is outside the sequence.")
    alternate = alt.upper()
    if alternate not in _DNA_BASES:
        raise RequestError("Alternate allele must be one of A/C/G/T.")
    if sequence[offset].islower():
        alternate = alternate.lower()
    return f"{sequence[:offset]}{alternate}{sequence[offset + 1 :]}"


def generate_null_snvs(request: NullSamplingRequest) -> list[Snv]:
    """Sample deterministic random in-window SNVs with replacement.

    Args:
        request: Deterministic null sampling configuration.

    Returns:
        Random SNVs of length null_n.

    Raises:
        RequestError: The sample count is invalid or no candidates exist.
    """
    if request.null_n <= 0:
        raise RequestError("null_n must be positive.")
    sequence = request.sequence.upper()
    encoded = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    eligible_offsets = np.flatnonzero(np.isin(encoded, list(map(ord, _DNA_BASES))))
    if eligible_offsets.size == 0:
        raise RequestError("The approved window has no eligible null SNVs.")
    generator = np.random.default_rng(request.seed)
    normalized_query = Snv(
        chrom=request.query.chrom,
        position=request.query.position,
        ref=request.query.ref.upper(),
        alt=request.query.alt.upper(),
    )
    variants: list[Snv] = []
    while len(variants) < request.null_n:
        offset = int(eligible_offsets[generator.integers(eligible_offsets.size)])
        ref = sequence[offset]
        alternatives = [base for base in _DNA_BASES if base != ref]
        alt = alternatives[int(generator.integers(len(alternatives)))]
        candidate = Snv(
            chrom=request.chrom,
            position=request.window_start + offset + 1,
            ref=ref,
            alt=alt,
        )
        if candidate != normalized_query:
            variants.append(candidate)
    return variants


def signed_absmax(delta: np.ndarray, axis: int = 0) -> np.ndarray:
    """Select values with maximum absolute magnitude while preserving sign.

    This matches the canonical VEP ``utils_scoring._absmax`` reduction.

    Args:
        delta: Signed alternate-minus-reference effects.
        axis: Axis to reduce.

    Returns:
        Signed effects with the selected axis removed.
    """
    indices = np.argmax(np.abs(delta), axis=axis)
    expanded = np.expand_dims(indices, axis=axis)
    return np.take_along_axis(delta, expanded, axis=axis).squeeze(axis=axis)


def empirical_absolute_percentile(
    query_effect: float,
    null_effects: np.ndarray,
) -> float:
    """Compute a finite-sample empirical absolute-effect rank percentile.

    The percentile is ``100 * (1 + count(|null| <= |query|)) / (null_n + 1)``.

    Args:
        query_effect: Signed query effect.
        null_effects: Signed null effects, shape (null_n,).

    Returns:
        Percentile in the interval (0, 100].
    """
    flattened = np.asarray(null_effects).reshape(-1)
    if flattened.size == 0:
        raise RequestError("At least one null effect is required.")
    rank = 1 + int(np.count_nonzero(np.abs(flattened) <= abs(query_effect)))
    return 100.0 * rank / (flattened.size + 1)


def require_disrupt_artifacts_absent(run_dir: str) -> Path:
    """Return the target directory after checking overwrite safety.

    Args:
        run_dir: Existing run directory.

    Returns:
        Target disrupt directory.

    Raises:
        ArtifactError: A target artifact already exists.
    """
    disrupt_dir = Path(run_dir) / "disrupt"
    existing = [name for name in _ARTIFACT_NAMES if (disrupt_dir / name).exists()]
    figures_dir = Path(run_dir) / "figures"
    existing.extend(
        f"figures/{name}"
        for name in ("disrupt_log2fc.png", "disrupt_log2fc.pdf")
        if (figures_dir / name).exists()
    )
    if existing:
        raise ArtifactError(
            f"Refusing to overwrite existing disrupt artifacts: {existing}."
        )
    return disrupt_dir


def write_disrupt_artifacts(
    run_dir: str,
    artifacts: DisruptArtifacts,
) -> Path:
    """Write all disruption artifacts without overwriting any existing file.

    Args:
        run_dir: Existing run directory.
        artifacts: Values to write.

    Returns:
        Path to the created disrupt directory.

    Raises:
        ArtifactError: A target artifact already exists.
    """
    disrupt_dir = require_disrupt_artifacts_absent(run_dir=run_dir)
    if not artifacts.log2fc_png or not artifacts.log2fc_pdf:
        raise ArtifactError("Disruption requires PNG and PDF log-fold-change plots.")
    disrupt_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = Path(run_dir) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    (disrupt_dir / "variant.json").write_text(
        json.dumps(artifacts.variant, indent=2) + "\n"
    )
    np.save(disrupt_dir / "ref_tracks.npy", artifacts.ref_tracks)
    np.save(disrupt_dir / "alt_tracks.npy", artifacts.alt_tracks)
    np.save(disrupt_dir / "delta_tracks.npy", artifacts.delta_tracks)
    fieldnames = list(artifacts.top_rows[0]) if artifacts.top_rows else []
    with (disrupt_dir / "disrupt_top.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(artifacts.top_rows)
    (disrupt_dir / "effect_percentile.json").write_text(
        json.dumps(artifacts.effect_percentile, indent=2) + "\n"
    )
    (figures_dir / "disrupt_log2fc.png").write_bytes(artifacts.log2fc_png)
    (figures_dir / "disrupt_log2fc.pdf").write_bytes(artifacts.log2fc_pdf)
    return disrupt_dir
