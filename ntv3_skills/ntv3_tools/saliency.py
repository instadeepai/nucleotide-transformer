"""Per-track saliency maps from the token-embedding gradient.

The public NTv3 tracks model takes tokens, so ``token_gradients`` gives the
gradient with respect to each position's one-hot token. Saliency reads it at
A, C, G, and T in the reference base's case (lowercase bases are repeats and
keep their soft-mask tokens), centres it across the four bases (their mean is
subtracted), and keeps it at the reference base: a centred gradient times input.
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ntv3_tools import plotting
from ntv3_tools.locus import LocusRequest
from ntv3_tools.model import (
    LoadedModel,
    bigwig_track_names,
    load_model,
    token_gradients,
)
from ntv3_tools.plan import SaliencyLogoConfig
from ntv3_tools.registry import resolve_track
from ntv3_tools.sequence import crop_to_multiple, fetch_ucsc_sequence

logger = logging.getLogger("ntv3_tools")

_BASES = "ACGT"


def _one_hot(sequence: str) -> np.ndarray:
    """Encode ACGT as a one-hot matrix.

    Args:
        sequence: DNA string.

    Returns:
        Array of shape (L, 4).
    """
    index: dict[str, int] = {_BASES[i]: i for i in range(4)}
    out = np.zeros((len(sequence), 4), dtype=np.float32)
    for pos, base in enumerate(sequence.upper()):
        if base in index:
            out[pos, index[base]] = 1.0
    return out


def _track_saliency(
    loaded: LoadedModel,
    sequence: str,
    species: str,
    track_index: int,
) -> np.ndarray:
    """Centred gradient of the mean track prediction times the reference one-hot.

    The gradient at A, C, G, and T is centred (their mean subtracted) before it
    is kept at the reference base.

    Args:
        loaded: Loaded JAX model.
        sequence: Cropped DNA string, soft mask kept.
        species: Species string.
        track_index: Bigwig track index.

    Returns:
        Signed per-base map, shape (L, 4). Non-reference bases and N are zero.
    """
    gradients = token_gradients(
        loaded=loaded, sequence=sequence, species=species, track_index=track_index
    )
    upper = [loaded.tokenizer.token_to_id(base) for base in _BASES]
    lower = [loaded.tokenizer.token_to_id(base.lower()) for base in _BASES]
    repeat = np.array([base.islower() for base in sequence])[:, None]
    per_base = np.where(repeat, gradients[:, lower], gradients[:, upper])
    centered = per_base - per_base.mean(axis=1, keepdims=True)
    return _one_hot(sequence=sequence) * centered


@dataclass(frozen=True)
class SaliencyRequest:
    """Inputs for per-track saliency generation."""

    model_name: str
    locus: LocusRequest
    track_ids: list[str]
    track_labels: list[str]
    out_dir: str
    logo: SaliencyLogoConfig
    device: str = "auto"
    loaded: LoadedModel | None = None
    sequence: str | None = None


def saliency_per_track(request: SaliencyRequest) -> list[Path]:
    """Write one saliency map per requested track.

    Args:
        request: Model, locus, tracks, labels, output, and runtime dependencies.

    Returns:
        Paths of written ``.npy`` files.
    """
    if len(request.track_ids) != len(request.track_labels):
        raise ValueError("Track identifiers and plot labels must have the same length.")
    loaded = request.loaded or load_model(
        model_name=request.model_name, device=request.device
    )
    sequence = request.sequence or crop_to_multiple(
        fetch_ucsc_sequence(
            assembly=request.locus.assembly,
            chrom=request.locus.chrom,
            start=request.locus.start,
            end=request.locus.end,
            keep_soft_mask=True,
        )
    )
    names = bigwig_track_names(loaded=loaded, species=request.locus.species)
    accessions = [resolve_track(name_or_accession=name) for name in request.track_ids]
    missing = [acc for acc in accessions if acc not in names]
    if missing:
        raise ValueError(f"Tracks not available for this model/species: {missing}.")
    root = Path(request.out_dir) / "interpret"
    root.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for accession, track_label in zip(accessions, request.track_labels, strict=True):
        track_index = names.index(accession)
        logger.info(f"Saliency for {accession} (index {track_index})")
        scores = _track_saliency(
            loaded=loaded,
            sequence=sequence,
            species=request.locus.species,
            track_index=track_index,
        )
        path = root / f"{accession}.npy"
        np.save(path, scores.astype(np.float16))
        figure_relative_path = Path("figures") / f"saliency_{accession}.png"
        logo_window = plotting.plot_saliency(
            plot=plotting.SaliencyPlot(
                scores=scores,
                sequence=sequence,
                genomic_start=request.locus.start,
                chrom=request.locus.chrom,
                assembly=request.locus.assembly,
                label=f"{track_label} saliency",
                out_path=str(Path(request.out_dir) / figure_relative_path),
                window_bp=request.logo.window_bp,
                center_bp=request.logo.center_bp,
            )
        )
        (root / f"{accession}.json").write_text(
            json.dumps(
                {
                    "accession": accession,
                    "track_index": track_index,
                    "chrom": request.locus.chrom,
                    "start": request.locus.start,
                    "end": request.locus.start + len(sequence),
                    "shape": list(scores.shape),
                    "method": "token_gradient_x_input",
                    "checkpoint": loaded.checkpoint,
                    "seq": sequence,
                    "logo": {
                        "selection": logo_window.selection,
                        "center_bp": logo_window.center_bp,
                        "start": logo_window.start,
                        "end": logo_window.end,
                        "window_bp": logo_window.end - logo_window.start,
                        "png": str(figure_relative_path),
                        "pdf": str(figure_relative_path.with_suffix(".pdf")),
                    },
                },
                indent=2,
            )
            + "\n"
        )
        written.append(path)
    return written
