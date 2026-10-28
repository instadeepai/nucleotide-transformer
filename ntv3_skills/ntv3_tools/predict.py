"""Track prediction at a locus.

Fetch the window sequence from UCSC (soft mask kept), run the NTv3 tracks
model on only the requested head columns, and read bigwig-track logits (shape
``(out_len, n_tracks)`` for one sequence).
"""

import logging
from dataclasses import dataclass

import numpy as np

from ntv3_tools.locus import LocusRequest
from ntv3_tools.model import (
    LoadedModel,
    as_numpy,
    bigwig_track_names,
    center_fraction,
    load_model,
    prediction_region,
    run_forward,
)
from ntv3_tools.registry import head_indices, resolve_track
from ntv3_tools.sequence import crop_to_multiple, fetch_ucsc_sequence

logger = logging.getLogger("ntv3_tools")


@dataclass
class TrackPrediction:
    """Predicted bigwig tracks for a set of selected accessions.

    Attributes:
        tracks: Selected track logits, shape (n_selected, out_len).
        accessions: Accession per selected track, length n_selected.
        prediction_start: Genomic start of the predicted (center 37.5%) region.
        prediction_end: Genomic end of the predicted region.
    """

    tracks: np.ndarray
    accessions: list[str]
    prediction_start: int
    prediction_end: int


@dataclass(frozen=True)
class TrackPredictionRequest:
    """Inputs for prediction with an optional pre-fetched reference."""

    model_name: str
    locus: LocusRequest
    track_ids: list[str]
    device: str = "auto"
    loaded: LoadedModel | None = None
    sequence: str | None = None


def _bigwig_logits(out) -> np.ndarray:
    """Extract single-sequence bigwig logits as a NumPy array of shape (out_len, n_tracks).

    Args:
        out: A model forward output mapping.

    Returns:
        Array of shape (out_len, n_tracks).
    """
    logits = as_numpy(out["bigwig_tracks_logits"])
    if logits.ndim == 3:
        return logits[0]
    return logits


def _validate_tracks(names: list[str], accessions: list[str]) -> list[int]:
    """Map requested head ids to output indices; a missing id raises and is named.

    Args:
        names: Ordered accession list from the model config.
        accessions: Requested model-head ids.

    Returns:
        Output track indices aligned with ``accessions``.
    """
    return head_indices(names=names, heads=accessions)


def predict_tracks(
    request: TrackPredictionRequest,
) -> TrackPrediction:
    """Predict tracks with optional model and sequence reuse.

    Args:
        request: Model, locus, tracks, runtime, and optional reference sequence.

    Returns:
        The selected track logits with their accessions and predicted genomic region.
    """
    loaded = request.loaded or load_model(
        model_name=request.model_name,
        device=request.device,
    )
    seq = request.sequence
    if seq is None:
        seq = crop_to_multiple(
            fetch_ucsc_sequence(
                assembly=request.locus.assembly,
                chrom=request.locus.chrom,
                start=request.locus.start,
                end=request.locus.end,
                keep_soft_mask=True,
            )
        )
    names = bigwig_track_names(loaded=loaded, species=request.locus.species)
    accessions = [
        resolve_track(name_or_accession=track_id) for track_id in request.track_ids
    ]
    indices = _validate_tracks(names=names, accessions=accessions)
    out = run_forward(
        loaded=loaded,
        sequence=seq,
        species=request.locus.species,
        track_indices=indices,
    )
    logits = _bigwig_logits(out)
    if logits.ndim != 2 or logits.shape[1] != len(indices):
        raise RuntimeError(
            f"Forward returned bigwig logits of shape {logits.shape}; expected "
            f"(positions, {len(indices)}) for the requested track_indices."
        )
    pred_start, pred_end = prediction_region(
        start=request.locus.start,
        end=request.locus.end,
        center_frac=center_fraction(loaded=loaded),
    )
    logger.info(
        f"Predicted {len(accessions)} tracks over "
        f"{request.locus.chrom}:{pred_start}-{pred_end}"
    )
    return TrackPrediction(
        tracks=logits.T,
        accessions=accessions,
        prediction_start=pred_start,
        prediction_end=pred_end,
    )


def predict_tracks_at_locus(
    model_name: str,
    locus: LocusRequest,
    track_ids: list[str],
    device: str = "auto",
    loaded: LoadedModel | None = None,
) -> TrackPrediction:
    """Predict tracks through the original caller-compatible interface.

    Args:
        model_name: ``NTv3_tracks``.
        locus: Genomic window.
        track_ids: Track aliases or model-head identifiers.
        device: JAX device selector.
        loaded: Optional model bundle to reuse.

    Returns:
        The selected track logits and predicted genomic region.
    """
    return predict_tracks(
        request=TrackPredictionRequest(
            model_name=model_name,
            locus=locus,
            track_ids=track_ids,
            device=device,
            loaded=loaded,
        )
    )
