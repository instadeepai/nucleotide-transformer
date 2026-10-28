"""Science-environment ENCODE BigWig comparison."""

import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import numpy as np

try:
    import pyBigWig
except ImportError:  # Installed by setup.py; remote reads need a libcurl build.
    pyBigWig = None

from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.gt_compare import finite_pair_correlations, read_binned_bigwig

logger = logging.getLogger("ntv3_tools")
_SCALE_CAVEAT = (
    "Model logits and experimental signal differ; correlations compare shape."
)


def bigwig_reader_problem() -> str | None:
    """Return why this Python cannot read remote ENCODE BigWigs, or None.

    The PyPI ``pyBigWig`` wheel imports but is built without libcurl, so it
    cannot open an HTTPS URL; ``pyBigWig.remote`` says which build this is.
    Reads also need its NumPy support (``pyBigWig.numpy``).
    """
    if pyBigWig is None:
        return "pyBigWig is not installed in this Python"
    if not getattr(pyBigWig, "remote", 0):
        return (
            "this pyBigWig build cannot open remote files (pyBigWig.remote is 0; "
            "the PyPI wheel is built without libcurl)"
        )
    if not getattr(pyBigWig, "numpy", 0):
        return "this pyBigWig build has no NumPy support (pyBigWig.numpy is 0)"
    return None


def can_read_bigwig() -> bool:
    """Return whether this Python can read remote ENCODE BigWigs."""
    return bigwig_reader_problem() is None


@dataclass(frozen=True)
class GroundTruthComparison:
    """Inputs for prediction-aligned experimental comparison.

    Attributes:
        run_dir: Approved run directory containing lookup artifacts.
        chrom: Prediction chromosome.
        start: Zero-based prediction start.
        end: Half-open prediction end.
        bin_size: Model output bin size in base pairs.
    """

    run_dir: Path
    chrom: str
    start: int
    end: int
    bin_size: int


def write_ground_truth_comparison(
    comparison: GroundTruthComparison,
    predicted: np.ndarray,
) -> tuple[Path, Path]:
    """Read approved ENCODE BigWigs and save aligned agreement artifacts.

    A track with no ENCODE experiment gets a NaN row in ``gt/observed.npy`` and
    no agreement row; the other tracks are still compared.

    Args:
        comparison: Prediction span and approved run directory.
        predicted: Model output values, shape (n_tracks, n_bins).

    Returns:
        Paths to ``gt/observed.npy`` and ``tables/track_agreement.csv``.
    """
    predictions = np.asarray(predicted, dtype=np.float64)
    if predictions.ndim != 2 or not np.isfinite(predictions).all():
        raise ValueError(
            f"Predictions must be finite with shape (n_tracks, n_bins); "
            f"got {predictions.shape}."
        )
    prediction_span = comparison.end - comparison.start
    if (
        comparison.bin_size <= 0
        or prediction_span != predictions.shape[1] * comparison.bin_size
    ):
        raise ValueError(
            f"Prediction span {prediction_span} does not match "
            f"{predictions.shape[1]} bins at {comparison.bin_size} bp."
        )
    observed_path = comparison.run_dir / "gt" / "observed.npy"
    agreement_path = comparison.run_dir / "tables" / "track_agreement.csv"
    for path in (observed_path, agreement_path):
        if path.exists():
            raise ArtifactError(f"Refusing to overwrite generated GT artifact: {path}.")

    resolve = _read_mapping(comparison.run_dir / "track_resolve.json")
    manifest = _read_list(comparison.run_dir / "gt" / "manifest.json")
    track_rows = resolve.get("tracks")
    if not isinstance(track_rows, list) or len(track_rows) != predictions.shape[0]:
        raise RequestError(
            "track_resolve.json track count does not match prediction rows."
        )
    sources_by_accession = {
        str(row.get("accession")): row
        for row in manifest
        if isinstance(row, dict) and row.get("accession")
    }
    observed_rows: list[np.ndarray] = []
    agreement_rows: list[dict[str, object]] = []
    for index, track_row in enumerate(track_rows):
        if not isinstance(track_row, dict):
            raise RequestError(f"Resolved track row {index} is not a mapping.")
        accession = str(track_row.get("gt_accession") or "")
        if not accession:
            observed_rows.append(np.full(predictions.shape[1], np.nan))
            continue
        source = sources_by_accession.get(accession)
        if source is None:
            raise RequestError(
                f"No approved gt/manifest.json source for GT accession "
                f"{accession or '<missing>'}."
            )
        url = str(source.get("url") or "")
        _validate_public_https(url=url, accession=accession)
        observed = _read_remote_bigwig(
            url=url,
            comparison=comparison,
            n_bins=predictions.shape[1],
        )
        agreement = finite_pair_correlations(
            predicted=predictions[index],
            observed=observed,
        )
        observed_rows.append(observed)
        caveats = [_SCALE_CAVEAT]
        if agreement.caveat:
            caveats.append(agreement.caveat)
        agreement_rows.append(
            {
                "track_query": str(track_row.get("query") or ""),
                "model_head": str(track_row.get("model_head") or ""),
                "gt_accession": accession,
                "file_accession": str(source.get("encff") or ""),
                "n_valid": agreement.n_valid,
                "pearson": agreement.pearson,
                "spearman": agreement.spearman,
                "caveat": " ".join(caveats),
            }
        )

    if not agreement_rows:
        raise RequestError("No requested track has an ENCODE experiment to compare.")
    observed_array = np.stack(observed_rows, axis=0)
    observed_path.parent.mkdir(parents=True, exist_ok=True)
    agreement_path.parent.mkdir(parents=True, exist_ok=True)
    with observed_path.open("xb") as handle:
        np.save(handle, observed_array)
    with agreement_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(agreement_rows[0]))
        writer.writeheader()
        writer.writerows(agreement_rows)
    logger.info(
        f"Wrote prediction-aligned GT for {len(agreement_rows)} tracks to "
        f"{observed_path}"
    )
    return observed_path, agreement_path


def _read_remote_bigwig(
    url: str,
    comparison: GroundTruthComparison,
    n_bins: int,
) -> np.ndarray:
    """Open one remote BigWig and read the exact prediction interval."""
    try:
        bigwig = pyBigWig.open(url)
        if bigwig is None:
            raise RuntimeError("pyBigWig returned no reader.")
        try:
            return read_binned_bigwig(
                bigwig=bigwig,
                chrom=comparison.chrom,
                start=comparison.start,
                end=comparison.end,
                n_bins=n_bins,
            )
        finally:
            bigwig.close()
    except (RuntimeError, OSError, ValueError) as exc:
        raise RequestError(
            f"Could not read public GT BigWig over "
            f"{comparison.chrom}:{comparison.start}-{comparison.end}: {url}. "
            "No prediction-vs-ground-truth comparison was claimed."
        ) from exc


def _validate_public_https(url: str, accession: str) -> None:
    """Require the approved GT source to be a public HTTPS URL."""
    parsed = urlparse(url)
    if parsed.scheme != "https" or not parsed.netloc:
        raise RequestError(
            f"GT accession {accession} does not have a public HTTPS BigWig URL "
            "in gt/manifest.json."
        )


def _read_mapping(path: Path) -> dict[str, Any]:
    """Read one required JSON mapping."""
    if not path.exists():
        raise RequestError(f"Required lookup artifact is missing: {path}.")
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise RequestError(f"Expected a JSON mapping at {path}.")
    return payload


def _read_list(path: Path) -> list[Any]:
    """Read one required JSON list."""
    if not path.exists():
        raise RequestError(f"Required GT artifact is missing: {path}.")
    payload = json.loads(path.read_text())
    if not isinstance(payload, list):
        raise RequestError(f"Expected a JSON list at {path}.")
    return payload
