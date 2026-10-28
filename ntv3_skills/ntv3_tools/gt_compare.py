"""Pure ground-truth bin alignment and correlation policy."""

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class BigWigValues(Protocol):
    """The pyBigWig value-reading boundary used by alignment tests."""

    def values(
        self,
        chrom: str,
        start: int,
        end: int,
        numpy: bool,
    ) -> np.ndarray:
        """Read base-resolution values over a half-open genomic interval."""


@dataclass(frozen=True)
class CorrelationResult:
    """Finite-pair track agreement.

    Attributes:
        n_valid: Number of bins finite in both arrays.
        pearson: Pearson correlation, or NaN when undefined.
        spearman: Spearman correlation, or NaN when undefined.
        caveat: Explicit reason when either correlation is undefined.
    """

    n_valid: int
    pearson: float
    spearman: float
    caveat: str


def aggregate_base_values(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Aggregate base values into equal model bins.

    Args:
        values: Base-resolution values, shape (prediction_span,).
        n_bins: Number of model output bins.

    Returns:
        Mean signal per model bin, shape (n_bins,). Bins without finite values
        are NaN.
    """
    base_values = np.asarray(values, dtype=np.float64)
    if base_values.ndim != 1 or n_bins <= 0:
        raise ValueError(
            f"Expected one-dimensional values and positive n_bins; got "
            f"{base_values.shape} and {n_bins}."
        )
    if base_values.size % n_bins != 0:
        raise ValueError(
            f"Prediction span {base_values.size} does not divide into {n_bins} bins."
        )
    reshaped = base_values.reshape(n_bins, base_values.size // n_bins)
    finite = np.isfinite(reshaped)
    counts = finite.sum(axis=1)
    sums = np.where(finite, reshaped, 0.0).sum(axis=1)
    means = np.full(n_bins, np.nan, dtype=np.float64)
    np.divide(sums, counts, out=means, where=counts > 0)
    return means


def read_binned_bigwig(
    bigwig: BigWigValues,
    chrom: str,
    start: int,
    end: int,
    n_bins: int,
) -> np.ndarray:
    """Read the exact prediction span and aggregate it to model bins.

    Args:
        bigwig: Open pyBigWig-compatible reader.
        chrom: Prediction chromosome.
        start: Zero-based prediction start.
        end: Half-open prediction end.
        n_bins: Number of model output bins.

    Returns:
        Experimental values, shape (n_bins,).
    """
    if start < 0 or end <= start:
        raise ValueError(f"Invalid prediction span {chrom}:{start}-{end}.")
    values = np.asarray(
        bigwig.values(chrom, start, end, numpy=True),
        dtype=np.float64,
    )
    expected = end - start
    if values.shape != (expected,):
        raise ValueError(
            f"BigWig returned shape {values.shape} for span length {expected}."
        )
    return aggregate_base_values(values=values, n_bins=n_bins)


def finite_pair_correlations(
    predicted: np.ndarray,
    observed: np.ndarray,
) -> CorrelationResult:
    """Compute stable Pearson and Spearman correlations on finite pairs.

    Args:
        predicted: Model values, shape (n_bins,).
        observed: Experimental values, shape (n_bins,).

    Returns:
        Correlations and an explicit undefined-value caveat.
    """
    predicted_values = np.asarray(predicted, dtype=np.float64)
    observed_values = np.asarray(observed, dtype=np.float64)
    if predicted_values.ndim != 1 or predicted_values.shape != observed_values.shape:
        raise ValueError(
            f"Correlation arrays must have one matching shape; got "
            f"{predicted_values.shape} and {observed_values.shape}."
        )
    finite = np.isfinite(predicted_values) & np.isfinite(observed_values)
    paired_predicted = predicted_values[finite]
    paired_observed = observed_values[finite]
    n_valid = int(finite.sum())
    if n_valid < 2:
        return CorrelationResult(
            n_valid=n_valid,
            pearson=np.nan,
            spearman=np.nan,
            caveat="Undefined correlation: fewer than two finite paired bins.",
        )
    constant_reasons: list[str] = []
    if np.ptp(paired_predicted) == 0:
        constant_reasons.append("constant predicted signal")
    if np.ptp(paired_observed) == 0:
        constant_reasons.append("constant observed signal")
    if constant_reasons:
        return CorrelationResult(
            n_valid=n_valid,
            pearson=np.nan,
            spearman=np.nan,
            caveat=f"Undefined correlation: {' and '.join(constant_reasons)}.",
        )
    pearson = _centered_correlation(
        left=paired_predicted,
        right=paired_observed,
    )
    spearman = _centered_correlation(
        left=_average_ranks(values=paired_predicted),
        right=_average_ranks(values=paired_observed),
    )
    return CorrelationResult(
        n_valid=n_valid,
        pearson=pearson,
        spearman=spearman,
        caveat="",
    )


def _centered_correlation(left: np.ndarray, right: np.ndarray) -> float:
    """Compute correlation from centered float64 vectors."""
    left_centered = left - left.mean()
    right_centered = right - right.mean()
    denominator = np.linalg.norm(left_centered) * np.linalg.norm(right_centered)
    if denominator == 0:
        return float("nan")
    value = np.dot(left_centered, right_centered) / denominator
    return float(np.clip(value, -1.0, 1.0))


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Assign stable average ranks, including tied values."""
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(values.size, dtype=np.float64)
    cursor = 0
    while cursor < values.size:
        stop = cursor + 1
        while stop < values.size and sorted_values[stop] == sorted_values[cursor]:
            stop += 1
        ranks[order[cursor:stop]] = (cursor + stop - 1) / 2.0 + 1.0
        cursor = stop
    return ranks
