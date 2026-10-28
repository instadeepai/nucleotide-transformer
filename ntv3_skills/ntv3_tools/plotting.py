"""Plotting helpers for NTv3 skill reports.

Track plotting mirrors tutorial 01's ``plot_tracks`` (light fill, despined). Styling follows the
impress.py conventions (Arial, thin spines, 600 dpi, embeddable fonts).
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import Affine2D

from ntv3_tools.colors import COLORS

logger = logging.getLogger("ntv3_tools")

LINE_WIDTH = 0.5
_FILL = COLORS["blue_3"]
_LINE = COLORS["blue_0"]
_NUCLEOTIDE_COLORS = {
    "A": COLORS["green_1"],
    "C": COLORS["blue_1"],
    "G": COLORS["gold_1"],
    "T": COLORS["red_1"],
}


def _format_megabases(position: float, _: int) -> str:
    """Format one genomic coordinate in megabases."""
    return f"{position / 1e6:.3f}"


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


def _style(ax: plt.Axes) -> None:
    """Apply the shared axis style (despine + thin spines/ticks).

    Args:
        ax: Axis to style in place.
    """
    sns.despine(ax=ax)
    for spine in ax.spines.values():
        spine.set_linewidth(LINE_WIDTH)
    ax.tick_params(width=LINE_WIDTH, labelsize=7)


def _save(fig: plt.Figure, out_path: str) -> str:
    """Save a figure at 600 dpi and close it.

    Args:
        fig: Figure to save.
        out_path: Destination path.

    Returns:
        The output path.
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    destination = Path(out_path)
    fig.savefig(destination, dpi=600, bbox_inches="tight")
    fig.savefig(destination.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out_path


@dataclass(frozen=True)
class TrackPlot:
    """Inputs for a publication-quality predicted-track figure."""

    tracks: np.ndarray
    labels: list[str]
    start: int
    end: int
    chrom: str
    assembly: str
    out_path: str


def plot_tracks(plot: TrackPlot) -> str:
    """Plot a stack of predicted functional tracks over the predicted region.

    Args:
        plot: Track values, labels, genomic span, and output path.

    Returns:
        The output path.
    """
    sns.set_style("white")
    tracks = plot.tracks
    if tracks.ndim != 2 or tracks.shape[0] != len(plot.labels):
        raise ValueError(
            f"Expected tracks shape (n_tracks, n_bins) aligned with labels; "
            f"got {tracks.shape} and {len(plot.labels)} labels."
        )
    n = tracks.shape[0]
    fig, axes = plt.subplots(n, 1, figsize=(8, 1.1 * n), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, row, label in zip(axes, tracks, plot.labels, strict=True):
        bin_width = (plot.end - plot.start) / row.shape[0]
        x = plot.start + (np.arange(row.shape[0]) + 0.5) * bin_width
        ax.fill_between(x, row, color=_FILL, linewidth=0)
        ax.plot(x, row, color=_LINE, linewidth=LINE_WIDTH)
        ax.set_ylabel(label, fontsize=8, rotation=0, ha="right", va="center")
        _style(ax)
    axes[-1].set_xlabel(
        f"Position on {plot.chrom} (Mb; {plot.assembly})",
        fontsize=9,
    )
    axes[-1].xaxis.set_major_formatter(FuncFormatter(_format_megabases))
    fig.align_ylabels(axes)
    return _save(fig, plot.out_path)


@dataclass(frozen=True)
class SaliencyWindow:
    """Selected sequence-logo interval.

    Attributes:
        left: Zero-based start offset in the saliency array.
        right: Exclusive end offset in the saliency array.
        start: Absolute genomic start.
        end: Absolute genomic end.
        center_bp: Genomic coordinate used to center the window.
        selection: ``strongest_peak`` or ``explicit_center``.
    """

    left: int
    right: int
    start: int
    end: int
    center_bp: int
    selection: str


def select_saliency_window(
    scores: np.ndarray,
    genomic_start: int,
    window_bp: int,
    center_bp: int | None = None,
) -> SaliencyWindow:
    """Select a fixed-width saliency window without sign cancellation.

    Args:
        scores: Signed attribution values, shape (sequence_length, 4).
        genomic_start: Absolute coordinate of the first score.
        window_bp: Requested number of bases.
        center_bp: Optional absolute genomic center.

    Returns:
        Selected array offsets and genomic coordinates.
    """
    if scores.ndim != 2 or scores.shape[1] != 4:
        raise ValueError(
            f"Expected saliency shape (sequence_length, 4), got {scores.shape}."
        )
    sequence_length = scores.shape[0]
    if sequence_length == 0:
        raise ValueError("Saliency scores must contain at least one base.")
    if window_bp <= 0:
        raise ValueError("Sequence-logo window_bp must be positive.")
    width = min(window_bp, sequence_length)
    if center_bp is None:
        peak_index = int(np.argmax(np.abs(scores).sum(axis=1)))
        selected_center = genomic_start + peak_index
        selection = "strongest_peak"
    else:
        peak_index = center_bp - genomic_start
        if not 0 <= peak_index < sequence_length:
            raise ValueError(
                f"center_bp={center_bp} is outside saliency interval "
                f"[{genomic_start}, {genomic_start + sequence_length})."
            )
        selected_center = center_bp
        selection = "explicit_center"
    left = max(0, min(peak_index - width // 2, sequence_length - width))
    right = left + width
    return SaliencyWindow(
        left=left,
        right=right,
        start=genomic_start + left,
        end=genomic_start + right,
        center_bp=selected_center,
        selection=selection,
    )


def _draw_logo(ax: plt.Axes, scores: np.ndarray) -> None:
    """Draw signed ACGT letters, each as tall as its score, stacked from zero.

    Args:
        ax: Axes to draw on; x is the position index.
        scores: Signed scores, shape (n_positions, 4), columns A, C, G, T.
    """
    font = FontProperties(family="sans-serif", weight="bold")
    low, high = 0.0, 0.0
    for position, row in enumerate(scores):
        floors = [0.0, 0.0]
        for column in np.argsort(np.abs(row)):
            value = float(row[column])
            if not value:
                continue
            base = "ACGT"[column]
            path = TextPath((0, 0), base, size=1, prop=font)
            box = path.get_extents()
            up = value > 0
            transform = (
                Affine2D()
                .translate(-box.x0, -box.y0)
                .scale(0.9 / box.width, abs(value) / box.height)
                .translate(position - 0.45, floors[0] if up else floors[1] + value)
            )
            ax.add_patch(
                PathPatch(
                    transform.transform_path(path),
                    color=_NUCLEOTIDE_COLORS[base],
                    linewidth=0,
                )
            )
            floors[0 if up else 1] += value
        low, high = min(low, floors[1]), max(high, floors[0])
    pad = 0.05 * ((high - low) or 1.0)
    ax.set_xlim(-0.5, len(scores) - 0.5)
    ax.set_ylim(low - pad, high + pad)
    ax.axhline(0, color="black", linewidth=LINE_WIDTH, alpha=0.5)


@dataclass(frozen=True)
class SaliencyPlot:
    """Inputs for a focused signed saliency figure."""

    scores: np.ndarray
    sequence: str
    genomic_start: int
    chrom: str
    assembly: str
    label: str
    out_path: str
    window_bp: int = 96
    center_bp: int | None = None


def plot_saliency(plot: SaliencyPlot) -> SaliencyWindow:
    """Plot the strongest local saliency window with paper nucleotide colors.

    Args:
        plot: Saliency matrix, sequence, coordinates, and output path.

    Returns:
        The selected sequence-logo interval. PNG and PDF files are written.
    """
    if plot.scores.ndim != 2 or plot.scores.shape[1] != 4:
        raise ValueError(
            f"Expected saliency shape (sequence_length, 4), got {plot.scores.shape}."
        )
    if len(plot.sequence) != plot.scores.shape[0]:
        raise ValueError(
            f"Sequence length {len(plot.sequence)} does not match saliency "
            f"length {plot.scores.shape[0]}."
        )
    window = select_saliency_window(
        scores=plot.scores,
        genomic_start=plot.genomic_start,
        window_bp=plot.window_bp,
        center_bp=plot.center_bp,
    )
    fig, ax = plt.subplots(figsize=(12, 3))
    _draw_logo(ax=ax, scores=plot.scores[window.left : window.right])
    ax.set_title(plot.label, fontsize=11)
    ax.set_ylabel("Attribution", fontsize=9)
    tick_step = max(1, (window.right - window.left) // 6)
    ticks = np.arange(0, window.right - window.left, tick_step)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{window.start + int(tick):,}" for tick in ticks],
        fontsize=7,
    )
    ax.set_xlabel(f"{plot.chrom} position ({plot.assembly})", fontsize=9)
    _style(ax)
    fig.tight_layout()
    _save(fig=fig, out_path=plot.out_path)
    return window
