"""Paper-style log-fold-change plots for saved disruption tracks."""

from dataclasses import dataclass
from io import BytesIO

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ntv3_tools.colors import COLORS

_EPSILON = 1e-8
_LINE_WIDTH = 0.5

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


@dataclass(frozen=True)
class DisruptionPlotData:
    """Inputs for per-bin alternate-versus-reference log2 fold change."""

    reference: np.ndarray
    alternate: np.ndarray
    labels: list[str]
    variant_bin: int
    bin_size_bp: int
    window_bp: int | None


@dataclass(frozen=True)
class DisruptionPlotArtifacts:
    """PNG and PDF bytes for one multi-track log2 fold-change figure."""

    png: bytes
    pdf: bytes


def _save_figure(figure: plt.Figure) -> DisruptionPlotArtifacts:
    """Serialize one figure as publication-ready PNG and PDF bytes."""
    png = BytesIO()
    pdf = BytesIO()
    figure.savefig(png, format="png", dpi=600, bbox_inches="tight")
    figure.savefig(pdf, format="pdf", dpi=600, bbox_inches="tight")
    plt.close(figure)
    return DisruptionPlotArtifacts(png=png.getvalue(), pdf=pdf.getvalue())


def plot_log2fc_tracks(data: DisruptionPlotData) -> DisruptionPlotArtifacts:
    """Plot paper-style per-bin log2(Alt/Ref) for approved tracks."""
    reference = np.asarray(data.reference, dtype=float)
    alternate = np.asarray(data.alternate, dtype=float)
    if (
        reference.ndim != 2
        or alternate.shape != reference.shape
        or reference.shape[0] != len(data.labels)
        or not reference.size
    ):
        raise ValueError(
            "Reference and alternate tracks must share shape "
            "(n_tracks, n_bins) and match labels."
        )
    if (
        not np.isfinite(reference).all()
        or not np.isfinite(alternate).all()
        or np.any(reference < 0)
        or np.any(alternate < 0)
    ):
        raise ValueError("Log-fold-change tracks must be finite and nonnegative.")
    if data.bin_size_bp <= 0:
        raise ValueError("Variant bin and bin size are invalid.")

    if data.window_bp is None:
        start, end = 0, reference.shape[1]
    else:
        if data.variant_bin < 0 or data.variant_bin >= reference.shape[1]:
            raise ValueError("A cropped plot requires the variant inside the output.")
        half_bins = max(1, data.window_bp // data.bin_size_bp)
        start = max(0, data.variant_bin - half_bins)
        end = min(reference.shape[1], data.variant_bin + half_bins + 1)
    variant_bin = data.variant_bin - start
    relative_bp = (np.arange(end - start) - variant_bin) * data.bin_size_bp
    log2fc = np.log2(
        (alternate[:, start:end] + _EPSILON) / (reference[:, start:end] + _EPSILON)
    )

    sns.set_style("white")
    figure, axes = plt.subplots(
        len(data.labels),
        1,
        figsize=(8, max(1.8, 1.6 * len(data.labels))),
        sharex=True,
        squeeze=False,
    )
    for axis, values, label in zip(axes[:, 0], log2fc, data.labels, strict=True):
        positive = np.where(values >= 0, values, 0)
        negative = np.where(values < 0, values, 0)
        axis.fill_between(relative_bp, positive, color=COLORS["red_2"])
        axis.plot(relative_bp, positive, color=COLORS["red_1"], linewidth=_LINE_WIDTH)
        axis.fill_between(relative_bp, negative, color=COLORS["blue_3"])
        axis.plot(relative_bp, negative, color=COLORS["blue_2"], linewidth=_LINE_WIDTH)
        axis.axhline(0, color="black", linewidth=_LINE_WIDTH, alpha=0.5)
        axis.axvline(0, color=COLORS["red_1"], linewidth=0.8)
        axis.set_ylabel(label, fontsize=8)
        axis.tick_params(width=_LINE_WIDTH, labelsize=7)
        sns.despine(ax=axis)
    axes[-1, 0].set_xlabel("Position relative to SNV (bp)", fontsize=9)
    figure.supylabel("Log2(Alt/Ref)", fontsize=9)
    figure.tight_layout()
    return _save_figure(figure=figure)
