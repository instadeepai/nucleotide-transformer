"""The Resolved block: each value a run will use and the rule that chose it."""

import os
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from ntv3_tools.browser import CATALOG_LABELS
from ntv3_tools.doctor import DoctorReport, jax_python_line
from ntv3_tools.errors import ArtifactError
from ntv3_tools.plan import AnalysisPlan, cpu_window_warning, gpu_saliency_warning
from ntv3_tools.registry import DEFAULT_CHECKPOINT
from ntv3_tools.resolve import TrackResolve
from ntv3_tools.sequence import CENTER_FRACTION, prediction_region

TRACKS_REPO = f"InstaDeepAI/{DEFAULT_CHECKPOINT}"
# ClinVar is human only. mm10 knownGene rows carry no gene symbol, so mouse
# uses ncbiRefSeq, whose rows carry it in name2.
_GENE_AND_CLINVAR_LANES: dict[str, tuple[str, ...]] = {
    "hg38": ("knownGene", "clinvarMain"),
    "mm10": ("ncbiRefSeq",),
}


@dataclass(frozen=True)
class CatalogChoice:
    """One UCSC browser lane the lane question names.

    Attributes:
        track: UCSC ``track=`` id.
        label: Label shown in the browser.
        rule: Why it was chosen.
    """

    track: str
    label: str
    rule: str


def catalog_label(track: str, assembly: str) -> str | None:
    """Return the browser lane label for a UCSC track, or None if unknown.

    Args:
        track: UCSC ``track=`` id.
        assembly: Genome assembly.

    Returns:
        The label, or None.
    """
    return CATALOG_LABELS.get(assembly, {}).get(track)


def _catalog_choice(track: str, assembly: str, rule: str) -> CatalogChoice:
    """Return one catalog choice, refusing an id the browser has no lane for."""
    label = catalog_label(track=track, assembly=assembly)
    if label is None:
        raise ArtifactError(
            f"Which browser lane? '{track}' has no lane label for {assembly}. Pass "
            "--ucsc-tracks with a lane the browser draws. No lane is assumed."
        )
    return CatalogChoice(track=track, label=label, rule=rule)


def gene_and_clinvar_lanes(assembly: str) -> list[CatalogChoice]:
    """Return the gene lane, and ClinVar on human, that the lane question names.

    None of them is added to a run: a lane is in the plan only when
    ``--ucsc-tracks`` names it, and ``--ucsc-tracks none`` gives no lane.

    Args:
        assembly: Genome assembly.

    Returns:
        Catalog choices in display order.

    Raises:
        ArtifactError: A chosen id has no browser lane label.
    """
    return [
        _catalog_choice(track=track, assembly=assembly, rule="named in the request")
        for track in _GENE_AND_CLINVAR_LANES.get(assembly, ())
    ]


def _span(chrom: str, start: int, end: int) -> str:
    """Format a zero-based half-open span with its length."""
    return f"{chrom}:{start:,}-{end:,} ({end - start:,} bp)"


def hub_weights_cached(repo: str, environ: Mapping[str, str] = os.environ) -> bool:
    """Return whether the Hugging Face cache already holds a repo's JAX weights.

    Reads only the cache directory listing; nothing is downloaded.

    Args:
        repo: Hub repo id, such as ``InstaDeepAI/NTv3_tracks``.
        environ: Process environment (``HF_HUB_CACHE`` / ``HF_HOME``).

    Returns:
        True when a snapshot of ``repo`` has a ``.msgpack``.
    """
    hub = environ.get("HF_HUB_CACHE") or str(
        Path(environ.get("HF_HOME") or Path.home() / ".cache" / "huggingface") / "hub"
    )
    snapshots = Path(hub) / f"models--{repo.replace('/', '--')}" / "snapshots"
    return snapshots.is_dir() and any(snapshots.rglob("*.msgpack"))


def tracks_weights_cached(environ: Mapping[str, str] = os.environ) -> bool:
    """Return whether the Hugging Face cache already holds NTv3_tracks weights."""
    return hub_weights_cached(repo=TRACKS_REPO, environ=environ)


def resolved_block(
    plan: AnalysisPlan,
    tracks: TrackResolve,
    doctor: DoctorReport,
    window_rule: str,
    run_dir: str = "results/<run>",
    weights_cached: bool | None = None,
    gene_from_lookup: bool = False,
) -> list[str]:
    """Return the Resolved block shown once before approval.

    Args:
        plan: The plan the approval would write.
        tracks: Resolved track ids.
        doctor: JAX Python discovery result.
        window_rule: Why the window is this span.
        run_dir: Run directory the approval would create.
        weights_cached: Whether NTv3_tracks is cached; None checks the cache.
        gene_from_lookup: Whether a UCSC gene lookup placed the window.

    Returns:
        Lines from the "Resolved" title to the approval question.
    """
    if weights_cached is None:
        weights_cached = tracks_weights_cached()
    crop_start, crop_end = prediction_region(start=plan.start, end=plan.end)
    device_warning = (
        cpu_window_warning(plan.window_bp, plan.analyses)
        if plan.device in {"cpu", "mps"}
        else gpu_saliency_warning(plan.window_bp, plan.analyses)
    )
    lines = [
        "Resolved",
        f"Assembly: {plan.assembly} ({plan.species}) — as you named it",
        (
            f"Locus: {plan.chrom} start {plan.start} end {plan.end} "
            "— 0-based, half-open; the window below"
        ),
        (
            f"Gene: {plan.gene} — found by UCSC gene lookup"
            if gene_from_lookup
            else f"Gene: none looked up — your coordinates (label '{plan.gene}')"
        ),
        f"Tissue: {plan.tissue} — as you named it",
        (
            f"Tracks: {', '.join(row.query for row in tracks.tracks)} "
            "— as you named them; each head below"
        ),
        (
            f"Window: {_span(plan.chrom, plan.start, plan.end)} on {plan.device} "
            f"— {window_rule}"
        ),
        (
            f"Crop: {_span(plan.chrom, crop_start, crop_end)} "
            f"— central {CENTER_FRACTION:.1%} of the window, the span the model "
            "predicts"
        ),
        *([device_warning] if device_warning else []),
        (
            f"Model: {plan.model} ({TRACKS_REPO}) — loads only after approve, "
            "through nucleotide_transformer_v3.pretrained.load_tracks_model; only "
            "the tracks below are computed"
        ),
        (
            f"Weights: {TRACKS_REPO} is in the Hugging Face cache — no download"
            if weights_cached
            else f"Weights: {TRACKS_REPO} is not in the Hugging Face cache — "
            "predict's first load downloads it (711M parameters)"
        ),
    ]
    for row in tracks.tracks:
        lines.append(
            f"{row.query} in {tracks.tissue}: head {row.model_head or 'none'}, "
            f"ground truth {row.gt_accession or 'none'} — {row.source}"
        )
    if plan.ucsc_tracks:
        labeled = "; ".join(
            f"{track} ({catalog_label(track, plan.assembly) or 'no browser label'})"
            for track in plan.ucsc_tracks
        )
        lines.append(
            f"Browser lanes: {labeled} — named in the request; lookup fetches "
            "these for the window; strike one to drop it"
        )
    else:
        lines.append("Browser lanes: none — you passed none; no gene or ClinVar lane")
    can_start = doctor.jax_python is not None
    lines.append(f"Analyses: {', '.join(plan.analyses)} — as requested")
    if "saliency" in plan.analyses:
        lines.append(
            f"Saliency: {', '.join(plan.attribution_tracks or plan.tracks)} — "
            "gradient of each track's mean over the crop at each base, centred "
            "across A, C, G and T and kept at the reference base; one forward "
            "and backward per track"
        )
    if plan.variant is not None:
        variant = plan.variant
        inside = crop_start < variant.position <= crop_end
        lines.append(
            f"Compare: {variant.chrom}:{variant.position:,} {variant.ref}>"
            f"{variant.alt} (1-based) — "
            + ("inside" if inside else "in the window but outside")
            + " the crop; one more forward for the alternate allele, then "
            "alt minus ref for each track above"
        )
    lines += [
        (
            "Report: ntv3_report.html, published as the NTv3 report artifact "
            "— the lanes above, the predicted tracks, stored ground truth when "
            "predict compared it, saliency and alt-minus-ref lanes when run, and "
            "an NTv3 annotation lane when the run has a GFF3"
        ),
        (
            "Critic: on — default; the report runs it and records the result. "
            "Approve with --no-critic to turn it off"
            if plan.want_critic
            else "Critic: off — turned off at approval; the report says so"
        ),
        f"Run directory: {run_dir} — named from assembly, chrom, start, end, analyses",
        jax_python_line(report=doctor),
    ]
    approve = f"`ntv3 approve --run-dir {run_dir}`"
    if not can_start:
        lines.append(
            "Predict cannot start from this install. Do not approve a forward."
        )
    elif weights_cached:
        lines.append(
            "Approve to run lookup, predict, and the report, or strike a lane or "
            f"change a line. Approve with {approve}; nothing loads before it."
        )
    else:
        lines.append(
            f"Approve to download {TRACKS_REPO} and run lookup, predict, and the "
            f"report, or strike a lane or change a line. Approve with {approve}; "
            "nothing loads before it."
        )
    return lines
