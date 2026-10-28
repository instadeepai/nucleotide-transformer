"""JAX-backed NTv3 prediction entry point.

Loads ``NTv3_tracks`` through this repo's ``load_tracks_model`` once, predicts
the approved tracks, then runs the approved saliency and Compare on the same
model and window. ``read_plan`` refuses any analysis outside the skill, so
nothing here is fetched or loaded for one.
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np

from ntv3_tools import plotting
from ntv3_tools.approval import require_approval
from ntv3_tools.disrupt import require_disrupt_artifacts_absent, validate_snv
from ntv3_tools.disrupt_science import CompareInputs, run_compare
from ntv3_tools.errors import Ntv3Error, RequestError
from ntv3_tools.gt_bigwig import (
    GroundTruthComparison,
    bigwig_reader_problem,
    write_ground_truth_comparison,
)
from ntv3_tools.locus import LocusRequest
from ntv3_tools.model import bigwig_track_names, load_model
from ntv3_tools.plan import cpu_window_warning, read_plan
from ntv3_tools.predict import TrackPredictionRequest, predict_tracks
from ntv3_tools.provenance import build_manifest, sha256_text
from ntv3_tools.registry import model_head_for, resolve_model
from ntv3_tools.report import ReportSections, save_report_bundle
from ntv3_tools.saliency import SaliencyRequest, saliency_per_track
from ntv3_tools.sequence import crop_to_multiple, fetch_ucsc_sequence

logger = logging.getLogger("ntv3_tools")


def _compare_ground_truth(
    run_dir: Path,
    track_rows: list[dict],
    comparison: GroundTruthComparison,
    predicted: np.ndarray,
) -> tuple[bool, str | None]:
    """Compare with ground truth where it exists, and say what was not compared.

    A track with no ENCODE experiment is left out; the other tracks are still
    compared. Any other ground-truth problem becomes the note, so the saved
    forward is never discarded.

    Args:
        run_dir: Approved run directory.
        track_rows: ``track_resolve.json`` rows, one per predicted track.
        comparison: Prediction span to compare.
        predicted: Model output values, shape (n_tracks, n_bins).

    Returns:
        Whether ``gt/observed.npy`` and the agreement table were written, and
        a note for ``run_manifest.json`` and the report, or None.
    """
    without_accession = [
        str(row.get("query") or "track")
        for row in track_rows
        if not row.get("gt_accession")
    ]
    skipped = (
        f"no ENCODE experiment for {', '.join(without_accession)}"
        if without_accession
        else None
    )
    if len(without_accession) == len(track_rows):
        return False, skipped
    if not (run_dir / "gt" / "manifest.json").is_file():
        return False, "lookup saved no gt/manifest.json"
    problem = bigwig_reader_problem()
    if problem is not None:
        return False, problem
    try:
        write_ground_truth_comparison(comparison=comparison, predicted=predicted)
    except (Ntv3Error, ValueError, OSError) as exc:
        # Includes an existing gt/observed.npy, non-finite values, a corrupt
        # manifest, and an unreadable BigWig.
        return False, str(exc).rstrip(".")
    return True, skipped


def run_prediction(plan_path: str) -> None:
    """Run track prediction from one approved plan.

    Args:
        plan_path: Path to the approved plan.yaml.
    """
    plan = read_plan(path=plan_path)
    run_dir = Path(plan_path).resolve().parent
    # Nothing is fetched or loaded for a plan the user has not approved.
    require_approval(run_dir=run_dir, plan=plan)
    resolved_path = run_dir / "track_resolve.json"
    if not resolved_path.exists():
        raise RequestError("Run lookup before prediction.")
    resolved = json.loads(resolved_path.read_text())
    track_rows = resolved.get("tracks") or []
    if len(track_rows) != len(plan.tracks):
        raise RequestError(
            f"Plan requests {len(plan.tracks)} tracks but lookup resolved "
            f"{len(track_rows)}."
        )
    try:
        model_name = resolve_model(name=plan.model)
    except KeyError as exc:
        raise RequestError(str(exc)) from exc

    locus = LocusRequest(
        chrom=plan.chrom,
        start=plan.start,
        end=plan.end,
        species=plan.species,
    )
    locus.validate()
    sequence = crop_to_multiple(
        fetch_ucsc_sequence(
            assembly=locus.assembly,
            chrom=locus.chrom,
            start=locus.start,
            end=locus.end,
            keep_soft_mask=True,
        )
    )
    if plan.variant is not None:
        require_disrupt_artifacts_absent(run_dir=str(run_dir))
        validate_snv(
            variant=plan.variant,
            window_chrom=plan.chrom,
            window_start=plan.start,
            sequence=sequence,
        )
    loaded = load_model(model_name=model_name, device=plan.device)
    if str(loaded.device.platform) == "cpu" and (
        warning := cpu_window_warning(plan.window_bp, plan.analyses)
    ):
        logger.warning(warning)
    available_heads = set(bigwig_track_names(loaded=loaded, species=plan.species))
    track_ids = [
        model_head_for(row=row, available=available_heads) for row in track_rows
    ]
    track_labels = [
        str(row.get("query") or model_head)
        for row, model_head in zip(track_rows, track_ids, strict=True)
    ]
    prediction = predict_tracks(
        request=TrackPredictionRequest(
            model_name=model_name,
            locus=locus,
            track_ids=track_ids,
            device=plan.device,
            loaded=loaded,
            sequence=sequence,
        )
    )
    prediction_span = prediction.prediction_end - prediction.prediction_start
    if (
        prediction.tracks.ndim != 2
        or prediction.tracks.shape[1] == 0
        or prediction_span % prediction.tracks.shape[1] != 0
    ):
        raise RuntimeError(
            f"Prediction shape {prediction.tracks.shape} does not divide genomic "
            f"span {prediction_span}."
        )
    bin_size = prediction_span // prediction.tracks.shape[1]

    # The forward is saved first, so no ground-truth problem can discard it.
    tracks_dir = run_dir / "tracks"
    tracks_dir.mkdir(parents=True, exist_ok=True)
    np.save(tracks_dir / "tracks.npy", prediction.tracks)
    (tracks_dir / "accessions.json").write_text(
        json.dumps(prediction.accessions, indent=2) + "\n"
    )
    (tracks_dir / "metadata.json").write_text(
        json.dumps(
            {
                "chrom": locus.chrom,
                "start": prediction.prediction_start,
                "end": prediction.prediction_end,
                "bin_size": bin_size,
            },
            indent=2,
        )
        + "\n"
    )
    ground_truth_compared, ground_truth_note = _compare_ground_truth(
        run_dir=run_dir,
        track_rows=track_rows,
        comparison=GroundTruthComparison(
            run_dir=run_dir,
            chrom=locus.chrom,
            start=prediction.prediction_start,
            end=prediction.prediction_end,
            bin_size=bin_size,
        ),
        predicted=prediction.tracks,
    )
    if not ground_truth_compared:
        logger.warning(f"Ground truth not compared: {ground_truth_note}.")
    elif ground_truth_note is not None:
        logger.warning(f"Ground truth compared, except: {ground_truth_note}.")
    figures = [
        plotting.plot_tracks(
            plot=plotting.TrackPlot(
                tracks=prediction.tracks,
                labels=track_labels,
                start=prediction.prediction_start,
                end=prediction.prediction_end,
                chrom=locus.chrom,
                assembly=locus.assembly,
                out_path=str(run_dir / "figures" / "tracks.png"),
            )
        )
    ]
    if "saliency" in plan.analyses:
        chosen = set(plan.attribution_tracks or track_labels)
        rows = [
            (head, label)
            for head, label in zip(track_ids, track_labels, strict=True)
            if label in chosen
        ]
        saliency_per_track(
            request=SaliencyRequest(
                model_name=model_name,
                locus=locus,
                track_ids=[head for head, _ in rows],
                track_labels=[label for _, label in rows],
                out_dir=str(run_dir),
                logo=plan.saliency_logo,
                loaded=loaded,
                sequence=sequence,
            )
        )
        figures += [
            str(run_dir / "figures" / f"saliency_{head}.png") for head, _ in rows
        ]
    if plan.variant is not None:
        run_compare(
            inputs=CompareInputs(
                plan=plan,
                run_dir=run_dir,
                loaded=loaded,
                sequence=sequence,
                reference=prediction.tracks,
            ),
            variant=plan.variant,
        )
        figures.append(str(run_dir / "figures" / "disrupt_log2fc.png"))
    manifest = build_manifest(
        model_name=model_name,
        genome_build=locus.assembly,
        objective="functional_track_prediction",
        extra={
            "window": f"{locus.chrom}:{locus.start}-{locus.end}",
            "predicted_region": (
                f"{locus.chrom}:{prediction.prediction_start}-"
                f"{prediction.prediction_end}"
            ),
            "tracks": prediction.accessions,
            "analyses": plan.analyses,
            "device": str(loaded.device.platform),
            "device_requested": plan.device,
            "ground_truth_compared": ground_truth_compared,
            "ground_truth_note": ground_truth_note,
            # The sequence the model saw, soft mask kept, so a changed input shows.
            "input_sha256": sha256_text(sequence),
        },
    )
    sections = ReportSections(
        title=f"NTv3 track prediction — {locus.chrom}:{locus.start}-{locus.end}",
        conclusion=(
            f"Predicted functional tracks {', '.join(prediction.accessions)} over "
            f"{locus.chrom}:{prediction.prediction_start}-"
            f"{prediction.prediction_end} ({locus.assembly})"
            + (
                f". Ground truth was not compared: {ground_truth_note}."
                if not ground_truth_compared
                else " and compared them with aligned public ENCODE measurements"
                + (f", except: {ground_truth_note}." if ground_truth_note else ".")
            )
        ),
        limitations=[
            "Model logits and experimental signal differ. Correlations compare shape."
        ],
    )
    save_report_bundle(
        out_dir=str(run_dir),
        manifest=manifest,
        sections=sections,
        figure_paths=figures,
    )


def main() -> None:
    """Run the JAX-backed command."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    parser = argparse.ArgumentParser(prog="ntv3-science")
    parser.add_argument("--plan", required=True)
    args = parser.parse_args()
    try:
        run_prediction(plan_path=args.plan)
    except (Ntv3Error, KeyError, RuntimeError, ValueError) as exc:
        logger.error(str(exc))
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
