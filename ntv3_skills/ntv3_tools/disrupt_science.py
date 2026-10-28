"""JAX single-SNV Compare: the approved tracks for the reference and alternate allele.

Two token forwards of ``NTv3_tracks`` over the approved window, one per allele,
on only the approved head columns; the reference forward is reused from
predict when it ran first. Seeded null SNVs in the same window, when asked
for, give each effect a percentile. Compare on the default path runs none.
"""

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ntv3_tools.approval import require_approval
from ntv3_tools.disrupt import (
    DisruptArtifacts,
    NullSamplingRequest,
    Snv,
    apply_snv,
    empirical_absolute_percentile,
    generate_null_snvs,
    require_disrupt_artifacts_absent,
    signed_absmax,
    validate_snv,
    write_disrupt_artifacts,
)
from ntv3_tools.disrupt_plotting import DisruptionPlotData, plot_log2fc_tracks
from ntv3_tools.errors import Ntv3Error, RequestError
from ntv3_tools.model import (
    LoadedModel,
    bigwig_track_names,
    center_fraction,
    load_model,
    prediction_region,
    run_forward,
)
from ntv3_tools.plan import AnalysisPlan, read_plan
from ntv3_tools.registry import head_indices, model_head_for, resolve_model
from ntv3_tools.sequence import fetch_ucsc_sequence

logger = logging.getLogger("ntv3_tools")


@dataclass(frozen=True)
class DisruptRequest:
    """Inputs for one approved single-SNV Compare run.

    Attributes:
        plan_path: Path to the approved plan.
        variant: Explicit user-approved SNV.
        null_n: Number of empirical null SNVs; 0 runs none.
        seed: Random seed for the null SNVs.
    """

    plan_path: str
    variant: Snv
    null_n: int = 0
    seed: int = 0


@dataclass(frozen=True)
class CompareInputs:
    """One approved plan, its loaded model, and its window sequence.

    Attributes:
        plan: Approved analysis plan.
        run_dir: Run directory holding track_resolve.json.
        loaded: Loaded ``NTv3_tracks``.
        sequence: Window sequence, soft mask kept, shape (window_bp,).
        reference: Reference predictions from predict, shape (n_tracks, n_bins),
            or None to run the reference forward here.
    """

    plan: AnalysisPlan
    run_dir: Path
    loaded: LoadedModel
    sequence: str
    reference: np.ndarray | None = None


def _approved_track_heads(
    plan: AnalysisPlan,
    run_dir: Path,
    loaded: LoadedModel,
) -> tuple[list[str], list[str], list[int]]:
    """Load the approved model heads and their output indices.

    Args:
        plan: Approved analysis plan.
        run_dir: Run directory containing track_resolve.json.
        loaded: Loaded checkpoint.

    Returns:
        Model heads, display labels, and exact output indices.
    """
    resolve_path = run_dir / "track_resolve.json"
    if not resolve_path.exists():
        raise RequestError("Run lookup before Compare.")
    rows = json.loads(resolve_path.read_text()).get("tracks") or []
    if len(rows) != len(plan.tracks):
        raise RequestError(
            f"Plan requests {len(plan.tracks)} tracks but lookup resolved {len(rows)}."
        )
    available = bigwig_track_names(loaded=loaded, species=plan.species)
    heads = [model_head_for(row=row, available=available) for row in rows]
    indices = head_indices(names=available, heads=heads)
    labels = [
        str(row.get("query") or head) for row, head in zip(rows, heads, strict=True)
    ]
    return heads, labels, indices


def _tracks(
    loaded: LoadedModel, sequence: str, species: str, indices: list[int]
) -> np.ndarray:
    """Return the selected tracks for one sequence, shape (n_tracks, n_bins)."""
    out = run_forward(
        loaded=loaded, sequence=sequence, species=species, track_indices=indices
    )
    return np.asarray(out["bigwig_tracks_logits"][0], dtype=np.float32).T


def run_compare(
    inputs: CompareInputs, variant: Snv, null_n: int = 0, seed: int = 0
) -> Path:
    """Score one SNV: alt-minus-ref of each approved track, and its null rank.

    Args:
        inputs: Approved plan, run directory, loaded model, and window sequence.
        variant: One-based SNV inside the window.
        null_n: Number of seeded null SNVs; 0 runs none and gives no percentile.
        seed: Seed for the null SNVs.

    Returns:
        Path to the created ``disrupt/`` directory.
    """
    plan, run_dir, loaded, sequence = (
        inputs.plan,
        inputs.run_dir,
        inputs.loaded,
        inputs.sequence,
    )
    require_disrupt_artifacts_absent(run_dir=str(run_dir))
    offset = validate_snv(
        variant=variant,
        window_chrom=plan.chrom,
        window_start=plan.start,
        sequence=sequence,
    )
    heads, labels, indices = _approved_track_heads(
        plan=plan, run_dir=run_dir, loaded=loaded
    )
    ref_tracks = inputs.reference
    if ref_tracks is None:
        ref_tracks = _tracks(loaded, sequence, plan.species, indices)
    alt_tracks = _tracks(
        loaded,
        apply_snv(sequence=sequence, offset=offset, alt=variant.alt),
        plan.species,
        indices,
    )
    if alt_tracks.shape != ref_tracks.shape:
        raise RuntimeError(
            f"Reference tracks {ref_tracks.shape} and alternate tracks "
            f"{alt_tracks.shape} differ."
        )
    delta_tracks = alt_tracks - ref_tracks
    track_effects = signed_absmax(delta_tracks, axis=1)
    effect = float(signed_absmax(track_effects, axis=0))

    percentile: float | None = None
    track_percentiles: list[float | None] = [None] * len(heads)
    if null_n > 0:
        null_track_effects = np.stack(
            [
                signed_absmax(
                    _tracks(
                        loaded,
                        apply_snv(
                            sequence=sequence,
                            offset=null.position - 1 - plan.start,
                            alt=null.alt,
                        ),
                        plan.species,
                        indices,
                    )
                    - ref_tracks,
                    axis=1,
                )
                for null in generate_null_snvs(
                    request=NullSamplingRequest(
                        sequence=sequence,
                        chrom=plan.chrom,
                        window_start=plan.start,
                        query=variant,
                        null_n=null_n,
                        seed=seed,
                    )
                )
            ]
        )
        percentile = empirical_absolute_percentile(
            query_effect=effect, null_effects=signed_absmax(null_track_effects, axis=1)
        )
        track_percentiles = [
            empirical_absolute_percentile(
                query_effect=float(track_effects[index]),
                null_effects=null_track_effects[:, index],
            )
            for index in range(len(heads))
        ]

    prediction_start, prediction_end = prediction_region(
        start=plan.start, end=plan.end, center_frac=center_fraction(loaded=loaded)
    )
    n_bins = ref_tracks.shape[1]
    if n_bins <= 0 or (prediction_end - prediction_start) % n_bins != 0:
        raise RuntimeError(
            f"Track shape {ref_tracks.shape} does not divide the prediction span "
            f"{prediction_end - prediction_start}."
        )
    bin_size = (prediction_end - prediction_start) // n_bins
    method = {
        "forward": "one NTv3_tracks token forward per allele, approved heads only",
        "track_reduction": "signed abs-max of alt-minus-ref across output bins",
        "overall_reduction": "signed abs-max across approved tracks",
        "null_sampling": (
            "uniform over eligible in-window A/C/G/T SNVs with replacement; "
            "query SNV excluded"
            if null_n
            else "none"
        ),
        "percentile": (
            "100 * (1 + count(abs(null) <= abs(query))) / (null_n + 1)"
            if null_n
            else "none"
        ),
    }
    per_track = [
        {
            "track": label,
            "model_head": head,
            "signed_effect": float(track_effect),
            "percentile": track_percentile,
        }
        for label, head, track_effect, track_percentile in zip(
            labels, heads, track_effects, track_percentiles, strict=True
        )
    ]
    plot = plot_log2fc_tracks(
        data=DisruptionPlotData(
            reference=ref_tracks,
            alternate=alt_tracks,
            labels=labels,
            variant_bin=(variant.position - 1 - prediction_start) // bin_size,
            bin_size_bp=bin_size,
            window_bp=None,
        )
    )
    output_dir = write_disrupt_artifacts(
        run_dir=str(run_dir),
        artifacts=DisruptArtifacts(
            variant={
                **asdict(variant),
                "position_0based": variant.position - 1,
                "window": {"chrom": plan.chrom, "start": plan.start, "end": plan.end},
                "checkpoint": loaded.checkpoint,
                "track_heads": heads,
                "prediction": {
                    "chrom": plan.chrom,
                    "start": prediction_start,
                    "end": prediction_end,
                    "bin_size": bin_size,
                },
                "in_prediction": prediction_start
                <= variant.position - 1
                < prediction_end,
                "null_n": null_n,
                "seed": seed,
                "method": method,
            },
            ref_tracks=ref_tracks,
            alt_tracks=alt_tracks,
            delta_tracks=delta_tracks,
            top_rows=[
                {"rank": rank, **row, "absolute_effect": abs(row["signed_effect"])}
                for rank, row in enumerate(
                    sorted(per_track, key=lambda row: -abs(row["signed_effect"])),
                    start=1,
                )
            ],
            effect_percentile={
                "signed_effect": effect,
                "absolute_effect": abs(effect),
                "percentile": percentile,
                "null_n": null_n,
                "seed": seed,
                "per_track": per_track,
                "method": method,
            },
            log2fc_png=plot.png,
            log2fc_pdf=plot.pdf,
        ),
    )
    logger.info(
        f"Compared {variant.chrom}:{variant.position} {variant.ref}>{variant.alt} "
        f"on {len(heads)} tracks into {output_dir}"
    )
    return output_dir


def run_disrupt(request: DisruptRequest) -> Path:
    """Load the model and score one SNV against an approved plan.

    Args:
        request: Approved plan, query variant, and null configuration.

    Returns:
        Path to the created ``disrupt/`` directory.
    """
    plan = read_plan(path=request.plan_path)
    run_dir = Path(request.plan_path).resolve().parent
    require_approval(run_dir=run_dir, plan=plan)
    require_disrupt_artifacts_absent(run_dir=str(run_dir))
    try:
        model_name = resolve_model(name=plan.model)
    except KeyError as exc:
        raise RequestError(str(exc)) from exc
    sequence = fetch_ucsc_sequence(
        assembly=plan.assembly,
        chrom=plan.chrom,
        start=plan.start,
        end=plan.end,
        keep_soft_mask=True,
    )
    if len(sequence) != plan.window_bp:
        raise RequestError(
            f"UCSC returned {len(sequence)} bp for the approved {plan.window_bp} bp window."
        )
    validate_snv(
        variant=request.variant,
        window_chrom=plan.chrom,
        window_start=plan.start,
        sequence=sequence,
    )
    return run_compare(
        inputs=CompareInputs(
            plan=plan,
            run_dir=run_dir,
            loaded=load_model(model_name=model_name, device=plan.device),
            sequence=sequence,
        ),
        variant=request.variant,
        null_n=request.null_n,
        seed=request.seed,
    )


def main() -> None:
    """Run the JAX-backed Compare command."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    parser = argparse.ArgumentParser(prog="ntv3-disrupt-science")
    parser.add_argument("--plan", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--alt", required=True)
    parser.add_argument("--null-n", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    try:
        run_disrupt(
            request=DisruptRequest(
                plan_path=args.plan,
                variant=Snv(
                    chrom=args.chrom,
                    position=args.position,
                    ref=args.ref.upper(),
                    alt=args.alt.upper(),
                ),
                null_n=args.null_n,
                seed=args.seed,
            )
        )
    except (Ntv3Error, KeyError, RuntimeError, ValueError) as exc:
        logger.error(str(exc))
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
