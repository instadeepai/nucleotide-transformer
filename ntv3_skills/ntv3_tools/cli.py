"""NTv3 skill CLI: Predict, Annotate, Compare, and Explain one approved run."""

import argparse
import importlib.util
import json
import logging
import os
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from ntv3_tools.annotate import (
    ANNOTATION_DECODERS,
    ANNOTATION_MODELS,
    CLADE_ALIASES,
    FORWARD_COMMAND,
    AnnotateRequest,
    annotation_resolved_block,
    plan_annotation,
    predict_run_window,
)
from ntv3_tools.annotation_report import ANNOTATION_PLAN
from ntv3_tools.approval import (
    RESOLVED_BLOCK,
    approve_annotation,
    approve_run,
    pending_annotation_plan,
    require_approval,
    write_resolved_block,
)
from ntv3_tools.critic import Claim, audit_claims, write_audit_log
from ntv3_tools.disrupt import Snv
from ntv3_tools.doctor import (
    INSTALL_COMMAND,
    DoctorReport,
    device_for,
    doctor_payload,
    gpu_line,
    jax_python_candidates,
    jax_python_line,
    run_doctor,
)
from ntv3_tools.errors import Ntv3Error, RequestError
from ntv3_tools.gt import search_encode_experiment
from ntv3_tools.literature import package_literature
from ntv3_tools.locus import resolve_gene_locus
from ntv3_tools.lookups import LookupRequest, run_lookup
from ntv3_tools.plan import (
    EXCLUDED_ANALYSES,
    AnalysisPlan,
    SaliencyLogoConfig,
    read_plan,
    refuse_excluded_analyses,
    require_window,
    run_dir_name,
    validate_plan,
    write_plan,
)
from ntv3_tools.registry import DEFAULT_CHECKPOINT, resolve_model
from ntv3_tools.report import REPORT_STAGE_LABELS, write_viewer_bundle
from ntv3_tools.report_app import REPORT_APP_HTML
from ntv3_tools.resolve import resolve_tracks
from ntv3_tools.resolved import resolved_block

logger = logging.getLogger("ntv3_tools")
_WORKTREE_ROOT = Path(__file__).resolve().parents[2]
_RESULTS_ROOT = (_WORKTREE_ROOT / "results").resolve()


def _require_run_directory(path: str) -> str:
    """Return an absolute run directory confined to worktree-root results/.

    Args:
        path: CLI output or run-directory value.

    Returns:
        Absolute validated run-directory path.

    Raises:
        RequestError: The path is outside results/ or names results/ itself.
    """
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = _WORKTREE_ROOT / candidate
    resolved = candidate.resolve()
    if resolved == _RESULTS_ROOT or _RESULTS_ROOT not in resolved.parents:
        raise RequestError(
            f"Which run directory? {resolved} is not a named directory under "
            f"{_RESULTS_ROOT}. Pass a results/<run> path. No directory is assumed."
        )
    return str(resolved)


def _require_checkpoint(model: str | None) -> str:
    """Return the predict model, ``NTv3_tracks``, refusing any other name.

    Args:
        model: CLI --model value.

    Returns:
        The resolved model name.

    Raises:
        RequestError: A name other than ``NTv3_tracks``.
    """
    try:
        return resolve_model(name=model or DEFAULT_CHECKPOINT)
    except KeyError as exc:
        raise RequestError(str(exc)) from exc


def _cmd_predict(args: argparse.Namespace) -> None:
    """Run the isolated JAX-backed command with the current Python."""
    plan_path = Path(args.plan).resolve()
    _require_run_directory(path=str(plan_path.parent))
    _require_approved_plan(plan_path=plan_path)
    if importlib.util.find_spec("jax") is None:
        raise RequestError(
            f"Predict requires the JAX Python from `ntv3 doctor`. Install it "
            f"from the repo root with `{INSTALL_COMMAND}`."
        )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ntv3_tools.science_cli",
            "--plan",
            str(plan_path),
        ],
        check=False,
    )
    if completed.returncode:
        raise SystemExit(completed.returncode)


def _variant(values: list[str] | None) -> Snv | None:
    """Return the ``--variant CHROM POS REF ALT`` SNV, or None."""
    if values is None:
        return None
    chrom, position, ref, alt = values
    if not position.isdigit():
        raise RequestError(
            f"Variant position must be a 1-based integer, not {position}."
        )
    return Snv(chrom=chrom, position=int(position), ref=ref.upper(), alt=alt.upper())


def _cmd_disrupt(args: argparse.Namespace) -> None:
    """Dispatch single-SNV Compare to the current JAX Python."""
    plan_path = Path(args.plan).resolve()
    _require_run_directory(path=str(plan_path.parent))
    _require_approved_plan(plan_path=plan_path)
    if importlib.util.find_spec("jax") is None:
        raise RequestError(
            f"Compare requires the JAX Python from `ntv3 doctor`. Install it "
            f"from the repo root with `{INSTALL_COMMAND}`."
        )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ntv3_tools.disrupt_science",
            "--plan",
            str(plan_path),
            "--chrom",
            args.chrom,
            "--position",
            str(args.position),
            "--ref",
            args.ref,
            "--alt",
            args.alt,
            "--null-n",
            str(args.null_n),
            "--seed",
            str(args.seed),
        ],
        check=False,
    )
    if completed.returncode:
        raise SystemExit(completed.returncode)


def _run_doctor() -> DoctorReport:
    """Find this repo's JAX Python and its backend."""
    return run_doctor(
        candidates=jax_python_candidates(environ=os.environ, repo_root=_WORKTREE_ROOT)
    )


def _cmd_annotate(args: argparse.Namespace) -> None:
    """Plan an Annotate ask, or print its Resolved block; never loads a model."""
    if (args.fasta is None) == (args.locus_run is None):
        raise RequestError(
            "What should be annotated? Pass --fasta with --out, or --locus-run with "
            "a Predict run. Not both. No input is assumed."
        )
    if args.fasta is not None and args.out is None:
        raise RequestError(
            "Which run directory? Pass --out results/<run> for a FASTA. No "
            "directory is assumed."
        )
    run_arg = args.out if args.fasta is not None else args.locus_run
    out_dir = _require_run_directory(path=run_arg)
    request = AnnotateRequest(
        clade=args.clade,
        decoder=args.decoder,
        out_dir=out_dir,
        fasta=args.fasta,
    )
    window = None if args.fasta else predict_run_window(run_dir=Path(out_dir))
    doctor = _run_doctor()
    lines = annotation_resolved_block(
        request=request, run_dir=run_arg, doctor=doctor, window=window
    )
    sys.stdout.write("\n".join(lines) + "\n")
    if args.resolve_only:
        return
    plan_annotation(request=request, block=lines, window=window)
    logger.info(
        f"Planned {run_arg}. Approve with `ntv3 approve --run-dir {run_arg}`, then "
        f"run `{FORWARD_COMMAND.format(plan=Path(run_arg) / ANNOTATION_PLAN)}`."
    )


def _cmd_annotate_forward(args: argparse.Namespace) -> None:
    """Run an approved annotation in this Python, which must be the JAX one."""
    plan_path = Path(args.plan).resolve()
    _require_run_directory(path=str(plan_path.parent))
    if importlib.util.find_spec("jax") is None:
        raise RequestError(
            f"The annotation forward requires the JAX Python from `ntv3 doctor`. "
            f"Install it from the repo root with `{INSTALL_COMMAND}`."
        )
    from ntv3_tools.annotate_science import run_annotation_forward

    run_annotation_forward(plan_path=plan_path)


# The Artifact tool publishes pages of at most 16 MB.
ARTIFACT_PAGE_LIMIT_BYTES = 16_000_000
SPECIES_QUESTION = "Which species? Pass --species. No species is assumed."
ASSEMBLY_QUESTION = "Which assembly? Pass --assembly. No assembly is assumed."
TISSUE_QUESTION = (
    "Which tissue or cell type should the tracks be predicted in? Pass --tissue. "
    "No tissue is assumed."
)
TRACKS_QUESTION = (
    "Which tracks should be predicted? Pass --tracks with assay names or exact "
    "model heads (ENCSR*/CNhs*). No track is chosen for you."
)
WINDOW_QUESTION = (
    "Which window length in bp? On CPU, 4096 is quick and 16384 or more runs with "
    "a warning (slow; Saliency can exhaust memory), up to 1048576. Pass "
    "--window-bp. No length is assumed."
)
GPU_WINDOW_QUESTION = (
    "Which window length in bp? This GPU runs any multiple of 128 up to 1048576, "
    "for example 131072 or 1048576. Pass --window-bp. No length is assumed."
)
LANES_QUESTION = (
    "Include gene and ClinVar browser lanes? Pass --ucsc-tracks with names, or "
    "pass none. No lane is added for you."
)


def _browser_lanes(values: list[str]) -> list[str]:
    """Return the named browser lanes; ``none`` alone (or nothing) means none."""
    if not values or [value.lower() for value in values] == ["none"]:
        return []
    if any(value.lower() == "none" for value in values):
        raise RequestError(
            "Which browser lanes? --ucsc-tracks has none and lane names together. "
            "Pass names, or pass none. No lane is added for you."
        )
    return list(values)


def _cmd_plan(args: argparse.Namespace) -> None:
    """Write plan.yaml and the Resolved block, or print the block and write nothing.

    The run directory is ``results/<assembly>_<chrom>_<start>_<end>_<analyses>``.
    Nothing here loads the model; ``ntv3 approve`` must follow.
    """
    refuse_excluded_analyses(analyses=args.analyses or [])
    # One question at a time, in this order, before any network call.
    if not (args.species or "").strip():
        raise RequestError(SPECIES_QUESTION)
    if not (args.assembly or "").strip():
        raise RequestError(ASSEMBLY_QUESTION)
    if not (args.tissue or "").strip():
        raise RequestError(TISSUE_QUESTION)
    if not args.tracks:
        raise RequestError(TRACKS_QUESTION)
    device = str(args.device).lower()
    doctor = None
    device_rule = ""
    if device == "auto":
        doctor = _run_doctor()
        device = device_for(report=doctor)
        device_rule = (
            "; cuda because this JAX python sees a GPU"
            if device == "cuda"
            else "; cpu because this JAX python sees no GPU"
        )
    on_cpu = device in {"cpu", "mps"}
    if args.window_bp is None:
        raise RequestError(WINDOW_QUESTION if on_cpu else GPU_WINDOW_QUESTION)
    require_window(window_bp=args.window_bp)
    if args.ucsc_tracks is None:
        raise RequestError(LANES_QUESTION)
    variant = _variant(values=args.variant)
    coordinates = (args.chrom, args.start, args.end)
    locus_metadata = None
    window_rule = f"your span equals the {device} {args.window_bp:,} bp window"
    if variant is not None and all(value is None for value in coordinates):
        chrom = variant.chrom
        start = max(0, variant.position - 1 - args.window_bp // 2)
        end = start + args.window_bp
        window_rule = f"centered on the variant {variant.chrom}:{variant.position:,}"
    elif all(value is None for value in coordinates):
        locus_metadata = resolve_gene_locus(
            gene=args.gene,
            assembly=args.assembly,
            window_bp=args.window_bp,
        )
        chrom = locus_metadata.chrom
        start = locus_metadata.start
        end = locus_metadata.end
        window_rule = (
            f"centered on the {args.gene} TSS ({locus_metadata.transcript}, "
            f"{locus_metadata.strand} strand) by {locus_metadata.source}"
        )
    elif any(value is None for value in coordinates):
        raise RequestError(
            "Which locus? Pass --chrom, --start, and --end together, or none of "
            "them to look up the gene. No coordinate is assumed."
        )
    else:
        chrom = args.chrom
        start = args.start
        end = args.end
    window_rule += device_rule
    ucsc_tracks = _browser_lanes(values=args.ucsc_tracks)
    plan = AnalysisPlan(
        species=args.species,
        assembly=args.assembly,
        gene=args.gene,
        chrom=chrom,
        start=start,
        end=end,
        tissue=args.tissue,
        tracks=args.tracks,
        device=device,
        window_bp=args.window_bp,
        want_report=args.want_report,
        want_notebook=args.want_notebook,
        # The critic is on until an approval turns it off.
        want_critic=True,
        ucsc_tracks=ucsc_tracks,
        analyses=args.analyses or [],
        model=args.model or DEFAULT_CHECKPOINT,
        attribution_tracks=args.attribution_tracks or [],
        saliency_logo=SaliencyLogoConfig(
            window_bp=args.saliency_logo_window_bp,
            center_bp=args.saliency_logo_center_bp,
        ),
        variant=variant,
    )
    validate_plan(plan=plan)
    run_dir = _RESULTS_ROOT / run_dir_name(plan=plan)
    run_label = f"results/{run_dir.name}"
    run_dir_question = (
        f"Which run directory? {run_label} already exists and is never "
        "overwritten. Change the locus or the analyses. No directory is reused."
    )
    if run_dir.exists():
        raise RequestError(run_dir_question)
    tracks = resolve_tracks(
        species=plan.species,
        tissue=plan.tissue,
        names=plan.tracks,
        encode_search=lambda tissue, assay: search_encode_experiment(
            tissue=tissue, assay=assay, species=plan.species
        ),
    )
    headless = [row for row in tracks.tracks if not row.model_head]
    if headless:
        raise RequestError(
            "Which model head should be predicted for "
            + "; ".join(
                f"'{row.query}'"
                + (
                    f" (ENCODE experiment {row.gt_accession} matched; an ENCODE "
                    "accession is not used as a model head)"
                    if row.gt_accession
                    else ""
                )
                for row in headless
            )
            + "? Pass an exact model head (ENCSR*/CNhs*) in --tracks. No head is "
            "picked for you."
        )
    doctor = doctor or _run_doctor()
    lines = resolved_block(
        plan=plan,
        tracks=tracks,
        doctor=doctor,
        window_rule=window_rule,
        run_dir=run_label,
        gene_from_lookup=locus_metadata is not None,
    )
    if not args.resolve_only:
        # Exclusive create: of two plans for one directory, only one writes.
        _RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
        try:
            run_dir.mkdir()
        except FileExistsError as exc:
            raise RequestError(run_dir_question) from exc
        write_plan(plan=plan, out_dir=str(run_dir))
        if locus_metadata is not None:
            Path(run_dir, "locus.json").write_text(
                json.dumps(asdict(locus_metadata), indent=2) + "\n"
            )
        write_resolved_block(run_dir=run_dir, lines=lines)
    sys.stdout.write("\n".join(lines) + "\n")
    if not args.resolve_only:
        logger.info(
            f"Wrote {run_label}/plan.yaml and {RESOLVED_BLOCK}. Nothing runs until "
            f"`ntv3 approve --run-dir {run_label}`."
        )


def _cmd_approve(args: argparse.Namespace) -> None:
    """Approve a planned run; the only way to set ``approved: true``.

    A pending annotation plan is approved first, so a locus annotation added to
    an approved Predict run is approved by the same command.
    """
    run_dir = _require_run_directory(path=args.run_dir)
    if pending_annotation_plan(run_dir=run_dir):
        approve_annotation(run_dir=run_dir, critic=not args.no_critic)
        logger.info(
            f"Approved the annotation in {run_dir}. Run `"
            f"{FORWARD_COMMAND.format(plan=Path(args.run_dir) / ANNOTATION_PLAN)}"
            "`, then report."
        )
        return
    plan = approve_run(run_dir=run_dir, critic=not args.no_critic)
    logger.info(
        f"Approved {run_dir} (critic {'on' if plan.want_critic else 'off'}). "
        "Run lookup, predict, and report."
    )


def _require_approved_plan(plan_path: Path) -> None:
    """Refuse an unapproved plan before the JAX Python is even started."""
    require_approval(run_dir=plan_path.parent, plan=read_plan(path=plan_path))


def _cmd_doctor(args: argparse.Namespace) -> None:
    """Report the first Python that imports JAX and this repo: paths and versions."""
    report = _run_doctor()
    sys.stdout.write(json.dumps(doctor_payload(report=report), indent=2) + "\n")
    logger.info(jax_python_line(report=report))
    gpu = gpu_line(report=report)
    if gpu is not None:
        logger.info(gpu)


def _cmd_lookup(args: argparse.Namespace) -> None:
    """Resolve dual ids, fetch GT BigWig URLs, optional UCSC overlays."""
    plan = read_plan(path=args.plan)
    run_dir = _require_run_directory(path=str(Path(args.plan).resolve().parent))
    resolved = run_lookup(
        request=LookupRequest(
            species=plan.species,
            tissue=plan.tissue,
            names=plan.tracks,
            out_dir=run_dir,
            chrom=plan.chrom,
            start=plan.start,
            end=plan.end,
            ucsc_tracks=plan.ucsc_tracks,
            gene=plan.gene,
        )
    )
    logger.info(
        f"Lookup wrote track_resolve.json "
        f"({len(resolved.tracks)} tracks) under {run_dir}"
    )


def _cmd_critic(args: argparse.Namespace) -> None:
    """KEEP / SUPPRESS claims against on-disk evidence."""
    run_dir = _require_run_directory(path=args.run_dir)
    claims_path = Path(args.claims) if args.claims else Path(run_dir) / "claims.json"
    if not claims_path.exists():
        raise RequestError(
            f"No claims file at {claims_path}. Write claims.json or pass --claims."
        )
    raw = json.loads(claims_path.read_text())
    claims = [
        Claim(
            claim_id=row["claim_id"],
            text=row["text"],
            evidence_paths=row.get("evidence_paths") or [],
            citekey=row.get("citekey"),
            kind=row.get("kind") or "other",
        )
        for row in raw
    ]
    bib = Path(run_dir) / "literature" / "citations.bib"
    decisions = audit_claims(
        claims=claims,
        run_dir=run_dir,
        bib_path=str(bib) if bib.exists() else None,
    )
    path = write_audit_log(decisions=decisions, out_dir=run_dir)
    logger.info(f"Wrote {path}")


def _cmd_literature(args: argparse.Namespace) -> None:
    """Validate and package agent-collected literature records."""
    run_dir = _require_run_directory(path=args.run_dir)
    bib_path, manifest_path = package_literature(
        run_dir=run_dir,
        records_path=args.records,
    )
    logger.info(f"Wrote {bib_path} and {manifest_path}")


def _cmd_report(args: argparse.Namespace) -> None:
    """Write report.json, the browser page to publish, and reproducibility files."""
    run_dir = _require_run_directory(path=args.run_dir)
    path = write_viewer_bundle(
        run_dir=run_dir,
        refresh=args.refresh,
        requested_stages=args.stages,
        contig=args.contig,
    )
    page_bytes = (path.parent / REPORT_APP_HTML).stat().st_size
    if page_bytes > ARTIFACT_PAGE_LIMIT_BYTES:
        logger.warning(
            f"{REPORT_APP_HTML} is {page_bytes / 1e6:.1f} MB, over the 16 MB "
            "artifact limit. Do not publish it; a shorter window or fewer tracks "
            "or analyses makes a smaller page."
        )
        return
    logger.info(
        f"Wrote {path} and {path.parent / REPORT_APP_HTML}. Publish "
        f"{REPORT_APP_HTML} with the Artifact tool; do not show its path."
    )


def main() -> None:
    """Entry point for the ``ntv3`` CLI."""
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s %(name)s: %(message)s"
    )
    parser = argparse.ArgumentParser(prog="ntv3", description="NTv3 skill CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    p_plan = sub.add_parser("plan", help="Write plan.yaml and refuse to guess")
    p_plan.add_argument(
        "--species", help="human or mouse. Required; asked for, never assumed."
    )
    p_plan.add_argument(
        "--assembly", help="hg38 or mm10. Required; asked for, never assumed."
    )
    p_plan.add_argument("--gene", required=True)
    p_plan.add_argument("--chrom")
    p_plan.add_argument("--start", type=int)
    p_plan.add_argument("--end", type=int)
    p_plan.add_argument(
        "--tissue", help="Tissue or cell type. Required; asked for, never assumed."
    )
    p_plan.add_argument(
        "--model",
        default=DEFAULT_CHECKPOINT,
        help=f"Predict model; only {DEFAULT_CHECKPOINT} loads",
    )
    p_plan.add_argument(
        "--tracks",
        nargs="+",
        help="Assay names or exact model heads. Required; asked for, never chosen.",
    )
    p_plan.add_argument(
        "--attribution-tracks",
        nargs="+",
        default=[],
        help=(
            "Selected prediction tracks that get saliency. Defaults to every "
            "selected track."
        ),
    )
    p_plan.add_argument(
        "--device",
        default="auto",
        type=str.lower,
        choices=["auto", "cpu", "cuda", "gpu"],
        help="default auto: cuda when the JAX python sees a GPU, else cpu",
    )
    p_plan.add_argument(
        "--window-bp",
        dest="window_bp",
        type=int,
        help="Window length in bp. Required; asked for, never assumed.",
    )
    p_plan.add_argument(
        "--saliency-logo-window-bp",
        type=int,
        default=96,
        help="Bases shown in each static saliency logo (default 96)",
    )
    p_plan.add_argument(
        "--saliency-logo-center-bp",
        type=int,
        help="Absolute genomic center; default is each track's strongest peak",
    )
    p_plan.add_argument("--want-report", dest="want_report", action="store_true")
    p_plan.add_argument("--want-notebook", dest="want_notebook", action="store_true")
    p_plan.add_argument(
        "--ucsc-tracks",
        dest="ucsc_tracks",
        nargs="*",
        default=None,
        help=(
            "UCSC browser lanes, e.g. knownGene clinvarMain (hg38) or "
            "ncbiRefSeq (mm10), or none. Required; asked for, never added."
        ),
    )
    p_plan.add_argument(
        "--analyses",
        nargs="*",
        default=[],
        help=(
            "Approved analyses: predict, saliency, compare. Saliency and compare "
            "require predict; compare requires --variant. The plan refuses any "
            "other analysis."
        ),
    )
    p_plan.add_argument(
        "--variant",
        nargs=4,
        metavar=("CHROM", "POS", "REF", "ALT"),
        help=(
            "The one SNV Compare scores, 1-based. Without chrom, start, and end "
            "the window is centered on it."
        ),
    )
    p_plan.add_argument(
        "--resolve-only",
        action="store_true",
        help="Print the Resolved block and write nothing",
    )
    p_plan.set_defaults(func=_cmd_plan)

    p_approve = sub.add_parser(
        "approve",
        help="Approve a planned run after its Resolved block; nothing loads before",
    )
    p_approve.add_argument("--run-dir", dest="run_dir", required=True)
    p_approve.add_argument(
        "--no-critic",
        action="store_true",
        help="Turn the critic off for this run (it is on by default)",
    )
    p_approve.set_defaults(func=_cmd_approve)

    p_doctor = sub.add_parser(
        "doctor",
        help=(
            "Check $NTV3_JAX_PYTHON and this repo's .venv for jax and "
            "nucleotide_transformer_v3. Prints paths and versions only."
        ),
    )
    p_doctor.set_defaults(func=_cmd_doctor)

    p_lookup = sub.add_parser(
        "lookup", help="Dual-id resolve + live ENCODE GT + optional UCSC overlays"
    )
    p_lookup.add_argument("--plan", required=True)
    p_lookup.set_defaults(func=_cmd_lookup)

    p_predict = sub.add_parser("predict", help="Predict tracks at a locus")
    p_predict.add_argument("--plan", required=True)
    p_predict.set_defaults(func=_cmd_predict)

    p_disrupt = sub.add_parser(
        "compare",
        aliases=["disrupt"],
        help=(
            "Optional null ranking only; disrupt is the same command under its "
            "old name, not a second Compare. The allele comparison runs in predict."
        ),
        description=(
            "Optional null ranking only. `ntv3 disrupt` is an alias of this "
            "command, kept for older scripts and replay; it is not a second "
            "Compare. The allele comparison itself runs in "
            "predict: approve `--analyses predict compare --variant CHROM POS REF "
            "ALT` and predict scores the alternate allele on its own load. Use this "
            "command only to rank one SNV against --null-n seeded random SNVs in "
            "the same window, on a plan without compare; each null SNV is one "
            "more forward."
        ),
    )
    p_disrupt.add_argument("--plan", required=True)
    p_disrupt.add_argument("--chrom", required=True)
    p_disrupt.add_argument("--position", required=True, type=int)
    p_disrupt.add_argument("--ref", required=True)
    p_disrupt.add_argument("--alt", required=True)
    p_disrupt.add_argument("--null-n", dest="null_n", type=int, default=0)
    p_disrupt.add_argument("--seed", type=int, default=0)
    p_disrupt.set_defaults(func=_cmd_disrupt)

    p_annotate = sub.add_parser(
        "annotate",
        help=(
            "Plan an annotation of a FASTA, or of the 1 Mb around a Predict run. "
            "mammals, plants, and insects are registered only at 1,048,576 bp "
            "windows. Writes annotation_plan.json and its Resolved block; the "
            "forward runs only after approve, with annotate-forward."
        ),
    )
    p_annotate.add_argument(
        "--clade",
        required=True,
        choices=sorted([*ANNOTATION_MODELS, *CLADE_ALIASES]),
        help="Model config name: mammals, plants, or insects (insecta is accepted)",
    )
    p_annotate.add_argument("--fasta", help="Plain or gzipped FASTA")
    p_annotate.add_argument(
        "--locus-run",
        dest="locus_run",
        help="Predict run whose window to annotate, 1 Mb around it",
    )
    p_annotate.add_argument(
        "--decoder",
        choices=ANNOTATION_DECODERS,
        default="hmm_constrained",
        help="State-path decoder (default hmm_constrained, as released)",
    )
    p_annotate.add_argument("--out", help="New results/<run> for a FASTA")
    p_annotate.add_argument(
        "--resolve-only",
        action="store_true",
        help="Print the Resolved block for approval and write nothing",
    )
    p_annotate.set_defaults(func=_cmd_annotate)

    p_forward = sub.add_parser(
        "annotate-forward",
        help="Run an approved annotation in the JAX Python and write its GFF3",
    )
    p_forward.add_argument("--plan", required=True)
    p_forward.set_defaults(func=_cmd_annotate_forward)

    p_crit = sub.add_parser("critic", help="KEEP / SUPPRESS claims vs evidence")
    p_crit.add_argument("--run-dir", dest="run_dir", required=True)
    p_crit.add_argument(
        "--claims", default=None, help="claims.json (default run-dir/claims.json)"
    )
    p_crit.set_defaults(func=_cmd_critic)

    p_literature = sub.add_parser(
        "literature",
        help="Validate records and package citations without network access",
    )
    p_literature.add_argument("--run-dir", dest="run_dir", required=True)
    p_literature.add_argument(
        "--records",
        required=True,
        help="JSON records file under the approved run directory",
    )
    p_literature.set_defaults(func=_cmd_literature)

    p_rep = sub.add_parser(
        "report",
        help="Write report.json, the browser page to publish, replay, and RO-Crate",
    )
    p_rep.add_argument("--run-dir", dest="run_dir", required=True)
    p_rep.add_argument(
        "--contig",
        help="Annotate run: the contig the page shows (default: the most genes)",
    )
    p_rep.add_argument(
        "--refresh",
        action="store_true",
        help="Refresh report metadata after a new approved analysis",
    )
    p_rep.add_argument(
        "--stages",
        nargs="*",
        choices=list(REPORT_STAGE_LABELS),
        default=None,
        help=(
            "Optional stages to show as pending until their artifacts exist. "
            "Omit to show only completed sections."
        ),
    )
    p_rep.set_defaults(func=_cmd_report)

    # Refused by name before parsing, so they are never offered as commands.
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        if command in EXCLUDED_ANALYSES:
            refuse_excluded_analyses(analyses=[command])
        args = parser.parse_args()
        args.func(args)
    except (Ntv3Error, KeyError) as exc:
        logger.error(str(exc))
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
