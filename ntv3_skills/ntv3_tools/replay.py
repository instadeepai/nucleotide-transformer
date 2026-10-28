"""Replay one approved NTv3 plan into a new local results directory."""

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from ntv3_tools.approval import APPROVAL, RESOLVED_BLOCK
from ntv3_tools.errors import RequestError


@dataclass(frozen=True)
class VariantReplay:
    """Explicit variant inputs for an optional disruption replay."""

    chrom: str
    position: int
    ref: str
    alt: str
    null_n: int
    seed: int


@dataclass(frozen=True)
class ReplayRequest:
    """Inputs for replaying an approved analysis."""

    source_run: Path
    output_run: str
    research_python: str
    claims: str | None
    literature_records: str | None
    variant: VariantReplay | None


def _reject_unsafe_input_path(path: str, label: str) -> Path:
    """Validate one explicit local input file.

    The returned path is absolute but not resolved: a venv's ``bin/python`` is a
    symlink to the base interpreter, and only the link runs inside the venv.
    """
    if "://" in path:
        raise RequestError(f"{label} must be a local file, not a cloud path.")
    candidate = Path(path).expanduser().absolute()
    resolved = candidate.resolve()
    for parts in (candidate.parts, resolved.parts):
        if ".env" in parts or "secrets" in parts:
            raise RequestError(f"{label} cannot come from secrets/ or .env.")
    if not resolved.is_file():
        raise RequestError(f"{label} does not exist: {candidate}.")
    return candidate


def _resolve_destination(source_run: Path, output_run: str) -> tuple[Path, Path]:
    """Resolve a new destination confined to the worktree results directory."""
    results_root = source_run.resolve().parent
    worktree_root = results_root.parent
    if results_root.name != "results":
        raise RequestError("The replay script must remain under results/<run>/code/.")
    candidate = Path(output_run)
    if not candidate.is_absolute():
        candidate = worktree_root / candidate
    destination = candidate.resolve()
    if destination == results_root or results_root not in destination.parents:
        raise RequestError(f"--output-run must name a directory under {results_root}.")
    if destination.exists():
        raise RequestError(f"Replay destination already exists: {destination}.")
    return destination, worktree_root


def _subprocess_environment(worktree_root: Path) -> dict[str, str]:
    """Build a local subprocess environment without credential variables."""
    blocked_fragments = (
        "AWS",
        "CREDENTIAL",
        "GCP",
        "GCS",
        "GOOGLE",
        "KEY",
        "MLFLOW",
        "PASSWORD",
        "S3",
        "SECRET",
        "TOKEN",
    )
    environment = {
        name: value
        for name, value in os.environ.items()
        if not any(fragment in name.upper() for fragment in blocked_fragments)
    }
    python_paths = [str(worktree_root / "ntv3_skills"), str(worktree_root)]
    existing = environment.get("PYTHONPATH")
    if existing:
        python_paths.append(existing)
    environment["PYTHONPATH"] = os.pathsep.join(python_paths)
    return environment


def _run_command(command: list[str], environment: dict[str, str]) -> None:
    """Run one replay command without a shell."""
    completed = subprocess.run(
        command,
        check=False,
        env=environment,
    )
    if completed.returncode:
        raise RequestError(
            f"Replay command failed with exit code {completed.returncode}: "
            f"{' '.join(command[:4])}."
        )


def _variant_from_args(args: argparse.Namespace) -> VariantReplay | None:
    """Validate the all-or-none optional variant argument group."""
    names = (
        "variant_chrom",
        "variant_position",
        "variant_ref",
        "variant_alt",
        "null_n",
        "seed",
    )
    values = [getattr(args, name) for name in names]
    if not any(value is not None for value in values):
        return None
    if any(value is None for value in values):
        raise RequestError(
            "Disruption replay requires variant chrom, position, ref, alt, "
            "null count, and seed."
        )
    return VariantReplay(
        chrom=args.variant_chrom,
        position=args.variant_position,
        ref=args.variant_ref,
        alt=args.variant_alt,
        null_n=args.null_n,
        seed=args.seed,
    )


def replay_analysis(request: ReplayRequest) -> Path:
    """Replay an approved plan into a new confined results directory.

    Args:
        request: Explicit replay inputs.

    Returns:
        New run directory.
    """
    source_run = request.source_run.resolve()
    source_plan = source_run / "plan.yaml"
    if not source_plan.is_file():
        raise RequestError(f"Approved plan does not exist: {source_plan}.")
    destination, worktree_root = _resolve_destination(
        source_run=source_run,
        output_run=request.output_run,
    )
    research_python = _reject_unsafe_input_path(
        path=request.research_python,
        label="--research-python",
    )
    claims = (
        _reject_unsafe_input_path(path=request.claims, label="--claims")
        if request.claims
        else None
    )
    records = (
        _reject_unsafe_input_path(
            path=request.literature_records,
            label="--literature-records",
        )
        if request.literature_records
        else None
    )

    destination.mkdir(parents=False)
    destination_plan = destination / "plan.yaml"
    shutil.copyfile(source_plan, destination_plan)
    if destination_plan.read_bytes() != source_plan.read_bytes():
        raise RequestError("Copied plan differs from the approved source plan.")
    # The source run's approval travels with its plan; predict checks both.
    for name in (RESOLVED_BLOCK, APPROVAL):
        if (source_run / name).is_file():
            shutil.copyfile(source_run / name, destination / name)

    environment = _subprocess_environment(worktree_root=worktree_root)
    cli_prefix = [sys.executable, "-m", "ntv3_tools.cli"]
    _run_command(
        command=[*cli_prefix, "lookup", "--plan", str(destination_plan)],
        environment=environment,
    )
    _run_command(
        command=[
            str(research_python),
            "-m",
            "ntv3_tools.cli",
            "predict",
            "--plan",
            str(destination_plan),
        ],
        environment=environment,
    )
    if request.variant is not None:
        variant = request.variant
        _run_command(
            command=[
                str(research_python),
                "-m",
                "ntv3_tools.cli",
                "disrupt",
                "--plan",
                str(destination_plan),
                "--chrom",
                variant.chrom,
                "--position",
                str(variant.position),
                "--ref",
                variant.ref,
                "--alt",
                variant.alt,
                "--null-n",
                str(variant.null_n),
                "--seed",
                str(variant.seed),
            ],
            environment=environment,
        )
    if records is not None:
        records_target = destination / "literature" / "records.json"
        records_target.parent.mkdir()
        shutil.copyfile(records, records_target)
        _run_command(
            command=[
                *cli_prefix,
                "literature",
                "--run-dir",
                str(destination),
                "--records",
                str(records_target),
            ],
            environment=environment,
        )
    if claims is not None:
        claims_target = destination / "claims.json"
        shutil.copyfile(claims, claims_target)
        _run_command(
            command=[
                *cli_prefix,
                "critic",
                "--run-dir",
                str(destination),
                "--claims",
                str(claims_target),
            ],
            environment=environment,
        )
    _run_command(
        command=[*cli_prefix, "report", "--run-dir", str(destination)],
        environment=environment,
    )
    return destination


def main(source_run: Path) -> None:
    """Parse generated-script arguments and run the replay."""
    parser = argparse.ArgumentParser(
        description="Replay an approved NTv3 plan into a new results directory."
    )
    parser.add_argument("--research-python", required=True)
    parser.add_argument("--output-run", required=True)
    parser.add_argument("--claims")
    parser.add_argument("--literature-records")
    parser.add_argument("--variant-chrom")
    parser.add_argument("--variant-position", type=int)
    parser.add_argument("--variant-ref")
    parser.add_argument("--variant-alt")
    parser.add_argument("--null-n", type=int)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    try:
        variant = _variant_from_args(args=args)
        replay_analysis(
            request=ReplayRequest(
                source_run=source_run,
                output_run=args.output_run,
                research_python=args.research_python,
                claims=args.claims,
                literature_records=args.literature_records,
                variant=variant,
            )
        )
    except RequestError as exc:
        parser.error(str(exc))
