"""Find the Python that runs Predict and Annotate. Paths and versions only.

Candidates, in order: ``$NTV3_JAX_PYTHON``, then this repo's ``.venv``, which
``pip install -e .`` makes. A candidate passes when ``import jax`` and
``import nucleotide_transformer_v3`` both work; predict loads its model through
that package. ``ntv3_skills/.venv`` runs the planner CLI and is never a
candidate. No checkpoint is loaded and no ``.env`` is read.
"""

import os
import subprocess
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

JAX_PYTHON_ENV = "NTV3_JAX_PYTHON"
PACKAGE = "nucleotide_transformer_v3"
INSTALL_COMMAND = "python3.11 -m venv .venv && .venv/bin/pip install -e ."
# PyPI's jaxlib is CPU-only, so a GPU host needs the extra for Annotate's forward.
GPU_INSTALL_COMMAND = ".venv/bin/pip install -e '.[gpu]'"
REPO_VENV_RULE = ".venv in this repo (pip install -e .)"
_PROBE = (
    "import sys\n"
    "print(sys.version.split()[0], flush=True)\n"
    "import jax\n"
    "print(jax.__version__, flush=True)\n"
    f"import {PACKAGE}\n"
    f"print({PACKAGE!r}, flush=True)\n"
    "print(jax.default_backend(), flush=True)\n"
)
# The probe inherits only these, so no secret reaches the child or its output. The
# last three decide which devices JAX sees, as they do for predict.
_PROBE_ENV_KEYS = (
    "PATH",
    "HOME",
    "TMPDIR",
    "CUDA_VISIBLE_DEVICES",
    "JAX_PLATFORMS",
    "LD_LIBRARY_PATH",
)

RunCommand = Callable[..., subprocess.CompletedProcess]


@dataclass(frozen=True)
class JaxCandidate:
    """One interpreter that was checked.

    Attributes:
        python: Interpreter path.
        rule: Why it is a candidate.
        passed: Whether ``import jax`` and the package import both worked.
        python_version: Interpreter version, when it ran.
        jax_version: JAX version, when ``import jax`` worked.
        imports_package: Whether ``import nucleotide_transformer_v3`` worked.
        detail: What passed, "not found", or the import failure.
        jax_backend: JAX's default backend ("gpu" or "cpu"), when it imported.
    """

    python: str
    rule: str
    passed: bool
    python_version: str | None
    jax_version: str | None
    imports_package: bool
    detail: str
    jax_backend: str | None = None


@dataclass(frozen=True)
class DoctorReport:
    """The selected Python and every candidate checked before it.

    Attributes:
        jax_python: First passing interpreter, or None when none passed.
        candidates: Checked candidates in order.
    """

    jax_python: str | None
    candidates: list[JaxCandidate]


def jax_python_candidates(
    environ: Mapping[str, str],
    repo_root: Path,
) -> list[tuple[str, str]]:
    """List interpreters to try, in order, with the rule for each.

    Args:
        environ: Process environment. Only ``NTV3_JAX_PYTHON`` is read.
        repo_root: Root of this repo.

    Returns:
        (interpreter path, rule) pairs.
    """
    candidates: list[tuple[str, str]] = []
    configured = environ.get(JAX_PYTHON_ENV)
    if configured:
        path = Path(configured).expanduser()
        if path.is_dir():
            path = path / "bin" / "python"
        candidates.append((str(path), f"${JAX_PYTHON_ENV}"))
    candidates.append((str(repo_root / ".venv" / "bin" / "python"), REPO_VENV_RULE))
    return candidates


def _probe_env() -> dict[str, str]:
    """Return the only environment a probe inherits."""
    return {key: os.environ[key] for key in _PROBE_ENV_KEYS if key in os.environ}


def check_jax_python(
    python: str,
    rule: str,
    run: RunCommand = subprocess.run,
) -> JaxCandidate:
    """Check whether one interpreter imports JAX and this repo's package.

    Args:
        python: Interpreter path.
        rule: Why it is a candidate.
        run: Subprocess runner (inject in tests).

    Returns:
        The checked candidate.
    """
    if not Path(python).is_file():
        return JaxCandidate(python, rule, False, None, None, False, "not found")
    try:
        completed = run(
            [python, "-c", _PROBE],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
            env=_probe_env(),
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return JaxCandidate(python, rule, False, None, None, False, type(exc).__name__)
    printed = [line.strip() for line in (completed.stdout or "").splitlines()]
    python_version, jax_version, package, backend = [*printed, *[None] * 4][:4]
    imports_package = package == PACKAGE
    if completed.returncode or not imports_package:
        lines = [line for line in (completed.stderr or "").splitlines() if line.strip()]
        detail = lines[-1].strip()[:200] if lines else f"exit {completed.returncode}"
        return JaxCandidate(
            python, rule, False, python_version, jax_version, False, detail
        )
    return JaxCandidate(
        python,
        rule,
        True,
        python_version,
        jax_version,
        True,
        f"imports jax and {PACKAGE}",
        backend,
    )


def run_doctor(
    candidates: list[tuple[str, str]],
    run: RunCommand = subprocess.run,
) -> DoctorReport:
    """Check candidates in order and stop at the first that passes.

    Args:
        candidates: (interpreter path, rule) pairs from ``jax_python_candidates``.
        run: Subprocess runner (inject in tests).

    Returns:
        The selected interpreter and the candidates checked.
    """
    checked: list[JaxCandidate] = []
    for python, rule in candidates:
        candidate = check_jax_python(python=python, rule=rule, run=run)
        checked.append(candidate)
        if candidate.passed:
            return DoctorReport(jax_python=python, candidates=checked)
    return DoctorReport(jax_python=None, candidates=checked)


def doctor_payload(report: DoctorReport) -> dict[str, Any]:
    """Return the JSON-ready report: paths, rules, and versions.

    Args:
        report: Doctor result.

    Returns:
        JSON-ready mapping. ``install`` is present only when nothing passed;
        ``gpu`` only when the selected Python is not on a GPU.
    """
    payload: dict[str, Any] = {
        "jax_python": report.jax_python,
        "candidates": [asdict(candidate) for candidate in report.candidates],
        "planner": "ntv3_skills/.venv runs plan, lookup, annotate, report, and "
        "critic. It is not a JAX environment.",
    }
    if report.jax_python is None:
        payload["install"] = INSTALL_COMMAND
    gpu = gpu_line(report=report)
    if gpu is not None:
        payload["gpu"] = gpu
    return payload


def device_for(report: DoctorReport) -> str:
    """Return the Predict device: ``cuda`` when the JAX Python sees a GPU, else cpu."""
    if report.jax_python is not None and report.candidates[-1].jax_backend == "gpu":
        return "cuda"
    return "cpu"


def gpu_line(report: DoctorReport) -> str | None:
    """Return what this JAX Python cannot run without a GPU, or None."""
    if report.jax_python is None:
        return None
    backend = report.candidates[-1].jax_backend
    if backend in (None, "gpu"):
        return None
    return (
        f"GPU: this JAX python runs on {backend}, so Annotate's forward refuses and "
        "Predict windows of 16,384 bp or more are slow. On a GPU host, run "
        f"`{GPU_INSTALL_COMMAND}`."
    )


def jax_python_line(report: DoctorReport) -> str:
    """Return the Resolved-block line for the JAX Python.

    Args:
        report: Doctor result.

    Returns:
        The selected interpreter with its rule, or the candidates that failed
        and the install command.
    """
    if report.jax_python is not None:
        chosen = report.candidates[-1]
        return (
            f"JAX python: {chosen.python} (Python {chosen.python_version}, "
            f"jax {chosen.jax_version}, imports {PACKAGE}) — first candidate "
            f"that passes, via {chosen.rule}"
        )
    failed = "; ".join(
        f"{candidate.python} ({candidate.detail})" for candidate in report.candidates
    )
    missing = (
        "this repo's .venv is missing"
        if all(candidate.detail == "not found" for candidate in report.candidates)
        else f"no candidate imports jax and {PACKAGE}"
    )
    return (
        f"JAX python: none — {failed}. Predict cannot start: {missing}. Install "
        f"this repo from its root with `{INSTALL_COMMAND}` (CPU wheels; not the "
        "gpu extra)."
    )
