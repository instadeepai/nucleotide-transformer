"""Reproducibility helpers: input hashing and run manifests.

Standalone (no repo dependency). Every skill run writes a manifest so results are reproducible from
the recorded model revision, window, seed, and input hash.
"""

import hashlib
import platform
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

_WORKTREE_ROOT = Path(__file__).resolve().parents[2]


def sha256_text(text: str) -> str:
    """Return the SHA-256 hex digest of a string.

    Args:
        text: Arbitrary input text (e.g. a serialized locus or an input sequence).

    Returns:
        Hex-encoded SHA-256 digest.
    """
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def build_manifest(
    model_name: str,
    genome_build: str,
    objective: str,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a reproducibility manifest for a run.

    Args:
        model_name: Internal JAX checkpoint name.
        genome_build: Genome assembly (e.g. "hg38").
        objective: The scalar objective / analysis type for the run.
        extra: Additional run-specific fields to merge in (windows, device, seed, input hash).

    Returns:
        Manifest dictionary ready to serialize to run_manifest.json.
    """
    manifest: dict[str, Any] = {
        "model": model_name,
        "genome_build": genome_build,
        "objective": objective,
        "platform": platform.platform(),
        "created_at": datetime.now(tz=UTC).isoformat(),
        **repository_state(),
    }
    if extra:
        manifest.update(extra)
    return manifest


def repository_state() -> dict[str, str | bool | None]:
    """Return the read-only Git provenance for the current worktree.

    Outside a git checkout (a pip install, a container) the commit is
    ``unknown`` and the dirty flag is None, so a finished forward is still saved.
    """
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=_WORKTREE_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=_WORKTREE_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return {"commit": "unknown", "worktree_dirty": None}
    return {
        "commit": commit,
        "worktree_dirty": bool(status.strip()),
    }
