"""The approval gate: no model loads for a plan the user has not approved.

``ntv3 plan`` writes ``plan.yaml`` (``approved: false``) and the Resolved
block, ``resolved.txt``. Only ``ntv3 approve`` sets ``approved: true``, and it
records digests of the plan and the block in ``approval.json``. Predict and
Compare call ``require_approval`` before any sequence fetch or model load, so
a plan edited after approval, or a hand-set ``approved: true``, is refused.

An Annotate plan, ``annotation_plan.json``, has its own block and record
(``annotation_resolved.txt``, ``annotation_approval.json``), so a locus run can
hold an approved Predict plan and an annotation plan side by side.
"""

import hashlib
import json
from dataclasses import asdict, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from ntv3_tools.annotation_report import ANNOTATION_PLAN
from ntv3_tools.errors import RequestError
from ntv3_tools.plan import AnalysisPlan, read_plan, validate_plan

RESOLVED_BLOCK = "resolved.txt"
APPROVAL = "approval.json"
ANNOTATION_RESOLVED_BLOCK = "annotation_resolved.txt"
ANNOTATION_APPROVAL = "annotation_approval.json"


def write_resolved_block(run_dir: str | Path, lines: list[str]) -> Path:
    """Write the Resolved block the user approves.

    Args:
        run_dir: Run directory holding ``plan.yaml``.
        lines: Resolved block lines.

    Returns:
        Path to ``resolved.txt``.
    """
    path = Path(run_dir) / RESOLVED_BLOCK
    path.write_text("\n".join(lines) + "\n")
    return path


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def plan_digest(plan: AnalysisPlan) -> str:
    """Return a digest of every plan field except ``approved`` itself."""
    payload = asdict(plan)
    payload.pop("approved")
    return _sha256(json.dumps(payload, sort_keys=True))


def _write_approval(record: Path, plan_sha256: str, block: Path, **extra: Any) -> None:
    """Record the digests of the approved plan and of the block shown for it."""
    payload = {
        "approved_at": datetime.now(UTC).isoformat(),
        "plan_sha256": plan_sha256,
        "resolved_sha256": _sha256(block.read_text()),
        **extra,
    }
    record.write_text(json.dumps(payload, indent=2) + "\n")


def _changed_since_approval(record: Path, plan_sha256: str, block: Path) -> bool:
    """Return whether the plan or its block differs from what was approved."""
    saved = json.loads(record.read_text())
    return saved.get("plan_sha256") != plan_sha256 or saved.get(
        "resolved_sha256"
    ) != _sha256(block.read_text())


def approve_run(run_dir: str | Path, critic: bool = True) -> AnalysisPlan:
    """Approve a planned run after its Resolved block exists.

    Args:
        run_dir: Run directory holding ``plan.yaml`` and ``resolved.txt``.
        critic: False turns the critic off; the default keeps the plan's.

    Returns:
        The approved plan.

    Raises:
        RequestError: No plan, no Resolved block, or already approved.
    """
    root = Path(run_dir)
    block = root / RESOLVED_BLOCK
    if not block.is_file():
        raise RequestError(
            f"Which run should be approved? No Resolved block at {block}; plan "
            "that run first. No run is approved for you."
        )
    plan = read_plan(path=root / "plan.yaml")
    if plan.approved:
        raise RequestError(
            f"Which run should be approved? {root} is already approved. Pass "
            "--run-dir of an unapproved run. No run is assumed."
        )
    approved = replace(plan, approved=True, want_critic=plan.want_critic and critic)
    validate_plan(plan=approved)
    (root / "plan.yaml").write_text(yaml.safe_dump(asdict(approved), sort_keys=False))
    _write_approval(
        root / APPROVAL, plan_digest(approved), block, critic=approved.want_critic
    )
    return approved


def require_approval(run_dir: str | Path, plan: AnalysisPlan) -> None:
    """Refuse a run the user has not approved, before any fetch or load.

    Args:
        run_dir: Run directory holding ``plan.yaml``.
        plan: The plan read from it.

    Raises:
        RequestError: Not approved, or the plan or block changed after approval.
    """
    root = Path(run_dir)
    command = f"ntv3 approve --run-dir {root}"
    if not plan.approved:
        raise RequestError(
            f"{root} is not approved. Which run should be predicted? Approve this "
            f"one with `{command}`. No run is approved for you; no model was loaded."
        )
    block = root / RESOLVED_BLOCK
    record_path = root / APPROVAL
    if not block.is_file() or not record_path.is_file():
        raise RequestError(
            f"{root} says approved but has no {RESOLVED_BLOCK} or {APPROVAL}; only "
            "`ntv3 approve` approves a run. No model was loaded."
        )
    if _changed_since_approval(record_path, plan_digest(plan), block):
        raise RequestError(
            f"{root}: the plan or its Resolved block changed after approval. Which "
            "run should be predicted? Plan and approve a new one. No model was "
            "loaded."
        )


def annotation_plan_digest(plan: dict[str, Any]) -> str:
    """Return a digest of every annotation plan field except ``approved``."""
    payload = {key: value for key, value in plan.items() if key != "approved"}
    return _sha256(json.dumps(payload, sort_keys=True))


def pending_annotation_plan(run_dir: str | Path) -> bool:
    """Return whether the run holds an annotation plan not yet approved."""
    path = Path(run_dir) / ANNOTATION_PLAN
    return path.is_file() and not json.loads(path.read_text()).get("approved")


def approve_annotation(run_dir: str | Path, critic: bool = True) -> dict[str, Any]:
    """Approve a planned annotation after its Resolved block exists.

    Args:
        run_dir: Run directory holding ``annotation_plan.json`` and its block.
        critic: False turns the critic off; the default keeps the plan's.

    Returns:
        The approved annotation plan.

    Raises:
        RequestError: No plan, no Resolved block, or already approved.
    """
    root = Path(run_dir)
    block = root / ANNOTATION_RESOLVED_BLOCK
    plan_path = root / ANNOTATION_PLAN
    if not plan_path.is_file() or not block.is_file():
        raise RequestError(
            f"Which annotation should be approved? {root} has no {ANNOTATION_PLAN} "
            f"and {ANNOTATION_RESOLVED_BLOCK}; plan it with `ntv3 annotate` first."
        )
    plan = json.loads(plan_path.read_text())
    if plan.get("approved"):
        raise RequestError(f"The annotation in {root} is already approved.")
    plan["approved"] = True
    plan["want_critic"] = bool(plan.get("want_critic", True)) and critic
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    _write_approval(root / ANNOTATION_APPROVAL, annotation_plan_digest(plan), block)
    return plan


def require_annotation_approval(run_dir: str | Path, plan: dict[str, Any]) -> None:
    """Refuse an annotation the user has not approved, before any load.

    Args:
        run_dir: Run directory holding ``annotation_plan.json``.
        plan: The plan read from it.

    Raises:
        RequestError: Not approved, or the plan or block changed after approval.
    """
    root = Path(run_dir)
    if not plan.get("approved"):
        raise RequestError(
            f"The annotation in {root} is not approved. Approve it with "
            f"`ntv3 approve --run-dir {root}`. No model was loaded."
        )
    block = root / ANNOTATION_RESOLVED_BLOCK
    record = root / ANNOTATION_APPROVAL
    if (
        not block.is_file()
        or not record.is_file()
        or _changed_since_approval(record, annotation_plan_digest(plan), block)
    ):
        raise RequestError(
            f"{root}: the annotation plan or its Resolved block changed after "
            "approval, or `ntv3 approve` did not approve it. Plan and approve a "
            "new one. No model was loaded."
        )
