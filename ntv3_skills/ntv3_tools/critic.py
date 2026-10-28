"""KEEP / SUPPRESS claims against on-disk evidence. No fake confidence %."""

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Claim:
    """One report sentence the critic can keep or drop.

    Attributes:
        claim_id: Stable id.
        text: The sentence.
        evidence_paths: Files that must exist to KEEP the claim, relative to the
            run directory and inside it. At least one is required.
        citekey: Optional bib key. Required for literature claims.
        kind: "model" | "literature" | "other".
    """

    claim_id: str
    text: str
    evidence_paths: list[str]
    citekey: str | None
    kind: str


@dataclass(frozen=True)
class CriticDecision:
    """One KEEP / SUPPRESS row.

    Attributes:
        claim_id: Claim id.
        text: Original sentence.
        verdict: KEEP or SUPPRESS.
        reason: Why.
    """

    claim_id: str
    text: str
    verdict: str
    reason: str


def _inside(root: Path, path: Path) -> bool:
    """Return whether ``path``, with symlinks followed, is below ``root``.

    The run directory itself ("" or ".") is not inside it.
    """
    return root in path.resolve().parents


def audit_claims(
    claims: list[Claim],
    run_dir: str,
    bib_path: str | None = None,
) -> list[CriticDecision]:
    """Drop unsupported claims. Keep rejects in the returned log.

    Args:
        claims: Candidate sentences.
        run_dir: Run directory (evidence paths are relative to this).
        bib_path: Optional citations.bib. Lit claims need a citekey in it.

    Returns:
        One decision per claim.
    """
    root = Path(run_dir).resolve()
    bib_text = (
        Path(bib_path).read_text() if bib_path and Path(bib_path).exists() else ""
    )
    decisions: list[CriticDecision] = []
    for claim in claims:
        if not claim.evidence_paths:
            decisions.append(
                CriticDecision(
                    claim_id=claim.claim_id,
                    text=claim.text,
                    verdict="SUPPRESS",
                    reason="No evidence path",
                )
            )
            continue
        outside = [
            path
            for path in claim.evidence_paths
            if Path(path).is_absolute() or not _inside(root=root, path=root / path)
        ]
        if outside:
            decisions.append(
                CriticDecision(
                    claim_id=claim.claim_id,
                    text=claim.text,
                    verdict="SUPPRESS",
                    reason=f"Evidence outside the run directory: {outside}",
                )
            )
            continue
        # A directory is not evidence; each path names one file.
        missing = [path for path in claim.evidence_paths if not (root / path).is_file()]
        if missing:
            decisions.append(
                CriticDecision(
                    claim_id=claim.claim_id,
                    text=claim.text,
                    verdict="SUPPRESS",
                    reason=f"Missing evidence: {missing}",
                )
            )
            continue
        if claim.kind == "literature" and (
            not claim.citekey or f"{{{claim.citekey}" not in bib_text
        ):
            decisions.append(
                CriticDecision(
                    claim_id=claim.claim_id,
                    text=claim.text,
                    verdict="SUPPRESS",
                    reason="Literature claim has no @citekey in citations.bib",
                )
            )
            continue
        lowered = claim.text.lower()
        if (
            claim.kind.lower().startswith("causal")
            or "proven" in lowered
            or "causal loop" in lowered
        ):
            decisions.append(
                CriticDecision(
                    claim_id=claim.claim_id,
                    text=claim.text,
                    verdict="SUPPRESS",
                    reason="Overclaim (proven / causal). Model output is not causality.",
                )
            )
            continue
        decisions.append(
            CriticDecision(
                claim_id=claim.claim_id,
                text=claim.text,
                verdict="KEEP",
                reason="Evidence files present",
            )
        )
    return decisions


def write_audit_log(decisions: list[CriticDecision], out_dir: str) -> Path:
    """Write critic/audit_log.json.

    Args:
        decisions: KEEP / SUPPRESS rows.
        out_dir: Run output directory.

    Returns:
        Path to the audit log.
    """
    root = Path(out_dir) / "critic"
    root.mkdir(parents=True, exist_ok=True)
    path = root / "audit_log.json"
    path.write_text(json.dumps([asdict(row) for row in decisions], indent=2) + "\n")
    return path
