"""Run one approved Annotate plan on a GPU and write its GFF3.

The forward is ``nucleotide_transformer_v3.gene_annotation``, the pipeline behind
``ntv3-annotate``, called with the plan's FASTA, clade, and decoder. A locus
plan's GFF3 is shifted to chromosome coordinates so the report draws it next to
the run's tracks.
"""

import json
import logging
import tempfile
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import jax

from ntv3_tools.annotate import ANNOTATION_OUTPUT, LocusWindow, shift_gff3
from ntv3_tools.annotation_report import ANNOTATION_MANIFEST, read_gff3
from ntv3_tools.approval import require_annotation_approval
from ntv3_tools.errors import ArtifactError, RequestError

logger = logging.getLogger("ntv3_tools")


def require_gpu() -> Any:
    """Return the first GPU JAX sees, or refuse.

    Raises:
        RequestError: JAX sees no GPU.
    """
    try:
        return jax.devices("gpu")[0]
    except RuntimeError as exc:
        raise RequestError(
            "The annotation forward needs a GPU and JAX sees none. On a GPU host, "
            "install the gpu extra: `.venv/bin/pip install -e '.[gpu]'`. No model "
            "was loaded."
        ) from exc


def run_annotation_forward(
    plan_path: str | Path,
    annotate_main: Callable[[list[str]], int] | None = None,
) -> Path:
    """Annotate an approved plan's input and write ``annotation.gff3``.

    Args:
        plan_path: The approved ``annotation_plan.json``.
        annotate_main: ``ntv3-annotate``'s entry point; the package's by default.

    Returns:
        Path to the run's ``annotation.gff3``.

    Raises:
        RequestError: The plan is not approved, or JAX sees no GPU.
        ArtifactError: The run already has a GFF3, or the pipeline wrote none.
    """
    plan_file = Path(plan_path).resolve()
    run_dir = plan_file.parent
    plan = json.loads(plan_file.read_text())
    require_annotation_approval(run_dir=run_dir, plan=plan)
    output = run_dir / ANNOTATION_OUTPUT
    if output.exists():
        raise ArtifactError(f"{output} exists; a run holds one annotation.")
    device = require_gpu()
    if annotate_main is None:
        from nucleotide_transformer_v3.gene_annotation.cli import main as annotate_main

    with tempfile.TemporaryDirectory() as tmp:
        # A locus window is annotated in window coordinates, then shifted.
        written = Path(tmp, ANNOTATION_OUTPUT) if plan["locus"] else output
        argv = [
            "--genome", str(plan["fasta"]),
            "--model_cfg", str(plan["clade"]),
            "--out", str(written),
            "--override", f"decoder={plan['decoder']}",
        ]  # fmt: skip
        logger.info(
            f"Annotating on {device.device_kind}: ntv3-annotate {' '.join(argv)}"
        )
        started = time.monotonic()
        status = annotate_main(argv)
        seconds = round(time.monotonic() - started, 1)
        if status or not written.is_file():
            raise ArtifactError(f"ntv3-annotate exited {status} and wrote no GFF3.")
        if plan["locus"]:
            shift_gff3(
                source=written, destination=output, window=LocusWindow(**plan["locus"])
            )

    genes = len(read_gff3(output).genes)
    manifest = {
        "model": plan["model"],
        "clade": plan["clade"],
        "decoder": plan["decoder"],
        "device": device.device_kind,
        "seconds": seconds,
        "genes": genes,
    }
    (run_dir / ANNOTATION_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    logger.info(f"Wrote {output}: {genes} genes in {seconds} s.")
    return output
