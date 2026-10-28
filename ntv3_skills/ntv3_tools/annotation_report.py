"""Read a finished NTv3 GFF3 into the report's annotation lane.

The GFF3 is what ``ntv3-annotate`` writes, through ``ntv3 annotate-forward``: ``gene``, ``CDS``, and ``intron`` rows, 1-based closed, with ``##ntv3-``
provenance comments. Each gene becomes one browser lane feature with its CDS runs
as exon blocks, in 0-based half-open coordinates. Explain values are read from
the header or ``annotation_plan.json`` and never recomputed. No per-base
probabilities are written, so none are reported.
"""

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ntv3_tools.errors import ArtifactError

ANNOTATION_LABEL = "NTv3 annotation"
ANNOTATION_PLAN = "annotation_plan.json"
ANNOTATION_MANIFEST = "annotation_manifest.json"
_GFF3_PATTERNS = ("*.gff3", "*.gff3.gz", "*.gff", "*.gff.gz")


@dataclass(frozen=True)
class Gff3Genes:
    """Genes and header facts from one GFF3.

    Attributes:
        provenance: ``##ntv3-<key> <value>`` comments, keyed without the prefix.
        regions: ``##sequence-region`` lengths as {"chrom", "length"}.
        genes: Report genes: chrom, name, start, end, strand, and CDS exons.
    """

    provenance: dict[str, str]
    regions: list[dict[str, Any]]
    genes: list[dict[str, Any]]


def _gff3_text(path: Path) -> str:
    """Return the text of a plain or gzipped GFF3."""
    if path.name.lower().endswith(".gz"):
        return gzip.decompress(path.read_bytes()).decode("utf-8")
    return path.read_text()


def read_gff3(path: str | Path) -> Gff3Genes:
    """Read genes and their CDS runs from an NTv3 GFF3.

    Args:
        path: Plain or gzipped GFF3.

    Returns:
        Provenance, sequence regions, and genes in file order.

    Raises:
        ArtifactError: A row is malformed or a CDS names a gene that is absent.
    """
    gff3_path = Path(path)
    provenance: dict[str, str] = {}
    regions: list[dict[str, Any]] = []
    genes: dict[str, dict[str, Any]] = {}
    for number, raw in enumerate(_gff3_text(gff3_path).splitlines(), start=1):
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("##ntv3-"):
            key, _, value = line[len("##ntv3-") :].partition(" ")
            provenance[key] = value.strip()
            continue
        if line.startswith("##sequence-region"):
            fields = line.split()
            if len(fields) == 4:
                regions.append({"chrom": fields[1], "length": int(fields[3])})
            continue
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 9:
            raise ArtifactError(f"{gff3_path} line {number} does not have 9 columns.")
        chrom, _, kind, start, end, _, strand, _, attributes = cols
        if kind not in ("gene", "CDS"):
            continue
        fields = dict(
            item.split("=", 1) for item in attributes.split(";") if "=" in item
        )
        try:
            first, last = int(start), int(end)
        except ValueError as exc:
            raise ArtifactError(
                f"{gff3_path} line {number} has non-integer coordinates."
            ) from exc
        if first < 1 or last < first or strand not in ("+", "-"):
            raise ArtifactError(
                f"{gff3_path} line {number} needs 1 <= start <= end and strand + or -."
            )
        span = {"start": first - 1, "end": last}
        if kind == "gene":
            name = fields.get("ID") or f"gene_{len(genes) + 1}"
            genes[name] = {
                "chrom": chrom,
                "name": name,
                **span,
                "strand": strand,
                "exons": [],
            }
            continue
        parent = genes.get(fields.get("Parent", ""))
        if parent is None:
            raise ArtifactError(
                f"{gff3_path} line {number}: CDS parent "
                f"'{fields.get('Parent', '')}' is not a gene above it."
            )
        parent["exons"].append(span)
    return Gff3Genes(provenance=provenance, regions=regions, genes=list(genes.values()))


def find_gff3(root: Path, plan: dict[str, Any] | None) -> Path | None:
    """Return the run's GFF3: the plan's output if written, else the only one.

    Args:
        root: Run directory.
        plan: ``annotation_plan.json`` payload, or None.

    Returns:
        The GFF3 path, or None when the run has none.

    Raises:
        ArtifactError: More than one GFF3 sits in the run directory.
    """
    named = (plan or {}).get("output")
    if named and (root / str(named)).is_file():
        return root / str(named)
    found = sorted({path for pattern in _GFF3_PATTERNS for path in root.glob(pattern)})
    if len(found) > 1:
        raise ArtifactError(
            f"More than one GFF3 in {root}: {', '.join(path.name for path in found)}. "
            "Keep one per run."
        )
    return found[0] if found else None


def _int_or_none(value: Any) -> int | None:
    """Return an int, or None when the value was not recorded."""
    return None if value in (None, "") else int(value)


def _float_or_none(value: Any) -> float | None:
    """Return a float, or None when the value was not recorded."""
    return None if value in (None, "") else float(value)


def _bool_or_none(value: Any) -> bool | None:
    """Return a bool from a header string or a JSON bool, or None."""
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def annotation_payload(root: Path) -> dict[str, Any] | None:
    """Return the report's annotation block, or None when the run has no GFF3.

    Header values win over ``annotation_plan.json`` because they record what
    produced the file. A value neither recorded stays None.

    Args:
        root: Run directory.

    Returns:
        JSON-ready annotation block, or None.

    Raises:
        ArtifactError: The GFF3 is malformed or ambiguous.
    """
    plan_path = root / ANNOTATION_PLAN
    plan = json.loads(plan_path.read_text()) if plan_path.is_file() else None
    if plan is not None and not isinstance(plan, dict):
        raise ArtifactError(f"{plan_path} must contain a JSON object.")
    gff3_path = find_gff3(root=root, plan=plan)
    if gff3_path is None:
        return None
    gff3 = read_gff3(gff3_path)
    header = gff3.provenance
    planned = plan or {}

    def recorded(key: str, plan_key: str) -> Any:
        return header.get(key) or planned.get(plan_key)

    regions = gff3.regions or [
        {"chrom": str(row["seqid"]), "length": int(row["length"])}
        for row in planned.get("records") or []
        if isinstance(row, dict)
    ]
    return {
        "label": ANNOTATION_LABEL,
        "model": planned.get("model"),
        "modelName": header.get("model-name"),
        "checkpoint": header.get("checkpoint"),
        "clade": planned.get("clade"),
        "decoder": recorded("decoder", "decoder"),
        "windowBp": _int_or_none(recorded("sequence-length", "window_bp")),
        "keepFraction": _float_or_none(
            recorded("keep-target-center-fraction", "keep_fraction")
        ),
        "minCdsLength": _int_or_none(recorded("min-cds-length", "min_cds_length")),
        "inFrameStopFilter": _bool_or_none(
            recorded("apply-biological-filters", "apply_biological_filters")
        ),
        "gff3": gff3_path.relative_to(root).as_posix(),
        "sources": [
            *([ANNOTATION_PLAN] if plan is not None else []),
            gff3_path.relative_to(root).as_posix(),
        ],
        "regions": regions,
        "genes": gff3.genes,
    }
