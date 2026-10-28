"""Build validated NTv3 report artifacts from one approved run directory."""

import json
import logging
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from ntv3_tools.annotation_report import (
    ANNOTATION_MANIFEST,
    ANNOTATION_PLAN,
    annotation_payload,
    read_gff3,
)
from ntv3_tools.browser import round_sig, saved_ucsc, write_browser_html
from ntv3_tools.critic import Claim, audit_claims, write_audit_log
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.plan import read_plan
from ntv3_tools.report_app import write_report_app
from ntv3_tools.reproducibility import (
    refuse_generated_overwrite,
    write_analysis_notebook,
    write_replay_script,
    write_ro_crate,
)

logger = logging.getLogger("ntv3_tools")

# Saved UCSC row fields the report app reads (report/src/browser/ucsc.ts).
UCSC_ROW_FIELDS = frozenset(
    {
        "blockSizes", "ccre", "cCRE_class", "cdsEnd", "cdsStart", "chrom",
        "chromEnd", "chromStart", "chromStarts", "description", "encodeLabel",
        "end", "exonEnds", "exonStarts", "geneName", "genoEnd", "genoStart",
        "name", "name2", "repClass", "repName", "start", "strand", "thickEnd",
        "thickStart", "tier", "txEnd", "txStart", "value",
    }
)  # fmt: skip


def slim_ucsc(payload: dict[str, Any]) -> dict[str, Any]:
    """Keep a saved UCSC payload's envelope and only the row fields the page draws."""
    return {
        key: [
            {field: v for field, v in row.items() if field in UCSC_ROW_FIELDS}
            for row in value
        ]
        if isinstance(value, list) and value and isinstance(value[0], dict)
        else value
        for key, value in payload.items()
    }


def _rounded(values: list[float | None]) -> list[float | None]:
    """Four significant digits for the page; the run directory keeps full precision."""
    return [None if value is None else round_sig(value) for value in values]


REPORT_STAGE_LABELS = {
    "disruption": "Compare (SNV)",
    "literature": "Literature",
    "claims": "Claim review",
}
# Optional analyses appear once their artifacts exist. None is shown as pending.
DEFAULT_REPORT_STAGES: tuple[str, ...] = ()
_EXPLICIT_ONLY_STEPS = frozenset({"literature"})


@dataclass
class ReportSections:
    """The human-facing content of a report.

    Attributes:
        title: Report title.
        conclusion: Plain-language, calibrated conclusion (first thing a biologist reads).
        limitations: Bullet points on caveats / what was not done.
        next_steps: Bullet points on recommended follow-up (may name deferred skills).
    """

    title: str
    conclusion: str
    limitations: list[str] = field(default_factory=list)
    next_steps: list[str] = field(default_factory=list)


def save_report_bundle(
    out_dir: str,
    manifest: dict[str, Any],
    sections: ReportSections,
    tables: dict[str, pd.DataFrame] | None = None,
    figure_paths: list[str] | None = None,
) -> str:
    """Write a complete results bundle to disk.

    Args:
        out_dir: Destination results directory (created if absent).
        manifest: Reproducibility manifest (serialized to run_manifest.json).
        sections: Report prose (title, conclusion, limitations, next steps).
        tables: Named tables written to tables/<name>.csv and linked in the report.
        figure_paths: Figure paths (already saved); embedded in the report as images.

    Returns:
        The output directory path.
    """
    root = Path(out_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "figures").mkdir(parents=True, exist_ok=True)

    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))

    table_links: list[str] = []
    for name, frame in (tables or {}).items():
        rel = f"tables/{name}.csv"
        frame.to_csv(root / rel, index=False)
        table_links.append(rel)

    report_md = _render_report_md(
        manifest=manifest,
        sections=sections,
        tables=tables or {},
        figure_paths=figure_paths or [],
    )
    (root / "report.md").write_text(report_md)
    logger.info(f"Wrote report bundle to {root}")
    return str(root)


def _render_report_md(
    manifest: dict[str, Any],
    sections: ReportSections,
    tables: dict[str, pd.DataFrame],
    figure_paths: list[str],
) -> str:
    """Render report.md content.

    Args:
        manifest: Reproducibility manifest.
        sections: Report prose.
        tables: Named tables (rendered inline as markdown).
        figure_paths: Figure paths to embed.

    Returns:
        The report markdown string.
    """
    lines: list[str] = [
        f"# {sections.title}",
        "",
        "## Conclusion",
        sections.conclusion,
        "",
    ]

    lines += ["## Inputs & context"]
    for key in ("model", "genome_build", "objective"):
        if key in manifest:
            lines.append(f"- {key}: {manifest[key]}")
    for key, value in manifest.items():
        if key in {"model", "genome_build", "objective", "platform"}:
            continue
        lines.append(f"- {key}: {value}")
    lines.append("")

    for name, frame in tables.items():
        rendered_table = frame.to_markdown(index=False)
        if not isinstance(rendered_table, str):
            raise ArtifactError(f"Could not render report table '{name}'.")
        lines += [
            f"## {name.replace('_', ' ').title()}",
            "",
            rendered_table,
            "",
        ]

    if figure_paths:
        lines.append("## Figures")
        for path in figure_paths:
            rel = path.split("figures/")[-1]
            lines.append(f"![{rel}](figures/{rel})")
        lines.append("")

    if sections.limitations:
        lines += ["## Limitations", *[f"- {item}" for item in sections.limitations], ""]
    if sections.next_steps:
        lines += [
            "## Suggested follow-up",
            *[f"- {item}" for item in sections.next_steps],
            "",
        ]

    lines += [
        "## Reproducibility",
        "See run_manifest.json (model, windows, device, seed, input hash).",
        "",
    ]
    return "\n".join(lines)


def write_viewer_bundle(
    run_dir: str,
    *,
    refresh: bool = False,
    requested_stages: list[str] | None = None,
    contig: str | None = None,
) -> Path:
    """Assemble report.json, ntv3_report.html, browser.html, notebook, replay, RO-Crate.

    ``ntv3_report.html`` is the page the Artifact tool publishes; nobody opens
    it from disk. ``browser.html`` is the older static page and is not published.

    A Predict run needs the approved plan and prediction artifacts; a GFF3 in
    it adds the annotation lane. An Annotate run needs ``annotation_plan.json``
    and its GFF3. Optional analyses appear once their artifacts exist; a stage
    is shown as pending only when ``requested_stages`` names it. The report
    never invents data.

    Args:
        run_dir: Run output directory.
        refresh: Replace only generated report metadata after a new analysis.
        requested_stages: Stages to show as pending until their artifacts exist.
            None keeps the previous report's selection, or shows none.
        contig: For an Annotate run, the contig the page shows; by default the
            one with the most genes.

    Returns:
        Path to ``report.json``.
    """
    root = Path(run_dir)
    if not root.is_dir():
        raise ArtifactError(f"Run directory does not exist: {root}.")
    report_path = root / "report.json"
    if refresh:
        if not report_path.is_file():
            raise ArtifactError(
                "Cannot refresh before the initial report exists. Run report once "
                "after prediction."
            )
    else:
        refuse_generated_overwrite(root=root)
    selected_stages = _resolve_requested_stages(
        existing_report=_read_json(report_path) if refresh else None,
        requested_stages=requested_stages,
    )
    if not (root / "plan.yaml").exists() and (root / ANNOTATION_PLAN).is_file():
        report_path.write_text(
            json.dumps(_annotation_report(root=root, contig=contig), indent=2) + "\n"
        )
        write_browser_html(run_dir=root)
        write_report_app(run_dir=root)
        write_ro_crate(root=root, main_entity=ANNOTATION_PLAN)
        logger.info(f"Wrote viewer bundle to {report_path}")
        return report_path
    if contig is not None:
        raise RequestError(
            "--contig chooses the contig of an Annotate run's page; a Predict run "
            "shows its own window."
        )
    try:
        plan = asdict(read_plan(path=root / "plan.yaml"))
    except RequestError as exc:
        raise ArtifactError(str(exc)) from exc
    resolve = _require_mapping(root / "track_resolve.json")
    critic_path = root / "critic" / "audit_log.json"
    critic_status = _run_critic(
        root=root,
        enabled=bool(plan.get("want_critic")),
        claims=_run_claims(root=root, plan=plan),
        log_path=critic_path,
    )
    critic_payload = _read_json(critic_path)
    if critic_path.exists() and (
        not isinstance(critic_payload, list) or not critic_payload
    ):
        raise ArtifactError(
            "The critic ran but critic/audit_log.json contains no decisions."
        )
    critic = critic_payload or []
    tracks_npy = root / "tracks" / "tracks.npy"
    accessions_path = root / "tracks" / "accessions.json"
    track_metadata_path = root / "tracks" / "metadata.json"
    if (
        not tracks_npy.exists()
        or not accessions_path.exists()
        or not track_metadata_path.exists()
    ):
        raise ArtifactError(
            "Report requires tracks/tracks.npy, tracks/accessions.json, and "
            "tracks/metadata.json. "
            "Run prediction before report generation."
        )
    gene = str(plan["gene"])
    cell = str(plan["tissue"])
    assembly = str(plan["assembly"])
    chrom, start, end = _locus_coords(plan=plan)
    gt_rows = _read_json(root / "gt" / "manifest.json") or []
    run_manifest = _read_json(root / "run_manifest.json") or {}
    # Predict records False when its Python could not read BigWigs.
    gt_compared = bool(gt_rows) and run_manifest.get("ground_truth_compared", True)
    assays = _assay_names(plan=plan, resolve=resolve)
    model_heads = _require_string_list(accessions_path)
    track_metadata = _require_mapping(track_metadata_path)
    try:
        track_chrom = str(track_metadata["chrom"])
        track_start = int(track_metadata["start"])
        track_end = int(track_metadata["end"])
        track_bin_size = int(track_metadata["bin_size"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ArtifactError(
            "tracks/metadata.json requires chrom, integer start, integer end, "
            "and integer bin_size."
        ) from exc
    if (
        track_chrom != chrom
        or track_start < start
        or track_end > end
        or track_end <= track_start
    ):
        raise ArtifactError(
            f"Prediction span {track_chrom}:{track_start}-{track_end} is not "
            f"inside the planned locus {chrom}:{start}-{end}."
        )
    predicted, observed = _load_track_arrays(
        tracks_npy=tracks_npy,
        observed_npy=root / "gt" / "observed.npy",
        n_assays=len(assays),
        require_observed=gt_compared,
    )
    if len(predicted) != len(assays) or len(model_heads) != len(assays):
        raise ArtifactError(
            f"Track artifact mismatch: {len(assays)} requested assays, "
            f"{len(model_heads)} model heads, and {len(predicted)} prediction rows."
        )
    prediction_span = track_end - track_start
    if track_bin_size <= 0 or any(
        len(values) * track_bin_size != prediction_span for values in predicted
    ):
        raise ArtifactError(
            f"Prediction arrays do not match span {track_start}-{track_end} "
            f"at bin_size={track_bin_size}."
        )
    annotation = annotation_payload(root=root)
    saliency_payload = _load_saliency(root=root / "interpret")
    disruption_payload = _load_disruption(root=root)
    literature_payload = _load_literature(root=root)
    labels = [
        lane_label(query=query, model_head=head)
        for query, head in zip(assays, model_heads, strict=True)
    ]
    labels_by_model_head = dict(zip(model_heads, labels, strict=True))
    for saliency_track in saliency_payload:
        model_head = str(saliency_track["modelHead"])
        saliency_track["name"] = labels_by_model_head.get(model_head, model_head)
    if not refresh:
        write_replay_script(root=root)
        if bool(plan.get("want_notebook")):
            write_analysis_notebook(root=root, gene=gene, cell=cell)
    claims = _claims_from_critic(critic=critic)
    # locus_metadata.json exists only with a knownGene lane; otherwise the TSS
    # and strand come from locus.json, whose source is not a {track, url} lane.
    locus_metadata = (
        _read_json(root / "locus_metadata.json")
        or _read_json(root / "locus.json")
        or {}
    )
    agreement_rows = _load_track_agreement(
        path=root / "tables" / "track_agreement.csv",
        required=gt_compared,
    )
    completed_steps = {"report", "ro-crate"}
    payload = {
        "tool": "ntv3 0.1.0",
        "commit": run_manifest.get("commit", "unrecorded"),
        "gene": gene,
        "cell": cell,
        "assembly": assembly,
        "modelShort": str(plan["model"]),
        "windowBp": max(end - start, 1),
        "locus": {
            "chrom": chrom,
            "start": start,
            "end": end,
            "tss": locus_metadata.get("tss"),
            "strand": locus_metadata.get("strand"),
            "focusStart": track_start,
            "focusEnd": track_end,
            "elements": list(locus_metadata.get("elements") or []),
            "source": _lane_source(locus_metadata.get("source")),
            "assays": labels,
        },
        "tracks": {
            "chrom": chrom,
            "start": track_start,
            "end": track_end,
            "assays": [
                _track_assay(
                    name=name,
                    model_head=model_heads[i],
                    predicted=predicted[i],
                    observed=observed[i],
                )
                for i, name in enumerate(labels)
            ],
        },
        "claims": claims,
        "trace": _trace_from_tree(root=root, completed_steps=completed_steps),
        "runParameters": {
            "run_dir": ".",
            "species": str(plan["species"]),
        },
        "critic": critic_status,
        "ucscTracks": list(plan.get("ucsc_tracks") or []),
        # The UCSC payloads lookup saved, so the report app draws them offline.
        "ucscData": {
            track: slim_ucsc(payload)
            for track, payload in saved_ucsc(root=root).items()
            if track in (plan.get("ucsc_tracks") or [])
        },
        "explain": _explain(
            root=root,
            plan=plan,
            run_manifest=run_manifest,
            prediction=(track_chrom, track_start, track_end, track_bin_size),
            track_rows=[
                {"query": query, "label": label, "modelHead": head}
                for query, label, head in zip(assays, labels, model_heads, strict=True)
            ],
            saliency_tracks=saliency_payload,
        ),
        "gtSources": _ground_truth_sources(
            gt_rows=gt_rows,
            resolve=resolve,
            model_heads_by_query=dict(zip(assays, model_heads, strict=True)),
        ),
        "trackAgreement": agreement_rows,
        "saliencyTracks": saliency_payload,
        "disruption": disruption_payload,
        "literature": literature_payload,
        "annotation": annotation,
        "analysisStages": _analysis_stages(
            requested_stages=selected_stages,
            completed={
                "disruption": disruption_payload is not None,
                "literature": literature_payload is not None,
                "claims": bool(critic),
            },
        ),
    }
    report_path.write_text(json.dumps(payload, indent=2) + "\n")
    write_browser_html(run_dir=root)
    write_report_app(run_dir=root)
    write_ro_crate(root=root)
    logger.info(f"Wrote viewer bundle to {report_path}")
    return report_path


def _annotation_report(root: Path, contig: str | None = None) -> dict[str, Any]:
    """Build report.json for an Annotate run: the browser on one contig's genes.

    Args:
        root: Run directory holding ``annotation_plan.json`` and its GFF3.
        contig: Contig to show; by default the one with the most genes.

    Returns:
        JSON-ready report with no predicted tracks.

    Raises:
        ArtifactError: The run has no GFF3 yet, or it has nothing to show.
        RequestError: ``contig`` is not in the GFF3.
    """
    annotation = annotation_payload(root=root)
    if annotation is None:
        raise ArtifactError(
            f"{root} has {ANNOTATION_PLAN} and no GFF3: the annotation forward "
            "has not run. Nothing to report yet."
        )
    lengths = {str(row["chrom"]): int(row["length"]) for row in annotation["regions"]}
    genes = annotation["genes"]
    counts = Counter(str(gene["chrom"]) for gene in genes)
    contigs = list(dict.fromkeys([*lengths, *counts]))
    if not contigs:
        raise ArtifactError(
            f"{annotation['gff3']} has no genes and no ##sequence-region to show."
        )
    if contig is None:
        chrom = max(contigs, key=lambda name: counts[name])
    elif contig in contigs:
        chrom = contig
    else:
        shown = ", ".join(contigs[:10]) + (", …" if len(contigs) > 10 else "")
        raise RequestError(
            f"Which contig? {contig} is not in {annotation['gff3']}, which has "
            f"{len(contigs)}: {shown}."
        )
    end = lengths.get(chrom) or max(
        int(gene["end"]) for gene in genes if gene["chrom"] == chrom
    )
    plan = _require_mapping(root / ANNOTATION_PLAN)
    critic_path = root / "critic" / "audit_log.json"
    critic_status = _run_critic(
        root=root,
        enabled=bool(plan.get("want_critic", True)),
        claims=_annotation_claims(root=root),
        log_path=critic_path,
    )
    model = str(annotation["model"] or annotation["modelName"] or "unrecorded")
    trace = [
        {"name": name, "artifact": rel, "status": "ok"}
        for name, rel in (
            ("annotate", ANNOTATION_PLAN),
            ("gff3", annotation["gff3"]),
            ("report", "report.json"),
            ("ro-crate", "ro-crate-metadata.json"),
        )
    ]
    return {
        "tool": "ntv3 0.1.0",
        "commit": "unrecorded",
        "gene": annotation["label"],
        "cell": Path(str(plan.get("fasta") or annotation["gff3"])).name,
        "assembly": "",
        "modelShort": model,
        "windowBp": int(annotation["windowBp"] or end),
        "locus": {
            "chrom": chrom,
            "start": 0,
            "end": end,
            "tss": None,
            "strand": None,
            "focusStart": 0,
            "focusEnd": end,
            "elements": [],
            "source": None,
            "assays": [],
        },
        "tracks": {"chrom": chrom, "start": 0, "end": end, "assays": []},
        "claims": _claims_from_critic(critic=_read_json(critic_path) or []),
        "trace": trace,
        "runParameters": {
            "run_dir": ".",
            "clade": str(annotation["clade"] or "unrecorded"),
            "contig": chrom,
            "genes_on_contig": counts[chrom],
            "contigs": len(contigs),
            "genes": len(genes),
        },
        "critic": critic_status,
        "ucscTracks": [],
        "ucscData": {},
        "explain": None,
        "gtSources": [],
        "trackAgreement": [],
        "saliencyTracks": [],
        "disruption": None,
        "literature": None,
        "annotation": annotation,
        "analysisStages": [],
    }


def _explain(
    root: Path,
    plan: dict[str, Any],
    run_manifest: dict[str, Any],
    prediction: tuple[str, int, int, int],
    track_rows: list[dict[str, str]],
    saliency_tracks: list[dict[str, Any]],
) -> dict[str, Any]:
    """Say which model, window, crop, and head produced each reported number.

    Reads only provenance already on disk. Nothing is recomputed.

    Args:
        root: Approved run directory.
        plan: Approved plan.
        run_manifest: Prediction manifest, or an empty mapping.
        prediction: Prediction chrom, start, end, and bin size.
        track_rows: Typed query, lane label, and model head for each predicted
            track.
        saliency_tracks: Loaded saliency letter tracks.

    Returns:
        JSON-ready explain summary.
    """
    chrom, start, end, bin_size = prediction
    window_bp = int(plan["end"]) - int(plan["start"])
    resolved = {
        str(row.get("query")): row
        for row in (_read_json(root / "track_resolve.json") or {}).get("tracks") or []
        if isinstance(row, dict)
    }
    # The ENCODE output type each ground-truth BigWig was saved with. Lookup
    # writes the accession as the label when ENCODE gave none.
    gt_units = {
        str(row["accession"]): str(row.get("label") or "")
        for row in _read_json(root / "gt" / "manifest.json") or []
        if isinstance(row, dict)
        and row.get("accession")
        and row.get("label")
        and row.get("label") != row.get("accession")
    }
    explain_tracks = []
    for row in track_rows:
        accession = str((resolved.get(row["query"]) or {}).get("gt_accession") or "")
        explain_tracks.append(
            {
                "label": row["label"],
                "modelHead": row["modelHead"],
                "gtAccession": accession,
                "gtUnit": gt_units.get(accession, "signal") if accession else "",
            }
        )
    variant = _read_json(root / "disrupt" / "variant.json")
    compare = None
    if isinstance(variant, dict):
        try:
            compare = {
                "label": (
                    f"{variant['chrom']}:{int(variant['position'])} "
                    f"{variant['ref']}>{variant['alt']}"
                ),
                "checkpoint": str(variant.get("checkpoint") or plan["model"]),
                "window": _span(variant["window"]),
                "crop": _span(variant["prediction"]),
                "binSize": int(variant["prediction"]["bin_size"]),
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise ArtifactError("disrupt/variant.json provenance is invalid.") from exc
    sources = [
        rel
        for rel in (
            "plan.yaml",
            "track_resolve.json",
            "tracks/metadata.json",
            "tracks/accessions.json",
            "run_manifest.json",
            "disrupt/variant.json",
        )
        if (root / rel).is_file()
    ]
    return {
        "model": str(plan["model"]),
        "manifestModel": run_manifest.get("model"),
        "species": str(plan["species"]),
        "assembly": str(plan["assembly"]),
        # The device predict ran on; manifests from before it was recorded
        # hold the approved value, which was then the same on CPU.
        "device": str(run_manifest.get("device") or plan["device"]),
        "deviceRequested": str(plan["device"]),
        "window": {
            "chrom": str(plan["chrom"]),
            "start": int(plan["start"]),
            "end": int(plan["end"]),
            "bp": window_bp,
        },
        "crop": {"chrom": chrom, "start": start, "end": end, "bp": end - start},
        "binSize": bin_size,
        "nBins": (end - start) // bin_size,
        "tracks": explain_tracks,
        "letters": [
            {
                "label": str(track["name"]),
                "modelHead": str(track["modelHead"]),
                "start": int(track["start"]),
                "end": int(track["end"]),
            }
            for track in saliency_tracks
        ],
        "compare": compare,
        "sources": sources,
    }


def _span(payload: Any) -> dict[str, Any]:
    """Return a chrom/start/end/bp span from a provenance mapping."""
    start = int(payload["start"])
    end = int(payload["end"])
    return {
        "chrom": str(payload["chrom"]),
        "start": start,
        "end": end,
        "bp": end - start,
    }


def _resolve_requested_stages(
    existing_report: Any,
    requested_stages: list[str] | None,
) -> set[str]:
    """Resolve and validate the downstream stages selected in setup QA."""
    if requested_stages is not None:
        unknown = sorted(set(requested_stages) - set(REPORT_STAGE_LABELS))
        if unknown:
            raise RequestError(
                f"Unknown report stages: {unknown}. "
                f"Allowed: {sorted(REPORT_STAGE_LABELS)}."
            )
        return set(requested_stages)
    if isinstance(existing_report, dict):
        stages = existing_report.get("analysisStages")
        if isinstance(stages, list):
            return {
                str(stage["id"])
                for stage in stages
                if isinstance(stage, dict)
                and stage.get("id") in REPORT_STAGE_LABELS
                and stage.get("status") != "not_requested"
            }
    return set(DEFAULT_REPORT_STAGES)


def _analysis_stages(
    requested_stages: set[str],
    completed: dict[str, bool],
) -> list[dict[str, str]]:
    """Build deterministic display states for all downstream analyses."""
    stages: list[dict[str, str]] = []
    for stage_id, label in REPORT_STAGE_LABELS.items():
        if completed.get(stage_id, False):
            status = "complete"
            message = "Completed. Results are available below."
        elif stage_id in requested_stages:
            status = "pending"
            message = "Pending until this analysis is run."
        else:
            status = "not_requested"
            message = "Not selected, or explicitly declined, in the agent workflow."
        stages.append(
            {
                "id": stage_id,
                "label": label,
                "status": status,
                "message": message,
            }
        )
    return stages


def _lane_source(value: Any) -> dict[str, str] | None:
    """Return the gene lane's ``{track, url}`` source, or None.

    The page takes only that mapping or null. ``locus.json`` holds a text
    source (the UCSC search and transcript used), so it gives null.
    """
    if (
        isinstance(value, dict)
        and isinstance(value.get("track"), str)
        and isinstance(value.get("url"), str)
    ):
        return {"track": value["track"], "url": value["url"]}
    return None


def _read_json(path: Path) -> Any:
    """Read JSON if the file exists.

    Args:
        path: File path.

    Returns:
        Parsed object or None.
    """
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _read_yaml(path: Path) -> dict[str, Any] | None:
    """Read YAML if the file exists.

    Args:
        path: File path.

    Returns:
        Mapping or None.
    """
    if not path.exists():
        return None
    payload = yaml.safe_load(path.read_text())
    return payload if isinstance(payload, dict) else None


def _require_mapping(path: Path) -> dict[str, Any]:
    """Read a required JSON or YAML mapping.

    Args:
        path: Required artifact path.

    Returns:
        Parsed mapping.

    Raises:
        ArtifactError: The artifact is absent or is not a mapping.
    """
    if not path.exists():
        raise ArtifactError(f"Required artifact is missing: {path}.")
    if path.suffix == ".yaml":
        payload = _read_yaml(path)
    else:
        payload = _read_json(path)
    if not isinstance(payload, dict):
        raise ArtifactError(f"Required artifact is not a mapping: {path}.")
    return payload


def _require_string_list(path: Path) -> list[str]:
    """Read a required JSON list of strings."""
    payload = _read_json(path)
    if not isinstance(payload, list) or not all(
        isinstance(value, str) for value in payload
    ):
        raise ArtifactError(f"Required artifact is not a string list: {path}.")
    return payload


def _locus_coords(plan: dict[str, Any]) -> tuple[str, int, int]:
    """Read validated chrom/start/end from the approved plan."""
    try:
        chrom = str(plan["chrom"])
        start = int(plan["start"])
        end = int(plan["end"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ArtifactError(
            "plan.yaml must contain valid chrom, start, and end values."
        ) from exc
    if not chrom.startswith("chr") or start < 0 or end <= start:
        raise ArtifactError(f"Invalid plan locus: {chrom}:{start}-{end}.")
    return chrom, start, end


def lane_label(query: str, model_head: str) -> str:
    """Return a track's lane label: the words the user typed, then the id.

    An id typed as the track (``ENCSR...``/``CNhs...``) is its own label. No
    friendlier name is made up.

    Args:
        query: The track as the user named it in ``--tracks``.
        model_head: The model head it resolved to, or "" if none.

    Returns:
        ``"<query> <model_head>"``, or the query when the two are the same.
    """
    if not model_head or query == model_head:
        return query
    return f"{query} {model_head}"


def _assay_names(
    plan: dict[str, Any] | None,
    resolve: dict[str, Any] | None,
) -> list[str]:
    """Collect assay labels.

    Args:
        plan: Parsed plan.yaml.
        resolve: Parsed track_resolve.json.

    Returns:
        Assay names.
    """
    if plan and plan.get("tracks"):
        return list(plan["tracks"])
    if resolve and resolve.get("tracks"):
        return [row.get("query") or row.get("accession") for row in resolve["tracks"]]
    return []


def _load_track_arrays(
    tracks_npy: Path,
    observed_npy: Path,
    n_assays: int,
    require_observed: bool,
) -> tuple[list[list[float]], list[list[float | None]]]:
    """Load prediction and aligned experimental track arrays.

    Args:
        tracks_npy: Path to ``tracks/tracks.npy``.
        observed_npy: Path to ``gt/observed.npy``.
        n_assays: Number of assays.
        require_observed: Whether GT sources make the observed array mandatory.

    Returns:
        (predicted, observed) lists of lists.
    """
    arr = np.load(tracks_npy)
    if arr.ndim != 2 or arr.shape[0] != n_assays or not np.isfinite(arr).all():
        raise ArtifactError(
            f"Prediction array must be finite with shape (n_assays, n_bins); "
            f"got {arr.shape} for {n_assays} assays."
        )
    predicted = [_rounded(row.astype(float).tolist()) for row in arr]
    if not observed_npy.exists():
        if require_observed:
            raise ArtifactError(
                "GT sources exist but gt/observed.npy is missing. "
                "Run prediction with the science environment."
            )
        observed: list[list[float | None]] = [[] for _ in predicted]
        return predicted, observed
    observed_array = np.load(observed_npy)
    if observed_array.shape != arr.shape or np.isinf(observed_array).any():
        raise ArtifactError(
            f"Observed array must match prediction shape {arr.shape} and contain "
            f"only finite values or NaN; got {observed_array.shape}."
        )
    observed = [
        _rounded([float(value) if np.isfinite(value) else None for value in row])
        for row in observed_array
    ]
    return predicted, observed


def _track_assay(
    name: str,
    model_head: str,
    predicted: list[float],
    observed: list[float | None],
) -> dict[str, Any]:
    """Build one report assay from validated prediction artifacts."""
    finite = [value for value in predicted if np.isfinite(value)]
    value_range = [min(finite), max(finite)] if finite else [0.0, 0.0]
    return {
        "name": name,
        "modelHead": model_head,
        "unit": "model logit",
        "range": value_range,
        "source": f"NTv3 prediction · {model_head}",
        "observed": observed,
        "predicted": predicted,
    }


def _load_track_agreement(path: Path, required: bool) -> list[dict[str, Any]]:
    """Load the prediction-vs-ground-truth agreement table for the viewer.

    Args:
        path: Path to ``tables/track_agreement.csv``.
        required: Whether GT sources make the comparison mandatory.

    Returns:
        JSON-ready agreement rows.
    """
    if not path.exists():
        if required:
            raise ArtifactError(
                "GT sources exist but tables/track_agreement.csv is missing. "
                "Run prediction with the science environment."
            )
        return []
    frame = pd.read_csv(path)
    required_columns = {
        "track_query",
        "model_head",
        "gt_accession",
        "file_accession",
        "n_valid",
        "pearson",
        "spearman",
        "caveat",
    }
    missing = sorted(required_columns - set(frame.columns))
    if missing:
        raise ArtifactError(f"Track agreement table is missing columns: {missing}.")
    rows: list[dict[str, Any]] = []
    for row in frame.to_dict(orient="records"):
        rows.append(
            {
                "trackQuery": str(row["track_query"]),
                "modelHead": str(row["model_head"]),
                "gtAccession": str(row["gt_accession"]),
                "fileAccession": str(row["file_accession"]),
                "nValid": int(row["n_valid"]),
                "pearson": _finite_or_none(row["pearson"]),
                "spearman": _finite_or_none(row["spearman"]),
                "caveat": str(row["caveat"]),
            }
        )
    return rows


def _finite_or_none(value: Any) -> float | None:
    """Return a finite float or JSON null."""
    number = float(value)
    return number if np.isfinite(number) else None


def _optional_percentile(value: Any) -> float | None:
    """Return a percentile in [0, 100], or None when no null ran."""
    if value is None:
        return None
    number = float(value)
    if not 0.0 <= number <= 100.0:
        raise ArtifactError(f"Percentile {number} is outside [0, 100].")
    return number


def _load_disruption(root: Path) -> dict[str, Any] | None:
    """Load one optional single-SNV Compare summary and its alt-minus-ref tracks."""
    effect_path = root / "disrupt" / "effect_percentile.json"
    variant_path = root / "disrupt" / "variant.json"
    png_path = root / "figures" / "disrupt_log2fc.png"
    pdf_path = root / "figures" / "disrupt_log2fc.pdf"
    if not effect_path.exists() and not variant_path.exists():
        return None
    if not all(
        path.is_file() for path in (effect_path, variant_path, png_path, pdf_path)
    ):
        raise ArtifactError("Disruption artifact bundle is incomplete.")
    effect = _require_mapping(effect_path)
    variant = _require_mapping(variant_path)
    try:
        chrom = str(variant["chrom"])
        position = int(variant["position"])
        ref = str(variant["ref"])
        alt = str(variant["alt"])
        signed_effect = float(effect["signed_effect"])
        absolute_effect = float(effect["absolute_effect"])
        percentile = _optional_percentile(effect["percentile"])
        null_n = int(effect["null_n"])
        seed = int(effect["seed"])
        prediction = _span(variant["prediction"])
        bin_size = int(variant["prediction"]["bin_size"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ArtifactError("Disruption summary is invalid or incomplete.") from exc
    if (
        len(ref) != 1
        or len(alt) != 1
        or not all(np.isfinite(value) for value in (signed_effect, absolute_effect))
        or null_n < 0
        or (percentile is None) != (null_n == 0)
    ):
        raise ArtifactError("Disruption summary contains invalid values.")
    per_track = effect.get("per_track")
    if not isinstance(per_track, list):
        raise ArtifactError("Disruption summary requires per_track results.")
    delta_path = root / "disrupt" / "delta_tracks.npy"
    delta = np.load(delta_path) if delta_path.is_file() else None
    if delta is not None and (
        delta.ndim != 2
        or delta.shape[0] != len(per_track)
        or delta.shape[1] * bin_size != prediction["bp"]
        or not np.isfinite(delta).all()
    ):
        raise ArtifactError(
            f"disrupt/delta_tracks.npy must be finite with shape "
            f"({len(per_track)}, {prediction['bp'] // max(bin_size, 1)}); got "
            f"{delta.shape}."
        )
    return {
        "variant": {
            "chrom": chrom,
            "position": position,
            "ref": ref,
            "alt": alt,
            "label": f"{chrom}:{position} {ref}>{alt}",
        },
        "signedEffect": signed_effect,
        "absoluteEffect": absolute_effect,
        "percentile": percentile,
        "nullN": null_n,
        "seed": seed,
        "prediction": {**prediction, "binSize": bin_size},
        "perTrack": [
            {
                "track": lane_label(
                    query=str(row["track"]), model_head=str(row["model_head"])
                ),
                "modelHead": str(row["model_head"]),
                "signedEffect": float(row["signed_effect"]),
                "percentile": _optional_percentile(row.get("percentile")),
                "delta": None if delta is None else _rounded(delta[index].tolist()),
            }
            for index, row in enumerate(per_track)
            if isinstance(row, dict)
        ],
        "method": str((effect.get("method") or {}).get("track_reduction") or ""),
        "figure": {
            "png": "figures/disrupt_log2fc.png",
            "pdf": "figures/disrupt_log2fc.pdf",
        },
        "caveat": (
            "Predicted model output, not experimental evidence"
            + (
                "; ranked against a seeded empirical null."
                if null_n
                else "; no null was run, so no percentile."
            )
        ),
    }


def _load_literature(root: Path) -> dict[str, Any] | None:
    """Load optional validated primary-source citation records."""
    manifest_path = root / "literature" / "search_manifest.json"
    bib_path = root / "literature" / "citations.bib"
    if not manifest_path.exists() and not bib_path.exists():
        return None
    if not manifest_path.is_file() or not bib_path.is_file():
        raise ArtifactError("Literature artifact bundle is incomplete.")
    manifest = _require_mapping(manifest_path)
    records = manifest.get("records")
    if not isinstance(records, list):
        raise ArtifactError("Literature search manifest requires records.")
    return {
        "records": [
            {
                "citekey": str(row["citekey"]),
                "title": str(row["title"]),
                "authors": [str(author) for author in row["authors"]],
                "year": int(row["year"]),
                "sourceUrl": str(row["source_url"]),
                "doi": str(row["doi"]) if row.get("doi") else None,
                "retrievalDate": str(row["retrieval_date"]),
                "searchQuery": str(row["search_query"]),
                "relevance": str(row["relevance"]),
            }
            for row in records
            if isinstance(row, dict)
        ],
        "bib": "literature/citations.bib",
    }


def _load_saliency(root: Path) -> list[dict[str, Any]]:
    """Load per-track saliency maps written by ``ntv3 predict --saliency``.

    Args:
        root: ``interpret/`` directory.

    Returns:
        One row per ``*.json`` sidecar next to a ``.npy`` map.
    """
    if not root.exists():
        return []
    rows: list[dict[str, Any]] = []
    for meta_path in sorted(root.glob("*.json")):
        if meta_path.name == "motif_match.json":
            continue
        npy_path = meta_path.with_suffix(".npy")
        if not npy_path.exists():
            continue
        meta = _require_mapping(meta_path)
        scores = np.load(npy_path)
        if scores.ndim != 2 or scores.shape[1] != 4 or not np.isfinite(scores).all():
            raise ArtifactError(
                f"Saliency array {npy_path} must be finite with shape "
                f"(sequence_length, 4); got {scores.shape}."
            )
        try:
            model_head = str(meta["accession"])
            start = int(meta["start"])
            end = int(meta["end"])
            sequence = str(meta["seq"])
            logo = meta["logo"]
            if not isinstance(logo, dict):
                raise TypeError("logo must be a mapping")
            logo_selection = str(logo["selection"])
            logo_center_bp = int(logo["center_bp"])
            logo_start = int(logo["start"])
            logo_end = int(logo["end"])
            logo_png = str(logo["png"])
            logo_pdf = str(logo["pdf"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ArtifactError(
                f"Saliency metadata {meta_path} requires accession, start, end, "
                "seq, and complete logo metadata."
            ) from exc
        if end - start != scores.shape[0] or len(sequence) != scores.shape[0]:
            raise ArtifactError(
                f"Saliency metadata {meta_path} does not match array shape "
                f"{scores.shape}."
            )
        if (
            logo_selection not in {"strongest_peak", "explicit_center"}
            or not start <= logo_start < logo_end <= end
            or not logo_start <= logo_center_bp < logo_end
        ):
            raise ArtifactError(f"Saliency logo interval in {meta_path} is invalid.")
        run_root = root.parent.resolve()
        for relative_path in (logo_png, logo_pdf):
            figure_path = (run_root / relative_path).resolve()
            if run_root not in figure_path.parents or not figure_path.is_file():
                raise ArtifactError(
                    f"Saliency logo artifact is missing or outside the run: "
                    f"{relative_path}."
                )
        # Saliency is kept at the reference base, so the sum is that base's value.
        projected = _rounded(scores.sum(axis=1).astype(float).tolist())
        rows.append(
            {
                "name": model_head,
                "modelHead": model_head,
                "values": projected,
                "start": start,
                "end": end,
                "seq": sequence,
                "logo": {
                    "selection": logo_selection,
                    "centerBp": logo_center_bp,
                    "start": logo_start,
                    "end": logo_end,
                    "png": logo_png,
                    "pdf": logo_pdf,
                },
            }
        )
    return rows


def _ground_truth_sources(
    gt_rows: list[Any],
    resolve: dict[str, Any],
    model_heads_by_query: dict[str, str],
) -> list[dict[str, str]]:
    """Attach each ENCODE source to its requested assay and model head."""
    resolved_rows = resolve.get("tracks") or []
    by_accession = {
        str(row.get("gt_accession")): row
        for row in resolved_rows
        if isinstance(row, dict) and row.get("gt_accession")
    }
    sources: list[dict[str, str]] = []
    for row in gt_rows:
        if not isinstance(row, dict) or not row.get("url"):
            continue
        accession = str(row.get("accession") or "")
        resolved = by_accession.get(accession, {})
        query = str(resolved.get("query") or row.get("label") or accession)
        model_head = str(
            resolved.get("model_head") or model_heads_by_query.get(query, "")
        )
        sources.append(
            {
                "name": lane_label(query=query, model_head=model_head),
                "accession": accession,
                "fileAccession": str(row.get("encff") or ""),
                "modelHead": model_head,
                "url": str(row["url"]),
                "origin": str(row.get("origin") or ""),
            }
        )
    return sources


def _run_claims(root: Path, plan: dict[str, Any]) -> list[Claim]:
    """Claims the run makes, one per approved analysis, with its evidence files.

    A claim whose artifacts are missing is suppressed by the critic, so an
    approved analysis that did not write its outputs cannot be reported.
    """
    heads = _read_json(root / "tracks" / "accessions.json") or []
    queries = list(plan.get("tracks") or [])
    claims = [
        Claim(
            claim_id="predict",
            text=(
                f"{plan['model']} predicted {', '.join(queries)} over the crop in "
                "tracks/metadata.json."
            ),
            evidence_paths=[
                "tracks/tracks.npy",
                "tracks/metadata.json",
                "tracks/accessions.json",
                "run_manifest.json",
            ],
            citekey=None,
            kind="model",
        )
    ]
    if "saliency" in (plan.get("analyses") or []):
        chosen = set(plan.get("attribution_tracks") or queries)
        for query, head in zip(queries, heads, strict=False):
            if query in chosen:
                claims.append(
                    Claim(
                        claim_id=f"saliency-{head}",
                        text=f"Saliency for {query} (head {head}) over the window.",
                        evidence_paths=[
                            f"interpret/{head}.npy",
                            f"interpret/{head}.json",
                        ],
                        citekey=None,
                        kind="model",
                    )
                )
    if plan.get("variant"):
        variant = plan["variant"]
        claims.append(
            Claim(
                claim_id="compare",
                text=(
                    f"Compare scored {variant['chrom']}:{variant['position']} "
                    f"{variant['ref']}>{variant['alt']} against the reference."
                ),
                evidence_paths=[
                    "disrupt/effect_percentile.json",
                    "disrupt/delta_tracks.npy",
                ],
                citekey=None,
                kind="model",
            )
        )
    return claims + _annotation_claims(root=root)


def _annotation_claims(root: Path) -> list[Claim]:
    """The claim an approved annotation makes once its GFF3 is in the run."""
    plan = _read_json(root / ANNOTATION_PLAN)
    if not plan or not (root / plan["output"]).is_file():
        return []
    genes = len(read_gff3(root / plan["output"]).genes)
    contigs = ", ".join(str(row["seqid"]) for row in plan["records"])
    return [
        Claim(
            claim_id="annotate",
            text=(
                f"{plan['model']} predicted {genes} genes with the "
                f"{plan['decoder']} decoder on {contigs}."
            ),
            evidence_paths=[plan["output"], ANNOTATION_MANIFEST, ANNOTATION_PLAN],
            citekey=None,
            kind="model",
        )
    ]


def _run_critic(
    root: Path, enabled: bool, claims: list[Claim], log_path: Path
) -> dict[str, Any]:
    """Run the critic when it is on, and say what happened either way.

    The agent's ``claims.json`` is audited when present; otherwise ``claims``,
    the claims the run makes from its outputs. An existing log is kept unless
    the run now makes a claim it has not audited, such as a locus annotation
    added to a Predict run.

    Returns:
        JSON-ready critic status for report.json.
    """
    if not enabled:
        return {
            "enabled": False,
            "note": "The critic was off for this run; it was turned off at approval.",
        }
    claims_path = root / "claims.json"
    if claims_path.is_file():
        claims = [
            Claim(
                claim_id=row["claim_id"],
                text=row["text"],
                evidence_paths=row.get("evidence_paths") or [],
                citekey=row.get("citekey"),
                kind=row.get("kind") or "other",
            )
            for row in json.loads(claims_path.read_text())
        ]
    audited = {row.get("claim_id") for row in _read_json(log_path) or []}
    unaudited = {claim.claim_id for claim in claims} - audited
    if not log_path.exists() or (audited and unaudited and not claims_path.is_file()):
        bib = root / "literature" / "citations.bib"
        write_audit_log(
            decisions=audit_claims(
                claims=claims,
                run_dir=str(root),
                bib_path=str(bib) if bib.exists() else None,
            ),
            out_dir=str(root),
        )
    rows = _read_json(log_path) or []
    return {
        "enabled": True,
        "log": "critic/audit_log.json",
        "claims": "claims.json"
        if claims_path.is_file()
        else "the run's approved analyses",
        "kept": sum(row.get("verdict") == "KEEP" for row in rows),
        "suppressed": sum(row.get("verdict") == "SUPPRESS" for row in rows),
    }


def _claims_from_critic(critic: list[Any]) -> list[dict[str, Any]]:
    """Map critic audit rows to report claims.

    Args:
        critic: Audit log rows.

    Returns:
        Claim objects for the viewer.
    """
    claims: list[dict[str, Any]] = []
    for row in critic:
        claims.append(
            {
                "id": row.get("claim_id", "claim"),
                "text": row.get("text", ""),
                "support": row.get("reason", ""),
                "command": None,
                "keep": row.get("verdict") == "KEEP",
            }
        )
    return claims


def _trace_from_tree(
    root: Path,
    completed_steps: set[str] | None = None,
) -> list[dict[str, Any]]:
    """List present artifacts as trace steps.

    Steps that run only on explicit request are listed only once their artifact
    exists, so a default run never shows them as skipped.

    Args:
        root: Run directory.
        completed_steps: Steps completed by the current operation before its
            artifacts are written.

    Returns:
        Trace rows.
    """
    checks = [
        ("plan", "plan.yaml"),
        ("lookup", "track_resolve.json"),
        ("predict", "tracks/tracks.npy"),
        ("disrupt", "disrupt/effect_percentile.json"),
        ("saliency", "interpret"),
        ("literature", "literature/search_manifest.json"),
        ("critic", "critic/audit_log.json"),
        ("replay", "code/run_analysis.py"),
        ("notebook", "analysis.ipynb"),
        ("report", "report.json"),
        ("ro-crate", "ro-crate-metadata.json"),
    ]
    completed = completed_steps or set()
    rows: list[dict[str, Any]] = []
    for name, rel in checks:
        path = root / rel
        is_complete = path.exists() or name in completed
        if not is_complete and name in _EXPLICIT_ONLY_STEPS:
            continue
        rows.append(
            {
                "name": name,
                "artifact": rel if is_complete else None,
                "status": "ok" if is_complete else "skipped",
            }
        )
    return rows
