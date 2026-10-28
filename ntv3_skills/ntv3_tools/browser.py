"""Fill the static browser page with one run.

``browser_template.html`` is a small page with inline CSS and JS and one empty
data tag. ``write_browser_html`` puts a run's lanes in that tag and writes
``results/<run>/browser.html``, an older static page; the page published as
the artifact is ``ntv3_report.html`` (see ``report_app``). With no run, the
same template is the empty state. The page fetches nothing.
"""

import json
import math
from pathlib import Path
from typing import Any

from ntv3_tools.errors import ArtifactError

TEMPLATE_PATH = Path(__file__).with_name("browser_template.html")
BROWSER_HTML = "browser.html"
RUN_DATA_TAG = '<script type="application/json" id="ntv3-run">null</script>'

# Lane label for each UCSC track= id the skill fetches, per assembly.
CATALOG_LABELS: dict[str, dict[str, str]] = {
    "hg38": {
        "knownGene": "GENCODE",
        "clinvarMain": "ClinVar",
        "cCREregistry": "ENCODE cCREs",
        "wgEncodeRegDnaseClustered": "DNase clusters",
        "phastCons100way": "phastCons",
        "wgEncodeRegMarkH3k27acK562": "H3K27ac K562",
        "wgEncodeRegMarkH3k4me3K562": "H3K4me3 K562",
        "wgEncodeRegMarkH3k27acGm12878": "H3K27ac GM12878",
        "wgEncodeBroadHistoneGm12878H3k4me3StdSig": "H3K4me3 GM12878",
        "wgEncodeRegMarkH3k27acHuvec": "H3K27ac HUVEC",
        "wgEncodeRegMarkH3k4me3Huvec": "H3K4me3 HUVEC",
        "wgEncodeRegMarkH3k27acH1hesc": "H3K27ac H1-hESC",
        "wgEncodeRegMarkH3k4me3H1hesc": "H3K4me3 H1-hESC",
        "wgEncodeRegMarkH3k27acHsmm": "H3K27ac HSMM",
        "wgEncodeRegMarkH3k4me3Hsmm": "H3K4me3 HSMM",
        "wgEncodeRegMarkH3k27acNhek": "H3K27ac NHEK",
        "wgEncodeRegMarkH3k4me3Nhek": "H3K4me3 NHEK",
        "wgEncodeRegMarkH3k27acNhlf": "H3K27ac NHLF",
        "wgEncodeRegMarkH3k4me3Nhlf": "H3K4me3 NHLF",
    },
    "mm10": {"knownGene": "UCSC Genes", "ncbiRefSeq": "NCBI RefSeq"},
}
# UCSC tracks drawn as gene models: knownGene (hg38 bigGenePred, geneName) and
# ncbiRefSeq (mm10 genePred, name2). A row shows a symbol only if it has one.
GENE_TRACKS = frozenset({"knownGene", "ncbiRefSeq"})
NO_PREDICTION_NOTE = "Predicted tracks appear when a run has them."
NO_ANNOTATION_NOTE = "An NTv3 annotation lane appears when the run has a GFF3."


def json_for_script(data: Any) -> str:
    """Serialize data for an inline JSON script tag that it cannot close."""
    return (
        json.dumps(data, separators=(",", ":"))
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def render_browser_html(run_data: dict[str, Any] | None) -> str:
    """Return the browser page, filled with one run or empty.

    Args:
        run_data: Output of ``browser_data``, or None for the empty state.

    Returns:
        The HTML document.

    Raises:
        ArtifactError: The template does not hold exactly one empty data tag.
    """
    template = TEMPLATE_PATH.read_text()
    if template.count(RUN_DATA_TAG) != 1:
        raise ArtifactError(
            f"{TEMPLATE_PATH.name} must hold exactly one empty ntv3-run data tag."
        )
    payload = json_for_script(run_data)
    return template.replace(
        RUN_DATA_TAG, RUN_DATA_TAG.replace(">null<", f">{payload}<")
    )


def round_sig(value: float) -> float:
    """Round to four significant digits to keep the page small."""
    if not value or not math.isfinite(value):
        return 0.0
    return round(value, 3 - math.floor(math.log10(abs(value))))


def _rgb(value: Any) -> str | None:
    """Return a UCSC itemRgb (``r,g,b`` or a packed integer) as CSS, or None."""
    text = str(value or "").strip()
    parts = text.split(",")
    if len(parts) == 3 and all(part.strip().isdigit() for part in parts):
        return f"rgb({','.join(part.strip() for part in parts)})"
    if text.isdigit() and int(text) > 0:
        packed = int(text)
        return f"rgb({packed >> 16 & 255},{packed >> 8 & 255},{packed & 255})"
    return None


def _gene_item(row: dict[str, Any]) -> dict[str, Any]:
    """Turn one bigGenePred (hg38) or genePred (mm10) row into a gene with exon blocks."""
    start = int(row.get("chromStart", row.get("txStart")))
    if "exonStarts" in row:
        starts = [int(x) for x in str(row["exonStarts"]).split(",") if x]
        ends = [int(x) for x in str(row["exonEnds"]).split(",") if x]
        blocks = [list(pair) for pair in zip(starts, ends, strict=False)]
    else:
        sizes = [int(x) for x in str(row.get("blockSizes") or "").split(",") if x]
        offsets = [int(x) for x in str(row.get("chromStarts") or "").split(",") if x]
        blocks = [
            [start + offset, start + offset + size]
            for offset, size in zip(offsets, sizes, strict=False)
        ]
    return {
        "s": start,
        "e": int(row.get("chromEnd", row.get("txEnd"))),
        # The gene symbol when the row carries one (bigGenePred geneName,
        # genePred name2 or geneSymbol), else the transcript id.
        "n": str(
            row.get("geneName")
            or row.get("name2")
            or row.get("geneSymbol")
            or row.get("name")
            or ""
        ),
        "d": str(row.get("strand") or ""),
        "b": blocks,
    }


def _ucsc_lane(track: str, label: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Turn one saved UCSC getData payload into a browser lane."""
    rows = [row for row in payload.get(track) or [] if isinstance(row, dict)]
    lane: dict[str, Any] = {"label": label, "source": f"UCSC {track}"}
    if rows and "value" in rows[0]:
        items = [
            [int(row["start"]), int(row["end"]), round_sig(float(row["value"]))]
            for row in rows
        ]
        values = [item[2] for item in items]
        return {
            **lane,
            "kind": "signal",
            "items": items,
            "min": min(values),
            "max": max(values),
        }
    if rows and track in GENE_TRACKS:
        return {**lane, "kind": "genes", "items": [_gene_item(row) for row in rows]}
    items = []
    for row in rows:
        item = {
            "s": int(row["chromStart"]),
            "e": int(row["chromEnd"]),
            "n": str(row.get("name") or ""),
        }
        color = _rgb(row.get("reserved") or row.get("itemRgb"))
        if color:
            item["c"] = color
        items.append(item)
    return {**lane, "kind": "features", "items": items}


def _signal_lane(
    label: str, source: str, start: int, bin_size: int, values: list[Any], color: str
) -> dict[str, Any] | None:
    """Turn evenly binned values from ``start`` into a signal lane."""
    finite = [
        (index, float(value))
        for index, value in enumerate(values)
        if value is not None and math.isfinite(float(value))
    ]
    if not finite:
        return None
    return {
        "label": label,
        "source": source,
        "kind": "signal",
        "color": color,
        "items": [
            [start + index * bin_size, start + (index + 1) * bin_size, round_sig(value)]
            for index, value in finite
        ],
        "min": min(value for _, value in finite),
        "max": max(value for _, value in finite),
    }


def _prediction_lanes(tracks: dict[str, Any]) -> list[dict[str, Any]]:
    """Return a predicted lane and, when ground truth was saved, a measured lane."""
    lanes: list[dict[str, Any]] = []
    start, end = int(tracks.get("start") or 0), int(tracks.get("end") or 0)
    for assay in tracks.get("assays") or []:
        for kind, color, label in (
            ("predicted", "signal", f"Predicted {assay['name']}"),
            ("observed", "measured", f"Measured {assay['name']}"),
        ):
            values = list(assay.get(kind) or [])
            if not values:
                continue
            bin_size = max(1, (end - start) // len(values))
            source = (
                f"NTv3 head {assay['modelHead']}" if kind == "predicted" else "ENCODE"
            )
            lane = _signal_lane(label, source, start, bin_size, values, color)
            if lane is not None:
                lanes.append(lane)
    return lanes


def _saliency_lanes(report: dict[str, Any]) -> list[dict[str, Any]]:
    """Return one per-base saliency lane per track: centred gradient times input."""
    lanes = []
    for track in report.get("saliencyTracks") or []:
        lane = _signal_lane(
            label=f"Saliency {track['name']}",
            source=f"centred gradient x input, head {track['modelHead']}",
            start=int(track["start"]),
            bin_size=1,
            values=list(track.get("values") or []),
            color="saliency",
        )
        if lane is not None:
            lanes.append(lane)
    return lanes


def _compare_lanes(disruption: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Return one alt-minus-ref lane per track Compare scored."""
    if not isinstance(disruption, dict):
        return []
    prediction = disruption["prediction"]
    lanes = []
    for row in disruption.get("perTrack") or []:
        lane = _signal_lane(
            label=f"Alt − ref {row['track']}",
            source=f"{disruption['variant']['label']}, head {row['modelHead']}",
            start=int(prediction["start"]),
            bin_size=int(prediction["binSize"]),
            values=list(row.get("delta") or []),
            color="delta",
        )
        if lane is not None:
            lanes.append(lane)
    return lanes


def _annotation_lane(annotation: dict[str, Any], chrom: str) -> dict[str, Any]:
    """Return the NTv3 annotation lane: genes with their CDS blocks."""
    return {
        "label": str(annotation.get("label") or "NTv3 annotation"),
        "source": str(annotation.get("gff3") or "GFF3"),
        "kind": "genes",
        "color": "annot",
        "items": [
            {
                "s": int(gene["start"]),
                "e": int(gene["end"]),
                "n": str(gene["name"]),
                "d": str(gene["strand"]),
                "b": [
                    [int(e["start"]), int(e["end"])] for e in gene.get("exons") or []
                ],
            }
            for gene in annotation.get("genes") or []
            if str(gene.get("chrom")) == chrom
        ],
        "empty": "No predicted genes on this chromosome.",
    }


def _clinvar_id(variant: dict[str, Any], ucsc: dict[str, Any]) -> str | None:
    """Return the saved ClinVar id whose own genomic HGVS is exactly this SNV.

    Matches ``g.<position><ref>><alt>`` in the record's HGVS table, so the
    position and both alleles must agree; nothing is inferred or fetched.
    """
    hgvs = f"g.{int(variant['position'])}{variant['ref']}>{variant['alt']}"
    for row in (ucsc.get("clinvarMain") or {}).get("clinvarMain") or []:
        if not isinstance(row, dict) or str(row.get("chrom")) != variant["chrom"]:
            continue
        try:
            table = json.loads(str(row.get("_jsonHgvsTable") or "[]"))
        except json.JSONDecodeError:
            continue
        if any(isinstance(pair, list) and pair and pair[0] == hgvs for pair in table):
            return str(row.get("name") or "") or None
    return None


def _explain_rows(
    report: dict[str, Any], clinvar: str | None = None
) -> list[list[str]]:
    """Return Explain rows (label, value) read from report.json provenance."""
    rows: list[list[str]] = []
    explain = report.get("explain")
    if isinstance(explain, dict):
        window, crop = explain["window"], explain["crop"]
        rows += [
            ["Model", str(explain["model"])],
            ["Species", f"{explain['species']} ({explain['assembly']})"],
            ["Device", str(explain["device"])],
            [
                "Window",
                (
                    f"{window['chrom']}:{window['start']:,}-{window['end']:,} "
                    f"({window['bp']:,} bp)"
                ),
            ],
            [
                "Crop",
                (
                    f"{crop['chrom']}:{crop['start']:,}-{crop['end']:,} "
                    f"({crop['bp']:,} bp, central share predicted)"
                ),
            ],
            ["Bins", f"{explain['nBins']:,} of {explain['binSize']:,} bp"],
        ]
        for track in explain.get("tracks") or []:
            rows.append(
                [
                    str(track["label"]),
                    (
                        f"head {track['modelHead']}, ground truth "
                        f"{track.get('gtAccession') or 'none'}"
                    ),
                ]
            )
        if explain.get("compare"):
            compare = explain["compare"]
            record = f" (ClinVar {clinvar})" if clinvar else ""
            rows.append(
                ["Compare", f"{compare['label']}{record} with {compare['checkpoint']}"]
            )
        rows.append(["Sources", ", ".join(explain.get("sources") or [])])
    annotation = report.get("annotation")
    if isinstance(annotation, dict):
        rows += [
            ["Annotation model", str(annotation.get("model") or "unrecorded")],
            ["Clade", str(annotation.get("clade") or "unrecorded")],
            ["Decoder", str(annotation.get("decoder") or "unrecorded")],
            [
                "Annotation window",
                (
                    f"{annotation.get('windowBp') or 'unrecorded'} bp, central "
                    f"{annotation.get('keepFraction') or 'unrecorded'} kept"
                ),
            ],
            [
                "Filters",
                f"CDS shorter than {annotation.get('minCdsLength') or '?'} bp "
                "dropped; internal in-frame stop "
                + ("dropped" if annotation.get("inFrameStopFilter") else "kept"),
            ],
            ["Annotation sources", ", ".join(annotation.get("sources") or [])],
        ]
    critic = report.get("critic")
    if isinstance(critic, dict):
        rows.append(
            [
                "Critic",
                (
                    f"ran on {critic.get('claims')}: {critic.get('kept')} kept, "
                    f"{critic.get('suppressed')} suppressed ({critic.get('log')})"
                    if critic.get("enabled")
                    else str(critic.get("note") or "off")
                ),
            ]
        )
    return rows


def _compare(report: dict[str, Any]) -> dict[str, Any] | None:
    """Return the Compare table, largest change first, when the run scored an SNV."""
    disruption = report.get("disruption")
    if not isinstance(disruption, dict):
        return None
    ranked = f"{disruption['percentile']:.1f}" if disruption["nullN"] else None
    return {
        "summary": (
            f"{disruption['variant']['label']}: largest alt − ref change "
            f"{disruption['signedEffect']:.3g}"
            + (
                f", percentile {ranked} against {disruption['nullN']} seeded null "
                f"variants (seed {disruption['seed']})"
                if ranked
                else ""
            )
            + f". {disruption['caveat']}"
        ),
        "rows": [
            [
                row["track"],
                row["modelHead"],
                f"{row['signedEffect']:.3g}",
                "—" if row["percentile"] is None else f"{row['percentile']:.1f}",
            ]
            for row in sorted(
                disruption.get("perTrack") or [],
                key=lambda row: -abs(row["signedEffect"]),
            )
        ],
    }


def browser_data(report: dict[str, Any], ucsc: dict[str, Any]) -> dict[str, Any]:
    """Build the page data for one run from report.json and saved UCSC payloads.

    Args:
        report: Parsed ``report.json``.
        ucsc: Saved UCSC payloads keyed by ``track=`` id.

    Returns:
        JSON-ready page data: title, view, lanes, variant marks, notes,
        Explain, and Compare.
    """
    locus = report["locus"]
    chrom = str(locus["chrom"])
    start, end = int(locus["start"]), int(locus["end"])
    assembly = str(report.get("assembly") or "")
    labels = CATALOG_LABELS.get(assembly, {})
    lanes = [
        _ucsc_lane(track=track, label=labels.get(track, track), payload=ucsc[track])
        for track in report.get("ucscTracks") or []
        if track in ucsc
    ]
    predicted = _prediction_lanes(report.get("tracks") or {})
    disruption = report.get("disruption")
    lanes += predicted + _saliency_lanes(report) + _compare_lanes(disruption)
    notes: list[str] = [] if predicted else [NO_PREDICTION_NOTE]
    measured = any(lane.get("color") == "measured" for lane in predicted)
    if report.get("gtSources") and not measured:
        notes.append("Ground truth was not compared in this run.")
    view = [start, end]
    annotation = report.get("annotation")
    if isinstance(annotation, dict):
        lane = _annotation_lane(annotation=annotation, chrom=chrom)
        lanes.append(lane)
        if not predicted and lane["items"]:
            first = min(item["s"] for item in lane["items"])
            last = max(item["e"] for item in lane["items"])
            pad = max(500, (last - first) // 20)
            view = [max(start, first - pad), min(end, last + pad)]
    else:
        notes.append(NO_ANNOTATION_NOTE)
    focus = None
    if predicted:
        focus = [int(locus["focusStart"]), int(locus["focusEnd"])]
    subtitle = " · ".join(
        part
        for part in (
            f"{chrom}:{start:,}-{end:,}",
            assembly,
            str(report.get("cell") or ""),
            str(report.get("modelShort") or ""),
        )
        if part
    )
    if focus:
        subtitle += " · shaded: the span the model predicts"
    marks = []
    clinvar = None
    if isinstance(disruption, dict):
        variant = disruption["variant"]
        clinvar = _clinvar_id(variant=variant, ucsc=ucsc)
        label = variant["label"] + (f" (ClinVar {clinvar})" if clinvar else "")
        marks.append({"pos": int(variant["position"]) - 1, "label": label})
    return {
        "title": f"{report.get('gene') or 'NTv3'} · NTv3 browser",
        "subtitle": subtitle,
        "chrom": chrom,
        "bounds": [start, end],
        "view": {"start": view[0], "end": view[1]},
        "focus": focus,
        "marks": marks,
        "lanes": lanes,
        "notes": notes,
        "explain": _explain_rows(report, clinvar=clinvar),
        "compare": _compare(report),
    }


def saved_ucsc(root: Path) -> dict[str, Any]:
    """Return every UCSC payload lookup saved, keyed by track."""
    ucsc_dir = root / "ucsc"
    if not ucsc_dir.is_dir():
        return {}
    return {
        path.stem: json.loads(path.read_text())
        for path in sorted(ucsc_dir.glob("*.json"))
    }


def write_browser_html(run_dir: str | Path) -> Path:
    """Write ``browser.html`` from ``report.json`` and the saved UCSC payloads.

    Args:
        run_dir: Run directory holding ``report.json``.

    Returns:
        Path to ``browser.html`` (not the published page; see ``report_app``).

    Raises:
        ArtifactError: ``report.json`` is missing.
    """
    root = Path(run_dir)
    report_path = root / "report.json"
    if not report_path.is_file():
        raise ArtifactError(f"{BROWSER_HTML} needs {report_path}. Run report first.")
    data = browser_data(
        report=json.loads(report_path.read_text()), ucsc=saved_ucsc(root=root)
    )
    path = root / BROWSER_HTML
    path.write_text(render_browser_html(run_data=data))
    return path
