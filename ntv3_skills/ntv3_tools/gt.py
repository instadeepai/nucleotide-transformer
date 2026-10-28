"""Ground-truth BigWig lookup on ENCODE. Do not invent ENCFF ids."""

import json
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import requests

from ntv3_tools.errors import ExternalServiceError, RequestError

GetJson = Callable[..., dict[str, Any]]

_ENCODE_SEARCH = "https://www.encodeproject.org/search/"
_ENCODE_BASE = "https://www.encodeproject.org"
_ASSEMBLY_ON_ENCODE = {"hg38": "GRCh38", "mm10": "mm10"}
_ENCODE_ORGANISM = {"human": "Homo sapiens", "mouse": "Mus musculus"}
_OUTPUT_PREF = (
    "fold change over control",
    "signal p-value",
    "plus strand signal of unique reads",
    "minus strand signal of unique reads",
    "signal of unique reads",
    "signal",
    "read-depth normalized signal",
    "RPKM",
)


@dataclass(frozen=True)
class GtSource:
    """One ground-truth BigWig.

    Attributes:
        accession: ENCSR experiment id.
        url: HTTP or s3 URI the browser can open.
        origin: "encode-public" or "in-house-s3".
        label: Display name.
        encff: Chosen file accession, if known.
    """

    accession: str
    url: str
    origin: str
    label: str
    encff: str = ""


def default_get_json(url: str, params: dict[str, str] | None = None) -> dict[str, Any]:
    """GET JSON from ENCODE.

    Args:
        url: Absolute URL.
        params: Query string.

    Returns:
        Parsed object. ENCODE answers an empty search with HTTP 404 and a JSON
        body that has ``@graph``; that body is returned, not raised.
    """
    try:
        response = requests.get(
            url,
            params=params,
            timeout=60,
            headers={"Accept": "application/json", "User-Agent": "ntv3-skills"},
        )
        if response.status_code == 404:
            try:
                body = response.json()
            except ValueError:
                body = None
            if isinstance(body, dict) and "@graph" in body:
                return body
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as exc:
        status = getattr(getattr(exc, "response", None), "status_code", "network")
        raise ExternalServiceError(
            f"ENCODE request failed ({status}) for {url}."
        ) from exc


def pick_encode_bigwig(
    files: list[dict[str, Any]],
    assembly: str,
) -> dict[str, Any]:
    """Choose one released BigWig from an ENCODE file list.

    Args:
        files: ENCODE File objects.
        assembly: Genome assembly (hg38 / mm10).

    Returns:
        One File object.

    Raises:
        KeyError: The list is empty.
    """
    if not files:
        raise KeyError("No ENCODE BigWig files to pick from.")
    encode_asm = _ASSEMBLY_ON_ENCODE.get(assembly, assembly)
    matching = [row for row in files if row.get("assembly") in {encode_asm, assembly}]
    if not matching:
        available = sorted(
            {str(row.get("assembly")) for row in files if row.get("assembly")}
        )
        raise KeyError(
            f"No released BigWig matches assembly '{assembly}'. "
            f"Available assemblies: {available}."
        )

    def rank(row: dict[str, Any]) -> tuple[int, str]:
        output_type = str(row.get("output_type") or "").lower()
        for index, preferred in enumerate(_OUTPUT_PREF):
            if preferred.lower() == output_type:
                return (index, str(row.get("accession") or ""))
        return (len(_OUTPUT_PREF), str(row.get("accession") or ""))

    return min(matching, key=rank)


def _file_url(row: dict[str, Any]) -> str:
    """Cloud or portal URL for one ENCODE file."""
    href = row.get("cloud_metadata", {}).get("url") or row.get("href")
    if href and str(href).startswith("/"):
        href = "https://www.encodeproject.org" + str(href)
    if not href:
        raise KeyError("ENCODE file hit had no URL.")
    return str(href)


def search_encode_files(
    accession: str,
    get_json: GetJson = default_get_json,
) -> list[dict[str, Any]]:
    """List usable BigWigs for one ENCSR experiment.

    Args:
        accession: ENCSR experiment id.
        get_json: HTTP GET (inject in tests).

    Returns:
        File objects.
    """
    try:
        payload = get_json(
            _ENCODE_SEARCH,
            {
                "type": "File",
                "dataset": f"/experiments/{accession}/",
                "file_format": "bigWig",
                "status": "released",
                "format": "json",
                "limit": "all",
            },
        )
    except ExternalServiceError:
        payload = {}
    rows = list(payload.get("@graph") or [])
    if rows:
        return rows

    experiment = get_json(
        f"{_ENCODE_BASE}/experiments/{accession}/",
        {"format": "json"},
    )
    released: list[dict[str, Any]] = []
    archived: list[dict[str, Any]] = []
    for file_ref in experiment.get("files") or []:
        if isinstance(file_ref, dict):
            row = file_ref
        else:
            file_path = str(file_ref)
            if not file_path.startswith("/"):
                file_path = f"/files/{file_path.strip('/')}/"
            row = get_json(
                f"{_ENCODE_BASE}{file_path}",
                {"format": "json"},
            )
        if str(row.get("file_format") or "").lower() != "bigwig":
            continue
        if row.get("status") == "released":
            released.append(row)
        elif row.get("status") == "archived":
            archived.append(row)
    return released or archived


def search_encode_experiment(
    tissue: str,
    assay: str,
    get_json: GetJson = default_get_json,
    species: str | None = None,
) -> str:
    """Find one released ENCSR for a tissue + assay. Ask if several match.

    Tries the exact ``assay_title`` first, then a ``searchTerm`` query on the
    assay words (filtered by organism when the species is known).

    Args:
        tissue: Biosample term (e.g. K562).
        assay: Assay words (e.g. DNase or CAGE).
        get_json: HTTP GET (inject in tests).
        species: Model species ("human" / "mouse"), if known.

    Returns:
        One ENCSR id.

    Raises:
        RequestError: Zero or many experiments.
    """
    assay_clean = assay.replace(tissue, "").strip() or assay
    base = {
        "type": "Experiment",
        "status": "released",
        "biosample_ontology.term_name": tissue,
    }
    exact = {**base, "assay_title": assay_clean, "format": "json", "limit": "10"}
    broad = {**base, "searchTerm": assay_clean, "format": "json", "limit": "10"}
    organism = _ENCODE_ORGANISM.get((species or "").lower())
    if organism:
        broad["replicates.library.biosample.donor.organism.scientific_name"] = organism
    for params in (exact, broad):
        try:
            payload = get_json(_ENCODE_SEARCH, params)
        except ExternalServiceError as exc:
            raise RequestError(
                f"ENCODE could not be reached for {assay_clean} in {tissue}. Which "
                f"ENCSR id should be predicted for {assay_clean} in {tissue}? Pass "
                "it in --tracks. No id is assumed."
            ) from exc
        hits = list(payload.get("@graph") or [])
        ids = [str(row.get("accession")) for row in hits if row.get("accession")]
        if len(ids) == 1:
            return ids[0]
        if ids:
            total = payload.get("total")
            shown = (
                f" {total} released; showing {len(ids)}."
                if isinstance(total, int) and total > len(ids)
                else ""
            )
            raise RequestError(
                f"ENCODE has several released {assay_clean} experiments for "
                f"{tissue}: {', '.join(ids)}.{shown} Which one should be used? "
                "Pass its ENCSR id in --tracks. None is picked for you."
            )
    raise RequestError(
        f"ENCODE has no released {assay_clean} experiment for {tissue}. Which ENCSR "
        f"id should be predicted for {assay_clean} in {tissue}? Pass it in --tracks. "
        "No id is assumed."
    )


def resolve_gt(
    assembly: str,
    accession: str,
    get_json: GetJson = default_get_json,
) -> GtSource:
    """Resolve one GT BigWig from the public ENCODE portal.

    Args:
        assembly: Genome assembly.
        accession: ENCSR experiment id. Not a CNhs model head.
        get_json: HTTP GET (inject in tests).

    Returns:
        Source the browser / CLI can fetch.
    """
    if accession.startswith("CNhs"):
        raise KeyError(
            f"{accession} is a FANTOM / model-head id, not an ENCODE experiment. "
            "Resolve the GT ENCSR separately. Do not send CNhs* to ENCODE."
        )
    if not accession.upper().startswith("ENCSR"):
        raise KeyError(f"GT lookup needs an ENCSR experiment, got {accession!r}.")
    files = search_encode_files(accession=accession, get_json=get_json)
    if not files:
        raise KeyError(
            f"ENCODE has no released BigWig for {accession}. Which ENCSR id should "
            "ground truth come from? Pass it in --tracks. No file is assumed."
        )
    chosen = pick_encode_bigwig(files=files, assembly=assembly)
    if not chosen.get("cloud_metadata") and chosen.get("accession"):
        detail = get_json(
            f"{_ENCODE_BASE}/files/{chosen['accession']}/",
            {"format": "json"},
        )
        if isinstance(detail, dict) and detail.get("cloud_metadata"):
            chosen = {**chosen, **detail}
    return GtSource(
        accession=accession,
        url=_file_url(chosen),
        origin=(
            "encode-public-archived"
            if chosen.get("status") == "archived"
            else "encode-public"
        ),
        label=str(chosen.get("output_type") or accession),
        encff=str(chosen.get("accession") or ""),
    )


def write_gt_manifest(
    assembly: str,
    accessions: list[str],
    out_dir: str,
    get_json: GetJson = default_get_json,
) -> Path:
    """Write gt/manifest.json with URLs (no S3 writes).

    Args:
        assembly: Genome assembly.
        accessions: ENCSR ids from track_resolve.json (gt_accession).
        out_dir: Run output directory.
        get_json: HTTP GET (inject in tests).

    Returns:
        Path to the manifest.
    """
    sources = [
        resolve_gt(assembly=assembly, accession=accession, get_json=get_json)
        for accession in accessions
        if accession
    ]
    return write_gt_sources(sources=sources, out_dir=out_dir)


def write_gt_sources(sources: list[GtSource], out_dir: str) -> Path:
    """Write already-resolved public ENCODE sources.

    Args:
        sources: Validated ground-truth sources.
        out_dir: Run output directory.

    Returns:
        Path to the manifest.
    """
    root = Path(out_dir) / "gt"
    root.mkdir(parents=True, exist_ok=True)
    path = root / "manifest.json"
    path.write_text(json.dumps([asdict(source) for source in sources], indent=2) + "\n")
    return path
