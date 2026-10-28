"""Validate and package agent-collected primary-source literature records."""

import json
import re
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from ntv3_tools.errors import RequestError

_DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[-._;()/:A-Z0-9]+$", re.IGNORECASE)
_CITEKEY_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_:-]*$")
_PLACEHOLDER_PATTERN = re.compile(
    r"(?:^(?:fake|placeholder|todo|tbd|unknown|xxx)|"
    r"(?:^|[-_.:/])(?:example|fake|placeholder|test|todo|tbd|unknown|xxx)"
    r"(?:$|[-_.:/]))",
    re.IGNORECASE,
)
_PRIMARY_HOSTS = {
    "academic.oup.com",
    "arxiv.org",
    "biorxiv.org",
    "cell.com",
    "doi.org",
    "elifesciences.org",
    "frontiersin.org",
    "journals.plos.org",
    "link.springer.com",
    "medrxiv.org",
    "nature.com",
    "ncbi.nlm.nih.gov",
    "onlinelibrary.wiley.com",
    "pnas.org",
    "pubmed.ncbi.nlm.nih.gov",
    "science.org",
    "sciencedirect.com",
    "www.biorxiv.org",
    "www.cell.com",
    "www.frontiersin.org",
    "www.nature.com",
    "www.ncbi.nlm.nih.gov",
    "www.pnas.org",
    "www.science.org",
    "www.sciencedirect.com",
}
_SEARCH_HOSTS = {
    "google.com",
    "scholar.google.com",
    "www.google.com",
}
_REQUIRED_FIELDS = (
    "citekey",
    "title",
    "authors",
    "year",
    "source_url",
    "retrieval_date",
    "search_query",
    "relevance",
)


@dataclass(frozen=True)
class LiteratureRecord:
    """One validated literature search result.

    Attributes:
        citekey: Stable BibTeX citation key.
        title: Primary-source title.
        authors: Ordered author names.
        year: Publication year.
        source_url: Direct HTTP URL to the primary source.
        retrieval_date: ISO retrieval date.
        search_query: Query used by the agent.
        relevance: Claim scope supported by the source.
        doi: Optional DOI without a URL prefix.
    """

    citekey: str
    title: str
    authors: list[str]
    year: int
    source_url: str
    retrieval_date: str
    search_query: str
    relevance: str
    doi: str | None = None


def package_literature(
    run_dir: str | Path,
    records_path: str | Path,
) -> tuple[Path, Path]:
    """Validate records and write deterministic literature artifacts.

    Args:
        run_dir: Approved run directory containing ``plan.yaml``.
        records_path: JSON records file located under the same run directory.

    Returns:
        Paths to ``citations.bib`` and ``search_manifest.json``.

    Raises:
        RequestError: An input is invalid, unsafe, or would overwrite output.
    """
    root = Path(run_dir).resolve()
    source = Path(records_path).resolve()
    _validate_paths(run_dir=root, records_path=source)
    records = _read_records(records_path=source)
    literature_dir = root / "literature"
    bib_path = literature_dir / "citations.bib"
    manifest_path = literature_dir / "search_manifest.json"
    existing = [path for path in (bib_path, manifest_path) if path.exists()]
    if existing:
        raise RequestError(f"Refusing to overwrite literature artifact: {existing[0]}.")

    bib_text = _render_bibtex(records=records)
    manifest_text = _render_manifest(records=records, records_path=source, root=root)
    literature_dir.mkdir(parents=True, exist_ok=True)
    try:
        with bib_path.open(mode="x", encoding="utf-8") as bib_file:
            bib_file.write(bib_text)
        with manifest_path.open(mode="x", encoding="utf-8") as manifest_file:
            manifest_file.write(manifest_text)
    except FileExistsError as exc:
        raise RequestError(
            f"Refusing to overwrite literature artifact: {exc.filename}."
        ) from exc
    return bib_path, manifest_path


def _validate_paths(run_dir: Path, records_path: Path) -> None:
    """Validate the approved run and confined input path."""
    if not run_dir.is_dir() or not (run_dir / "plan.yaml").is_file():
        raise RequestError(
            f"Literature packaging requires an approved run with plan.yaml: {run_dir}."
        )
    if not records_path.is_file():
        raise RequestError(f"Literature records file does not exist: {records_path}.")
    if run_dir not in records_path.parents:
        raise RequestError(
            f"Literature records must be located under the approved run: {run_dir}."
        )


def _read_records(records_path: Path) -> list[LiteratureRecord]:
    """Read and validate a JSON list of literature records."""
    try:
        payload = json.loads(records_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise RequestError(
            f"Could not read literature records: {records_path}."
        ) from exc
    if not isinstance(payload, list) or not payload:
        raise RequestError("Literature records JSON must be a non-empty list.")

    records: list[LiteratureRecord] = []
    citekeys: set[str] = set()
    for index, raw_record in enumerate(payload):
        record = _validate_record(raw_record=raw_record, index=index)
        if record.citekey in citekeys:
            raise RequestError(f"Duplicate literature citekey: {record.citekey}.")
        citekeys.add(record.citekey)
        records.append(record)
    return records


def _validate_record(raw_record: Any, index: int) -> LiteratureRecord:
    """Validate and normalize one raw JSON record."""
    if not isinstance(raw_record, dict):
        raise RequestError(f"Literature record {index} must be a JSON object.")
    missing = [field for field in _REQUIRED_FIELDS if field not in raw_record]
    if missing:
        raise RequestError(f"Literature record {index} is missing fields: {missing}.")

    values = {
        field: _require_text(value=raw_record[field], field=field, index=index)
        for field in _REQUIRED_FIELDS
        if field not in {"authors", "year"}
    }
    citekey = values["citekey"]
    if not _CITEKEY_PATTERN.fullmatch(citekey) or _looks_placeholder(citekey):
        raise RequestError(f"Literature record {index} has an invalid citekey.")

    authors = _validate_authors(value=raw_record["authors"], index=index)
    year = raw_record["year"]
    if isinstance(year, bool) or not isinstance(year, int) or not 1500 <= year <= 2100:
        raise RequestError(f"Literature record {index} has an invalid year.")

    source_url = values["source_url"]
    _validate_source_url(source_url=source_url, index=index)
    retrieval_date = values["retrieval_date"]
    try:
        parsed_date = date.fromisoformat(retrieval_date)
    except ValueError as exc:
        raise RequestError(
            f"Literature record {index} retrieval_date must be an ISO date."
        ) from exc
    if parsed_date.isoformat() != retrieval_date:
        raise RequestError(
            f"Literature record {index} retrieval_date must be an ISO date."
        )

    raw_doi = raw_record.get("doi")
    doi = None
    if raw_doi is not None:
        doi = _require_text(value=raw_doi, field="doi", index=index)
        if not _DOI_PATTERN.fullmatch(doi) or _looks_placeholder(doi):
            raise RequestError(f"Literature record {index} has an invalid DOI.")
    if doi is None and not _is_primary_source_url(source_url=source_url):
        raise RequestError(
            f"Literature record {index} needs a DOI or recognizable primary-source URL."
        )

    return LiteratureRecord(
        citekey=citekey,
        title=values["title"],
        authors=authors,
        year=year,
        source_url=source_url,
        retrieval_date=retrieval_date,
        search_query=values["search_query"],
        relevance=values["relevance"],
        doi=doi,
    )


def _require_text(value: Any, field: str, index: int) -> str:
    """Return a required non-blank string."""
    if not isinstance(value, str) or not value.strip():
        raise RequestError(
            f"Literature record {index} field '{field}' must be non-blank text."
        )
    return value.strip()


def _validate_authors(value: Any, index: int) -> list[str]:
    """Return a validated non-empty author list."""
    if isinstance(value, str):
        authors = [author.strip() for author in value.split(" and ")]
    elif isinstance(value, list):
        authors = [
            _require_text(value=author, field="authors", index=index)
            for author in value
        ]
    else:
        raise RequestError(
            f"Literature record {index} field 'authors' must be text or a list."
        )
    if not authors or any(not author for author in authors):
        raise RequestError(f"Literature record {index} has invalid authors.")
    return authors


def _validate_source_url(source_url: str, index: int) -> None:
    """Validate one direct HTTP source URL."""
    parsed = urlparse(source_url)
    hostname = (parsed.hostname or "").lower()
    if (
        parsed.scheme not in {"http", "https"}
        or not hostname
        or hostname in _SEARCH_HOSTS
        or _looks_placeholder(source_url)
        or "/search" in parsed.path.lower()
    ):
        raise RequestError(f"Literature record {index} has an invalid source URL.")


def _is_primary_source_url(source_url: str) -> bool:
    """Return whether a URL points to a recognizable primary-source host."""
    parsed = urlparse(source_url)
    hostname = (parsed.hostname or "").lower()
    if hostname not in _PRIMARY_HOSTS:
        return False
    path = parsed.path.lower()
    if hostname == "doi.org":
        return bool(_DOI_PATTERN.fullmatch(path.lstrip("/")))
    if hostname == "pubmed.ncbi.nlm.nih.gov":
        return bool(re.fullmatch(r"/\d+/?", path))
    return path not in {"", "/"}


def _looks_placeholder(value: str) -> bool:
    """Return whether a citation identifier or URL looks fabricated."""
    return bool(_PLACEHOLDER_PATTERN.search(value))


def _render_bibtex(records: list[LiteratureRecord]) -> str:
    """Render deterministic minimal BibTeX accepted by the critic."""
    entries: list[str] = []
    for record in records:
        fields = [
            ("title", record.title),
            ("author", " and ".join(record.authors)),
            ("year", str(record.year)),
            ("url", record.source_url),
            ("urldate", record.retrieval_date),
        ]
        if record.doi is not None:
            fields.append(("doi", record.doi))
        lines = [f"@article{{{record.citekey},"]
        lines.extend(
            f"  {name} = {{{_escape_bibtex(value=value)}}}," for name, value in fields
        )
        lines.append("}")
        entries.append("\n".join(lines))
    return "\n\n".join(entries) + "\n"


def _escape_bibtex(value: str) -> str:
    """Escape BibTeX special characters without changing Unicode text."""
    escapes = {
        "\\": r"\textbackslash{}",
        "{": r"\{",
        "}": r"\}",
        "#": r"\#",
        "%": r"\%",
        "&": r"\&",
        "_": r"\_",
        "^": r"\^{}",
        "~": r"\~{}",
    }
    return "".join(escapes.get(character, character) for character in value).replace(
        "\n", " "
    )


def _render_manifest(
    records: list[LiteratureRecord],
    records_path: Path,
    root: Path,
) -> str:
    """Render the deterministic structured search manifest."""
    payload = {
        "schema_version": 1,
        "input_records": str(records_path.relative_to(root)),
        "records": [asdict(record) for record in records],
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


__all__ = ["LiteratureRecord", "package_literature"]
