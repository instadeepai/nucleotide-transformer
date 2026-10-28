"""Validated locus values and UCSC gene-symbol resolution."""

import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import requests

from ntv3_tools.errors import ExternalServiceError, RequestError
from ntv3_tools.registry import assembly_for_species
from ntv3_tools.ucsc import fetch_track

GetJson = Callable[..., dict[str, Any]]
FetchTrack = Callable[..., dict[str, Any]]
_UCSC_SEARCH_URL = "https://api.genome.ucsc.edu/search"
_POSITION_PATTERN = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
_TRACK_PRIORITY = ("mane", "hgnc", "knownGene", "refGene")
# The UCSC gene-model track that places the TSS, and the row field holding the
# gene symbol. mm10 knownGene rows carry no symbol; ncbiRefSeq has it in name2.
_GENE_MODEL = {"hg38": ("knownGene", "geneName"), "mm10": ("ncbiRefSeq", "name2")}


@dataclass(frozen=True)
class LocusRequest:
    """A genomic window request for the NTv3 model.

    Attributes:
        chrom: Chromosome name (e.g. "chr11").
        start: 0-based half-open window start.
        end: 0-based half-open window end.
        species: Species string accepted by the model (e.g. "human").
    """

    chrom: str
    start: int
    end: int
    species: str

    @property
    def assembly(self) -> str:
        """Genome assembly for this locus (e.g. "hg38")."""
        return assembly_for_species(self.species)

    @property
    def length_bp(self) -> int:
        """Window length in base pairs."""
        return self.end - self.start

    def validate(self) -> None:
        """Validate the genomic window.

        Raises:
            RequestError: The coordinates are invalid.
        """
        if not self.chrom.startswith("chr"):
            raise RequestError(
                f"Chromosome '{self.chrom}' must use a UCSC chr-prefixed name."
            )
        if self.start < 0 or self.end <= self.start:
            raise RequestError(
                f"Invalid 0-based half-open coordinates: "
                f"{self.chrom}:{self.start}-{self.end}."
            )


@dataclass(frozen=True)
class GeneLocus:
    """A gene match and the model window centred on its TSS.

    Attributes:
        gene: Exact gene symbol requested by the user.
        assembly: UCSC genome assembly.
        chrom: Chromosome for the exact match.
        gene_start: Zero-based gene start.
        gene_end: Half-open gene end.
        start: Zero-based model-window start.
        end: Half-open model-window end.
        source: UCSC search track and gene-model transcript used.
        tss: Zero-based transcription start site the window is centred on.
        strand: Strand of that transcript, "+" or "-".
        transcript: Transcript id whose TSS was used.
    """

    gene: str
    assembly: str
    chrom: str
    gene_start: int
    gene_end: int
    start: int
    end: int
    source: str
    tss: int
    strand: str
    transcript: str


def default_get_json(url: str, params: dict[str, str] | None = None) -> dict[str, Any]:
    """Get JSON from the UCSC API."""
    try:
        response = requests.get(
            url,
            params=params,
            timeout=60,
            headers={"Accept": "application/json", "User-Agent": "ntv3-skills"},
        )
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as exc:
        status = getattr(getattr(exc, "response", None), "status_code", "network")
        raise ExternalServiceError(
            f"UCSC gene search failed ({status}) for {url}."
        ) from exc


def resolve_gene_locus(
    gene: str,
    assembly: str,
    window_bp: int,
    get_json: GetJson = default_get_json,
    fetch_gene_rows: FetchTrack = fetch_track,
) -> GeneLocus:
    """Resolve an exact gene symbol and centre the window on its TSS.

    The UCSC search API finds the locus. Matches that overlap are one locus;
    two that do not, on any chromosome, are a question. The TSS comes from one
    transcript of that gene in the assembly's UCSC gene-model track: the
    canonical one on hg38 (knownGene). mm10 ncbiRefSeq rows carry no canonical
    tag, so there it is the shortest curated NM_ transcript, never an XM_ model
    when an NM_ exists.

    Args:
        gene: Gene symbol.
        assembly: UCSC assembly.
        window_bp: Length of the model window.
        get_json: Injectable UCSC JSON request.
        fetch_gene_rows: Injectable UCSC track fetch for the gene model.

    Returns:
        Exact gene match and the window centred on its TSS.

    Raises:
        RequestError: No exact match, several distinct loci, or no transcript
            of that name to place a TSS.
    """
    symbol = gene.strip()
    if not symbol:
        raise RequestError(
            "Which gene or locus? Pass --gene, or --chrom, --start, and --end. "
            "No locus is assumed."
        )
    if window_bp <= 0:
        raise RequestError(f"window_bp must be positive, got {window_bp}.")

    payload = get_json(
        _UCSC_SEARCH_URL,
        {"search": symbol, "genome": assembly},
    )
    matches: list[tuple[int, str, str, int, int]] = []
    for group in payload.get("positionMatches") or []:
        track = str(group.get("trackName") or group.get("name") or "")
        priority = (
            _TRACK_PRIORITY.index(track)
            if track in _TRACK_PRIORITY
            else len(_TRACK_PRIORITY)
        )
        for row in group.get("matches") or []:
            if not _is_exact_symbol(symbol=symbol, row=row):
                continue
            parsed = _parse_position(position=str(row.get("position") or ""))
            if parsed is None:
                continue
            chrom, gene_start, gene_end = parsed
            matches.append((priority, track, chrom, gene_start, gene_end))

    if not matches:
        raise RequestError(
            f"UCSC found no exact genomic match for '{symbol}' on {assembly}. "
            "Which locus should be used? Pass --chrom, --start, and --end. "
            "No locus is assumed."
        )
    best_priority = min(priority for priority, _, _, _, _ in matches)
    preferred = [match for match in matches if match[0] == best_priority]
    loci = _overlap_clusters(
        spans=[(chrom, start, end) for _, _, chrom, start, end in preferred]
    )
    if len(loci) > 1:
        choices = [f"{chrom}:{start}-{end}" for chrom, start, end in loci]
        raise RequestError(
            f"UCSC returned several exact loci for '{symbol}': {choices}. "
            "Which one should be used? Pass its --chrom, --start, and --end. "
            "None is picked for you."
        )
    chrom, gene_start, gene_end = loci[0]
    search_track = preferred[0][1]

    if assembly not in _GENE_MODEL:
        raise RequestError(
            f"No UCSC gene model is set for {assembly} to place a TSS. Which locus "
            "should be used? Pass --chrom, --start, and --end. No locus is assumed."
        )
    track, name_field = _GENE_MODEL[assembly]
    payload = fetch_gene_rows(
        genome=assembly, track=track, chrom=chrom, start=gene_start, end=gene_end
    )
    rows = [
        row
        for row in _track_rows(payload=payload, track=track)
        if str(row.get("chrom") or chrom) == chrom
        and str(row.get(name_field) or "").strip().casefold() == symbol.casefold()
        and _valid_known_gene_row(row=row)
    ]
    if not rows:
        raise RequestError(
            f"UCSC {track} has no transcript named '{symbol}' on {assembly} to place "
            "a TSS. Which locus should be used? Pass --chrom, --start, and --end. "
            "No locus is assumed."
        )
    selected = min(rows, key=_known_gene_rank)
    strand = str(selected["strand"])
    transcript_start, transcript_end = _transcript_coords(row=selected)
    tss = transcript_start if strand == "+" else transcript_end - 1
    transcript = str(selected.get("name") or "")
    start = max(0, tss - window_bp // 2)
    end = start + window_bp
    return GeneLocus(
        gene=symbol,
        assembly=assembly,
        chrom=chrom,
        gene_start=gene_start,
        gene_end=gene_end,
        start=start,
        end=end,
        source=f"UCSC /search ({search_track}) and {track} ({transcript})",
        tss=tss,
        strand=strand,
        transcript=transcript,
    )


def _overlap_clusters(spans: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge overlapping spans; each result is one distinct locus."""
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(spans):
        if merged and merged[-1][0] == chrom and start < merged[-1][2]:
            last_chrom, last_start, last_end = merged[-1]
            merged[-1] = (last_chrom, last_start, max(last_end, end))
        else:
            merged.append((chrom, start, end))
    return merged


def _track_rows(payload: dict[str, Any], track: str) -> list[dict[str, Any]]:
    """Return a UCSC getData payload's rows for one track, as a flat list."""
    rows = payload.get(track) or []
    if isinstance(rows, dict):
        rows = [row for values in rows.values() for row in values]
    return [row for row in rows if isinstance(row, dict)]


def _is_exact_symbol(symbol: str, row: dict[str, Any]) -> bool:
    """Return whether a UCSC search row names the requested symbol exactly."""
    target = symbol.casefold()
    pos_name = str(row.get("posName") or "").strip()
    description = str(row.get("description") or "").strip()
    first_name = pos_name.split(" ", 1)[0].split("(", 1)[0].strip()
    return first_name.casefold() == target or description.casefold() == target


def _parse_position(position: str) -> tuple[str, int, int] | None:
    """Parse a UCSC display position into zero-based half-open coordinates."""
    match = _POSITION_PATTERN.match(position.replace(",", ""))
    if match is None:
        return None
    chrom, start_text, end_text = match.groups()
    start = int(start_text) - 1
    end = int(end_text)
    if start < 0 or end <= start:
        return None
    return chrom, start, end


def known_gene_metadata(
    gene: str,
    chrom: str,
    start: int,
    end: int,
    payload: dict[str, Any],
) -> dict[str, Any]:
    """Derive sourced gene context from saved UCSC knownGene rows.

    Args:
        gene: Requested gene symbol. It can be absent for coordinate runs.
        chrom: Requested chromosome.
        start: Zero-based locus start.
        end: Half-open locus end.
        payload: Saved UCSC knownGene response.

    Returns:
        JSON-ready locus metadata with TSS, strand, named elements, and source.
    """
    rows = [
        row
        for row in payload.get("knownGene", [])
        if isinstance(row, dict)
        and str(row.get("chrom") or chrom) == chrom
        and _valid_known_gene_row(row=row)
    ]
    target = gene.strip().casefold()
    exact = [
        row
        for row in rows
        if str(row.get("geneName") or "").strip().casefold() == target
    ]
    selected = min(exact, key=_known_gene_rank) if exact else None
    named_genes = sorted(
        {
            str(row.get("geneName")).strip()
            for row in rows
            if str(row.get("geneName") or "").strip()
            and _transcript_coords(row=row)[0] < end
            and _transcript_coords(row=row)[1] > start
        }
    )
    exact_name = str(selected.get("geneName")).strip() if selected is not None else None
    elements = [
        f"{name} · {'gene' if name == exact_name else 'nearby gene'}"
        for name in named_genes
    ]
    strand = str(selected["strand"]) if selected is not None else None
    tss = None
    transcript = None
    if selected is not None:
        transcript = str(selected.get("name") or "")
        transcript_start, transcript_end = _transcript_coords(row=selected)
        tss = transcript_start if strand == "+" else transcript_end - 1
    return {
        "gene": gene,
        "chrom": chrom,
        "start": start,
        "end": end,
        "tss": tss,
        "strand": strand,
        "transcript": transcript,
        "elements": elements,
        "source": {
            "track": str(payload.get("_label") or "knownGene"),
            "url": str(payload.get("_request_url") or ""),
        },
    }


def _valid_known_gene_row(row: dict[str, Any]) -> bool:
    """Return whether a knownGene row has valid transcript coordinates."""
    try:
        transcript_start, transcript_end = _transcript_coords(row=row)
    except (KeyError, TypeError, ValueError):
        return False
    return (
        transcript_start >= 0
        and transcript_end > transcript_start
        and str(row.get("strand") or "") in {"+", "-"}
    )


def _transcript_coords(row: dict[str, Any]) -> tuple[int, int]:
    """Return UCSC transcript coordinates across knownGene response versions."""
    start = row.get("chromStart", row.get("txStart"))
    end = row.get("chromEnd", row.get("txEnd"))
    return int(start), int(end)


def _known_gene_rank(row: dict[str, Any]) -> tuple[int, int, int, int, int, str]:
    """Rank canonical tags, explicit rank, RefSeq curation, then shortest."""
    tag_text = " ".join(
        str(row.get(key) or "") for key in ("tag", "tags", "transcriptClass", "status")
    ).casefold()
    is_preferred = bool(row.get("canonical")) or any(
        marker in tag_text for marker in ("canonical", "gencode primary", "mane")
    )
    try:
        rank = int(row.get("rank", row.get("transcriptRank", 2**31 - 1)))
    except (TypeError, ValueError):
        rank = 2**31 - 1
    # RefSeq: curated NM_ before NR_, then predicted XM_/XR_ models. knownGene
    # ENST ids all share one tier, so hg38 ranking is unchanged.
    name = str(row.get("name") or "")
    curation = next(
        (tier for tier, prefix in enumerate(("NM_", "NR_")) if name.startswith(prefix)),
        2,
    )
    transcript_start, transcript_end = _transcript_coords(row=row)
    return (
        0 if is_preferred else 1,
        rank,
        curation,
        transcript_end - transcript_start,
        transcript_start,
        str(row.get("name") or ""),
    )
