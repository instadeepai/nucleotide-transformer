"""One lookup pass: dual-id resolve, live ENCODE GT, UCSC overlays."""

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from ntv3_tools.errors import RequestError
from ntv3_tools.gt import (
    GetJson,
    default_get_json,
    resolve_gt,
    search_encode_experiment,
    write_gt_sources,
)
from ntv3_tools.locus import known_gene_metadata
from ntv3_tools.resolve import TrackResolve, resolve_tracks, write_track_resolve
from ntv3_tools.ucsc import fetch_track, write_ucsc_payloads

FetchUcsc = Callable[..., dict[str, object]]


@dataclass(frozen=True)
class LookupRequest:
    """Inputs for one deterministic metadata lookup.

    Attributes:
        species: Model species.
        tissue: Tissue or cell type.
        names: Requested assay or track names.
        out_dir: Run output directory.
        chrom: Optional chromosome for UCSC overlays.
        start: Optional zero-based window start.
        end: Optional half-open window end.
        ucsc_tracks: UCSC track names to fetch.
        gene: Requested gene symbol for knownGene context.
    """

    species: str
    tissue: str
    names: list[str]
    out_dir: str
    chrom: str | None = None
    start: int | None = None
    end: int | None = None
    ucsc_tracks: list[str] | None = None
    gene: str = ""


def run_lookup(
    request: LookupRequest,
    get_json: GetJson = default_get_json,
    fetch_ucsc: FetchUcsc = fetch_track,
) -> TrackResolve:
    """Resolve ids, write GT URLs, optionally fetch UCSC overlays.

    Args:
        request: Validated lookup inputs.
        get_json: ENCODE GET (inject in tests).
        fetch_ucsc: UCSC track fetch (inject in tests).

    Returns:
        The resolved track set (also written to ``track_resolve.json``).
    """
    ucsc_tracks = request.ucsc_tracks or []
    coordinates = (request.chrom, request.start, request.end)
    if any(value is None for value in coordinates) and any(
        value is not None for value in coordinates
    ):
        raise RequestError("UCSC context requires chrom, start, and end together.")
    chrom = cast(str, request.chrom)
    start = cast(int, request.start)
    end = cast(int, request.end)
    resolved = resolve_tracks(
        species=request.species,
        tissue=request.tissue,
        names=request.names,
        encode_search=lambda tissue_name, assay: search_encode_experiment(
            tissue=tissue_name,
            assay=assay,
            get_json=get_json,
            species=request.species,
        ),
    )
    gt_accessions = [row.gt_accession for row in resolved.tracks if row.gt_accession]
    gt_sources = [
        resolve_gt(
            assembly=resolved.assembly,
            accession=accession,
            get_json=get_json,
        )
        for accession in gt_accessions
    ]
    # Only the lanes the user approved; with no gene lane, the report reads the
    # TSS and strand from locus.json.
    context_tracks = list(dict.fromkeys(ucsc_tracks))
    ucsc_payloads = {
        track: fetch_ucsc(
            genome=resolved.assembly,
            track=track,
            chrom=chrom,
            start=start,
            end=end,
        )
        for track in context_tracks
        if all(value is not None for value in coordinates)
    }

    write_track_resolve(resolved=resolved, out_dir=request.out_dir)
    if gt_sources:
        write_gt_sources(sources=gt_sources, out_dir=request.out_dir)
    if ucsc_payloads:
        write_ucsc_payloads(payloads=ucsc_payloads, out_dir=request.out_dir)
    known_gene_payload = ucsc_payloads.get("knownGene")
    if known_gene_payload is not None:
        metadata = known_gene_metadata(
            gene=request.gene,
            chrom=chrom,
            start=start,
            end=end,
            payload=known_gene_payload,
        )
        metadata_path = Path(request.out_dir) / "locus_metadata.json"
        with metadata_path.open("x") as handle:
            json.dump(metadata, handle, indent=2)
            handle.write("\n")
    return resolved
