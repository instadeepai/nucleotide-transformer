"""Resolve tissue / track names to GT experiment + model head. Never invent an id."""

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from ntv3_tools.errors import RequestError
from ntv3_tools.registry import (
    SPECIES_TO_BUILD,
    TrackIds,
    assembly_for_species,
    find_track_hint,
    resolve_track_ids,
)

EncodeSearch = Callable[[str, str], str]


@dataclass(frozen=True)
class ResolvedTrack:
    """One requested track after alias / live lookup.

    Attributes:
        query: What the user typed.
        gt_accession: ENCODE experiment for GT, or empty.
        model_head: Id in the model cfg, or empty.
        species: Species used for the lookup.
        source: Provenance for the dual-id mapping.
        encode_hits: Experiments a live ENCODE search returned, or None when
            no search ran.
    """

    query: str
    gt_accession: str
    model_head: str
    species: str
    source: str
    encode_hits: tuple[str, ...] | None = None


@dataclass(frozen=True)
class TrackResolve:
    """Resolved track set for one locus.

    Attributes:
        species: Confirmed species.
        assembly: Confirmed assembly.
        tissue: Confirmed tissue / cell type.
        tracks: Resolved rows.
    """

    species: str
    assembly: str
    tissue: str
    tracks: list[ResolvedTrack]


def resolve_tracks(
    species: str,
    tissue: str,
    names: list[str],
    encode_search: EncodeSearch | None = None,
) -> TrackResolve:
    """Map friendly names to GT experiments and model heads.

    A hint matches the name alone or the tissue plus the name ("CAGE" in
    K562), and wins without a network call. Only names with no hint, ENCSR,
    or CNhs form reach the live ENCODE search.

    Args:
        species: Model species. Must already be unambiguous.
        tissue: Tissue or cell type. Must already be unambiguous.
        names: Assay words, aliases, ENCSR experiments, or CNhs heads.
        encode_search: Optional live ENCODE experiment lookup for unknown names.

    Returns:
        Resolved track set.

    Raises:
        RequestError: Species, tissue, or a track name is unset.
        KeyError: A name is not a hint, ENCSR, CNhs, or ENCODE hit.
    """
    if not species or species.lower() in {"unknown", "guess"}:
        raise RequestError("Which species? Pass --species. No species is assumed.")
    if species not in SPECIES_TO_BUILD:
        raise RequestError(
            f"Which species? This tool predicts {', '.join(SPECIES_TO_BUILD)}, not "
            f"'{species}'. Pass --species. No species is assumed."
        )
    if not tissue:
        raise RequestError(
            "Which tissue or cell type should the tracks be predicted in? Pass "
            "--tissue. No tissue is assumed."
        )
    if not names:
        raise RequestError(
            "Which tracks should be predicted? Pass --tracks with assay names or "
            "exact model heads (ENCSR*/CNhs*). No track is chosen for you."
        )
    rows: list[ResolvedTrack] = []
    for name in names:
        hinted = find_track_hint(name=name, tissue=tissue)
        if hinted is not None:
            key, hint = hinted
            if (
                hint.species.casefold() != species.casefold()
                or hint.tissue.casefold() != tissue.casefold()
            ):
                raise RequestError(
                    f"Track mapping '{name}' is validated for "
                    f"{hint.species}/{hint.tissue}, not {species}/{tissue}."
                )
            rows.append(
                ResolvedTrack(
                    query=name,
                    gt_accession=hint.gt or "",
                    model_head=hint.model or "",
                    species=species,
                    source=f"Track hint '{key}': {hint.source}",
                )
            )
            continue
        ids, encode_hits = _resolve_one(
            name=name, tissue=tissue, encode_search=encode_search
        )
        rows.append(
            ResolvedTrack(
                query=name,
                gt_accession=ids.gt_accession or "",
                model_head=ids.model_head or "",
                species=species,
                source=ids.source,
                encode_hits=encode_hits,
            )
        )
    return TrackResolve(
        species=species,
        assembly=assembly_for_species(species=species),
        tissue=tissue,
        tracks=rows,
    )


def _resolve_one(
    name: str,
    tissue: str,
    encode_search: EncodeSearch | None,
) -> tuple[TrackIds, tuple[str, ...] | None]:
    """Resolve one unhinted name, then the live ENCODE search if provided.

    Returns:
        The ids, and the search hits when the live search ran.
    """
    try:
        return resolve_track_ids(name_or_accession=name), None
    except KeyError:
        if encode_search is None:
            raise
        experiment = encode_search(tissue, name)
        ids = TrackIds(
            query=name,
            gt_accession=experiment,
            model_head=None,
            source=(
                "Live ENCODE experiment search; no model head was inferred. "
                "Validate a model head before prediction."
            ),
        )
        return ids, (experiment,)


def write_track_resolve(resolved: TrackResolve, out_dir: str) -> Path:
    """Write track_resolve.json.

    Args:
        resolved: Resolved track set.
        out_dir: Run output directory.

    Returns:
        Path to track_resolve.json.
    """
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "track_resolve.json"
    payload = {
        "species": resolved.species,
        "assembly": resolved.assembly,
        "tissue": resolved.tissue,
        "tracks": [
            {
                "query": row.query,
                "gt_accession": row.gt_accession,
                "model_head": row.model_head,
                "source": row.source,
                **(
                    {"encode_hits": list(row.encode_hits)}
                    if row.encode_hits is not None
                    else {}
                ),
            }
            for row in resolved.tracks
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path
