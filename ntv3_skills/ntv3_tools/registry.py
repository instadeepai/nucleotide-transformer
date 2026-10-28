"""NTv3 checkpoint names and track-id rules.

Three ID systems. Do not collapse them.

1. UCSC ``track=`` names — genes, ClinVar, H3K27ac. Never ``ENCSR*``.
2. ENCODE experiments (``ENCSR*``) — ground-truth BigWigs.
3. Model heads in ``cfg.bigwigs_per_file_assembly`` — often the same ENCSR,
   sometimes FANTOM ``CNhs*``.

Hints map assay language to both a GT experiment and a model head. They are
not the full catalog. Unknown names pass through as ENCSR / CNhs or go to a
live ENCODE search.
"""

from collections.abc import Collection
from dataclasses import dataclass
from typing import Any

from ntv3_tools.errors import RequestError

# The released functional-tracks model, loaded by
# ``nucleotide_transformer_v3.pretrained.load_tracks_model``.
DEFAULT_CHECKPOINT = "NTv3_tracks"
_HUB_ORG = "InstaDeepAI/"


@dataclass(frozen=True)
class TrackIds:
    """One requested track in both ID systems.

    Attributes:
        query: What the user typed.
        gt_accession: ENCODE experiment for GT, or None (FANTOM-only heads).
        model_head: Name in ``cfg.bigwigs_per_file_assembly``, or None if unknown.
        source: How the mapping was established.
    """

    query: str
    gt_accession: str | None
    model_head: str | None
    source: str


@dataclass(frozen=True)
class TrackHint:
    """Validated assay-language shortcut. Not an exhaustive catalog.

    Attributes:
        species: Model species for this mapping.
        tissue: Biosample or cell type for this mapping.
        gt: ENCODE experiment for ground truth.
        model: Model-head id (ENCSR or CNhs).
        source: Provenance for the mapping.
    """

    species: str
    tissue: str
    gt: str | None
    model: str | None
    source: str


SPECIES_TO_BUILD: dict[str, str] = {"human": "hg38", "mouse": "mm10"}

# Convenience only. Live ENCODE / cfg lookup is the real path.
TRACK_HINTS: dict[str, TrackHint] = {
    "K562 RNA-seq": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR056HPM",
        model="ENCSR056HPM",
        source="NTv3 checkpoint track metadata and ENCODE experiment",
    ),
    "K562 DNase": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR921NMD",
        model="ENCSR921NMD",
        source="NTv3 checkpoint track metadata and ENCODE experiment",
    ),
    "K562 H3K4me3": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR000DWD",
        model="ENCSR000DWD",
        source="NTv3 checkpoint track metadata and ENCODE experiment",
    ),
    "K562 CTCF": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR000AKO",
        model="ENCSR000AKO",
        source="NTv3 checkpoint track metadata and ENCODE experiment",
    ),
    "K562 GATA1": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR000EFT",
        model="ENCSR000EFT",
        source="NTv3 checkpoint track metadata and ENCODE experiment",
    ),
    # Paper / cfg head 5981 is FANTOM. ENCSR000CIO is the ENCODE GT file.
    "K562 CAGE": TrackHint(
        species="human",
        tissue="K562",
        gt="ENCSR000CIO",
        model="CNhs11250_P",
        source="NTv3 paper head 5981 and ENCODE experiment ENCSR000CIO",
    ),
}


def resolve_model(name: str) -> str:
    """Resolve a predict model name to the released tracks model.

    Args:
        name: ``NTv3_tracks`` or ``InstaDeepAI/NTv3_tracks``.

    Returns:
        ``NTv3_tracks``, the name ``load_tracks_model`` takes.

    Raises:
        KeyError: Any other name. Old internal checkpoint keys are not a load
            path and are never sent to the Hub.
    """
    if name.strip().removeprefix(_HUB_ORG) == DEFAULT_CHECKPOINT:
        return DEFAULT_CHECKPOINT
    raise KeyError(
        f"Predict loads only {DEFAULT_CHECKPOINT}; got '{name}'. Omit --model."
    )


def find_track_hint(name: str, tissue: str) -> tuple[str, TrackHint] | None:
    """Find the hint for a track name, alone or after the tissue.

    "K562 CAGE" matches its hint directly. "CAGE" with tissue K562 matches the
    same hint. Case is ignored. No network call is made.

    Args:
        name: Track name as the user typed it.
        tissue: Tissue or cell type from the same request.

    Returns:
        The hint key and hint, or None when no hint matches.
    """
    keys = {key.casefold(): key for key in TRACK_HINTS}
    query = " ".join(name.split()).casefold()
    with_tissue = " ".join(f"{tissue} {name}".split()).casefold()
    for candidate in (query, with_tissue):
        if candidate in keys:
            key = keys[candidate]
            return key, TRACK_HINTS[key]
    return None


def resolve_track_ids(name_or_accession: str) -> TrackIds:
    """Split a user track string into GT experiment vs model head.

    Args:
        name_or_accession: Alias, ENCSR experiment, or CNhs model head.

    Returns:
        Dual ids. ``model_head`` is what predict uses. ``gt_accession`` is
        what ``ntv3 gt`` uses.

    Raises:
        KeyError: Not a hint, not ENCSR, not CNhs.
    """
    query = name_or_accession.strip()
    if query in TRACK_HINTS:
        hint = TRACK_HINTS[query]
        return TrackIds(
            query=query,
            gt_accession=hint.gt,
            model_head=hint.model,
            source=hint.source,
        )
    if query.upper().startswith("ENCSR"):
        return TrackIds(
            query=query,
            gt_accession=query,
            model_head=query,
            source="Explicit ENCODE experiment supplied by the user",
        )
    if query.startswith("CNhs"):
        return TrackIds(
            query=query,
            gt_accession=None,
            model_head=query,
            source="Explicit FANTOM model head supplied by the user",
        )
    raise KeyError(
        f"Which model head should be predicted for '{query}'? It is not an ENCSR "
        "experiment, a CNhs model head, or a known track name. Pass it in "
        "--tracks. No head is assumed."
    )


def resolve_track(name_or_accession: str) -> str:
    """Return the model-head id for a user track string.

    Args:
        name_or_accession: Alias, ENCSR experiment, or CNhs model head.

    Returns:
        The id to look up in ``cfg.bigwigs_per_file_assembly``.
    """
    ids = resolve_track_ids(name_or_accession=name_or_accession)
    if ids.model_head:
        return ids.model_head
    if ids.gt_accession:
        return ids.gt_accession
    raise KeyError(f"Track '{name_or_accession}' has no model head.")


def head_indices(names: list[str], heads: list[str]) -> list[int]:
    """Return each requested head's column, refusing any id that is not an exact head.

    No other head is ever substituted for a missing one.

    Args:
        names: ``cfg.bigwigs_per_file_assembly[assembly]``.
        heads: Requested model-head ids.

    Returns:
        Output track indices aligned with ``heads``.

    Raises:
        RequestError: A requested id is not one of ``names``; the message names it.
    """
    missing = [head for head in heads if head not in names]
    if missing:
        raise RequestError(
            f"Track(s) {missing} are not model heads of this model "
            f"({len(names)} heads for this assembly). No other head is substituted. "
            "Use the exact model_head from track_resolve.json."
        )
    return [names.index(head) for head in heads]


def model_head_for(row: dict[str, Any], available: Collection[str]) -> str:
    """Return a resolved track's model head; never its ENCODE accession in its place.

    Args:
        row: One ``track_resolve.json`` track row.
        available: Model heads of the loaded model for this assembly.

    Returns:
        The row's ``model_head``.

    Raises:
        RequestError: The row has no model head. When its ENCODE accession is
            also a head id, the message names the substitution it refused.
    """
    model_head = str(row.get("model_head") or "")
    if model_head:
        return model_head
    query = str(row.get("query") or "track")
    accession = str(row.get("gt_accession") or "")
    if accession and accession in available:
        raise RequestError(
            f"Track '{query}' has no model head. Its ENCODE accession {accession} "
            f"would have been substituted as model head {accession}; an ENCODE "
            "accession is not a model head. Pass the model head explicitly."
        )
    raise RequestError(
        f"Track '{query}' has no validated model head (ENCODE accession "
        f"{accession or 'none'}). Pass the model head explicitly."
    )


def assembly_for_species(species: str) -> str:
    """Return the genome assembly for a species, refusing to guess unknown species.

    Args:
        species: Species string accepted by the model (e.g. "human").

    Returns:
        The genome assembly (e.g. "hg38").
    """
    if species not in SPECIES_TO_BUILD:
        raise KeyError(
            f"No genome build registered for species '{species}'. Known: {list(SPECIES_TO_BUILD)}. "
            "Never guess the build — add it explicitly."
        )
    return SPECIES_TO_BUILD[species]
