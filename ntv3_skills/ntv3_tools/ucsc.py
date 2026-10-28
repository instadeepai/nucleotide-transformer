"""UCSC REST fetch by documented ``track=`` names only."""

import json
from pathlib import Path

import requests

from ntv3_tools.errors import ExternalServiceError, RequestError

_UCSC_URL = "https://api.genome.ucsc.edu/getData/track"


def fetch_track(
    genome: str,
    track: str,
    chrom: str,
    start: int,
    end: int,
) -> dict[str, object]:
    """Fetch one UCSC track. Reject ENCSR* names.

    Args:
        genome: Assembly (e.g. "hg38").
        track: UCSC ``track=`` name.
        chrom: Chromosome.
        start: 0-based start.
        end: 0-based end.

    Returns:
        Parsed JSON plus the request URL.

    Raises:
        RequestError: The name is an ENCODE experiment id.
        requests.HTTPError: UCSC does not recognize the track or locus.
    """
    if track.upper().startswith("ENCSR"):
        raise RequestError(
            f"'{track}' is an ENCODE experiment id. UCSC /getData/track "
            "returns 400 for ENCSR*. Use a UCSC track name from "
            "the UCSC browser catalog or fetch GT through ENCODE."
        )
    url = (
        f"{_UCSC_URL}?genome={genome};track={track};"
        f"chrom={chrom};start={start};end={end}"
    )
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        status = getattr(getattr(exc, "response", None), "status_code", "network")
        raise ExternalServiceError(
            f"UCSC track request failed ({status}) for track '{track}'."
        ) from exc
    payload["_request_url"] = url
    payload["_label"] = track
    return payload


def write_ucsc_overlay(
    genome: str,
    chrom: str,
    start: int,
    end: int,
    tracks: list[str],
    out_dir: str,
) -> Path:
    """Fetch each track and write ucsc/<track>.json.

    Args:
        genome: Assembly.
        chrom: Chromosome.
        start: 0-based start.
        end: 0-based end.
        tracks: UCSC track names.
        out_dir: Run output directory.

    Returns:
        Directory containing the overlay JSON files.
    """
    payloads: dict[str, dict[str, object]] = {}
    for track in tracks:
        payloads[track] = fetch_track(
            genome=genome, track=track, chrom=chrom, start=start, end=end
        )
    return write_ucsc_payloads(payloads=payloads, out_dir=out_dir)


def write_ucsc_payloads(payloads: dict[str, dict[str, object]], out_dir: str) -> Path:
    """Write already-fetched UCSC overlays.

    Args:
        payloads: Mapping from UCSC track name to validated response.
        out_dir: Run output directory.

    Returns:
        Directory containing the overlay JSON files.
    """
    root = Path(out_dir) / "ucsc"
    root.mkdir(parents=True, exist_ok=True)
    for track, payload in payloads.items():
        (root / f"{track}.json").write_text(json.dumps(payload, indent=2) + "\n")
    return root
