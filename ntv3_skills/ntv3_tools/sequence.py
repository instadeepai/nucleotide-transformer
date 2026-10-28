"""Reference-sequence retrieval via the UCSC REST API.

Lifted from tutorial 01: sequences are fetched from ``api.genome.ucsc.edu`` rather than a local
FASTA, so the skill has no local-data or repo dependency. The model downsamples by 128, so sequence
length must be cropped to a multiple of 128.
"""

import logging

import requests

from ntv3_tools.errors import ExternalServiceError, RequestError

logger = logging.getLogger("ntv3_tools")

_UCSC_URL = "https://api.genome.ucsc.edu/getData/sequence"
_UCSC_CHROMOSOMES_URL = "https://api.genome.ucsc.edu/list/chromosomes"
_DOWNSAMPLE = 128
# keep_target_center_fraction of NTv3_tracks.
CENTER_FRACTION = 0.375


def fetch_ucsc_sequence(
    assembly: str,
    chrom: str,
    start: int,
    end: int,
    keep_soft_mask: bool = False,
) -> str:
    """Fetch the reference sequence for a window from the UCSC API.

    Args:
        assembly: Genome assembly (e.g. "hg38").
        chrom: Chromosome name (e.g. "chr11").
        start: 0-based window start.
        end: 0-based window end.
        keep_soft_mask: Keep UCSC's lowercase repeats. The tracks model
            tokenizes them as their own ids, so predict keeps them.

    Returns:
        Reference sequence spanning [start, end), uppercase unless
        ``keep_soft_mask``.
    """
    url = f"{_UCSC_URL}?genome={assembly};chrom={chrom};start={start};end={end}"
    try:
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError) as exc:
        status = getattr(getattr(exc, "response", None), "status_code", "network")
        raise ExternalServiceError(
            f"UCSC sequence request failed ({status}) for "
            f"{chrom}:{start}-{end} ({assembly})."
        ) from exc
    if "dna" not in data:
        raise ValueError(
            f"UCSC returned no 'dna' for {chrom}:{start}-{end} ({assembly}): {data}"
        )
    seq = data["dna"] if keep_soft_mask else data["dna"].upper()
    logger.info(
        f"Fetched {len(seq)} bp from UCSC for {chrom}:{start}-{end} ({assembly})"
    )
    return seq


def fetch_chrom_size(assembly: str, chrom: str) -> int:
    """Return a chromosome's length in bp from the UCSC API.

    Args:
        assembly: Genome assembly (e.g. "hg38").
        chrom: Chromosome name (e.g. "chr11").

    Returns:
        Chromosome length.

    Raises:
        ExternalServiceError: The request failed.
        RequestError: The assembly has no such chromosome.
    """
    try:
        resp = requests.get(
            _UCSC_CHROMOSOMES_URL, params={"genome": assembly}, timeout=60
        )
        resp.raise_for_status()
        sizes = resp.json()["chromosomes"]
    except (requests.RequestException, ValueError, KeyError) as exc:
        raise ExternalServiceError(
            f"UCSC chromosome list request failed for {assembly}."
        ) from exc
    if chrom not in sizes:
        raise RequestError(f"{assembly} has no chromosome {chrom}.")
    return int(sizes[chrom])


def crop_to_multiple(sequence: str, multiple: int = _DOWNSAMPLE) -> str:
    """Crop a sequence to a multiple of the model downsampling factor.

    Args:
        sequence: Input DNA sequence.
        multiple: Downsampling factor the length must divide (128 for NTv3).

    Returns:
        The sequence truncated to ``len // multiple * multiple``.
    """
    return sequence[: (len(sequence) // multiple) * multiple]


def prediction_region(
    start: int, end: int, center_frac: float = CENTER_FRACTION
) -> tuple[int, int]:
    """Return the genomic span the model predicts (center crop).

    Args:
        start: 0-based window start.
        end: 0-based window end.
        center_frac: Fraction of the window kept (config default 0.375).

    Returns:
        The ``(prediction_start, prediction_end)`` genomic coordinates.
    """
    window = end - start
    offset_frac = (1.0 - center_frac) / 2.0
    pred_start = start + int(window * offset_frac)
    pred_end = pred_start + int(window * center_frac)
    return pred_start, pred_end
