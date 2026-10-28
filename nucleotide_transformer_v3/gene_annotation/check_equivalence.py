r"""Compare two annotation GFF3s, and say whether a difference looks structural.

Tier 3 of the production pipeline's regression contract, and the only tier that
needs a real checkpoint and a GPU, so it is a command you run rather than a test
CI runs.

Its job is to answer one question about a set of differing genes: is this a
stitching, cropping or padding error, or is it numerical noise? A stitching
error clusters at window seams; numerical noise does not. So the tool reports
the distance from each differing gene to the nearest seam alongside what a
uniform distribution would predict, which is the comparison that makes the
answer readable rather than a matter of opinion.

Comment lines are ignored, so production's provenance header does not register
as a difference against a GFF3 from another tool that has none.

Run both pipelines with the filters *off* and bf16 *off*, or the comparison is
not measuring what it appears to (bf16 on a GPU is not bit-reproducible):

    ntv3-annotate --genome <contig.fa> --out /tmp/ported.gff3 \\
      --override weights_dir=<dir> --seq_len <L> \\
      --override use_bf16_compute=false \\
      --override apply_biological_filters=false --override min_cds_length=0

    python -m nucleotide_transformer_v3.gene_annotation.check_equivalence \\
      /tmp/ported.gff3 /tmp/reference.gff3

Exits non-zero on any difference, so it can gate a job.
"""

import argparse
import gzip
import sys
from collections import Counter



from pathlib import Path

# How close to a seam a differing gene has to be to count as "at" one. A
# stitching error lands within a few bp; 2 kb is generous on purpose, so the
# diagnostic errs towards flagging a structural problem rather than missing one.
SEAM_WINDOW_BP = 2_000


def report(message: str) -> None:
    """Write one line of the comparison report to stdout.

    Deliberately ``print`` rather than ``logger.info``. This tool exists to say
    what it found, and on a real GPU run every ``logger.info`` here vanished:
    ``get_logger`` configures a level but the handler that actually renders it
    is set up by the Hydra entry points, which a bare script never goes through.
    The job still exited 0, so the run looked like a success and produced no
    verdict at all. A reporting tool must not depend on someone else's logging
    configuration.
    """
    print(message, flush=True)


def read_feature_lines(path: str | Path) -> list[str]:
    """Read a GFF3's feature lines, ignoring comments and blanks.

    Args:
        path: GFF3 to read, gzipped or plain.

    Returns:
        The feature lines, in file order.
    """
    return [ln for ln in _read_text(path).splitlines() if ln and not ln.startswith("#")]


# 0 identical, 1 a real DIFFERENT verdict, 2 no verdict reached. Distinct because the
# shell jobs write `exit=$?` into their report, and with both failures exiting 1 that
# line could not tell "the pipelines differ" from "the comparison never ran" - which is
# why the guard has to read the report text at all.
UNREADABLE_EXIT_CODE = 2


def _explain_an_unreadable_path(path: str, error: Exception) -> SystemExit:
    """Turn any failure to read a GFF3 into an actionable failure.

    Deliberately broad: this catches every read failure, a local typo included,
    not only the remote-path case. That is right for a tool whose job is to say
    what went wrong rather than die deep inside a dependency, and the message
    hedges accordingly with "if this is a remote path".

    Remote paths are the motivating story. This package reads with ``pathlib``
    and imports no cloud SDK, so ``Path("s3://bucket/x.gff3")`` constructs
    happily and only ``read_bytes()`` fails, several frames down, with no hint
    of what to do. A plant equivalence job once spent both of its GPU arms and
    produced no verdict that way: the research pipeline writes to S3, the
    comparison was pointed straight at it, and the failure arrived only after
    all the compute was gone.

    The read is ATTEMPTED rather than refused by scheme, so that a path this
    process happens to be able to read is not rejected on the strength of how
    it is spelled.
    """
    # Printed rather than carried as the SystemExit's argument, because the
    # argument is now the exit CODE: SystemExit("text") exits 1, which is also
    # a real DIFFERENT verdict's code.
    print(
        f"{path}: could not be read ({type(error).__name__}: {error}). If this "
        "is a remote path, stage it locally first: this package reads local "
        "files only and imports no cloud SDK by design.",
        file=sys.stderr,
    )
    return SystemExit(2)


def _read_text(path: str | Path) -> str:
    """Read a GFF3 as text, transparently decompressing ``.gz``."""
    resolved = Path(path)
    try:
        raw = resolved.read_bytes()
    except Exception as error:
        raise _explain_an_unreadable_path(str(path), error) from error
    # Case-insensitive, so a file this pipeline is willing to write is a file
    # this tool is willing to read.
    if str(resolved).lower().endswith(".gz"):
        return gzip.decompress(raw).decode("utf-8")
    return raw.decode("utf-8")


def read_provenance(path: str | Path) -> dict[str, str]:
    """The ``##ntv3-*`` provenance a production GFF3 carries, if any.

    A GFF3 from another tool has none, which is why the seam stride falls back
    to a flag rather than requiring it.

    Args:
        path: GFF3 to read, gzipped or plain.

    Returns:
        Provenance keys without their ``##ntv3-`` prefix, mapped to their raw
        values. Empty for a file that carries none.
    """
    out: dict[str, str] = {}
    for line in _read_text(path).splitlines():
        if not line.startswith("#"):
            break
        if line.startswith("##ntv3-"):
            key, _, value = line[len("##ntv3-") :].partition(" ")
            out[key] = value.strip()
    return out


def stride_from_provenance(path: str | Path) -> int | None:
    """The decoded stride a production GFF3 implies, or None if it cannot say.

    Computed through ``compute_window_geometry`` rather than as
    ``sequence_length * crop``, because the two differ by one whenever the
    uncropped remainder is odd, and that error accumulates across windows.

    Args:
        path: GFF3 to read.

    Returns:
        Stride in bp, or None when the file carries no usable geometry.
    """
    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        compute_window_geometry,
    )

    provenance = read_provenance(path)
    try:
        length = int(provenance["sequence-length"])
        crop_fraction = float(provenance["keep-target-center-fraction"])
    except (KeyError, ValueError):
        return None
    if length <= 0 or not 0.0 < crop_fraction <= 1.0:
        return None
    return length - 2 * compute_window_geometry(length, crop_fraction)


def gene_starts(lines: list[str]) -> list[tuple[str, int]]:
    """``(seqid, start)`` of every gene row."""
    return [
        (cols[0], int(cols[3]))
        for cols in (line.split("\t") for line in lines)
        if cols[2] == "gene"
    ]


def coordinate_key(line: str) -> str:
    """A feature's identity without its attributes column.

    Gene ids are assigned sequentially, so one extra gene early in a contig
    renumbers every gene after it and an attribute-inclusive comparison reports
    a near-total difference that says nothing about whether the two pipelines
    found the same genes. Comparing on seqid, type, coordinates, strand and
    phase answers the question actually being asked; the id difference is
    reported separately.
    """
    return "\t".join(line.split("\t")[:8])


def distance_to_seam(start: int, stride: int) -> int:
    """Distance in bp from a 1-based coordinate to the nearest window seam.

    Seams sit at multiples of the decoded stride, which is where two windows'
    retained regions meet and where a crop or padding error would surface. That
    includes each contig's last seam: ``plan_contig_windows`` drops the shifted
    tail window's overlap, so its retained slice also starts on the grid.
    """
    offset = (start - 1) % stride
    return min(offset, stride - offset)


def summarise_seam_proximity(
    starts: list[tuple[str, int]], stride: int
) -> tuple[int, float, list[int]]:
    """How many differing genes sit near a seam, and how many would by chance.

    Args:
        starts: ``(seqid, start)`` of the differing genes.
        stride: Decoded bp per window.

    Returns:
        The count within ``SEAM_WINDOW_BP`` of a seam, the count a uniform
        distribution would predict, and every distance, ascending.
    """
    distances = sorted(distance_to_seam(start, stride) for _, start in starts)
    near = sum(1 for d in distances if d <= SEAM_WINDOW_BP)
    # A uniform start lands within SEAM_WINDOW_BP of a seam with probability
    # 2 * SEAM_WINDOW_BP / stride, capped at 1.
    expected = len(distances) * min(1.0, 2 * SEAM_WINDOW_BP / stride)
    return near, expected, distances


def _report_side(label: str, only: list[str], limit: int) -> None:
    """Log the features present on one side only."""
    report(f"  only in {label}: {len(only)}")
    for line in only[:limit]:
        report(f"    {line}")
    if len(only) > limit:
        report(f"    ... and {len(only) - limit} more")


def compare(left_path: str, right_path: str, stride: int | None, limit: int) -> bool:
    """Compare two GFF3s and log a verdict.

    Args:
        left_path: First GFF3, conventionally production's.
        right_path: Second GFF3, conventionally the reference one's.
        stride: Decoded bp per window, enabling the seam diagnostic. None skips
            it, since without the geometry a seam distance cannot be computed.
        limit: Maximum differing features to print per side.

    Returns:
        Whether the two are identical.
    """
    left = read_feature_lines(left_path)
    right = read_feature_lines(right_path)

    # Coordinates are the primary comparison; whole lines are reported too, so
    # an id-only difference is visible but does not swamp the verdict.
    left_keys = [coordinate_key(ln) for ln in left]
    right_keys = [coordinate_key(ln) for ln in right]
    left_set, right_set = set(left_keys), set(right_keys)
    only_left = [
        ln for ln, k in zip(left, left_keys, strict=True) if k not in right_set
    ]
    only_right = [
        ln for ln, k in zip(right, right_keys, strict=True) if k not in left_set
    ]

    left_genes, right_genes = gene_starts(left), gene_starts(right)
    report(f"{left_path}: {len(left)} features, {len(left_genes)} genes")
    report(f"{right_path}: {len(right)} features, {len(right_genes)} genes")

    if len(left_keys) != len(left_set) or len(right_keys) != len(right_set):
        # Through report(), not logger: this is one of the two most consequential
        # things the tool can say and it printed nothing as a bare script, which
        # is the exact failure this module's docstring was written about.
        report(
            "  NOTE: duplicate coordinate keys present "
            f"({len(left_keys) - len(left_set)} left, "
            f"{len(right_keys) - len(right_set)} right). Counted with "
            "multiplicity below, so the verdict accounts for them."
        )

    # Multiplicity, not set membership. Sets collapse duplicates, so two files
    # with different feature counts could satisfy `not only_left and not
    # only_right` and be reported IDENTICAL with exit 0, against this module's
    # promise to exit non-zero on any difference. Production output cannot
    # contain duplicate keys (genes are maximal non-intergenic runs), but this
    # tool takes arbitrary GFF3 on both sides, and a reference annotation whose
    # isoforms share CDS intervals has duplicates routinely.
    counts_match = Counter(left_keys) == Counter(right_keys)

    identical_lines = left == right
    if not only_left and not only_right and counts_match:
        if identical_lines:
            report("IDENTICAL: every feature matches, attributes included.")
            return True
        report(
            "IDENTICAL COORDINATES: every feature matches on seqid, type, "
            "coordinates, strand and phase, but the attribute columns differ. "
            "That is gene-id renumbering, not a difference in what was found."
        )
        return True

    if not only_left and not only_right and not counts_match:
        left_counts, right_counts = Counter(left_keys), Counter(right_keys)
        repeated = [
            (key, left_counts[key], right_counts[key])
            for key in left_counts | right_counts
            if left_counts[key] != right_counts[key]
        ]
        report(
            f"DIFFERENT: both files use the same {len(left_set)} coordinate "
            f"key(s), but {len(repeated)} of them appear a different number of "
            f"times ({len(left)} vs {len(right)} features in total)."
        )
        for key, n_left, n_right in sorted(repeated, key=str)[:limit]:
            report(f"  {key}: {n_left} left, {n_right} right")
        return False

    shared = len(left_set & right_set)
    total = max(len(left_set), len(right_set))
    report(
        f"DIFFERENT: {shared}/{total} features share coordinates "
        f"({100 * shared / total if total else 0:.4f}%)."
    )
    _report_side("left", only_left, limit)
    _report_side("right", only_right, limit)

    differing_genes = [
        (cols[0], int(cols[3]))
        for cols in (line.split("\t") for line in only_left + only_right)
        if cols[2] == "gene"
    ]
    per_contig = Counter(seqid for seqid, _ in differing_genes)
    if per_contig:
        report(f"  differing genes per contig: {dict(per_contig)}")

    if stride is None:
        report(
            "  seam diagnostic skipped: neither file carries ##ntv3- "
            "provenance to read the geometry from. Pass --stride "
            "<decoded bp per window> to tell a stitching error from noise."
        )
        return False
    if not differing_genes:
        report("  no differing gene rows, so there is no seam test to run.")
        return False

    if 2 * SEAM_WINDOW_BP >= stride:
        # Every position is then "near" a seam, so the test cannot discriminate.
        # Saying "31/31 within a seam, 31.0 expected" would read as a finding.
        report(
            f"  NOTE: seam diagnostic carries NO information at stride "
            f"{stride:,} bp: "
            f"the {SEAM_WINDOW_BP:,} bp window covers the whole stride, so every "
            "position counts as near a seam. Real inference strides (tens of kb "
            "upwards) are far larger; this verdict is only meaningful there."
        )
        return False

    near, expected, distances = summarise_seam_proximity(differing_genes, stride)
    report(
        f"  seam proximity (stride {stride:,} bp, window {SEAM_WINDOW_BP:,} bp): "
        f"{near}/{len(distances)} within a seam, {expected:.1f} expected by chance"
    )
    report(
        f"  nearest {distances[0]:,} bp, median {distances[len(distances) // 2]:,} bp"
    )
    if near > expected:
        report(
            "  WARNING: differing genes cluster AT the seams, more than chance "
            "predicts. "
            "That is the signature of a stitching, cropping or padding error, "
            "not of numerical noise. Investigate the window bookkeeping."
        )
    else:
        report(
            "  differing genes do NOT cluster at the seams, which is consistent "
            "with numerical noise rather than a stitching error. Confirm by "
            "re-running one pipeline twice in separate processes: if that also "
            "differs, these are noise."
        )
    return False


def main() -> None:
    """Compare two GFF3s and exit non-zero if they differ."""
    parser = argparse.ArgumentParser(
        description="Compare two annotation GFF3s feature by feature.",
        epilog="Run both pipelines with filters off and bf16 off.",
    )
    parser.add_argument("left", help="First GFF3 (conventionally production's).")
    parser.add_argument("right", help="Second GFF3 (the reference).")
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help=(
            "Decoded bp per window, for the seam diagnostic. Usually "
            "unnecessary: it is read from the GFF3's own ##ntv3- provenance. "
            "Pass it for a file that carries none, or to override."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Maximum differing features to print per side (default 10).",
    )
    args = parser.parse_args()

    if args.stride is not None and args.stride <= 0:
        parser.error("--stride must be positive.")

    stride = args.stride
    if stride is None:
        for candidate in (args.left, args.right):
            stride = stride_from_provenance(candidate)
            if stride is not None:
                report(
                    f"Seam stride {stride:,} bp, read from {candidate}'s own "
                    "provenance."
                )
                break

    identical = compare(args.left, args.right, stride, args.limit)
    sys.exit(0 if identical else 1)


if __name__ == "__main__":
    main()
