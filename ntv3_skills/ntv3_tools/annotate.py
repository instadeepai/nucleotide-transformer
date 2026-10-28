"""Annotate ask: input, plan, and Resolved block for the NTv3 annotation models.

The ``InstaDeepAI/NTv3_annot_<clade>`` models predict one of 15 gene-structure
states for every base on each strand. ``ntv3 annotate`` writes
``annotation_plan.json`` and the Resolved block; after ``ntv3 approve``,
``ntv3 annotate-forward`` runs the plan on a GPU (:mod:`ntv3_tools.annotate_science`).

The input is a FASTA, or the 1 Mb centred on an existing Predict run's window, so
that the predicted genes are drawn next to its tracks.
"""

import gzip
import json
import math
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from ntv3_tools.annotation_report import ANNOTATION_PLAN
from ntv3_tools.approval import ANNOTATION_RESOLVED_BLOCK
from ntv3_tools.doctor import DoctorReport, jax_python_line
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.plan import read_plan
from ntv3_tools.resolved import hub_weights_cached
from ntv3_tools.sequence import fetch_chrom_size, fetch_ucsc_sequence

# Keyed by the model config names in nucleotide_transformer_v3/gene_annotation/model_cfg.
ANNOTATION_MODELS: dict[str, str] = {
    "mammals": "InstaDeepAI/NTv3_annot_mammals",
    "plants": "InstaDeepAI/NTv3_annot_plants",
    "insects": "InstaDeepAI/NTv3_annot_insects",
}
CLADE_ALIASES: dict[str, str] = {"insecta": "insects"}
ANNOTATION_DECODERS = ("argmax", "hmm", "hmm_constrained")
# The released models read 1 Mb windows and keep the central 75% of each.
ANNOTATION_WINDOW_BP = 1_048_576
ANNOTATION_KEEP_FRACTION = 0.75
MIN_CDS_LENGTH = 201
ANNOTATION_INPUT = "annotation_input.fa"
ANNOTATION_OUTPUT = "annotation.gff3"
FORWARD_COMMAND = (
    "PYTHONPATH=ntv3_skills <jax-python> -m ntv3_tools.cli annotate-forward "
    "--plan {plan}"
)
_IUPAC = frozenset("ACGTNRYKMSWBDHV")


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record.

    Attributes:
        seqid: Record name, the first word of the header.
        sequence: Bases as written. Lowercase marks soft-masked repeats.
    """

    seqid: str
    sequence: str


@dataclass(frozen=True)
class LocusWindow:
    """The genomic window a locus ask annotates, 0-based half-open.

    Attributes:
        assembly: UCSC assembly.
        chrom: Chromosome.
        start: Window start.
        end: Window end.
    """

    assembly: str
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class AnnotateRequest:
    """An Annotate ask, collected before any run.

    Attributes:
        clade: "mammals", "plants", or "insects" ("insecta" is accepted).
        decoder: Decoder name from ``ANNOTATION_DECODERS``.
        out_dir: Run directory: a new ``results/<run>`` for a FASTA, or the
            Predict run for a locus.
        fasta: Plain or gzipped FASTA; None for a locus ask.
    """

    clade: str
    decoder: str
    out_dir: str
    fasta: str | None = None


def clade_config(clade: str) -> str:
    """Return the model config name for a clade.

    Args:
        clade: Clade name as asked; ``insecta`` is accepted for ``insects``.

    Returns:
        "mammals", "plants", or "insects".

    Raises:
        RequestError: The clade has no released model.
    """
    name = clade.strip().lower()
    name = CLADE_ALIASES.get(name, name)
    if name not in ANNOTATION_MODELS:
        raise RequestError(
            f"Which clade? There is no NTv3 annotation model for '{clade}'; the "
            f"models are {', '.join(ANNOTATION_MODELS)}. Pass --clade. No clade is "
            "assumed."
        )
    return name


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read and validate a plain or gzipped FASTA.

    Args:
        path: FASTA path.

    Returns:
        Records in file order.

    Raises:
        RequestError: The file is missing, empty, or malformed.
    """
    fasta_path = Path(path).expanduser()
    if not fasta_path.is_file():
        raise RequestError(
            f"Which FASTA? {fasta_path} does not exist. Pass --fasta with an absolute "
            "path. No file is assumed."
        )
    opener = gzip.open if fasta_path.suffix == ".gz" else open
    names: list[str] = []
    chunks: dict[str, list[str]] = {}
    with opener(fasta_path, "rt") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                seqid = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                if not seqid or seqid in chunks:
                    raise RequestError(
                        f"FASTA header on line {line_number} is empty or repeated."
                    )
                names.append(seqid)
                chunks[seqid] = []
                continue
            if not names:
                raise RequestError("FASTA must start with a '>' header line.")
            invalid = set(line.upper()) - _IUPAC
            if invalid:
                raise RequestError(
                    f"FASTA line {line_number} has non-IUPAC characters "
                    f"{sorted(invalid)}."
                )
            chunks[names[-1]].append(line)
    records = [FastaRecord(name, "".join(chunks[name])) for name in names]
    if not records:
        raise RequestError(f"FASTA has no records: {fasta_path}.")
    empty = [record.seqid for record in records if not record.sequence]
    if empty:
        raise RequestError(f"FASTA records have no bases: {empty}.")
    return records


def count_tiles(length: int) -> int:
    """Return how many model windows a contig of ``length`` bp is read in."""
    return math.ceil(length / (ANNOTATION_WINDOW_BP * ANNOTATION_KEEP_FRACTION))


def locus_window(center: int, chrom_size: int) -> tuple[int, int]:
    """Return the 1 Mb centred on ``center``, shifted to fit the chromosome.

    Args:
        center: 0-based centre.
        chrom_size: Chromosome length.

    Returns:
        (start, end), 0-based half-open; the whole chromosome when it is shorter.
    """
    start = max(0, center - ANNOTATION_WINDOW_BP // 2)
    end = min(chrom_size, start + ANNOTATION_WINDOW_BP)
    return max(0, end - ANNOTATION_WINDOW_BP), end


def predict_run_window(
    run_dir: Path, chrom_size: Callable[[str, str], int] = fetch_chrom_size
) -> LocusWindow:
    """Return the 1 Mb centred on a Predict run's window.

    Args:
        run_dir: Predict run holding ``plan.yaml``.
        chrom_size: ``(assembly, chrom) -> length``.

    Returns:
        The window to annotate.

    Raises:
        RequestError: The directory is not a Predict run, or already has an
            annotation plan.
    """
    plan_path = run_dir / "plan.yaml"
    if not plan_path.is_file():
        raise RequestError(
            f"Which Predict run? {run_dir} has no plan.yaml. Pass --locus-run with a "
            "Predict run directory, or --fasta. No locus is assumed."
        )
    if (run_dir / ANNOTATION_PLAN).exists():
        raise RequestError(
            f"{run_dir} already has {ANNOTATION_PLAN}. A run holds one annotation."
        )
    plan = read_plan(path=plan_path)
    size = chrom_size(plan.assembly, plan.chrom)
    start, end = locus_window(center=(plan.start + plan.end) // 2, chrom_size=size)
    return LocusWindow(plan.assembly, plan.chrom, start, end)


def _annotation_plan(
    request: AnnotateRequest, window: LocusWindow | None
) -> dict[str, Any]:
    """Validate the ask and return its ``annotation_plan.json`` payload.

    Args:
        request: The ask.
        window: The window, for a locus ask.

    Returns:
        The plan, unapproved.

    Raises:
        RequestError: The clade, decoder, or FASTA is invalid.
    """
    clade = clade_config(request.clade)
    if request.decoder not in ANNOTATION_DECODERS:
        raise RequestError(
            f"Unknown decoder '{request.decoder}'. Allowed: {ANNOTATION_DECODERS}."
        )
    if window is None:
        fasta = Path(str(request.fasta)).expanduser().resolve()
        records = [
            {"seqid": record.seqid, "length": len(record.sequence)}
            for record in read_fasta(fasta)
        ]
    else:
        fasta = Path(request.out_dir).resolve() / ANNOTATION_INPUT
        records = [{"seqid": window.chrom, "length": window.end - window.start}]
    return {
        "clade": clade,
        "model": ANNOTATION_MODELS[clade],
        "decoder": request.decoder,
        "fasta": str(fasta),
        "locus": None if window is None else asdict(window),
        "records": records,
        "tiles": sum(count_tiles(int(row["length"])) for row in records),
        "window_bp": ANNOTATION_WINDOW_BP,
        "keep_fraction": ANNOTATION_KEEP_FRACTION,
        "min_cds_length": MIN_CDS_LENGTH,
        "apply_biological_filters": True,
        "output": ANNOTATION_OUTPUT,
        "want_critic": True,
        "approved": False,
    }


def annotation_resolved_block(
    request: AnnotateRequest,
    run_dir: str,
    doctor: DoctorReport,
    window: LocusWindow | None = None,
    weights_cached: bool | None = None,
) -> list[str]:
    """Return the Resolved block for an Annotate ask. Writes nothing.

    Args:
        request: The ask.
        run_dir: Run directory as the user sees it.
        doctor: JAX Python discovery result.
        window: The window, for a locus ask.
        weights_cached: Whether the clade's weights are cached; None checks.

    Returns:
        Lines from the "Resolved" title to the approval question.

    Raises:
        RequestError: The clade, decoder, or FASTA is invalid.
    """
    plan = _annotation_plan(request=request, window=window)
    model = plan["model"]
    if weights_cached is None:
        weights_cached = hub_weights_cached(repo=model)
    total_bp = sum(int(row["length"]) for row in plan["records"])
    if window is None:
        input_line = (
            f"FASTA: {plan['fasta']} — {len(plan['records'])} record(s), "
            f"{total_bp:,} bp, {plan['tiles']} window(s), read from the file"
        )
        lane = "the NTv3 annotation lane, after `ntv3 report`"
    else:
        input_line = (
            f"Locus: {window.chrom}:{window.start:,}-{window.end:,} ({window.assembly}, "
            f"{plan['tiles']} window(s)) — the 1 Mb centred on the Predict run's "
            "window, fetched from UCSC with its soft mask"
        )
        lane = "the NTv3 annotation lane next to its tracks, after `report --refresh`"
    clade_rule = (
        "the clade you named"
        if plan["clade"] == request.clade.strip().lower()
        else f"the clade you named, as its model config '{plan['clade']}'"
    )
    decoder_rule = (
        "released default; say argmax or hmm to change it"
        if request.decoder == "hmm_constrained"
        else "named in your request"
    )
    if weights_cached:
        weights_line = f"Weights: {model} is in the Hugging Face cache — no download"
    else:
        weights_line = (
            f"Weights: {model} is not in the Hugging Face cache — the forward "
            "downloads 2.8 GB; the repo is gated, so it needs a Hugging Face "
            "login with approved access"
        )
    backend = doctor.candidates[-1].jax_backend if doctor.jax_python else None
    gpu_note = (
        f"; this JAX python runs on {backend}, so the forward will refuse"
        if backend not in (None, "gpu")
        else ""
    )
    return [
        "Resolved",
        input_line,
        f"Clade: {plan['clade']} — {clade_rule}",
        f"Model: {model} — the Hub model for {plan['clade']}",
        weights_line,
        f"Decoder: {plan['decoder']} — {decoder_rule}",
        (
            f"Window: {ANNOTATION_WINDOW_BP:,} bp, central "
            f"{ANNOTATION_KEEP_FRACTION:.0%} kept — fixed by the model config"
        ),
        (
            f"Filters: CDS shorter than {MIN_CDS_LENGTH} bp dropped; internal "
            "in-frame stop dropped — released defaults"
        ),
        (
            "Forward: on a GPU, after approve — about 25 GB of GPU memory per window "
            f"and about 2 minutes on one H100{gpu_note}"
        ),
        f"Output: {ANNOTATION_OUTPUT} in {run_dir} — drawn as {lane}",
        (
            "Critic: on — default; the report audits the run's gene claim. "
            "Approve with --no-critic to turn it off"
        ),
        jax_python_line(report=doctor),
        (
            "Approve to "
            + ("" if weights_cached else f"download {model} and ")
            + "run the forward, or change a line."
        ),
    ]


def plan_annotation(
    request: AnnotateRequest,
    block: list[str],
    window: LocusWindow | None = None,
) -> Path:
    """Write ``annotation_plan.json`` and its Resolved block. Loads no model.

    A locus ask also writes the window's sequence to ``annotation_input.fa``.

    Args:
        request: The ask.
        block: The Resolved block lines the user is shown.
        window: The window, for a locus ask.

    Returns:
        Path to ``annotation_plan.json``.

    Raises:
        RequestError: The clade, decoder, or FASTA is invalid.
        ArtifactError: A FASTA run directory already holds files.
    """
    plan = _annotation_plan(request=request, window=window)
    out_dir = Path(request.out_dir)
    if window is None:
        if out_dir.exists() and any(out_dir.iterdir()):
            raise ArtifactError(
                f"Refusing to use nonempty run directory: {out_dir}. "
                "Choose a new results/<run> name."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        sequence = fetch_ucsc_sequence(**asdict(window), keep_soft_mask=True)
        lines = [sequence[i : i + 80] for i in range(0, len(sequence), 80)]
        Path(plan["fasta"]).write_text(f">{window.chrom}\n" + "\n".join(lines) + "\n")
    plan_path = out_dir / ANNOTATION_PLAN
    with plan_path.open("x") as handle:
        json.dump(plan, handle, indent=2)
        handle.write("\n")
    (out_dir / ANNOTATION_RESOLVED_BLOCK).write_text("\n".join(block) + "\n")
    return plan_path


def shift_gff3(source: Path, destination: Path, window: LocusWindow) -> None:
    """Copy a window's GFF3 into chromosome coordinates.

    Args:
        source: GFF3 whose one seqid is the window's chromosome, with
            coordinates relative to the window.
        destination: GFF3 to write.
        window: The window ``source`` was annotated on.

    Raises:
        ArtifactError: A feature names another seqid.
    """
    lines: list[str] = []
    for line in source.read_text().splitlines():
        if line.startswith("##sequence-region"):
            _, seqid, first, last = line.split()
            first, last = int(first) + window.start, int(last) + window.start
            lines.append(f"##sequence-region {seqid} {first} {last}")
        elif line.startswith("#") or not line.strip():
            lines.append(line)
        else:
            cols = line.split("\t")
            if cols[0] != window.chrom:
                raise ArtifactError(
                    f"{source} names seqid {cols[0]}, not {window.chrom}."
                )
            cols[3] = str(int(cols[3]) + window.start)
            cols[4] = str(int(cols[4]) + window.start)
            lines.append("\t".join(cols))
    destination.write_text("\n".join(lines) + "\n")
