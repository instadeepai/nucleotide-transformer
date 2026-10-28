"""Write the planning questionnaire to disk and refuse to guess."""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path

import yaml

from ntv3_tools.disrupt import Snv
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.registry import DEFAULT_CHECKPOINT, assembly_for_species

REQUIRED_FIELDS = (
    "species",
    "assembly",
    "gene",
    "chrom",
    "start",
    "end",
    "tissue",
    "tracks",
    "device",
    "window_bp",
    "analyses",
    "want_report",
    "want_notebook",
)
# From this length on, a CPU run is slow and Saliency's memory can run out.
CPU_WARN_WINDOW_BP = 16_384
# NTv3_tracks is released at this length.
MAX_WINDOW_BP = 1_048_576
# Above this length Saliency gets a memory warning on either device: it ran at
# 131,072 bp on an 80 GB GPU and ran out of memory at 1,048,576 bp.
SALIENCY_WARN_BP = 131_072
# What a plan may ask for. Everything else is outside the skill. Annotate has
# its own command, `ntv3 annotate`, and is not a plan analysis.
ALLOWED_ANALYSES = ("predict", "saliency", "compare")
EXCLUDED_ANALYSES = (
    "integrated_gradients",
    "smoothgrad",
    "ism",
    "attention",
    "knockout",
    "design",
    "ledidi",
    "motif",
)


@dataclass(frozen=True)
class SaliencyLogoConfig:
    """Configuration for each static saliency sequence logo.

    Attributes:
        window_bp: Number of bases shown in the sequence logo.
        center_bp: Optional absolute genomic coordinate. A null value selects
            each track's strongest attribution peak.
    """

    window_bp: int = 96
    center_bp: int | None = None


@dataclass(frozen=True)
class IntegratedGradientsConfig:
    """Bounded Integrated Gradients runtime configuration."""

    n_steps: int = 16


@dataclass(frozen=True)
class SmoothGradConfig:
    """Bounded SmoothGrad runtime configuration."""

    n_samples: int = 8
    noise_sigma: float = 0.15
    seed: int = 0


@dataclass(frozen=True)
class IsmConfig:
    """Bounded ISM runtime configuration using genomic coordinates."""

    chunk_size: int = 512
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class AttributionConfig:
    """Runtime settings for the public attribution methods."""

    integrated_gradients: IntegratedGradientsConfig = field(
        default_factory=IntegratedGradientsConfig
    )
    smoothgrad: SmoothGradConfig = field(default_factory=SmoothGradConfig)
    ism: IsmConfig = field(default_factory=IsmConfig)


@dataclass(frozen=True)
class AnalysisPlan:
    """Answers collected before any model run.

    Attributes:
        species: Model species string (e.g. "human"). Never defaulted silently.
        assembly: Genome assembly (e.g. "hg38").
        gene: Gene symbol or explicit locus label.
        chrom: Chromosome (e.g. "chr11"). Required; do not invent it.
        start: 0-based window start.
        end: 0-based window end.
        tissue: Tissue or cell type (e.g. "K562").
        tracks: Requested track aliases or ENCODE accessions.
        device: "cpu", "cuda", "mps", or "auto".
        window_bp: Input window length the user approved.
        want_report: Offer a collaborator report.
        want_notebook: Offer an editable notebook.
        want_critic: Launch a critic subagent before the report.
        ucsc_tracks: UCSC ``track=`` names to overlay.
        analyses: Named analyses the user approved.
        model: Internal JAX checkpoint key.
        attribution_tracks: Requested tracks that receive attribution analysis.
        saliency_logo: Static sequence-logo crop configuration.
        attribution: Settings for the refused attribution methods.
        variant: The one SNV that Compare scores, 1-based, or None.
        approved: Set only by ``ntv3 approve`` after the Resolved block exists.
    """

    species: str
    assembly: str
    gene: str
    chrom: str
    start: int
    end: int
    tissue: str
    tracks: list[str]
    device: str
    window_bp: int
    want_report: bool
    want_notebook: bool
    want_critic: bool
    ucsc_tracks: list[str]
    analyses: list[str]
    model: str = DEFAULT_CHECKPOINT
    attribution_tracks: list[str] = field(default_factory=list)
    saliency_logo: SaliencyLogoConfig = field(default_factory=SaliencyLogoConfig)
    attribution: AttributionConfig = field(default_factory=AttributionConfig)
    variant: Snv | None = None
    approved: bool = False


def run_dir_name(plan: AnalysisPlan) -> str:
    """Return the run directory name, built only from resolved fields.

    ``<assembly>_<chrom>_<start>_<end>_<analyses>``, analyses joined with ``+``
    in plan order. No gene, tissue, or free text.

    Args:
        plan: A validated plan.

    Returns:
        The directory name under ``results/``.
    """
    analyses = "+".join(plan.analyses)
    return f"{plan.assembly}_{plan.chrom}_{plan.start}_{plan.end}_{analyses}"


def refuse_excluded_analyses(analyses: Iterable[str]) -> None:
    """Refuse any analysis outside the skill, before any fetch or model load.

    ``validate_plan`` calls it on a plan's analyses, and the CLI on the
    knockout, design, and attention commands.

    Args:
        analyses: Analysis or command names.

    Raises:
        RequestError: An excluded or unknown analysis was named.
    """
    names = list(analyses)
    excluded = [name for name in EXCLUDED_ANALYSES if name in names]
    if excluded:
        raise RequestError(
            f"{', '.join(excluded)} is outside this skill. Drop it from --analyses. "
            f"Allowed: {', '.join(ALLOWED_ANALYSES)}."
        )
    unknown = sorted(set(names) - set(ALLOWED_ANALYSES))
    if unknown:
        hint = (
            " Annotate has its own command: `ntv3 annotate`."
            if "annotate" in unknown
            else ""
        )
        raise RequestError(
            f"Unsupported analyses: {unknown}. Allowed: {', '.join(ALLOWED_ANALYSES)}."
            + hint
        )


def require_window(window_bp: int) -> None:
    """Refuse a window that is not a multiple of 128 or exceeds ``MAX_WINDOW_BP``."""
    if window_bp % 128 != 0:
        raise RequestError(f"window_bp={window_bp} must be divisible by 128.")
    if window_bp > MAX_WINDOW_BP:
        raise RequestError(
            f"A {window_bp:,} bp window exceeds the released {MAX_WINDOW_BP:,} bp "
            "of NTv3_tracks. Pass a shorter --window-bp. No length is assumed."
        )


def cpu_window_warning(window_bp: int, analyses: list[str]) -> str | None:
    """Return the warning for a CPU window of ``CPU_WARN_WINDOW_BP`` or more, or None."""
    if window_bp < CPU_WARN_WINDOW_BP:
        return None
    saliency = (
        "; Saliency needs several times more memory and can exhaust this machine's RAM"
        if "saliency" in analyses and window_bp > SALIENCY_WARN_BP
        else ""
    )
    return (
        f"CPU warning: a {window_bp:,} bp window on CPU is slow and needs a lot of "
        f"memory, more as the window grows{saliency}. A GPU host runs it faster."
    )


def gpu_saliency_warning(window_bp: int, analyses: list[str]) -> str | None:
    """Return the warning for GPU Saliency above ``SALIENCY_WARN_BP``, or None."""
    if "saliency" not in analyses or window_bp <= SALIENCY_WARN_BP:
        return None
    return (
        f"GPU warning: Saliency over a {window_bp:,} bp window may run out of GPU "
        "memory."
    )


def _validate_variant(plan: AnalysisPlan) -> None:
    """Raise unless Compare and its one in-window SNV come together."""
    if ("compare" in plan.analyses) != (plan.variant is not None):
        raise RequestError(
            "Compare needs exactly one variant (--variant CHROM POS REF ALT), and "
            "a variant needs compare in --analyses."
        )
    variant = plan.variant
    if variant is None:
        return
    if "predict" not in plan.analyses:
        raise RequestError("Compare requires predict.")
    if variant.chrom != plan.chrom or not plan.start < variant.position <= plan.end:
        raise RequestError(
            f"Variant {variant.chrom}:{variant.position} (1-based) is outside the "
            f"window {plan.chrom}:{plan.start}-{plan.end}."
        )
    ref, alt = variant.ref.upper(), variant.alt.upper()
    if len(ref) != 1 or len(alt) != 1 or not {ref, alt} <= set("ACGT") or ref == alt:
        raise RequestError("Compare needs one SNV: two different A/C/G/T alleles.")


def validate_plan(plan: AnalysisPlan) -> None:
    """Raise if any required field is empty.

    Args:
        plan: Collected planning answers.

    Raises:
        RequestError: A required field is missing or inconsistent.
    """
    payload = asdict(plan)
    missing = [
        name
        for name in REQUIRED_FIELDS
        if payload.get(name) is None
        or payload.get(name) == ""
        or payload.get(name) == []
    ]
    if missing:
        raise RequestError(f"Which {', '.join(missing)}? Pass them. None is assumed.")
    if plan.species.lower() in {"", "unknown", "guess"}:
        raise RequestError("Which species? Pass --species. No species is assumed.")
    try:
        expected_assembly = assembly_for_species(species=plan.species)
    except KeyError as exc:
        raise RequestError(
            f"Which species? This tool predicts human (hg38) or mouse (mm10), not "
            f"'{plan.species}'. Pass --species. No species is assumed."
        ) from exc
    if plan.assembly != expected_assembly:
        raise RequestError(
            f"Which assembly? {plan.species} uses {expected_assembly}, not "
            f"{plan.assembly}. Pass --assembly. No assembly is assumed."
        )
    if plan.start < 0 or plan.end <= plan.start:
        raise RequestError(
            f"Invalid 0-based half-open coordinates: "
            f"{plan.chrom}:{plan.start}-{plan.end}."
        )
    coordinate_length = plan.end - plan.start
    if coordinate_length != plan.window_bp:
        raise RequestError(
            f"window_bp={plan.window_bp} does not match the coordinate span "
            f"({coordinate_length} bp)."
        )
    require_window(window_bp=plan.window_bp)
    refuse_excluded_analyses(analyses=plan.analyses)
    unknown_attribution_tracks = sorted(set(plan.attribution_tracks) - set(plan.tracks))
    if unknown_attribution_tracks:
        raise RequestError(
            "Attribution tracks must be selected prediction tracks. "
            f"Unknown: {unknown_attribution_tracks}."
        )
    if plan.attribution_tracks and "saliency" not in plan.analyses:
        raise RequestError("Attribution tracks require saliency.")
    if "saliency" in plan.analyses and "predict" not in plan.analyses:
        raise RequestError("Saliency requires predict.")
    _validate_variant(plan=plan)
    if not 2 <= plan.attribution.integrated_gradients.n_steps <= 256:
        raise RequestError("attribution.integrated_gradients.n_steps must be 2-256.")
    smoothgrad = plan.attribution.smoothgrad
    if not 1 <= smoothgrad.n_samples <= 128:
        raise RequestError("attribution.smoothgrad.n_samples must be 1-128.")
    if not 0.0 < smoothgrad.noise_sigma <= 1.0:
        raise RequestError("attribution.smoothgrad.noise_sigma must be in (0, 1].")
    ism = plan.attribution.ism
    if not 1 <= ism.chunk_size <= 4096:
        raise RequestError("attribution.ism.chunk_size must be 1-4096.")
    if plan.saliency_logo.window_bp <= 0:
        raise RequestError("saliency_logo.window_bp must be positive.")
    if plan.saliency_logo.window_bp > plan.window_bp:
        raise RequestError(
            f"saliency_logo.window_bp={plan.saliency_logo.window_bp} exceeds "
            f"the approved {plan.window_bp} bp locus window."
        )
    center_bp = plan.saliency_logo.center_bp
    if center_bp is not None and not plan.start <= center_bp < plan.end:
        raise RequestError(
            f"saliency_logo.center_bp={center_bp} is outside the approved "
            f"half-open locus [{plan.start}, {plan.end})."
        )
    if not plan.want_report and not plan.want_notebook:
        raise RequestError("Select a report, a notebook, or both.")


def write_plan(plan: AnalysisPlan, out_dir: str) -> Path:
    """Write plan.yaml under the run directory.

    Args:
        plan: Validated planning answers.
        out_dir: Run output directory.

    Returns:
        Path to plan.yaml.
    """
    validate_plan(plan=plan)
    root = Path(out_dir)
    if root.exists() and any(root.iterdir()):
        raise ArtifactError(
            f"Refusing to use nonempty run directory: {root}. "
            "Choose a new results/<run> name."
        )
    root.mkdir(parents=True, exist_ok=True)
    path = root / "plan.yaml"
    path.write_text(yaml.safe_dump(asdict(plan), sort_keys=False))
    return path


def read_plan(path: str | Path) -> AnalysisPlan:
    """Read and validate an approved analysis plan.

    Args:
        path: Path to plan.yaml.

    Returns:
        Validated immutable plan.

    Raises:
        RequestError: The file is missing, malformed, or invalid.
    """
    plan_path = Path(path)
    if not plan_path.exists():
        raise RequestError(f"Plan does not exist: {plan_path}.")
    payload = yaml.safe_load(plan_path.read_text())
    if not isinstance(payload, dict):
        raise RequestError(f"Plan is not a YAML mapping: {plan_path}.")
    logo_payload = payload.pop("saliency_logo", None)
    attribution_payload = payload.pop("attribution", None)
    variant_payload = payload.pop("variant", None)
    try:
        if logo_payload is None:
            saliency_logo = SaliencyLogoConfig()
        elif isinstance(logo_payload, dict):
            saliency_logo = SaliencyLogoConfig(**logo_payload)
        else:
            raise TypeError("saliency_logo must be a YAML mapping")
        if attribution_payload is None:
            attribution = AttributionConfig()
        elif isinstance(attribution_payload, dict):
            attribution = AttributionConfig(
                integrated_gradients=IntegratedGradientsConfig(
                    **(attribution_payload.get("integrated_gradients") or {})
                ),
                smoothgrad=SmoothGradConfig(
                    **(attribution_payload.get("smoothgrad") or {})
                ),
                ism=IsmConfig(**(attribution_payload.get("ism") or {})),
            )
        else:
            raise TypeError("attribution must be a YAML mapping")
        plan = AnalysisPlan(
            **payload,
            saliency_logo=saliency_logo,
            attribution=attribution,
            variant=None if variant_payload is None else Snv(**variant_payload),
        )
    except TypeError as exc:
        raise RequestError(f"Plan fields are invalid: {exc}") from exc
    validate_plan(plan=plan)
    return plan
