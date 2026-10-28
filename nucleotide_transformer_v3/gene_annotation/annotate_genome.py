# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Annotate a genome FASTA with a trained NTv3 annotation head, and write a GFF3
or a GTF.

The format follows ``output_path``'s extension, as Tiberius' does, and
``codingseq_path`` / ``protseq_path`` add the spliced-CDS and protein FASTAs its
``--codingseq`` and ``--protseq`` produce. See ``formats.py`` for the GTF dialect
and ``sequences.py`` for the splicing, phase and translation rules.

The production counterpart of ``annotate/annotate.py``. Same model, same window
geometry and the same decoder — that code is imported, not reimplemented — but
none of the research scaffolding: no MLflow, no ground-truth targets, no confusion
matrices, no Tiberius evaluation, and no grain dataloader. The caller hands over
FASTA files, so each contig is tokenized once in memory and its windows are slices
of that array.

Loads a model from a directory of converted weights, always the NTv3 BED-head
family. The research pipeline dispatched between three checkpoint families on
the head recorded in the saved config; only one of them was ever released, so
that dispatch does not travel.

``run`` takes a :class:`ProductionAnnotationConfig` directly:

    from nucleotide_transformer_v3.gene_annotation.annotate_genome import run
    from nucleotide_transformer_v3.gene_annotation.config import (
        ProductionAnnotationConfig,
    )

    run(ProductionAnnotationConfig(
        fasta_paths=["/data/cow/genome"],
        model_cfg="mammals",
        output_path="/data/cow.gff3.gz",
    ))

``model_cfg`` names a registered model (see ``model_cfg/``) and carries the
weights together with the window geometry they must be run at. An unregistered
model still works, with the geometry given explicitly:

    run(ProductionAnnotationConfig(
        fasta_paths=["/data/genome.fa"],
        weights_dir="/data/my_model",
        output_path="/data/out.gff3",
        sequence_length=1048576,
    ))

``cli.py`` is the same pipeline behind a Tiberius-style flag interface.
"""

import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from jax.sharding import Mesh, NamedSharding, PartitionSpec

from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
    TOKEN_TO_BASE_CODE,
    build_annotation_inference_fn,
    build_rc_token_lut,
    compute_window_geometry,
    pin_keep_target_center_fraction,
    resolve_condition_inputs,
)
from nucleotide_transformer_v3.gene_annotation.config import ProductionAnnotationConfig
from nucleotide_transformer_v3.gene_annotation.formats import resolve_output_format
from nucleotide_transformer_v3.gene_annotation.gff_utils import find_genes
from nucleotide_transformer_v3.gene_annotation.model_registry import (
    apply_model_config,
    load_model_config,
)
from nucleotide_transformer_v3.gene_annotation.sequences import (
    gene_sequences,
    merge_fasta_parts,
    render_fasta_records,
    write_fasta_part,
)
from nucleotide_transformer_v3.gene_annotation.utils import (
    assemble_runs,
    build_part_identity,
    build_provenance_lines,
    build_settings_digest,
    build_token_lut,
    build_window_row,
    describe_part_mismatch,
    is_fasta_name,
    iter_fasta_records,
    max_gene_id_in_parts,
    merge_parts,
    part_path_for,
    plan_contig_windows,
    read_part,
    resolve_fasta_paths,
    tokenize_contig,
    write_record_part,
)
from nucleotide_transformer_v3.logging import get_logger
from nucleotide_transformer_v3.pretrained import (
    load_annotation_model as load_weights,
)
from nucleotide_transformer_v3.pretrained import read_provenance
from nucleotide_transformer_v3.tokenizers import NucleotideTokenizer

logger = get_logger()


NUM_ANNOTATION_CLASSES = 15
VALID_DECODERS = ("hmm_constrained", "hmm", "argmax")


def resolve_devices(backend: str) -> list[jax.Device]:
    """Every device the pass runs on, for the configured backend.

    Args:
        backend: ``"gpu"``, ``"cpu"``, or empty to accept whatever backend JAX
            picks by default. Empty is what lets a CPU-only host and a GPU host
            share one config.

    Returns:
        The backend's devices in JAX's own order, which becomes the ``dp`` mesh
        axis, of which the first is used.

    Raises:
        ValueError: If the backend is unknown or has no visible device, which is
            otherwise a bare RuntimeError from deep inside JAX.
    """
    try:
        devices = jax.devices(backend or None)
    except RuntimeError as error:
        raise ValueError(
            f"backend={backend!r} has no visible device ({error}). "
            "backend=cpu runs but is not a practical alternative: one 1 Mb "
            "window takes 188s on 16 CPU cores against 0.8s on an H100, so a "
            "2.7 Gb genome is weeks rather than the ~45 min it takes on a GPU. "
            "It is there for tests, not for annotating a genome. See the "
            "hardware requirements in the README."
        ) from error
    if not devices:
        raise ValueError(f"backend={backend!r} has no visible device.")
    _warn_if_silently_on_cpu(backend, devices)
    return devices


def _warn_if_silently_on_cpu(backend: str, devices: list[jax.Device]) -> None:
    """Say so when a genome is about to be annotated on the CPU by accident.

    A plain ``pip install`` pulls the PyPI jaxlib, which carries no CUDA, so
    JAX falls back to the CPU on a machine with a perfectly good GPU. That
    turns ~45 minutes into ~7.5 days, and the only sign is one line JAX prints
    among a dozen other backend probes. Nothing fails; the run just looks slow.

    The driver is read from the OS rather than from JAX, because a CPU-only
    jaxlib cannot see the GPU either -- asking JAX would answer "no CUDA
    device" on exactly the machine this exists to catch.
    """
    if backend or devices[0].platform != "cpu":
        return
    driver_present = Path("/proc/driver/nvidia/version").exists()
    remedy = (
        "An NVIDIA driver is present, so this is the CPU-only jaxlib from "
        "PyPI. Reinstall with the GPU extra: pip install "
        "'nucleotide_transformer[gpu]'."
        if driver_present
        else "No NVIDIA driver is present. Pass backend=cpu to say so "
        "deliberately and silence this."
    )
    logger.warning(
        "Annotating on the CPU: one 1 Mb window takes ~188 s here against "
        f"~0.8 s on an H100, so a genome is weeks, not ~45 minutes. {remedy}"
    )


def load_annotation_model(
    config: ProductionAnnotationConfig,
    sharding: NamedSharding,
) -> tuple[Any, NucleotideTokenizer, list[str], float]:
    """Load the model from a directory of converted weights.

    One way in, unlike the research pipeline's two. The tokenizer is serialized
    beside the weights, so it is restored rather than rebuilt from a
    ``use_softmasked_vocab`` flag, and the geometry the model is *served* at
    comes from the registry rather than from the weights, which record the
    geometry they were *trained* at.

    Args:
        config: The resolved annotation configuration.
        sharding: Unused; the loader replicates onto one device. Accepted so the
            signature matches what ``run`` passes.

    Returns:
        The model in eval mode, its sequence tokenizer, the condition names its
        config declares, and the resolved keep-target-center fraction.
    """
    del sharding
    model, sequence_tokenizer = load_weights(
        config.weights_dir, use_bf16_compute=config.use_bf16_compute
    )
    provenance = read_provenance(config.weights_dir)
    _reject_weights_for_another_model(config, provenance)

    keep_target_center_fraction = pin_keep_target_center_fraction(
        model.config,
        config.keep_target_center_fraction,
        # The export substitutes the served crop into the config and records
        # the trained one here, so the config alone cannot tell them apart.
        trained=provenance.get("trained_keep_target_center_fraction"),
    )
    model.eval()
    disable_remat_for_inference(model)

    if model.config.bed_num_classes != NUM_ANNOTATION_CLASSES:
        raise ValueError(
            f"Expected bed_num_classes={NUM_ANNOTATION_CLASSES} "
            f"(IR/I0-I2/E0-E2/ASS0-ASS2/DSS0-DSS2/START/STOP), got "
            f"{model.config.bed_num_classes}. This model's head is not a "
            "genome-annotation head."
        )

    conditions_names = list(model.config.conditions_names or [])
    # The trunk builds one embedding table per conditions_vocab_size entry and
    # zips it against the conditions passed in, so a config naming fewer
    # conditions than it sizes dies inside the forward pass with a bare
    # "zip() argument 2 is longer than argument 1".
    if len(conditions_names) != len(model.config.conditions_vocab_size or []):
        raise ValueError(
            f"This model declares conditions_names={conditions_names} but "
            f"conditions_vocab_size={list(model.config.conditions_vocab_size or [])}; "
            "the two must have the same length. Its hyperparameters.json is "
            "inconsistent."
        )
    _reject_registry_conditioning_drift(
        config.model_cfg, species_conditioning_is_consumed(model.config)
    )
    return model, sequence_tokenizer, conditions_names, keep_target_center_fraction


def _reject_remote_paths(config: ProductionAnnotationConfig) -> None:
    """Refuse a remote URI on any path this package will open with ``pathlib``.

    ``resolve_fasta_paths`` already rejects one on the read side, where the
    failure is at least loud. The write side is worse: nothing opens the path
    until the end, ``pathlib`` turns ``s3://bucket/x.gff3`` into the relative
    directory ``s3:/bucket/``, and the run creates it, writes there and reports
    success -- after a genome of GPU time, having logged a destination that
    reads almost right.

    Checked here, beside the extension, so a mistyped destination costs a
    second rather than a run.

    Raises:
        ValueError: If any output path or the weights directory is a URI.
    """
    fields = (
        "output_path",
        "output_parts_dir",
        "codingseq_path",
        "protseq_path",
        "weights_dir",
    )
    for field in fields:
        value = getattr(config, field)
        if value and "://" in str(value):
            raise ValueError(
                f"Remote path {value} is not supported for {field}: this "
                "package reads and writes with pathlib and imports no cloud "
                "SDK. Write locally and copy afterwards."
            )


def _reject_weights_for_another_model(
    config: ProductionAnnotationConfig, provenance: dict[str, Any]
) -> None:
    """Fail when the weights are an export of a model other than the one asked for.

    A released name is unambiguous, but ``NTV3_ANNOTATE_WEIGHTS_DIR`` points
    at a directory the user staged, and nothing about a directory proves which
    clade is inside it. Put the plant export in a folder called ``mammals``
    and every other check passes: all three clades share their 15 classes,
    their 1 Mb geometry, their alphabet and their inert conditioning. The run
    would annotate a mammalian genome with the plant model and report success,
    which for an annotation tool is the worst available outcome -- a plausible
    GFF3 that is wrong throughout.

    The registry entry names the model, the export records which model it came
    from, and this is where the two meet.

    Args:
        config: The resolved configuration, after ``apply_model_config``.
        provenance: The export's provenance block.

    Raises:
        ValueError: When the two names disagree.
    """
    expected = config.model_name
    found = provenance.get("source_model_name", "")
    if not expected or not found or expected == found:
        return
    raise ValueError(
        f"model_cfg={config.model_cfg!r} expects {expected!r}, but the weights "
        f"in {config.weights_dir} record source_model_name={found!r}. These "
        "are an export of a different model. Unset NTV3_ANNOTATE_WEIGHTS_DIR "
        "to take the released one, point it at a directory holding the right "
        "export, or pass weights_dir directly to use these deliberately."
    )


def _reject_registry_conditioning_drift(model_cfg: str, consumed: bool) -> None:
    """Fail when a registry entry's ``species_conditioned`` contradicts its weights.

    The registry is hand-written, the checkpoint is authoritative, and the two
    disagreed on the very first entries: both declared ``species_conditioned:
    true`` for weights that consume no species. A YAML comment saying "read
    from the checkpoint" is not a check, so this is the check. It costs one
    comparison per run and turns a silent wrong required flag into a startup
    error naming the file to edit.

    Args:
        model_cfg: The registered name in use, or ``""`` when the caller passed
            a bare ``weights_dir`` and there is nothing to cross-check.
        consumed: What the checkpoint's own config implies.

    Raises:
        ValueError: When the entry and the checkpoint disagree.
    """
    if not model_cfg:
        return
    declared = load_model_config(model_cfg).species_conditioned
    if declared == consumed:
        return
    raise ValueError(
        f"model_cfg={model_cfg!r} declares species_conditioned={declared}, but "
        f"its checkpoint's own config implies {consumed}. The checkpoint is "
        f"authoritative: edit species_conditioned in model_cfg/{model_cfg}.yaml "
        "to match. Declaring conditioning that does not exist makes species a "
        "required flag that cannot change the output; declaring none where it "
        "does exist silently drops a real input."
    )


def species_conditioning_is_consumed(model_config: Any) -> bool:
    """Whether the species condition can reach this model's annotation output.

    A checkpoint that saves a species tokenizer is not necessarily conditioned.
    The mammal and plant checkpoints save ``conditions_names: ['species']`` and
    build a species embedding table, yet set ``enable_trunk_conditioning``,
    ``condition_annotation_head`` and ``condition_bigwig_head`` all ``False``;
    the insect one saves no conditions at all. The table's only reader is then
    the auxiliary ``species_logits`` head, which predicts species *from* the
    trunk embedding and which annotation never reads, so all three released
    exports drop the conditioning outright, leaving ``conditions_vocab_size``
    empty and recording what was dropped under
    ``provenance.dropped_conditioning``. Deriving "is it conditioned?" from the
    tokenizer's presence, as the registry first did, makes ``--species`` a
    required flag that changes nothing, and invites attributing output
    differences to conditioning.

    The two families gate this differently:

    * NTv3 post-training (``bed_num_classes``): per-site flags. The annotation
      output is reachable from the trunk FiLM and from the BED head, so either
      one being on counts. ``condition_bigwig_head`` is deliberately not
      consulted; the bigwig head does not feed the annotation logits.
    * LongRangeNT genome-annotation (``num_classes``): no gate. Its conditions
      flow through the trunk's adaptive norms whenever they exist.

    Args:
        model_config: The restored model's config.

    Returns:
        True when a species input can change the annotation logits.
    """
    if not getattr(model_config, "conditions_vocab_size", []):
        return False
    if hasattr(model_config, "condition_annotation_head"):
        return bool(getattr(model_config, "enable_trunk_conditioning", False)) or bool(
            model_config.condition_annotation_head
        )
    return True


def disable_remat_for_inference(model: Any) -> bool:
    """Turn rematerialisation off. Inference has no backward pass to trade for.

    ``use_remat_in_transformer`` / ``use_remat_in_convs`` are ``True`` in the
    shipped checkpoints because they were set for *training*, where remat trades
    recomputation for activation memory during backpropagation. There is no
    backpropagation here, so at best the setting is inert.

    Deliberately not a flag. Both fields are read inside ``NTv3.__call__``, so
    flipping them is enough, and the output is bit-identical: measured on a real
    1 Mb window of cow chr25, the decoded state path is unchanged in every one of
    3,072 positions. A knob whose only correct setting is off is a knob a user
    can only get wrong.

    What it buys is the compile, not the throughput. Measured on one H100-80GB
    (``9ab6998d``): steady state 0.7085s -> 0.7060s per window, which is inside
    the run-to-run noise, but **compile 67.7s -> 28.2s**. That is ~0.8% of a
    whole-genome run and most of a short one, and it is paid once per process, so
    it matters most exactly where the pipeline is used in parallel across contigs.

    Args:
        model: A loaded model whose config may carry the remat fields.

    Returns:
        Whether anything was changed, so a config family without these fields is
        distinguishable from a no-op.
    """
    config = getattr(model, "config", None)
    fields = ("use_remat_in_transformer", "use_remat_in_convs")
    present = [name for name in fields if hasattr(config, name)]
    if not present:
        return False
    was_on = [name for name in present if getattr(config, name)]
    for name in present:
        setattr(config, name, False)
    if was_on:
        logger.info(
            f"Disabled rematerialisation for inference ({', '.join(was_on)}): "
            "output is unchanged and the XLA compile is ~2.4x faster."
        )
    return bool(was_on)


def resolve_conditions(
    config: ProductionAnnotationConfig,
    conditions_names: list[str],
    data_sharding: NamedSharding,
    *,
    species_consumed: bool,
) -> tuple[list[jax.Array], list[jax.Array] | None]:
    """Encode the species condition, which no released model reads.

    Two cases. A model with no condition tokenizer needs nothing, so a supplied
    ``species`` is warned about rather than silently ignored. A model that
    declares a species condition but reads nothing from it gets a masked
    placeholder: it cannot get empty lists, because this model family embeds
    conditions with ``zip(strict=True)``.

    There is deliberately no third case. A genuinely conditioned checkpoint
    would need its species tokenizer restored and its out-of-vocabulary
    behaviour reproduced; none of the released models is one
    (``conditions_dims=[]`` and ``condition_annotation_head=False`` in all
    three), so that branch would be untestable code on a path nothing reaches.
    :func:`load_annotation_model` raises before here if one ever appears.

    Args:
        config: The run's settings.
        conditions_names: Condition names the model's config declares.
        data_sharding: Mesh to place the condition tokens on.
        species_consumed: Whether the model reads the species condition, from
            :func:`species_conditioning_is_consumed`. Not inferable from
            ``conditions_names``, which only says a tokenizer exists.

    Returns:
        The condition tokens and optional masks to call the model with.

    Raises:
        NotImplementedError: If the model declares a condition other than
            species, or genuinely consumes one.
    """
    if not conditions_names:
        if config.species:
            logger.warning(
                f"species={config.species!r} is ignored: this model does not "
                "read a species condition. It annotates any organism from "
                "sequence alone, and the output is identical whatever you pass."
            )
        return [], None

    for name in conditions_names:
        if name != "species":
            raise NotImplementedError(
                f"Unsupported condition {name!r}; only 'species' is supported."
            )

    if species_consumed:
        raise NotImplementedError(
            "This model genuinely consumes its species condition, which the "
            "released models do not and which this pipeline therefore does not "
            "implement. Annotating it here would silently drop a real input."
        )

    # A model can carry a species tokenizer and still read nothing from it. All
    # three released models do exactly that, so requiring a species would ask
    # the user to choose a value that cannot change a single base of the output.
    # The conditions cannot simply be dropped either: this model family embeds
    # them with zip(strict=True), so an empty list against a non-empty
    # conditions_vocab_size raises. A masked placeholder is the honest
    # substitute, and is the same input the masked path has always produced.
    if config.species:
        logger.warning(
            f"species={config.species!r} is ignored: this model declares a "
            "species condition, but its config disables every path that would "
            "read it. Annotating with the condition masked."
        )
    return resolve_condition_inputs(
        conditions_names,
        {},
        "",
        mask_species_condition=True,
        data_sharding=data_sharding,
    )


def annotate_record(
    *,
    record_id: str,
    sequence: str,
    config: ProductionAnnotationConfig,
    keep_target_center_fraction: float,
    inference_fn: Callable[..., dict[str, jax.Array]],
    conditions_tokens_1d: list[jax.Array],
    conditions_masks_1d: list[jax.Array] | None,
    token_lut: np.ndarray,
    n_token_id: int,
    sharding: NamedSharding,
    gene_id_offset: int,
    is_first_record: bool = True,
) -> tuple[list[str], int]:
    """Tile, decode and stitch one FASTA record into GFF3 lines.

    Returns the lines and the last gene id assigned, which the caller carries into
    the next record so ids stay unique across the whole input.

    ``sharding`` places each batch of token rows on the single mesh device;

    ``is_first_record`` only affects how the first batch's timing is reported: the
    graph is traced and compiled once per process, so labelling every record's
    first batch as a cold start makes a warm forward look like one.
    """
    contig_length = len(sequence)
    tokens, num_folded = tokenize_contig(sequence, token_lut)
    if num_folded:
        logger.info(
            f"{record_id}: {num_folded:,} base(s) "
            f"({100 * num_folded / max(contig_length, 1):.4f}%) outside the model's "
            "alphabet folded to N."
        )

    windows, center_offset = plan_contig_windows(
        sequence_length=config.sequence_length,
        keep_target_center_fraction=keep_target_center_fraction,
        contig_length=contig_length,
        chrom=record_id,
    )
    if not windows:
        logger.warning(f"{record_id}: no windows under the current geometry; skipping.")
        return [], gene_id_offset

    crop = center_offset
    states_plus: dict[int, np.ndarray] = {}
    states_minus: dict[int, np.ndarray] = {}
    base_codes_by_start: dict[int, np.ndarray] | None = (
        {} if config.apply_biological_filters else None
    )
    expected_starts = [w.keep_start for w in windows]

    logger.info(
        f"{record_id}: {contig_length:,} bp -> {len(windows)} window(s) of "
        f"{config.sequence_length:,} bp, decoding..."
    )
    start_time = time.time()
    first_batch_time: float | None = None
    batch_size = config.batch_size
    index = 0
    cursor = 0
    while cursor < len(windows):
        batch = windows[cursor : cursor + batch_size]
        cursor += len(batch)
        rows = np.stack([build_window_row(tokens, w, n_token_id) for w in batch])
        batch_tokens = jax.device_put(jnp.asarray(rows, dtype=jnp.int32), sharding)
        batch_start = time.time()
        try:
            decoded = np.asarray(
                inference_fn(batch_tokens, conditions_tokens_1d, conditions_masks_1d)[
                    config.decoder
                ]
            )
        # Matched on the exception TYPE and then on XLA's own status code, not
        # on a substring of any exception. `"out of memory" in str(error)` reads
        # an unrelated failure whose message happens to carry that text as
        # "lower batch_size", which sends a user to the one knob that will not
        # help. A device OOM arrives as JaxRuntimeError carrying
        # RESOURCE_EXHAUSTED; anything else is not this.
        except jax.errors.JaxRuntimeError as error:
            if "RESOURCE_EXHAUSTED" not in str(error):
                raise
            raise RuntimeError(
                f"Ran out of device memory on a batch of {len(batch)} window(s) of "
                f"{config.sequence_length:,} bp. batch_size is the only thing to "
                f"change: it is {config.batch_size} now, and 1 is the safe value. "
                "The output is identical at any batch size, so lowering it costs "
                "only speed."
            ) from error
        if first_batch_time is None:
            first_batch_time = time.time() - batch_start
            detail = (
                " (includes JIT compile / cache load)"
                if is_first_record
                else " (graph already compiled)"
            )
            logger.info(
                f"{record_id}: first forward pass + decode took "
                f"{first_batch_time:.1f}s{detail}."
            )
        batch_codes = (
            TOKEN_TO_BASE_CODE[rows][:, crop : rows.shape[1] - crop]
            if base_codes_by_start is not None
            else None
        )
        for i, window in enumerate(batch):
            kept = slice(window.slice_start, window.slice_stop)
            states_plus[window.keep_start] = decoded[i, kept, 0]
            states_minus[window.keep_start] = decoded[i, kept, 1]
            if base_codes_by_start is not None and batch_codes is not None:
                base_codes_by_start[window.keep_start] = batch_codes[i][kept]

        index += 1
        done = cursor
        if index % config.log_every_n_batches == 0 and done < len(windows):
            elapsed = time.time() - start_time
            rate = done / elapsed if elapsed else 0.0
            remaining = (len(windows) - done) / rate if rate else float("nan")
            logger.info(
                f"{record_id}: {min(done, len(windows))}/{len(windows)} windows, "
                f"{elapsed:.0f}s elapsed, {rate:.2f} window/s, "
                f"~{remaining:.0f}s remaining."
            )

    lines: list[str] = []
    gene_id = gene_id_offset
    for plus, minus, run_codes, run_start in assemble_runs(
        states_plus, states_minus, base_codes_by_start, expected_starts
    ):
        for states, strand in ((plus, "plus"), (minus, "minus")):
            strand_lines, gene_id = find_genes(
                states,
                strand,
                gene_id,
                chrom=record_id,
                genomic_start=run_start,
                base_codes=run_codes,
                min_cds_length=config.min_cds_length,
            )
            lines.extend(strand_lines)

    duration = time.time() - start_time
    logger.info(
        f"{record_id}: {gene_id - gene_id_offset} gene(s) from {len(windows)} "
        f"window(s) in {duration:.1f}s "
        f"({len(windows) / duration if duration else 0:.2f} window/s)."
    )
    return lines, gene_id


@dataclass(frozen=True)
class SequenceOutput:
    """One configured FASTA output and where its per-record parts are staged.

    Attributes:
        path: The FASTA to write; gzipped when it ends in ``.gz``.
        parts_dir: Staging directory, ``<path>.parts``, so a killed job keeps
            the records it finished exactly as the annotation parts do.
        translate_to_protein: Whether this output holds the translation rather
            than the spliced CDS.
    """

    path: str
    parts_dir: Path
    translate_to_protein: bool


def _same_file(a: Path, b: Path) -> bool:
    """Whether two paths name one file, as spelled or by identity."""
    return a == b or (a.exists() and b.exists() and a.samefile(b))


def resolve_sequence_outputs(
    config: ProductionAnnotationConfig,
) -> list[SequenceOutput]:
    """The FASTA outputs this run was asked for, coding first.

    Args:
        config: The run's configuration.

    Returns:
        One entry per configured path; empty when neither was configured, which
        is what makes the whole sequence stage cost nothing by default.

    Raises:
        ValueError: If a FASTA output resolves to another output or an input
            FASTA, or sits in an input directory, which the next run would
            read back as genome.
    """
    # Paths are compared resolved, so ``./out.gff3``, ``a//out.gff3`` and a
    # symlink count as the same file. Each of these collisions is silent: the
    # FASTAs are merged last, so the write lands on top, exit code 0.
    resolved = {
        label: Path(path).resolve()
        for label, path in (
            ("output_path", config.output_path),
            ("codingseq_path", config.codingseq_path),
            ("protseq_path", config.protseq_path),
        )
        if path
    }
    inputs = [Path(entry).resolve() for entry in config.fasta_paths]
    for label in ("codingseq_path", "protseq_path"):
        target = resolved.get(label)
        if target is None:
            continue
        for other in ("output_path", "codingseq_path"):
            if other != label and resolved.get(other) == target:
                raise ValueError(
                    f"{label} and {other} are both {str(target)!r}, so the "
                    "sequence write would silently overwrite it. Give them "
                    "different paths."
                )
        for source in inputs:
            # samefile also catches hard links and case-insensitive filesystems.
            if _same_file(target, source):
                raise ValueError(
                    f"{label} {str(target)!r} is an input FASTA and would be "
                    "overwritten. Write it elsewhere."
                )
            if (
                source.is_dir()
                and _same_file(target.parent, source)
                and is_fasta_name(target.name)
            ):
                raise ValueError(
                    f"{label} {str(target)!r} is a FASTA inside the input "
                    f"directory {str(source)!r}, so the next run would read it "
                    "as genome. Write it elsewhere."
                )
    return [
        SequenceOutput(
            path=path,
            parts_dir=Path(f"{path}.parts"),
            translate_to_protein=to_protein,
        )
        for path, to_protein in (
            (config.codingseq_path, False),
            (config.protseq_path, True),
        )
        if path
    ]


def write_sequence_parts(
    outputs: list[SequenceOutput],
    record_id: str,
    lines: list[str],
    sequence: str,
) -> None:
    """Stage one record's coding and protein FASTA parts.

    The spliced CDS is extracted once and reused for both outputs, so asking for
    the protein as well as the coding sequence costs one translation rather than
    a second pass over the record.

    Args:
        outputs: The configured FASTA outputs; nothing happens when empty.
        record_id: The record's seqid.
        lines: The record's GFF3 feature lines.
        sequence: The record's sequence, exactly as it was annotated.
    """
    if not outputs:
        return
    genes = gene_sequences(lines, sequence)
    for output in outputs:
        write_fasta_part(
            output.parts_dir,
            record_id,
            render_fasta_records(genes, output.translate_to_protein),
        )


def run(config: ProductionAnnotationConfig) -> None:
    """Load the checkpoint, annotate every input record, write the outputs.

    The annotation is one GFF3 or one GTF, chosen by ``output_path``'s
    extension, plus the optional spliced-CDS and protein FASTAs.
    """
    # Named models resolve to a checkpoint plus the geometry it must be run at,
    # before anything reads either. Both surfaces (this script and cli.py) come
    # through here, so a name means the same thing whichever was used.
    config = apply_model_config(config)
    if config.decoder not in VALID_DECODERS:
        raise ValueError(f"decoder={config.decoder!r} is not one of {VALID_DECODERS}.")
    if not config.weights_dir:
        raise ValueError(
            "weights_dir is required: pass --model_cfg <name>, or point it at a "
            "directory of converted weights."
        )
    if not config.output_path:
        raise ValueError("output_path is required.")
    # Beside the checks above for the same reason the extension is resolved
    # early: both of these otherwise die on the FIRST BATCH, which on a
    # genome-scale job is after the model load and the JIT compile have already
    # been paid for. batch_size=0 reaches np.stack([]) ("need at least one array
    # to stack") and log_every_n_batches=0 reaches a modulo by zero.
    if config.batch_size < 1:
        raise ValueError(
            f"batch_size={config.batch_size} must be at least 1. The output is "
            "invariant to it, so it trades only memory against speed."
        )
    if config.log_every_n_batches < 1:
        raise ValueError(
            f"log_every_n_batches={config.log_every_n_batches} must be at least "
            "1. Pass a large value to log rarely rather than 0 to log never."
        )
    _reject_remote_paths(config)
    # Resolved before anything expensive runs, so a mistyped extension fails in
    # a second rather than after a genome's worth of forward passes.
    output_format = resolve_output_format(config.output_path)

    fasta_paths = resolve_fasta_paths(list(config.fasta_paths))
    sequence_outputs = resolve_sequence_outputs(config)
    logger.info(f"Annotating {len(fasta_paths)} FASTA file(s) -> {config.output_path}")

    # One mesh, one sharding.
    # split along its leading axis across every device, while the params are
    # replicated so each device holds the whole model.
    # One device, deliberately. Sharding the batch across several was
    # implemented and then removed: it bought no measured speedup, it could
    # only ever be verified on simulated CPU devices here, and it required
    # padding a ragged final batch with rows whose decode had to be proven
    # unreachable. Annotating a genome is embarrassingly parallel across
    # contigs, so anyone with several accelerators is better served by several
    # jobs than by a fan-out inside one.
    devices = resolve_devices(config.backend)[:1]
    mesh = Mesh(np.array(devices), ("dp",))
    sharding = NamedSharding(mesh, PartitionSpec())

    (
        model,
        sequence_tokenizer,
        conditions_names,
        keep_target_center_fraction,
    ) = load_annotation_model(config, sharding)

    # Whether the species input can reach the output at all. Derived from the
    # model's own config rather than from the presence of a tokenizer, which is
    # the distinction the registry originally got wrong.
    # getattr because the golden tests drive run() with a stub head that has no
    # .config; every real model has one, and the helper reads None as "nothing
    # to consume", which is right for a stub that takes no conditions either.
    species_consumed = species_conditioning_is_consumed(getattr(model, "config", None))

    conditions_tokens_1d, conditions_masks_1d = resolve_conditions(
        config,
        conditions_names,
        sharding,
        species_consumed=species_consumed,
    )

    center_offset = compute_window_geometry(
        config.sequence_length, keep_target_center_fraction
    )
    # Always single-strand: the BED head predicts one strand and the minus
    # strand comes from a second pass on the reverse complement.
    rc_token_lut = jnp.asarray(build_rc_token_lut(sequence_tokenizer))
    inference_fn = build_annotation_inference_fn(
        model,
        decoders=[config.decoder],
        center_offset=center_offset,
        rc_token_lut=rc_token_lut,
        decode_chunk_len=config.decode_chunk_len,
    )

    token_lut, n_token_id = build_token_lut(
        sequence_tokenizer, sequence_tokenizer.include_softmasked_tokens
    )
    logger.info(
        f"Window {config.sequence_length:,} bp, keep fraction "
        f"{keep_target_center_fraction}, center offset {center_offset:,}, "
        f"batch size {config.batch_size}, "
        f"{devices[0].platform} device(s), "
        f"single-strand, decoder {config.decoder}, "
        f"bf16={config.use_bf16_compute}."
    )

    # Built before the records rather than after them, because a part's identity
    # is keyed on it: what the header tells a recipient the file means is exactly
    # what decides whether a part from an earlier job may be reused here.
    provenance = build_provenance_lines(
        # The released name, not the directory: the header exists to tell a
        # recipient which model produced the file, and it is hashed into the
        # resume digest, so keying it on a local path both says nothing and
        # re-annotates a whole genome whenever the weights are moved or
        # bind-mounted elsewhere.
        checkpoint_path=config.model_name or config.weights_dir,
        model_name=config.model_name or config.model_cfg,
        sequence_length=config.sequence_length,
        keep_target_center_fraction=keep_target_center_fraction,
        decoder=config.decoder,
        apply_biological_filters=config.apply_biological_filters,
        min_cds_length=config.min_cds_length,
        use_bf16_compute=config.use_bf16_compute,
        species=config.species,
        mask_species_condition=config.mask_species_condition,
        species_conditioning_consumed=species_consumed,
        # The resolved platform, not config.backend: that is "" when the caller
        # accepts whatever JAX picks, which is the case where recording the
        # request instead of the outcome would be most misleading.
        backend=devices[0].platform,
    )
    # Fold species into the digest only when the model reads it. Otherwise two
    # runs differing only in --species produce byte-identical output, and
    # digesting the flag would make resume refuse to reuse a perfectly valid
    # part over a setting that cannot have changed it.
    settings_digest = build_settings_digest(
        provenance,
        species=config.species or "" if species_consumed else "",
        mask_species_condition=(
            config.mask_species_condition if species_consumed else False
        ),
    )

    # Each record's lines go straight to their own part file, so a whole genome
    # is never held in memory and a killed job keeps everything it finished.
    parts_dir = Path(config.output_parts_dir or f"{config.output_path}.parts")
    # Every record this job covers, annotated or resumed. Only these are merged,
    # so a leftover part from an earlier job cannot leak into the output.
    covered_records: list[str] = []
    # Length per covered record, for the GFF3's ``##sequence-region`` header.
    record_lengths: dict[str, int] = {}
    # Seeded from whatever is already on disk, on EVERY run and not only a
    # resumed one. The ids a job hands out have to avoid colliding with the ids
    # in parts it is not touching, and `--resume` is not what decides whether
    # such parts exist: nothing ever clears this directory, and a non-resume run
    # rewrites only the records it covers. So
    #
    #   1. --fasta chrA chrB          -> parts chrA (genes 1-10), chrB (11-20)
    #   2. --fasta chrB (no resume)   -> chrB rewritten as 1-12; chrA still 1-10
    #   3. --fasta chrA chrB --resume -> both parts validate and are skipped,
    #                                    and the merged GFF3 carries ID=gene_1_p
    #                                    twice
    #
    # which is exactly the duplication this seeding exists to prevent, and its
    # consequence is the one `max_gene_id_in_parts` documents: our evaluator
    # drops the second `gene` row and appends its CDS to the first, so two genes
    # on two contigs become one gene spanning both.
    #
    # Seeding unconditionally is free. Ids are opaque, only `covered_records` are
    # merged, and a gap in the numbering costs nothing.
    gene_id = max_gene_id_in_parts(parts_dir)
    first_gene_id = gene_id
    num_records = 0
    num_skipped = 0
    num_invalidated = 0
    total_bp = 0
    job_start = time.time()
    for record_id, sequence in iter_fasta_records(fasta_paths):
        covered_records.append(record_id)
        # Recorded for every covered record, resumed ones included: the FASTA is
        # read either way, so a resumed run still knows the contig's extent.
        record_lengths[record_id] = len(sequence)
        identity = build_part_identity(sequence, settings_digest)
        lines: list[str] = []
        resumed = False
        if config.resume and part_path_for(parts_dir, record_id).exists():
            stored, part_lines = read_part(parts_dir, record_id)
            mismatch = describe_part_mismatch(stored, identity)
            if mismatch:
                # Reusing it would put features describing one input under a
                # provenance header describing another, so the part is not
                # finished work for this job and the record is done again.
                logger.warning(
                    f"{record_id}: the part staged in {parts_dir} cannot be "
                    f"resumed because {mismatch}; re-annotating this record and "
                    "replacing the part."
                )
                num_invalidated += 1
            else:
                resumed = True
                lines = part_lines
        if resumed:
            logger.info(f"{record_id}: part already present, skipping (resume).")
            num_skipped += 1
        else:
            lines, gene_id = annotate_record(
                record_id=record_id,
                sequence=sequence,
                config=config,
                keep_target_center_fraction=keep_target_center_fraction,
                inference_fn=inference_fn,
                conditions_tokens_1d=conditions_tokens_1d,
                conditions_masks_1d=conditions_masks_1d,
                token_lut=token_lut,
                n_token_id=n_token_id,
                sharding=sharding,
                gene_id_offset=gene_id,
                is_first_record=num_records == 0,
            )
            write_record_part(parts_dir, record_id, lines, identity=identity)
            num_records += 1
            total_bp += len(sequence)
        # Derived from the record's own features and its own sequence, so they
        # describe exactly the genes in the annotation. A resumed record reads
        # its features back from the part it skipped, which is also what lets a
        # run that only now asks for them fill in the records it already did.
        #
        # Sorted here so the FASTAs do not depend on WHETHER a record was
        # resumed. A freshly annotated record arrives in decode order (per run,
        # plus strand then minus); a resumed one arrives from `write_record_part`
        # already sorted by start. Without this the same job, run straight
        # through or resumed, emits byte-different codingseq/protseq FASTAs over
        # the same records with the same content in a different order. Same key
        # as `write_record_part`, so the sequence outputs get the ordering
        # guarantee the merged GFF3 already has.
        write_sequence_parts(
            sequence_outputs,
            record_id,
            sorted(lines, key=lambda line: int(line.split("\t")[3])),
            sequence,
        )

    num_features = merge_parts(
        config.output_path,
        parts_dir,
        provenance,
        covered_records,
        output_format=output_format,
        record_lengths=record_lengths,
    )
    for output in sequence_outputs:
        merge_fasta_parts(output.path, output.parts_dir, covered_records)

    duration = time.time() - job_start
    skipped = f", {num_skipped} skipped (resume)" if num_skipped else ""
    invalidated = (
        f", {num_invalidated} re-annotated (stale part)" if num_invalidated else ""
    )
    logger.info(
        f"Annotated {num_records} record(s){skipped}{invalidated} / {total_bp:,} bp "
        f"-> {gene_id - first_gene_id} new gene(s), {num_features} "
        f"{output_format.upper()} features in {duration:.1f}s "
        f"({total_bp / duration / 1e6 if duration else 0:.2f} Mbp/s)."
    )


if __name__ == "__main__":
    # `run` takes a config, and nothing injects one now that Hydra is gone, so
    # this module has no entry point of its own. Defer to the CLI, which builds
    # the config and then calls the same `run`.
    from nucleotide_transformer_v3.gene_annotation import cli

    raise SystemExit(cli.main())
