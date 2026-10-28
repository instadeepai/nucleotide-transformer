"""JAX NTv3 tracks-model load and forward through this repo's public loader.

``nucleotide_transformer_v3.pretrained.load_tracks_model("NTv3_tracks")`` is
the only load path, called as in ``docs/nucleotide_transformer_v3.md`` and
``notebooks/nucleotide_transformer_v3/inference_tracks.ipynb``. Every forward
passes ``track_indices``, so the 7,362-track head is never materialised.

The public model takes tokens, not ``inputs_embeds``. Saliency differentiates
a zero shift on the token embedding instead (``token_gradients``). The other
one-hot forwards (Integrated Gradients, SmoothGrad, ISM, design, knockout)
raise ``NotImplementedError``.

This module is imported only by the JAX science entry points.
"""

import logging
from dataclasses import dataclass
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx

from ntv3_tools.errors import RequestError
from ntv3_tools.registry import (
    DEFAULT_CHECKPOINT,
    assembly_for_species,
    resolve_model,
)
from ntv3_tools.sequence import CENTER_FRACTION, prediction_region
from nucleotide_transformer_v3.pretrained import load_tracks_model

logger = logging.getLogger("ntv3_tools")

_NO_ONE_HOT_FORWARD = (
    "The public NTv3 tracks model has no inputs_embeds path, so a one-hot "
    "forward (Integrated Gradients, SmoothGrad, ISM, Ledidi design, knockout) "
    "cannot run on it."
)


@dataclass
class LoadedModel:
    """A loaded NTv3 tracks model.

    Attributes:
        model: Flax NNX module on ``device``.
        tokenizer: Sequence tokenizer (soft-mask aware).
        condition_tokenizers: Condition tokenizers (species is index 0).
        config: Model config (source of track names).
        device: JAX device the weights sit on.
        checkpoint: Released model name passed to ``load_tracks_model``.
    """

    model: Any
    tokenizer: Any
    condition_tokenizers: list[Any]
    config: Any
    device: Any
    checkpoint: str


def as_numpy(value: Any) -> np.ndarray:
    """Convert a JAX / NumPy array-like value to ``np.ndarray``.

    Args:
        value: Array, scalar, or nested array-like.

    Returns:
        A NumPy array.
    """
    return np.asarray(value)


def resolve_device(device: str) -> Any:
    """Resolve ``auto`` / ``cuda`` / ``gpu`` / ``cpu`` to a JAX device.

    Args:
        device: ``auto`` | ``cpu`` | ``cuda`` | ``gpu`` | ``mps``.

    Returns:
        A JAX device. ``mps`` maps to CPU (no JAX MPS backend here).

    Raises:
        RequestError: ``cuda`` or ``gpu`` was approved and JAX sees no CUDA
            device. It never falls back to the CPU.
    """
    requested = device.lower()
    if requested == "mps":
        logger.info("JAX has no MPS path. Using CPU.")
        requested = "cpu"
    if requested in {"cuda", "gpu"}:
        try:
            return jax.devices("cuda")[0]
        except RuntimeError as exc:
            raise RequestError(
                f"--device {device} was approved, but JAX in this Python sees no "
                "CUDA device. Re-plan with --device cpu, or run on a GPU host with "
                "the gpu extra. Nothing was loaded."
            ) from exc
    if requested == "cpu":
        return jax.devices("cpu")[0]
    try:
        return jax.devices("cuda")[0]
    except RuntimeError:
        return jax.devices("cpu")[0]


def load_model(
    model_name: str = DEFAULT_CHECKPOINT, device: str = "cpu"
) -> LoadedModel:
    """Load the released NTv3 tracks model.

    Args:
        model_name: ``NTv3_tracks`` (``InstaDeepAI/NTv3_tracks`` is accepted).
        device: ``auto`` | ``cpu`` | ``cuda`` | ``gpu``.

    Returns:
        The loaded model bundle on the chosen device.
    """
    checkpoint = resolve_model(name=model_name)
    jax_device = resolve_device(device=device)
    logger.info(f"Loading {checkpoint} on {jax_device}")
    model, tokenizer, condition_tokenizers = load_tracks_model(checkpoint)
    nnx.update(model, jax.device_put(nnx.state(model), jax_device))
    return LoadedModel(
        model=model,
        tokenizer=tokenizer,
        condition_tokenizers=list(condition_tokenizers),
        config=model.config,
        device=jax_device,
        checkpoint=checkpoint,
    )


def bigwig_track_names(loaded: LoadedModel, species: str) -> list[str]:
    """Return the ordered bigwig-track accessions for a species assembly.

    Args:
        loaded: A loaded model bundle.
        species: Species string (e.g. ``human``).

    Returns:
        Ordered list of track accessions; index i is head column i.
    """
    assembly = assembly_for_species(species=species)
    names_by_assembly = loaded.config.bigwigs_per_file_assembly
    if assembly not in names_by_assembly:
        raise KeyError(
            f"Assembly '{assembly}' is not in cfg.bigwigs_per_file_assembly. "
            f"Available: {list(names_by_assembly)}."
        )
    return list(names_by_assembly[assembly])


def bed_element_names(loaded: LoadedModel) -> list[str]:
    """Return the ordered BED genomic-element names.

    Args:
        loaded: A loaded model bundle.

    Returns:
        Ordered list of element names (e.g. ``protein_coding_gene``).
    """
    return list(loaded.config.bed_elements_names)


def center_fraction(loaded: LoadedModel) -> float:
    """Return the fraction of the window the track head predicts.

    Args:
        loaded: A loaded model bundle.

    Returns:
        ``keep_target_center_fraction`` from the config, default 0.375.
    """
    return float(getattr(loaded.config, "keep_target_center_fraction", CENTER_FRACTION))


def species_ids(loaded: LoadedModel, species: str, batch_size: int = 1) -> Any:
    """Return the species condition ids, one per sequence.

    Args:
        loaded: A loaded model bundle.
        species: Species string (e.g. ``human``). ``token_to_id`` raises on a
            name the model does not know.
        batch_size: Number of sequences in the batch.

    Returns:
        Integer array of shape (batch_size,) on ``loaded.device``.
    """
    species_id = loaded.condition_tokenizers[0].token_to_id(species)
    return jax.device_put(
        jnp.asarray([species_id] * batch_size, dtype=jnp.int32), loaded.device
    )


def run_forward(
    loaded: LoadedModel,
    sequence: str,
    species: str,
    track_indices: list[int] | None = None,
    output_attentions: bool = False,
) -> dict[str, np.ndarray]:
    """Tokenize one sequence, condition on species, and run the tracks forward.

    Args:
        loaded: A loaded model bundle.
        sequence: DNA sequence, soft mask kept, length a multiple of 128.
        species: Species string to condition on (e.g. ``human``).
        track_indices: Head columns to compute, positions in
            ``config.bigwigs_per_file_assembly[assembly]``, shape (n_tracks,).
        output_attentions: Not supported by the public model.

    Returns:
        Output mapping with NumPy arrays; ``bigwig_tracks_logits`` has shape
        (1, int(0.375 * len(sequence)), n_tracks).

    Raises:
        ValueError: No track index was given.
        NotImplementedError: ``output_attentions`` was requested.
    """
    if output_attentions:
        raise NotImplementedError(
            "The public NTv3 tracks model does not return attention maps."
        )
    if not track_indices:
        raise ValueError("At least one approved track index is required.")
    tokens = jax.device_put(
        jnp.asarray(loaded.tokenizer.batch_np_tokenize([sequence]), dtype=jnp.int32),
        loaded.device,
    )
    out = loaded.model(
        tokens,
        [species_ids(loaded=loaded, species=species)],
        track_indices=jax.device_put(
            jnp.asarray(track_indices, dtype=jnp.int32), loaded.device
        ),
        skip_bed_head=True,
        skip_lm_head=True,
    )
    return {key: as_numpy(value) for key, value in out.items()}


def token_gradients(
    loaded: LoadedModel,
    sequence: str,
    species: str,
    track_index: int,
) -> np.ndarray:
    """Gradient of one track's mean prediction with respect to each input token.

    The model embeds token ids with ``embed_layer``. Differentiating a zero
    shift added to that embedding gives the gradient at each position's
    embedding, ``g`` of shape (L, D). A one-hot token ``x`` embeds as
    ``x @ E``, so the gradient with respect to the one-hot is ``g @ E.T``.
    The gradient runs on a merged copy of the model: the stem norm writes its
    running statistic on every call, and a traced write must not reach the
    loaded model.

    Args:
        loaded: A loaded model bundle.
        sequence: DNA sequence, soft mask kept, length a multiple of 128.
        species: Species string to condition on (e.g. ``human``).
        track_index: Head column, a position in
            ``config.bigwigs_per_file_assembly[assembly]``.

    Returns:
        Array of shape (L, alphabet_size); column v is token id v.
    """
    tokens = jax.device_put(
        jnp.asarray(loaded.tokenizer.batch_np_tokenize([sequence]), dtype=jnp.int32),
        loaded.device,
    )
    conditions = [species_ids(loaded=loaded, species=species)]
    indices = jax.device_put(jnp.asarray([track_index], dtype=jnp.int32), loaded.device)
    graphdef, state = nnx.split(loaded.model)

    def mean_prediction(shift: jax.Array) -> jax.Array:
        model = nnx.merge(graphdef, state)
        embed = model.embed_layer
        model.embed_layer = lambda ids: embed(ids) + shift
        out = model(
            tokens,
            conditions,
            track_indices=indices,
            skip_bed_head=True,
            skip_lm_head=True,
        )
        return out["bigwig_tracks_logits"][0, :, 0].astype(jnp.float32).mean()

    shape = jax.eval_shape(loaded.model.embed_layer, tokens)
    grad = jax.grad(mean_prediction)(jnp.zeros(shape.shape, shape.dtype))
    table = loaded.model.embed_layer.embedding.value.astype(jnp.float32)
    return as_numpy(grad[0].astype(jnp.float32) @ table.T)


def acgt_embedding(loaded: LoadedModel) -> Any:
    """Not available: the public model has no ``inputs_embeds`` path."""
    raise NotImplementedError(_NO_ONE_HOT_FORWARD)


def make_one_hot_forward(loaded: LoadedModel, species: str) -> Any:
    """Not available: the public model has no ``inputs_embeds`` path."""
    raise NotImplementedError(_NO_ONE_HOT_FORWARD)


def make_selected_tracks_one_hot_forward(
    loaded: LoadedModel,
    species: str,
    track_indices: list[int],
) -> Any:
    """Not available: the public model has no ``inputs_embeds`` path."""
    raise NotImplementedError(_NO_ONE_HOT_FORWARD)


def run_one_hot_batch(
    loaded: LoadedModel,
    one_hot_acgt: Any,
    species: str,
    track_indices: list[int],
) -> dict[str, Any]:
    """Not available: the public model has no ``inputs_embeds`` path."""
    raise NotImplementedError(_NO_ONE_HOT_FORWARD)


__all__ = [
    "DEFAULT_CHECKPOINT",
    "LoadedModel",
    "acgt_embedding",
    "as_numpy",
    "bed_element_names",
    "bigwig_track_names",
    "center_fraction",
    "load_model",
    "make_one_hot_forward",
    "make_selected_tracks_one_hot_forward",
    "prediction_region",
    "resolve_device",
    "run_forward",
    "run_one_hot_batch",
    "species_ids",
    "token_gradients",
]
