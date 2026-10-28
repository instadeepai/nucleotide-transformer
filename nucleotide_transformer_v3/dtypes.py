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

"""Switch a loaded NTv3 model's compute dtypes (fp32 training -> bf16 inference).

Submodules read ``dtype`` once in their own ``__init__``, so mutating
``model.config`` after a load is a no-op; the switch has to rebuild the graph and
transplant the loaded weights. Only the compute dtypes move, so the parameter
tree is unchanged and an fp32 checkpoint fits as-is.

Typical inference use, right after the checkpoint load::

    model = load_model_params(...)  # or load_ema_params(...)
    if use_bf16_compute:
        model = build_bf16_compute_model(model)
    model.eval()
"""

import dataclasses
from collections.abc import Callable
from typing import TypeVar, cast

import jax
import jax.numpy as jnp
from flax import nnx
from flax.typing import Dtype

from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.logging import get_logger

logger = get_logger()

# The trunk knobs build_bf16_compute_model flips: every *_compute_dtype
# NTv3Config exposes. Listed rather than scanned off the config so that a future
# knob has to be opted in deliberately; pinned against the config by
# test_compute_dtype_fields_membership_is_pinned.
COMPUTE_DTYPE_FIELDS: tuple[str, ...] = (
    "stem_compute_dtype",
    "down_convolution_compute_dtype",
    "up_convolution_compute_dtype",
    "transformer_qkvo_compute_dtype",
    "transformer_ffn_compute_dtype",
)

# Flipped by cast_head=True: the bigwig head only. Its projection is
# number_bigwig_tracks wide (~0.93x the whole trunk's FLOPs at 7362 tracks).
# The annotation head is ~0.19% of trunk FLOPs -- too small a share to be worth
# lowering, while argmax reads its logits directly, so bf16 there only shifts
# almost-tied positions.
HEAD_COMPUTE_DTYPE_FIELDS: tuple[str, ...] = ("bigwig_head_compute_dtype",)

ConfigT = TypeVar("ConfigT", bound=NTv3Config)
ModelT = TypeVar("ModelT", bound=nnx.Module)


def _canonical_dtype_name(dtype: Dtype) -> str:
    """Normalize a dtype to its canonical name, e.g. ``jnp.bfloat16`` -> "bfloat16".

    Names serialize into a saved config and are what ``_convert_config_bf16_to_f32``
    string-matches on; a type object would work for Flax but slip past it.
    """
    return str(jnp.dtype(dtype).name)


def _fields_to_flip(config: NTv3Config, *, cast_head: bool) -> list[str]:
    """Names of the compute-dtype fields to move on ``config``."""
    names = list(COMPUTE_DTYPE_FIELDS)
    if cast_head:
        names += [name for name in HEAD_COMPUTE_DTYPE_FIELDS if hasattr(config, name)]
    return names


def with_compute_dtype(
    config: ConfigT,
    compute_dtype: Dtype,
    *,
    cast_head: bool = False,
) -> ConfigT:
    """Return a copy of ``config`` with the compute dtypes set to ``compute_dtype``.

    Every ``*_param_dtype`` is left alone, so parameter storage is unaffected.

    Args:
        config: Source config (any :class:`NTv3Config` subclass).
        compute_dtype: Target compute dtype, e.g. ``jnp.bfloat16``.
        cast_head: Also set ``HEAD_COMPUTE_DTYPE_FIELDS``, where the config
            exposes them.

    Returns:
        A new config of the same type.
    """
    dtype_name = _canonical_dtype_name(compute_dtype)
    names = _fields_to_flip(config, cast_head=cast_head)
    return dataclasses.replace(config, **dict.fromkeys(names, dtype_name))


def _leaf_paths(state: nnx.State) -> set[str]:
    """Flattened variable paths of an nnx state, as ``"a/b/c"`` strings."""
    return {
        "/".join(str(key) for key in path)
        for path, _ in jax.tree_util.tree_flatten_with_path(state)[0]
    }


def build_bf16_compute_model(
    model: ModelT,
    *,
    compute_dtype: Dtype = jnp.bfloat16,
    cast_head: bool = False,
) -> ModelT:
    """Rebuild ``model`` with bfloat16 compute dtypes, keeping its loaded weights.

    Rebuilds the graph from a patched ``model.config`` with :func:`nnx.eval_shape`
    (so no second set of random parameters is materialized), then transplants the
    full variable state: ``nnx.Param``, ``nnx.BatchStat``
    (``RMSBatchNorm.var_ema``) and the train/eval mode flag alike. Dropping the
    BatchStats would leave ``var_ema`` at its ones-init and mis-scale every
    stem/encoder/decoder activation in eval mode.

    Parameter storage is untouched, which is what lets an fp32 checkpoint's
    weights land on the rebuilt graph unchanged. Norms, embeddings and the LM head
    stay float32.

    Args:
        model: A loaded model whose ``config`` is an :class:`NTv3Config` (or a
            subclass, e.g. ``DiscreteConditionedNTv3WithTracksHeadsConfig``).
        compute_dtype: Dtype for the matmul/conv paths. Defaults to
            ``jnp.bfloat16``; pass ``jnp.float32`` to rebuild an fp32 graph.
        cast_head: Also lower the bigwig head's projection. Defaults to False; it
            only pays on a wide head. A no-op on configs without the field.

    Returns:
        A new model of the same type carrying ``model``'s weights.

    Raises:
        TypeError: If ``model.config`` is not a dataclass exposing the
            compute-dtype fields (a model family without these knobs).
        ValueError: If the rebuilt graph's variable structure does not match the
            source model's, which would leave part of it uninitialized.
    """
    config = getattr(model, "config", None)
    if config is None or not dataclasses.is_dataclass(config):
        raise TypeError(
            f"{type(model).__name__} has no dataclass `config`; cannot patch "
            "compute dtypes."
        )
    # Only the trunk fields are required: the head fields live on the
    # tracks-heads config and are legitimately absent from a plain NTv3.
    missing = [name for name in COMPUTE_DTYPE_FIELDS if not hasattr(config, name)]
    if missing:
        raise TypeError(
            f"{type(config).__name__} does not expose compute-dtype fields "
            f"{missing}; this model family has no bf16 compute knobs."
        )
    patched = with_compute_dtype(config, compute_dtype, cast_head=cast_head)
    state = nnx.state(model)
    model_class = type(model)
    # Constructed dynamically from the runtime class; every NTv3 variant takes
    # ``(config, rngs=...)`` but that signature is not on the ``nnx.Module`` bound.
    build = cast(Callable[..., ModelT], model_class)
    rebuilt = nnx.eval_shape(lambda: build(patched, rngs=nnx.Rngs(0)))

    source_paths, rebuilt_paths = _leaf_paths(state), _leaf_paths(nnx.state(rebuilt))
    if source_paths != rebuilt_paths:
        raise ValueError(
            f"Rebuilding {model_class.__name__} with compute_dtype="
            f"{_canonical_dtype_name(compute_dtype)} changed the variable "
            f"structure: {len(rebuilt_paths - source_paths)} new and "
            f"{len(source_paths - rebuilt_paths)} dropped leaves. Refusing to "
            "return a partially-initialised model."
        )

    nnx.update(rebuilt, state)
    flipped = ", ".join(_fields_to_flip(config, cast_head=cast_head))
    logger.info(
        f"Rebuilt {model_class.__name__} with "
        f"{_canonical_dtype_name(compute_dtype)} compute dtypes "
        f"({flipped}); parameter dtypes unchanged."
    )
    return rebuilt
