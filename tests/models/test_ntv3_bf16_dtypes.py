# This source code is licensed under CC BY-NC-SA 4.0. See the LICENSE file.
"""The bf16 compute path, which is what ``ntv3-annotate`` runs by default.

``ProductionAnnotationConfig.use_bf16_compute`` defaults to ``True``, so this
is not an optional mode -- it is the one every user gets, and the acceptance
run on the GPU deliberately pinned fp32 because bf16 is not bit-reproducible
and would have measured the dtype rather than the port. That leaves the
default unguarded unless it is pinned here.
"""

import dataclasses

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx

from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.dtypes import (
    COMPUTE_DTYPE_FIELDS,
    HEAD_COMPUTE_DTYPE_FIELDS,
    build_bf16_compute_model,
)
from nucleotide_transformer_v3.model import NTv3

# Tiny but structurally complete: stem, 2 encoder/decoder stages, 2 trunk layers.
SMALL_CONFIG_KWARGS = {
    "alphabet_size": 11,
    "token_embed_dim": 8,
    "conv_init_embed_dim": 16,
    "embed_dim": 32,
    "num_downsamples": 2,
    "num_layers": 2,
    "attention_heads": 4,
    "num_kv_heads": 2,
    "key_size": 8,
    "value_size": 8,
}

CONFIG_COMPUTE_DTYPE_FIELDS = tuple(
    field.name
    for field in dataclasses.fields(NTv3Config)
    if field.name.endswith("_compute_dtype")
)
PARAM_DTYPE_FIELDS = tuple(
    field.name
    for field in dataclasses.fields(NTv3Config)
    if field.name.endswith("_param_dtype")
)


def _tokens() -> jax.Array:
    return jnp.asarray(
        np.random.default_rng(0).integers(3, 11, size=(2, 64)), dtype=jnp.int32
    )


def _param_dtypes(model: nnx.Module) -> dict[str, str]:
    """Flattened param path -> dtype name, for structural comparisons."""
    flat = jax.tree_util.tree_flatten_with_path(nnx.state(model, nnx.Param))[0]
    return {"/".join(str(key) for key in path): str(leaf.dtype) for path, leaf in flat}


def _fp32_model() -> NTv3:
    return NTv3(NTv3Config(**SMALL_CONFIG_KWARGS), rngs=nnx.Rngs(params=0))


def test_compute_dtype_fields_membership_is_pinned() -> None:
    """``COMPUTE_DTYPE_FIELDS`` is hand-written, so pin it against the config.

    ``dtypes.py`` names this test in a comment as the thing that keeps the two
    in step, and it was not ported with the module. A sixth
    ``*_compute_dtype`` knob would otherwise be missed silently, leaving that
    submodule in fp32 while the rest of the trunk ran bf16.
    """
    assert set(COMPUTE_DTYPE_FIELDS) == {
        "stem_compute_dtype",
        "down_convolution_compute_dtype",
        "up_convolution_compute_dtype",
        "transformer_qkvo_compute_dtype",
        "transformer_ffn_compute_dtype",
    }
    assert set(COMPUTE_DTYPE_FIELDS) == set(CONFIG_COMPUTE_DTYPE_FIELDS)


def test_bf16_compute_leaves_param_dtypes_untouched() -> None:
    """The guarantee that lets an fp32 export load onto a bf16 graph."""
    fp32 = _fp32_model()
    bf16 = build_bf16_compute_model(fp32)

    assert _param_dtypes(bf16) == _param_dtypes(fp32)
    assert set(_param_dtypes(bf16).values()) == {"float32"}
    for name in PARAM_DTYPE_FIELDS:
        assert getattr(bf16.config, name) == "float32", name


def test_build_bf16_compute_model_preserves_batch_stats_and_eval_mode() -> None:
    """The 34 ``var_ema`` buffers have to ride along or eval mis-scales.

    They are ``BatchStat``, not ``Param``, and the rebuild splits and restores
    the state; dropping them leaves the norm dividing by its ones-init, which
    produces garbage with no error.
    """
    model = _fp32_model()
    model.train()
    model(_tokens())  # move var_ema off its jnp.ones init
    model.eval()
    var_ema = model.stem.res_norm.var_ema.value
    assert not jnp.array_equal(var_ema, jnp.ones_like(var_ema))

    rebuilt = build_bf16_compute_model(model)

    assert bool(rebuilt.stem.res_norm.use_running_average)
    assert jnp.array_equal(rebuilt.stem.res_norm.var_ema.value, var_ema)


def test_bf16_forward_tracks_fp32_and_leaves_the_source_usable() -> None:
    """bf16 rounding only, and the rebuild does not consume the fp32 model."""
    tokens = _tokens()
    fp32 = _fp32_model()
    fp32.eval()
    expected = fp32(tokens)

    bf16 = build_bf16_compute_model(fp32)
    actual = bf16(tokens)

    assert actual.keys() == expected.keys()
    for key, reference in expected.items():
        scale = float(jnp.max(jnp.abs(reference))) + 1e-9
        error = float(jnp.max(jnp.abs(actual[key] - reference))) / scale
        assert error < 0.05, f"{key}: relative error {error}"
    # The source model still runs, and still in fp32.
    for name in COMPUTE_DTYPE_FIELDS:
        assert getattr(fp32.config, name) == "float32", name
    assert fp32(tokens).keys() == expected.keys()


def test_the_annotation_head_is_never_lowered_to_bf16() -> None:
    """argmax and Viterbi read the annotation logits, so the head stays fp32.

    It is a rounding error away from a different decoded class, and it is a
    fraction of a percent of trunk FLOPs, so there is nothing to win.
    """
    assert "annotation_head_compute_dtype" not in HEAD_COMPUTE_DTYPE_FIELDS
