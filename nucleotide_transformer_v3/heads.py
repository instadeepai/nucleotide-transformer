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

"""Prediction heads for the NTv3 U-Net/GQA architecture.

``AdaptiveLayerNorm`` and ``LMHead`` are carried here rather than imported so
that this package stands alone; both are copies, not reimplementations.
"""

from collections.abc import Sequence

import jax
import jax.numpy as jnp
from flax import nnx
from flax.typing import Dtype


class AdaptiveLayerNorm(nnx.LayerNorm):
    """Adaptive layer norm. Rescale a layer norm with one or several conditions."""

    def __init__(
        self,
        num_features: int,
        conditions_dims: list[int],
        epsilon: float = 1e-5,
        *,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """Initialize the AdaptiveLayerNorm layer.

        Args:
            num_features: input dimension.
            conditions_dims: conditions dimensions.
            epsilon: epsilon value.
            rngs: random number generator.
        """
        super().__init__(
            num_features=num_features,
            epsilon=epsilon,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )
        self.modulation_layers = [
            nnx.Linear(
                in_features=cond_dim,
                out_features=2 * num_features,
                kernel_init=nnx.initializers.constant(0.0),
                bias_init=nnx.initializers.constant(0.0),
                dtype=modulation_dtype,
                param_dtype=modulation_param_dtype,
                rngs=rngs,
            )
            for cond_dim in conditions_dims
        ]
        self._num_conditions = len(conditions_dims)
        self._dim = num_features

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = super().__call__(x=x)

        assert len(conditions) == self._num_conditions, "num conditions mismatch"

        if conditions_masks is None:
            conditions_masks = [
                jnp.ones(shape=(x.shape[0]), dtype=x.dtype) for _ in conditions
            ]

        assert len(conditions_masks) == self._num_conditions, (
            "num conditions masks mismatch"
        )

        scale = 1.0
        shift = 0.0
        for i, x_c_i in enumerate(conditions):
            tmp = self.modulation_layers[i](x_c_i)[:, None, :]
            shift_i, scale_i = jnp.split(tmp, indices_or_sections=2, axis=-1)
            condition_mask = conditions_masks[i]
            shift_i = jnp.where(condition_mask[:, None, None], shift_i, 0.0)
            scale_i = jnp.where(condition_mask[:, None, None], scale_i, 0.0)
            scale = scale * (1.0 + scale_i)
            shift = shift + shift_i
        x = x * scale + shift
        return x


class LMHead(nnx.Module):
    def __init__(
        self,
        embed_dim: int,
        num_hidden_layers: int,
        alphabet_size: int,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        self.hidden_layers = [
            nnx.Linear(
                in_features=embed_dim,
                out_features=embed_dim,
                rngs=rngs,
                dtype=dtype,
                param_dtype=param_dtype,
            )
            for _ in range(num_hidden_layers)
        ]
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=alphabet_size,
            rngs=rngs,
            dtype=dtype,
            param_dtype=param_dtype,
        )

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        x = jax.nn.gelu(x)
        for layer in self.hidden_layers:
            x = layer(x)
            x = jax.nn.gelu(x)
        x = self.head(x)
        return x

def _sliced_linear(
    linear: nnx.Linear, x: jax.Array, track_indices: jax.Array
) -> jax.Array:
    """Apply ``linear`` restricted to the ``track_indices`` output columns.

    Equivalent to ``linear(x)[..., track_indices]`` but never materialises the
    full ``(..., out_features)`` activation: the kernel and bias are sliced
    first, so the matmul is over ``len(track_indices)`` columns.

    It reuses ``linear``'s own ``promote_dtype``, ``dot_general`` and
    ``precision`` rather than re-implementing them, so it cannot silently
    desynchronise from :meth:`nnx.Linear.__call__` if either is later
    configured - a hand-rolled ``jnp.matmul`` here would match only while both
    stay at their defaults.
    """
    kernel = linear.kernel.value[:, track_indices]
    bias = None if linear.bias is None else linear.bias.value[track_indices]
    x, kernel, bias = nnx.nn.dtypes.promote_dtype((x, kernel, bias), dtype=linear.dtype)
    y = linear.dot_general(
        x,
        kernel,
        (((x.ndim - 1,), (0,)), ((), ())),
        precision=linear.precision,
    )
    if bias is not None:
        y += jnp.reshape(bias, (1,) * (y.ndim - 1) + (-1,))
    return y


class LinearHead(nnx.Module):
    """A linear head that predicts one scalar value per track."""

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        *,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        self.layer_norm = nnx.LayerNorm(embed_dim, rngs=rngs)
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels,
            rngs=rngs,
        )

    def __call__(
        self, x: jax.Array, *, track_indices: jax.Array | None = None
    ) -> jax.Array:
        x = self.layer_norm(x)
        if track_indices is not None:
            x = _sliced_linear(self.head, x, track_indices)
        else:
            x = self.head(x)
        x = jax.nn.softplus(x)
        return x


class ClassificationHead(nnx.Module):
    """A linear head that predicts num_classes values per label."""

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        num_classes: int,
        layer_norm: bool = True,
        *,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        if layer_norm:
            self.layer_norm = nnx.LayerNorm(embed_dim, rngs=rngs)
        else:
            self.layer_norm = None
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels * num_classes,
            rngs=rngs,
        )
        self.num_labels = num_labels
        self.num_classes = num_classes

    def __call__(self, x: jax.Array) -> jax.Array:
        if self.layer_norm is not None:
            x = self.layer_norm(x)
        x = self.head(x)
        x = jnp.reshape(x, (*x.shape[:-1], self.num_labels, self.num_classes))
        return x


def _make_output_bias_layers(
    conditions_dims: Sequence[int],
    out_features: int,
    rngs: nnx.Rngs,
) -> list[nnx.Linear]:
    """Zero-init Linear layers projecting each condition to ``out_features``.

    Used for the output-side (post-linear, pre-activation) additive bias of the
    adaptive heads. Zero init keeps the bias at 0 at step 0, so the head is
    identity to its unconditioned counterpart at init.
    """
    return [
        nnx.Linear(
            in_features=cond_dim,
            out_features=out_features,
            kernel_init=nnx.initializers.constant(0.0),
            bias_init=nnx.initializers.constant(0.0),
            rngs=rngs,
        )
        for cond_dim in conditions_dims
    ]


def _apply_output_bias(
    logits: jax.Array,
    bias_layers: list[nnx.Linear],
    conditions: list[jax.Array],
    conditions_masks: list[jax.Array] | None,
    *,
    track_indices: jax.Array | None = None,
) -> jax.Array:
    """Add a masked, zero-init per-condition additive bias to ``logits``.

    ``logits`` has shape ``(B, L, out_features)``; each bias is ``(B,
    out_features)`` broadcast over ``L``. When ``track_indices`` is given,
    ``logits`` carries only those labels, so each bias is restricted to the same
    columns before being added.
    """
    if conditions_masks is None:
        conditions_masks = [
            jnp.ones((logits.shape[0],), dtype=logits.dtype) for _ in conditions
        ]
    for cond, mask, bias_layer in zip(
        conditions, conditions_masks, bias_layers, strict=True
    ):
        if track_indices is None:
            bias = bias_layer(cond).astype(logits.dtype)  # (B, out_features)
        else:
            # Slice the projection rather than computing all out_features and
            # discarding most of them, so both per-track paths in this file
            # read the same way.
            bias = _sliced_linear(bias_layer, cond, track_indices).astype(logits.dtype)
        bias = bias * mask.reshape(-1, 1).astype(logits.dtype)
        logits = logits + bias[:, None, :]
    return logits


class AdaptiveLinearHead(nnx.Module):
    """FiLM-conditioned :class:`LinearHead`.

    Modulates the input embedding with an :class:`AdaptiveLayerNorm` (input-side
    FiLM) and, optionally, adds a zero-init per-condition additive bias on the
    raw logits (output-side FiLM) before the ``softplus``. Both FiLM paths are
    zero-init, so the head is identity to a plain :class:`LinearHead` at init.
    """

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        conditions_dims: Sequence[int],
        *,
        use_output_bias: bool = True,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        self._num_labels = num_labels
        self._use_output_bias = use_output_bias
        # `dtype` applies to the projection only; norm and output-bias stay fp32.
        self.layer_norm = AdaptiveLayerNorm(
            num_features=embed_dim,
            conditions_dims=list(conditions_dims),
            rngs=rngs,
        )
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        if use_output_bias:
            self.track_bias_layers = _make_output_bias_layers(
                conditions_dims, num_labels, rngs
            )

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        *,
        track_indices: jax.Array | None = None,
    ) -> jax.Array:
        """Predict per-position track values.

        Args:
            x: Input embedding of shape ``(B, L, embed_dim)``.
            conditions: Per-condition embeddings for the FiLM paths.
            conditions_masks: Optional per-example condition masks.
            track_indices: When given, compute only these tracks by slicing the
                projection kernel, as :class:`LinearHead` does. The full head
                materialises a full ``(B, L, num_labels)`` activation, which at 1 bp
                resolution over a long window is tens of GB for a wide head
                (projection, output bias, softplus), which at 1 bp resolution
                over a long window is tens of GB for a wide head; slicing keeps
                that proportional to the tracks actually wanted.

        Returns:
            Track predictions, ``(B, L, num_labels)`` or ``(B, L,
            len(track_indices))`` when ``track_indices`` is given.
        """
        x = self.layer_norm(x, conditions, conditions_masks)
        # `dtype` lowers the projection matmul only; logits go back to fp32.
        if track_indices is None:
            logits = self.head(x).astype(jnp.float32)
        else:
            logits = _sliced_linear(self.head, x, track_indices).astype(jnp.float32)
        if self._use_output_bias:
            logits = _apply_output_bias(
                logits,
                self.track_bias_layers,
                conditions,
                conditions_masks,
                track_indices=track_indices,
            )
        return jax.nn.softplus(logits)


class AdaptiveClassificationHead(nnx.Module):
    """FiLM-conditioned :class:`ClassificationHead`.

    Same FiLM mechanism as :class:`AdaptiveLinearHead` (input-side
    ``AdaptiveLayerNorm`` plus optional zero-init output bias), predicting
    ``num_classes`` logits per label. Identity to :class:`ClassificationHead`
    at init.
    """

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        num_classes: int,
        conditions_dims: Sequence[int],
        layer_norm: bool = True,
        *,
        use_output_bias: bool = True,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        self.num_labels = num_labels
        self.num_classes = num_classes
        self._use_output_bias = use_output_bias
        # As in AdaptiveLinearHead: `dtype` applies to the projection only.
        if layer_norm:
            self.layer_norm = AdaptiveLayerNorm(
                num_features=embed_dim,
                conditions_dims=list(conditions_dims),
                rngs=rngs,
            )
        else:
            self.layer_norm = None
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels * num_classes,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        if use_output_bias:
            self.track_bias_layers = _make_output_bias_layers(
                conditions_dims, num_labels * num_classes, rngs
            )

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        if self.layer_norm is not None:
            x = self.layer_norm(x, conditions, conditions_masks)
        logits = self.head(x).astype(jnp.float32)
        if self._use_output_bias:
            logits = _apply_output_bias(
                logits, self.track_bias_layers, conditions, conditions_masks
            )
        return jnp.reshape(
            logits, (*logits.shape[:-1], self.num_labels, self.num_classes)
        )
