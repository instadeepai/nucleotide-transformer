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

"""NTv3 building blocks: primitives, conv blocks, transformer blocks."""

from collections.abc import Sequence
from typing import Literal

import jax
import jax.numpy as jnp
from flax import nnx
from flax.typing import Dtype

_MAX_POSITION = 8192  # 1 Mb / 128 bp resolution at trunk


# ---------------------------------------------------------------------------
# Low-level primitives
# ---------------------------------------------------------------------------


class RMSBatchNorm(nnx.Module):
    r"""Batch normalization without mean centering (RMS-only).

    Normalizes ``x`` by the running RMS (root-mean-square) of its values along
    the batch and spatial axes, then applies a learnable per-feature scale and
    offset:

    .. math:: y = \text{scale} \cdot \frac{x}{\sqrt{\overline{x^2} + \epsilon}} + \text{offset}

    The ``\overline{x^2}`` quantity is the running EMA in eval mode and the
    current batch's mean of squares in train mode (which also updates the EMA).
    Mode is toggled by ``model.train()`` / ``model.eval()``.

    The ``use_running_average`` flag is stored as an :class:`nnx.Variable` of a
    JAX bool, not as a plain Python attribute, so toggling it (e.g. via
    ``model.eval()`` -> ``set_attributes(use_running_average=True, ...)``)
    becomes a *value* update rather than a *GraphDef* mutation. That makes it
    safe to flip inside an ``nnx.jit``-compiled function (the standard
    train/validation step pattern), which would otherwise hit
    ``ValueError: The graph structure of a node added to cached_partial was
    mutated inside the transformation``.

    Args:
        num_features: size of the last (feature) axis.
        ndim: total expected rank of the input (used to shape the per-feature
            params as ``(1, ..., 1, num_features)``).
        use_running_average: initial mode. ``False`` (default) puts the layer
            in train mode (use batch RMS, update EMA). ``True`` puts it in
            eval mode (use the EMA, skip the EMA update). Live toggling via
            ``model.train()`` / ``model.eval()`` works at any time, including
            inside ``nnx.jit``.
        ema_momentum: decay factor for ``var_ema`` updates. A value close to 1
            (default 0.999) tracks a slow EMA appropriate for long sequences
            and small batch sizes.
        epsilon: numerical stability constant added under the rsqrt.
    """

    def __init__(
        self,
        num_features: int,
        *,
        ndim: int,
        use_running_average: bool = False,
        ema_momentum: float = 0.999,
        epsilon: float = 1e-5,
    ) -> None:
        self.epsilon = epsilon
        self.ema_momentum = ema_momentum
        # Stored as a Variable so that flipping it via set_attributes inside
        # an nnx.jit'd function is a *value* update, not a GraphDef change.
        # The public ``use_running_average`` property below reads/writes it.
        self._use_running_average = nnx.Variable(jnp.asarray(use_running_average))

        broadcast_shape = (1,) * (ndim - 1) + (num_features,)
        self.scale = nnx.Param(jnp.ones(broadcast_shape))
        self.offset = nnx.Param(jnp.zeros(broadcast_shape))
        # ``BatchStat`` separates running stats from trainable params in the
        # state dict; this is what ``model.eval()`` toggles read.
        self.var_ema = nnx.BatchStat(jnp.ones(broadcast_shape))

    @property
    def use_running_average(self) -> jax.Array:
        """Current mode as a JAX scalar bool (True = eval, False = train)."""
        return self._use_running_average.value

    @use_running_average.setter
    def use_running_average(self, value: bool | jax.Array) -> None:
        # set_attributes(use_running_average=True/False) lands here. Update
        # the Variable's value in place; the Variable object identity (and
        # therefore the GraphDef) is preserved.
        self._use_running_average.value = jnp.asarray(value)

    def __call__(self, x: jax.Array) -> jax.Array:
        # Explicit fp32 floor: the casts below follow `x.dtype` instead of
        # promoting `x` up to the params' dtype, so a bf16 input would otherwise
        # carry bf16 through the norm.
        x = jnp.asarray(x, jnp.promote_types(x.dtype, jnp.float32))
        dtype = x.dtype
        reduce_axes = tuple(range(x.ndim - 1))

        # E[x^2] over batch + spatial axes. This is the "RMS" quantity (no
        # centering, so it is the second raw moment, not the variance).
        # Computed unconditionally (cheap); both modes need a value of the
        # right shape/dtype to feed into the jnp.where calls below.
        batch_second_moment = jnp.mean(
            jnp.square(x.astype(jnp.float32)), axis=reduce_axes, keepdims=True
        )

        # jnp.where, not zero-weighting: 0.0 * inf is NaN, so one overflowing batch
        # would poison var_ema forever, and only eval reads it so training looks fine.
        use_running = self.use_running_average  # JAX scalar bool
        ema_update = self.ema_momentum * self.var_ema.value + (
            1.0 - self.ema_momentum
        ) * batch_second_moment.squeeze(reduce_axes)
        next_var_ema = jnp.where(use_running, self.var_ema.value, ema_update)
        self.var_ema.value = jnp.where(
            jnp.isfinite(next_var_ema), next_var_ema, self.var_ema.value
        )

        # In eval mode normalize by the EMA; in train mode normalize by the
        # current batch's E[x^2].
        variance = jnp.where(
            use_running,
            jnp.asarray(self.var_ema.value, dtype=dtype),
            batch_second_moment.astype(dtype),
        )

        scale = jnp.asarray(self.scale.value, dtype=dtype)
        offset = jnp.asarray(self.offset.value, dtype=dtype)
        inv = scale * jax.lax.rsqrt(variance + self.epsilon).astype(dtype)
        return x * inv + offset


# ---------------------------------------------------------------------------
# FiLM (feature-wise linear modulation) conditioning
# ---------------------------------------------------------------------------


def _make_film_modulation_layers(
    conditions_dims: Sequence[int],
    num_features: int,
    rngs: nnx.Rngs,
) -> list[nnx.Linear]:
    """Zero-init Linear layers projecting each condition to ``2 * num_features``.

    One layer per condition. The output is split into ``(shift_i, scale_i)``.
    Zero init (kernel and bias) makes every FiLM norm an identity at step 0:
    ``scale_i = 0`` and ``shift_i = 0`` so the modulation reduces to the plain
    base norm output. Training moves the modulation away from identity.
    """
    return [
        nnx.Linear(
            in_features=cond_dim,
            out_features=2 * num_features,
            kernel_init=nnx.initializers.constant(0.0),
            bias_init=nnx.initializers.constant(0.0),
            rngs=rngs,
        )
        for cond_dim in conditions_dims
    ]


def _apply_film(
    x: jax.Array,
    num_features: int,
    modulation_layers: list[nnx.Linear],
    conditions: list[jax.Array] | None,
    conditions_masks: list[jax.Array] | None,
) -> jax.Array:
    """Apply FiLM modulation to an already-normalized tensor.

    Computes ``x * (1 + sum_i scale_i(cond_i)) + sum_i shift_i(cond_i)`` where
    each ``(shift_i, scale_i)`` is produced by a per-condition Linear and
    broadcast over the leading (spatial) axes. A per-condition mask (``True`` =
    apply) zeroes out both terms for masked examples, enabling condition
    dropout / unconditional inference without changing the graph.

    Args:
        x: normalized input of shape ``(B, ..., num_features)``.
        num_features: size of the last (feature) axis.
        modulation_layers: zero-init Linear layers, one per condition.
        conditions: per-example condition embeddings ``[(B, cond_dim_i)]`` or
            ``None`` (returns ``x`` unchanged).
        conditions_masks: per-example masks ``[(B,)]`` or ``None`` (all ones).

    Returns:
        Modulated tensor with the same shape as ``x``.
    """
    if conditions is None or len(modulation_layers) == 0:
        return x

    if len(conditions) != len(modulation_layers):
        raise ValueError(
            f"num conditions ({len(conditions)}) must match num modulation "
            f"layers ({len(modulation_layers)})"
        )

    if conditions_masks is None:
        conditions_masks = [jnp.ones((x.shape[0],), dtype=x.dtype) for _ in conditions]
    if len(conditions_masks) != len(conditions):
        raise ValueError("num conditions masks must match num conditions")

    # Broadcast (B, num_features) -> (B, 1, ..., 1, num_features).
    feature_shape = (x.shape[0],) + (1,) * (x.ndim - 2) + (num_features,)
    mask_shape = (x.shape[0],) + (1,) * (x.ndim - 1)

    scale = jnp.zeros((), dtype=x.dtype)
    shift = jnp.zeros((), dtype=x.dtype)
    for cond, mask, mod in zip(
        conditions, conditions_masks, modulation_layers, strict=True
    ):
        params = mod(cond).astype(x.dtype)  # (B, 2 * num_features)
        shift_i, scale_i = jnp.split(params, 2, axis=-1)
        shift_i = shift_i.reshape(feature_shape)
        scale_i = scale_i.reshape(feature_shape)
        mask_b = mask.reshape(mask_shape).astype(x.dtype)
        shift_i = shift_i * mask_b
        scale_i = scale_i * mask_b
        # Additive composition across conditions, matching the docstring:
        # x * (1 + sum_i scale_i) + sum_i shift_i. Zero-init => identity.
        scale = scale + scale_i
        shift = shift + shift_i
    return x * (1.0 + scale) + shift


class AdaptiveRMSBatchNorm(RMSBatchNorm):
    r""":class:`RMSBatchNorm` with FiLM conditioning (for conv blocks).

    Applies the plain RMS batch norm, then a zero-init FiLM modulation:

    .. math:: y = \text{RMSBatchNorm}(x) \cdot (1 + \sum_i s_i(c_i))
        + \sum_i t_i(c_i)

    At init the modulation layers are zero, so ``y == RMSBatchNorm(x)`` (the
    layer is an identity extension of its parent). Passing ``conditions=None``
    also reduces to the plain norm, so the same module can run unconditioned.
    """

    def __init__(
        self,
        num_features: int,
        conditions_dims: Sequence[int],
        *,
        ndim: int,
        use_running_average: bool = False,
        ema_momentum: float = 0.999,
        epsilon: float = 1e-5,
        rngs: nnx.Rngs,
    ) -> None:
        super().__init__(
            num_features,
            ndim=ndim,
            use_running_average=use_running_average,
            ema_momentum=ema_momentum,
            epsilon=epsilon,
        )
        self._num_features = num_features
        self.modulation_layers = _make_film_modulation_layers(
            conditions_dims, num_features, rngs
        )

    def __call__(  # type: ignore[override]
        self,
        x: jax.Array,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = super().__call__(x)
        return _apply_film(
            x, self._num_features, self.modulation_layers, conditions, conditions_masks
        )


class AdaptiveRMSNorm(nnx.RMSNorm):
    """:class:`nnx.RMSNorm` with FiLM conditioning (for transformer pre-norms).

    Applies the plain per-token RMS norm, then a zero-init FiLM modulation
    (identical formula to :class:`AdaptiveRMSBatchNorm`). Identity at init.
    """

    def __init__(
        self,
        num_features: int,
        conditions_dims: Sequence[int],
        *,
        epsilon: float = 1e-5,
        rngs: nnx.Rngs,
    ) -> None:
        super().__init__(num_features=num_features, epsilon=epsilon, rngs=rngs)
        self._num_features = num_features
        self.modulation_layers = _make_film_modulation_layers(
            conditions_dims, num_features, rngs
        )

    def __call__(  # type: ignore[override]
        self,
        x: jax.Array,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = super().__call__(x)
        return _apply_film(
            x, self._num_features, self.modulation_layers, conditions, conditions_masks
        )


def _make_conv_norm(
    num_features: int,
    conditions_dims: Sequence[int],
    *,
    epsilon: float = 1e-5,
    rngs: nnx.Rngs,
) -> RMSBatchNorm:
    """Build a conv-block norm: FiLM ``AdaptiveRMSBatchNorm`` if conditioned.

    When ``conditions_dims`` is empty, returns exactly the plain
    :class:`RMSBatchNorm` the unconditioned model builds (same args, same
    param names, no extra rng draw), preserving byte-identical parameters.
    """
    if conditions_dims:
        return AdaptiveRMSBatchNorm(
            num_features,
            conditions_dims,
            ndim=3,
            epsilon=epsilon,
            rngs=rngs,
        )
    return RMSBatchNorm(num_features=num_features, ndim=3, epsilon=epsilon)


def _make_transformer_norm(
    num_features: int,
    conditions_dims: Sequence[int],
    *,
    epsilon: float = 1e-5,
    rngs: nnx.Rngs,
) -> nnx.RMSNorm:
    """Build a transformer pre-norm: FiLM ``AdaptiveRMSNorm`` if conditioned.

    When ``conditions_dims`` is empty, returns exactly the plain
    :class:`nnx.RMSNorm` the unconditioned model builds (same args, same rng
    draw), preserving byte-identical parameters.
    """
    if conditions_dims:
        return AdaptiveRMSNorm(
            num_features,
            conditions_dims,
            epsilon=epsilon,
            rngs=rngs,
        )
    return nnx.RMSNorm(num_features=num_features, epsilon=epsilon, rngs=rngs)


def _apply_norm(
    norm: nnx.Module,
    x: jax.Array,
    conditions: list[jax.Array] | None,
    conditions_masks: list[jax.Array] | None,
) -> jax.Array:
    """Call ``norm`` with conditions iff it is an adaptive (FiLM) norm."""
    if isinstance(norm, AdaptiveRMSBatchNorm | AdaptiveRMSNorm):
        return norm(x, conditions, conditions_masks)
    return norm(x)


def approximate_sigmoid_gelu(x: jax.Array) -> jax.Array:
    """Sigmoid-approximated GELU: ``sigmoid(1.702 * x) * x``.

    Used by the conv encoder/decoder where the original architecture used this
    cheaper approximation. The transformer FFN uses the exact ``jax.nn.gelu``.
    """
    coef = jax.lax.convert_element_type(1.702, x.dtype)
    return jax.nn.sigmoid(coef * x) * x


class StandardizedConv1D(nnx.Module):
    """1D convolution with weight standardization.

    Each forward pass standardizes the kernel weights (zero mean, unit
    variance over the spatial x input-channel axes) before the convolution,
    then multiplies by a learnable per-output-channel scale. This keeps
    activations well-conditioned without depending on a separate norm layer.

    Args:
        num_channels: number of output channels.
        width: spatial kernel width.
        in_channels: number of input channels.
        dtype: compute dtype for the convolution. ``None`` (default) follows
            the input's dtype; set explicitly (e.g. ``jnp.bfloat16``) to force
            it regardless of the input, mirroring ``nnx.Conv``'s ``dtype``.
        param_dtype: storage dtype for ``w``/``scale``/``bias``.
    """

    def __init__(
        self,
        num_channels: int,
        width: int,
        *,
        in_channels: int,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
    ) -> None:
        self.num_channels = num_channels
        self.width = width
        self.in_channels = in_channels
        self.dtype = dtype
        self.param_dtype = param_dtype

        # Zero-init kernel: every call site adds this conv inside a residual
        # branch (TokenDnaEmbedder, ConvBlock in DownResBlock/UpResBlock), so
        # identity-at-init keeps the trunk untouched at step 0.
        self.w = nnx.Param(
            jnp.zeros((width, in_channels, num_channels), dtype=param_dtype)
        )
        self.scale = nnx.Param(jnp.ones((1, 1, num_channels), dtype=param_dtype))
        self.bias = nnx.Param(jnp.zeros((num_channels,), dtype=param_dtype))

    def __call__(self, x: jax.Array) -> jax.Array:
        dtype = self.dtype if self.dtype is not None else x.dtype

        # Standardize at >= fp32 whatever the conv computes at: `scale` multiplies
        # a whole output channel, so a bf16 mean/var/rsqrt puts a systematic
        # per-channel gain error on it, and the reduction is only kernel-sized.
        stat_dtype = jnp.result_type(self.param_dtype, jnp.float32)
        w = jnp.asarray(self.w.value, dtype=stat_dtype)
        scale_param = jnp.asarray(self.scale.value, dtype=stat_dtype)

        fan_in = self.width * x.shape[-1]

        w_centered = w - jnp.mean(w, axis=(0, 1), keepdims=True)
        # jnp.var already centers internally, so compute it on the raw w to
        # avoid the redundant subtraction.
        var_w = jnp.var(w, axis=(0, 1), keepdims=True)
        # 1e-4 floor on fan_in * var_w is a divide-by-zero guard for the
        # zero-init kernel (var_w = 0 at step 0); for any non-trivial kernel
        # fan_in * var_w >> 1e-4 and the floor is inert.
        scale = scale_param * jax.lax.rsqrt(jnp.maximum(fan_in * var_w, 1e-4))
        w_standardized = jnp.asarray(w_centered * scale, dtype=dtype)

        # `x` too, not just the kernel: conv_general_dilated requires lhs/rhs to
        # share a dtype, and with an explicit `self.dtype` the input may arrive
        # in a different one.
        x = jnp.asarray(x, dtype=dtype)
        bias = jnp.asarray(self.bias.value, dtype=dtype)

        out = jax.lax.conv_general_dilated(
            lhs=x,
            rhs=w_standardized,
            window_strides=(1,),
            padding="SAME",
            dimension_numbers=jax.lax.ConvDimensionNumbers(
                lhs_spec=(0, 2, 1),  # (B, S, C) → (B, C, S)
                rhs_spec=(2, 1, 0),  # (W, I, O) → (O, I, W)
                out_spec=(0, 2, 1),  # back to (B, S, C)
            ),
        )
        return out + jnp.broadcast_to(bias, out.shape)


def apply_rope(
    x: jax.Array,  # (B, S, H, C)
    positions: jax.Array | None,  # (B, S)
    max_position: int,
) -> jax.Array:
    """Rotary position embeddings with geometrically-spaced frequencies.

    Each pair of even/odd channels of ``x`` is rotated by an angle that is a
    function of the token's position. Frequencies are spaced geometrically up
    to ``max_position``. If ``positions`` is None, the standard arange of the
    sequence axis is used.

    Note: this is the genomic RoPE recipe and differs from the LLM-style
    ``base ** (-2 i / d)`` spacing used in LLaMA / Hugging Face RoPE
    implementations. Do not diff-replace against another RoPE module
    expecting the two to be identical.

    Args:
        x: input tensor of shape ``(B, S, H, C)`` with ``C`` even.
        positions: optional explicit positions of shape ``(B, S)``. If None,
            the integer positions ``[0, 1, ..., S-1]`` are used per example.
        max_position: largest position for which a unique frequency is
            allocated. For genomic data we use 8192 (1 Mb input at 128 bp
            trunk resolution).

    Returns:
        Tensor of the same shape as ``x``.
    """
    if positions is None:
        positions = jnp.arange(x.shape[1]).astype(x.dtype).reshape(1, x.shape[1])
    num_freq = x.shape[-1] // 2
    # Genomic RoPE recipe: ``inv_freq = 1 / (arange(num_freq) +
    # geomspace(1, max_position - num_freq + 1, num_freq))``. The geomspace
    # term shifts the low end of the frequency spectrum from 1 (high
    # frequency, ~ standard RoPE) up to ``max_position`` (lowest frequency,
    # period covers the full input). Adding ``arange`` puts the high-frequency
    # bands on a near-integer spacing while leaving the low-frequency tail
    # dominated by ``geomspace``, giving a smooth interpolation between
    # local-token periodicity and long-range positional structure.
    geom_offset = jnp.geomspace(1, max_position - num_freq + 1, num_freq)
    inv_freq = (1.0 / (jnp.arange(num_freq) + geom_offset)).astype(x.dtype)
    theta = jnp.einsum("bs,f->bsf", positions, inv_freq)
    theta = jnp.repeat(theta, 2, axis=-1)[..., None, :]  # (B, S, 1, C)
    x_rotated = jnp.stack([-x[..., 1::2], x[..., ::2]], axis=-1).reshape(x.shape)
    return x * jnp.cos(theta) + x_rotated * jnp.sin(theta)


def pool(x: jax.Array, window: int = 2, reduce: str = "max") -> jax.Array:
    """Down-pool along the sequence axis (axis=-2) with stride and window
    size equal to ``window`` and SAME padding.

    Input shape ``(..., S, D)`` becomes ``(..., ceil(S/window), D)``.

    For ``reduce="max"`` the pad value is ``-inf``; for ``"mean"``/``"avg"``
    padding is zero and the divisor is the *valid* count per window so the
    last (potentially partial) window averages over its real elements only.
    """
    if window == 1:
        return x
    if window <= 0:
        raise ValueError(f"`window` must be >= 1, got {window}")

    s, d = x.shape[-2], x.shape[-1]
    out_s = (s + window - 1) // window
    pad_total = out_s * window - s

    if pad_total > 0:
        if reduce == "max":
            pad_value = jnp.array(-jnp.inf, dtype=x.dtype)
        elif reduce in ("avg", "mean"):
            pad_value = jnp.array(0, dtype=x.dtype)
        else:
            raise NotImplementedError(f"Unknown reduce={reduce}")

        pad_width = [(0, 0)] * x.ndim
        pad_width[-2] = (0, pad_total)
        x = jnp.pad(x, pad_width, mode="constant", constant_values=pad_value)

    new_shape = (*x.shape[:-2], out_s, window, d)
    xw = jnp.reshape(x, new_shape)

    if reduce == "max":
        return jnp.max(xw, axis=-2)

    if reduce in ("avg", "mean"):
        summed = jnp.sum(xw, axis=-2)
        if pad_total == 0:
            return summed / jnp.array(window, dtype=x.dtype)
        counts = jnp.full((out_s,), window, dtype=x.dtype)
        counts = counts.at[-1].set(jnp.array(window - pad_total, dtype=x.dtype))
        reshape_counts = (1,) * (x.ndim - 2) + (out_s, 1)
        return summed / jnp.reshape(counts, reshape_counts)

    raise NotImplementedError(f"Unknown reduce={reduce}")


# ---------------------------------------------------------------------------
# Sequence stem
# ---------------------------------------------------------------------------


class TokenDnaEmbedder(nnx.Module):
    """Project token embeddings to ``conv_init_embed_dim`` channels.

    Takes per-token embeddings of dimension ``token_embed_dim`` and produces a
    ``conv_init_embed_dim``-channel feature stream over the same length
    sequence. A small residual conv branch with weight standardization adds
    higher-order interactions on top of the initial linear projection.
    """

    def __init__(
        self,
        token_embed_dim: int,
        conv_init_embed_dim: int,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ) -> None:
        self.conv1 = nnx.Conv(
            in_features=token_embed_dim,
            out_features=conv_init_embed_dim,
            kernel_size=(15,),
            padding="SAME",
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        # res_norm keeps fp32 via RMSBatchNorm's own floor, even though conv1's
        # output feeds straight into it.
        self.res_norm = RMSBatchNorm(
            num_features=conv_init_embed_dim,
            ndim=3,
        )
        self.res_conv = StandardizedConv1D(
            num_channels=conv_init_embed_dim,
            width=5,
            in_channels=conv_init_embed_dim,
            dtype=dtype,
            param_dtype=param_dtype,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.conv1(x)
        return x + self.res_conv(approximate_sigmoid_gelu(self.res_norm(x)))


# ---------------------------------------------------------------------------
# Convolutional blocks
# ---------------------------------------------------------------------------


class ConvBlock(nnx.Module):
    """Pre-norm convolutional block with ``gelu`` activation.

    Pattern: ``conv(approximate_sigmoid_gelu(rmsbatchnorm(x)))``. Used by both the encoder
    (``DownResBlock``) and decoder (``UpResBlock``).
    """

    def __init__(
        self,
        num_channels: int,
        width: int,
        *,
        in_channels: int,
        conditions_dims: Sequence[int] = (),
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ) -> None:
        self.num_channels = num_channels
        self.width = width

        # The pre-norm stays fp32 whatever the conv/proj compute dtype is.
        self.norm = _make_conv_norm(
            in_channels,
            conditions_dims,
            epsilon=1e-5,
            rngs=rngs,
        )

        if width == 1:
            self.proj = nnx.Linear(
                in_features=in_channels,
                out_features=num_channels,
                dtype=dtype,
                param_dtype=param_dtype,
                rngs=rngs,
            )
            self.conv = None
        else:
            self.conv = StandardizedConv1D(
                num_channels=num_channels,
                width=width,
                in_channels=in_channels,
                dtype=dtype,
                param_dtype=param_dtype,
            )
            self.proj = None

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = approximate_sigmoid_gelu(
            _apply_norm(self.norm, x, conditions, conditions_masks)
        )
        if self.width == 1:
            assert self.proj is not None
            return self.proj(x)
        assert self.conv is not None
        return self.conv(x)


class DownResBlock(nnx.Module):
    """Encoder residual block: two ``ConvBlock``s with a zero-padded skip.

    Channels can grow across the block; when ``out_channels > in_channels``
    the residual path is zero-padded along the feature axis to match. Pooling
    is *not* done here - the encoder loop pools externally.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ) -> None:
        if out_channels < in_channels:
            # Residual path is zero-padded along the feature axis to match the
            # block output, so this only handles channel growth (the encoder's
            # monotonic-widen invariant).
            raise ValueError(
                f"DownResBlock requires out_channels ({out_channels}) "
                f">= in_channels ({in_channels}); the residual path "
                "cannot shrink channels."
            )
        self._pad = out_channels - in_channels
        self.block1 = ConvBlock(
            out_channels,
            width=5,
            in_channels=in_channels,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.block2 = ConvBlock(
            out_channels,
            width=5,
            in_channels=out_channels,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        out = self.block1(x)
        if self._pad > 0:
            out = out + jnp.pad(x, [(0, 0), (0, 0), (0, self._pad)])
        else:
            out = out + x
        return out + self.block2(out)


class UpResBlock(nnx.Module):
    """Decoder upsampling block with U-Net skip and learnable residual scale.

    Sequence:

    1. ``conv_in`` over the trunk-side input ``x``, plus a channel-truncated
       residual.
    2. 2x nearest-neighbour upsampling along the sequence axis.
    3. Multiply by a learnable scalar ``residual_scale`` (init 1).
    4. Add a 1x1-projected U-Net skip from the matching encoder resolution.
    5. ``conv_out`` residual.

    Args:
        x_channels: Number of channels on the decoder-side input ``x`` (the
            trunk path being upsampled).
        unet_skip_channels: Number of channels on the U-Net skip connection
            coming from the matching encoder resolution. Also fixes the
            output channel count of this block (the trunk is widened or
            truncated to match this number before adding the skip).
        rngs: NNX random number generators for parameter initialisation.
    """

    def __init__(
        self,
        x_channels: int,
        unet_skip_channels: int,
        *,
        conditions_dims: Sequence[int] = (),
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ) -> None:
        num_channels = unet_skip_channels
        self.conv_in = ConvBlock(
            num_channels,
            width=5,
            in_channels=x_channels,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.pointwise_unet_skip = ConvBlock(
            num_channels,
            width=1,
            in_channels=unet_skip_channels,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.conv_out = ConvBlock(
            num_channels,
            width=5,
            in_channels=num_channels,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.residual_scale = nnx.Param(jnp.ones(()))

    def __call__(
        self,
        x: jax.Array,
        unet_skip: jax.Array,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        num_channels = unet_skip.shape[-1]
        out = self.conv_in(x, conditions, conditions_masks) + x[:, :, :num_channels]
        out = jnp.repeat(out, 2, axis=1)  # 2x upsample along sequence axis
        out = out * jnp.asarray(self.residual_scale.value, dtype=out.dtype)
        out = out + self.pointwise_unet_skip(unet_skip, conditions, conditions_masks)
        return out + self.conv_out(out, conditions, conditions_masks)


# ---------------------------------------------------------------------------
# Transformer blocks
# ---------------------------------------------------------------------------


class GQABlock(nnx.Module):
    """Grouped-query attention block.

    Features:

    * ``attention_heads`` query heads and ``num_kv_heads`` K/V heads
      (``attention_heads`` must be a multiple of ``num_kv_heads``). Default is
      GQA-2 (8 Q heads, 2 K/V heads, group size 4). Setting
      ``num_kv_heads=1`` recovers multi-query attention.
    * Per-head LayerNorm on each of Q, K, V.
    * Genomic rotary position embeddings (``apply_rope``) with frequencies
      tuned for ``max_position = _MAX_POSITION``.
    * Pre-norm and post-norm with ``nnx.RMSNorm`` (the post-norm is on the
      residual delta, before the add).
    * Attention is computed via ``jax.nn.dot_product_attention``. The
      ``attention_implementation`` argument selects the backend:

      - ``"xla"`` (default): the portable XLA implementation. Runs in the
        model dtype (fp32) and materializes the full S x S score matrix,
        so attention memory is O(S^2). This is what the model uses today.
      - ``"cudnn"``: the cuDNN FlashAttention kernel (memory ~O(S)). cuDNN
        only supports fp16/bf16, so Q/K/V are cast to bfloat16 for the
        kernel and the output is cast back to fp32. GPU-only and
        experimental here; validate numerics + the (B, 1, S, S) mask path
        on device before relying on it.

      Both paths require K and V head dims to be equal and broadcast K/V
      across Q head groups natively (so explicit head-tiling is unnecessary).

      Note: passing ``implementation=None`` to ``jax.nn.dot_product_attention``
      does NOT auto-select cuDNN in jax 0.5.3 (it falls back to XLA), so the
      backend is always set explicitly here.

    No tanh logit soft-cap: the per-head Q/K LayerNorms already keep logits
    bounded, and the cap is incompatible with the cuDNN kernel.
    """

    def __init__(
        self,
        *,
        embed_dim: int,
        attention_heads: int = 8,
        num_kv_heads: int = 2,
        key_size: int = 192,
        value_size: int = 192,
        max_position: int = _MAX_POSITION,
        attention_implementation: Literal["xla", "cudnn"] = "xla",
        conditions_dims: Sequence[int] = (),
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        epsilon: float = 1e-5,
    ) -> None:
        if attention_heads % num_kv_heads != 0:
            raise ValueError(
                f"attention_heads ({attention_heads}) must be divisible by "
                f"num_kv_heads ({num_kv_heads})"
            )
        if key_size != value_size:
            raise ValueError(
                "jax.nn.dot_product_attention requires key_size == "
                f"value_size; got {key_size} vs {value_size}."
            )
        if attention_implementation not in ("xla", "cudnn"):
            raise ValueError(
                "attention_implementation must be 'xla' or 'cudnn'; got "
                f"{attention_implementation!r}."
            )

        self.embed_dim = embed_dim
        self.attention_heads = attention_heads
        self.num_kv_heads = num_kv_heads
        self.key_size = key_size
        self.value_size = value_size
        self.max_position = max_position
        self.attention_implementation = attention_implementation

        self.norm_in = _make_transformer_norm(
            embed_dim,
            conditions_dims,
            epsilon=epsilon,
            rngs=rngs,
        )

        # QKVO projections: the matmul-heavy path `dtype` controls. The norms
        # below keep dtype=None, so Flax's promote_dtype upcasts them back to
        # fp32 from their own fp32 scale/bias.
        self.q_layer = nnx.Linear(
            embed_dim,
            attention_heads * key_size,
            use_bias=False,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.k_layer = nnx.Linear(
            embed_dim,
            num_kv_heads * key_size,
            use_bias=False,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.v_layer = nnx.Linear(
            embed_dim,
            num_kv_heads * value_size,
            use_bias=False,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )

        self.norm_q = nnx.LayerNorm(num_features=key_size, epsilon=epsilon, rngs=rngs)
        self.norm_k = nnx.LayerNorm(num_features=key_size, epsilon=epsilon, rngs=rngs)
        self.norm_v = nnx.LayerNorm(num_features=value_size, epsilon=epsilon, rngs=rngs)

        self.linear_embedding = nnx.Linear(
            attention_heads * value_size,
            embed_dim,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )

        self.norm_out = nnx.RMSNorm(num_features=embed_dim, epsilon=epsilon, rngs=rngs)

    def __call__(
        self,
        x: jax.Array,
        *,
        attention_mask: jax.Array | None = None,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        """Forward pass.

        Args:
            x: input of shape ``(B, S, embed_dim)``.
            attention_mask: optional boolean mask broadcastable to
                ``(B, 1, S, S)`` where ``True`` means "attend" and ``False``
                means "mask out". Passed directly to
                ``jax.nn.dot_product_attention``.
            conditions: optional per-example condition embeddings
                ``[(B, cond_dim_i)]`` consumed by the FiLM pre-norm.
            conditions_masks: optional per-example masks ``[(B,)]``.

        Returns:
            Output of shape ``(B, S, embed_dim)``.
        """
        b, s, _ = x.shape

        h = _apply_norm(self.norm_in, x, conditions, conditions_masks)

        q = self.q_layer(h).reshape(b, s, self.attention_heads, self.key_size)
        q = self.norm_q(q)
        k = self.k_layer(h).reshape(b, s, self.num_kv_heads, self.key_size)
        k = self.norm_k(k)
        v = self.v_layer(h).reshape(b, s, self.num_kv_heads, self.value_size)
        v = self.norm_v(v)

        q = apply_rope(q, None, max_position=self.max_position)
        k = apply_rope(k, None, max_position=self.max_position)

        # Default scale is 1/sqrt(key_size). When num_heads != num_kv_heads the
        # op broadcasts K/V across Q groups natively (ratio must be an integer;
        # enforced in __init__). ``attention_mask`` follows the
        # jax.nn.dot_product_attention convention: True = attend, False = mask
        # out. Cast defensively in case a caller passes int/bool-ish input.
        if attention_mask is not None and attention_mask.dtype != jnp.bool_:
            attention_mask = attention_mask.astype(bool)

        out_dtype = q.dtype
        if self.attention_implementation == "cudnn":
            # cuDNN FlashAttention only supports fp16/bf16; the rest of the
            # model is fp32, so cast Q/K/V into bf16 just for the kernel and
            # cast the result back. This is the only place lower precision is
            # used; softmax accumulation stays in fp32 inside the kernel.
            q = q.astype(jnp.bfloat16)
            k = k.astype(jnp.bfloat16)
            v = v.astype(jnp.bfloat16)

        y = jax.nn.dot_product_attention(
            query=q,
            key=k,
            value=v,
            mask=attention_mask,
            is_causal=False,
            implementation=self.attention_implementation,
        )
        y = y.astype(out_dtype)
        y = self.linear_embedding(y.reshape(b, s, -1))

        return self.norm_out(y)


class MLPBlock(nnx.Module):
    """Position-wise feed-forward block with GeGLU gating.

    Pattern: ``norm_in -> [gelu(W_g x)] * (W_v x) -> W_out -> norm_out``.

    Gated activation (GeGLU): the ``fc1_gate`` branch goes through ``gelu``
    and multiplies the ``fc1_value`` branch element-wise. Hidden width is
    ``4 * embed_dim / 3`` so total FFN params match a vanilla
    ``norm -> Linear(d, 2d) -> relu -> Linear(2d, d) -> norm`` block.

    Uses ``nnx.RMSNorm`` per token for the pre/post norms.
    """

    def __init__(
        self,
        embed_dim: int,
        *,
        conditions_dims: Sequence[int] = (),
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        epsilon: float = 1e-5,
    ) -> None:
        hidden = (embed_dim * 4) // 3
        self.norm1 = _make_transformer_norm(
            embed_dim,
            conditions_dims,
            epsilon=epsilon,
            rngs=rngs,
        )
        # fc1_gate/fc1_value/fc2 are the FFN's matmul path; norm2 is left at
        # dtype=None so it upcasts back to fp32 automatically (see GQABlock).
        self.fc1_gate = nnx.Linear(
            embed_dim, hidden, dtype=dtype, param_dtype=param_dtype, rngs=rngs
        )
        self.fc1_value = nnx.Linear(
            embed_dim, hidden, dtype=dtype, param_dtype=param_dtype, rngs=rngs
        )
        self.fc2 = nnx.Linear(
            hidden, embed_dim, dtype=dtype, param_dtype=param_dtype, rngs=rngs
        )
        self.norm2 = nnx.RMSNorm(num_features=embed_dim, epsilon=epsilon, rngs=rngs)

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        h = _apply_norm(self.norm1, x, conditions, conditions_masks)
        h = jax.nn.gelu(self.fc1_gate(h)) * self.fc1_value(h)
        h = self.fc2(h)
        return self.norm2(h)


def crop_sequence_center(
    sequences_array: jnp.ndarray, keep_target_center_fraction: float
) -> jnp.ndarray:
    """Crops an array to keep only the center portion specified by keep_center_fraction.

    Args:
        sequences_array: Input array of shape (batch_size, sequence_length, ...).
        keep_target_center_fraction: Fraction of the sequence to keep in the center
            (0.0 to 1.0). For example, 0.375 means keeping the middle 37.5% of the seq.

    Returns:
        Cropped array with the same number of dimensions but reduced sequence length.
    """
    assert 0.0 < keep_target_center_fraction <= 1.0, (
        f"{keep_target_center_fraction=} must be in (0, 1]"
    )
    sequence_length = sequences_array.shape[1]
    crop_length = int((1.0 - keep_target_center_fraction) * sequence_length)
    start_idx = crop_length // 2
    end_idx = sequence_length - crop_length // 2

    return sequences_array[:, start_idx:end_idx, ...]
