import math
from typing import List, Optional

import haiku as hk
import jax
from einops import rearrange
from jax import numpy as jnp

from nucleotide_transformer.borzoi.layers import get_positional_embed_borzoi


def get_positional_features_exponential(
    positions: jnp.ndarray, features: int, seq_len: int, min_half_life: float = 3.0
) -> jnp.ndarray:
    """
    Exponential features positional embeddings.

    Args:
        positions: positions.
        features: number of features.
        seq_len: sequence length.
        min_half_life: minimum half life.

    Returns:
        Positional embeddings.
    """
    max_range = math.log(seq_len) / math.log(2.0)
    half_life = 2 ** jnp.linspace(min_half_life, max_range, features)
    half_life = half_life[None, ...]
    positions = jnp.abs(positions)[..., None]
    return jnp.exp(-math.log(2.0) / half_life * positions)


def get_positional_features_central_mask(
    positions: jnp.ndarray, features: int
) -> jnp.ndarray:
    """
    Exponential features central mask.

    Args:
        positions: positions.
        features: number of features.

    Returns:
        Positional embeddings.
    """
    center_widths = 2 ** jnp.arange(1, features + 1).astype(jnp.float32)
    center_widths = center_widths - 1
    return (center_widths[None, ...] > jnp.abs(positions)[..., None]).astype(
        jnp.float32
    )


def gamma_pdf(
    x: jnp.ndarray, concentration: jnp.ndarray, rate: jnp.ndarray
) -> jnp.ndarray:
    """
    Gamma law PDF function.

    Args:
        x: input tensor.
        concentration: gamma concentration.
        rate: gamma rate.

    Returns:
        gamma pdf value.
    """
    log_unnormalized_prob = jax.scipy.special.xlogy(concentration - 1.0, x) - rate * x
    log_normalization = jax.lax.lgamma(concentration) - concentration * jnp.log(rate)
    return jnp.exp(log_unnormalized_prob - log_normalization)


def get_positional_features_gamma(
    positions: jnp.ndarray,
    features: int,
    seq_len: int,
    stddev: Optional[float] = None,
    start_mean: Optional[float] = None,
    eps: float = 1e-8,
) -> jnp.ndarray:
    """
    Get Gamma positional features.
    """
    if stddev is None:
        stddev = seq_len / (2 * features)

    if start_mean is None:
        start_mean = seq_len / features

    mean = jnp.linspace(start_mean, seq_len, features)
    mean = mean[None, ...]
    concentration = (mean / stddev) ** 2
    rate = mean / stddev**2
    probabilities = gamma_pdf(
        jnp.abs(positions.astype(jnp.float32))[..., None], concentration, rate
    )
    probabilities = probabilities + eps
    outputs = probabilities / jnp.amax(probabilities, axis=-1, keepdims=True)
    return outputs


def get_positional_embed_enformer(seq_len: int, feature_size: int) -> jnp.ndarray:
    """
    Compute positional embedding.
    """
    distances = jnp.arange(-seq_len + 1, seq_len)

    feature_functions = [
        get_positional_features_exponential,
        get_positional_features_central_mask,
        get_positional_features_gamma,
    ]

    num_components = len(feature_functions) * 2

    if (feature_size % num_components) != 0:
        raise ValueError(
            f"feature size is not divisible by number of components ({num_components})"
        )

    num_basis_per_class = feature_size // num_components

    embeddings = []
    embeddings.append(
        get_positional_features_exponential(distances, num_basis_per_class, seq_len)
    )
    embeddings.append(
        get_positional_features_central_mask(distances, num_basis_per_class)
    )
    embeddings.append(
        get_positional_features_gamma(distances, num_basis_per_class, seq_len)
    )

    embeddings = jnp.concatenate(embeddings, axis=-1)
    embeddings = jnp.concatenate(
        (embeddings, jnp.sign(distances)[..., None] * embeddings), axis=-1
    )
    return embeddings


def relative_shift(x: jnp.ndarray) -> jnp.ndarray:
    """
    Apply relative shift.
    """
    to_pad = jnp.zeros_like(x[..., :1])
    x = jnp.concatenate((to_pad, x), axis=-1)
    _, h, t1, t2 = x.shape
    x = x.reshape(-1, h, t2, t1)
    x = x[:, :, 1:, :]
    x = x.reshape(-1, h, t1, t2 - 1)
    return x[..., : ((t2 + 1) // 2)]


def exponential_linspace_int(
    start: int, end: int, num: int, divisible_by: int = 1
) -> List[int]:
    """
    Create list of dimensions to construct the Enformer model.
    """

    def _round(x: float) -> int:
        return int(round(x / divisible_by) * divisible_by)

    base = math.exp(math.log(end / start) / (num - 1))
    return [_round(start * base**i) for i in range(num)]


def gelu_fn(x: jnp.ndarray) -> jnp.ndarray:
    """
    Custom GELU activation function.
    """
    return jax.nn.sigmoid(1.702 * x) * x


class Attention(hk.Module):
    """
    Enformer Attention Layer.
    """

    def __init__(
        self,
        dim: int,
        *,
        num_rel_pos_features: int,
        heads: int = 8,
        dim_key: int = 64,
        dim_value: int = 64,
        positional_encoding_type: str = "enformer",
        name: Optional[str] = None,
    ):
        super().__init__(name=name)
        self.scale = dim_key**-0.5
        self.heads = heads

        self.to_q = hk.Linear(
            output_size=dim_key * heads, with_bias=False, name="attn_q"
        )
        self.to_k = hk.Linear(
            output_size=dim_key * heads, with_bias=False, name="attn_k"
        )
        self.to_v = hk.Linear(
            output_size=dim_value * heads, with_bias=False, name="attn_v"
        )

        self.to_out = hk.Linear(output_size=dim, name="attn_o")

        # relative positional encoding
        self.num_rel_pos_features = num_rel_pos_features

        self.to_rel_k = hk.Linear(
            output_size=dim_key * heads, with_bias=False, name="attn_to_rel_k"
        )

        w_init = hk.initializers.RandomNormal()
        self.rel_content_bias = hk.get_parameter(
            "rel_content_bias",
            shape=(1, heads, 1, dim_key),
            init=w_init,
        )
        w_init = hk.initializers.RandomNormal()
        self.rel_pos_bias = hk.get_parameter(
            "rel_pos_bias",
            shape=(1, heads, 1, dim_key),
            init=w_init,
        )

        # Set the way the position encodings are computed
        if positional_encoding_type == "enformer":
            self.get_positional_embed = get_positional_embed_enformer
        elif positional_encoding_type == "borzoi":
            self.get_positional_embed = get_positional_embed_borzoi

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        n = x.shape[-2]

        q = self.to_q(x)
        k = self.to_k(x)
        v = self.to_v(x)

        def _rearrange(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "b n (h d) -> b h n d", h=self.heads)

        q = _rearrange(q)
        k = _rearrange(k)
        v = _rearrange(v)

        q = q * self.scale

        content_logits = jnp.einsum(
            "b h i d, b h j d -> b h i j", q + self.rel_content_bias, k
        )

        positions = self.get_positional_embed(n, self.num_rel_pos_features)
        rel_k = self.to_rel_k(positions)

        def _rearrange_1(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "n (h d) -> h n d", h=self.heads)

        rel_k = _rearrange_1(rel_k)
        rel_logits = jnp.einsum(
            "b h i d, h j d -> b h i j", q + self.rel_pos_bias, rel_k
        )
        rel_logits = relative_shift(rel_logits)

        logits = content_logits + rel_logits
        attn = jax.nn.softmax(logits, axis=-1)

        out = jnp.einsum("b h i j, b h j d -> b h i d", attn, v)

        def _rearrange_2(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "b h n d -> b n (h d)")

        out = _rearrange_2(out)

        return self.to_out(out)


class AttentionPool(hk.Module):
    """Enformer Attention Pooling layer."""

    def __init__(
        self,
        dim: int,
        pool_size: int = 2,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            pool_size: pooling size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._pool_size = pool_size

        self._to_attn_logits = hk.Conv2D(
            output_channels=dim, kernel_shape=1, with_bias=False, data_format="NCHW"
        )

    def _pool_fn(self, x: jnp.ndarray) -> jnp.ndarray:
        b, d, n = x.shape
        x = jnp.reshape(x, (b, d, n // self._pool_size, self._pool_size))
        return x

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:

        b, d, n = x.shape
        remainder = n % self._pool_size
        needs_padding = remainder > 0

        if needs_padding:
            x = jnp.concatenate(
                [x, jnp.zeros(shape=(b, d, remainder), dtype=x.dtype)], axis=-1
            )
            mask = jnp.zeros((b, 1, n), dtype=jnp.bool_)
            mask = jnp.concatenate(
                [mask, jnp.ones(shape=(b, 1, remainder), dtype=jnp.bool_)], axis=-1
            )

        x = self._pool_fn(x)
        logits = self._to_attn_logits(x)

        if needs_padding:
            mask_value = -jnp.inf
            logits = jnp.where(
                self._pool_fn(mask),
                mask_value * jnp.ones_like(logits, dtype=logits.dtype),
                logits,
            )

        attn = jax.nn.softmax(logits, axis=-1)
        return jnp.sum((x * attn), axis=-1)


class ResidualConvBlock(hk.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._dim_out = dim_out
        self._kernel_size = kernel_size

    def __call__(self, x: jnp.ndarray, is_training: bool = False) -> jnp.ndarray:
        conv_block = ConvBlock(
            dim=self._dim, dim_out=self._dim_out, kernel_size=self._kernel_size
        )
        return x + conv_block(x, is_training)


class ConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._dim_out = dim_out
        self._kernel_size = kernel_size

    def __call__(self, x: jnp.ndarray, is_training: bool = False) -> jnp.ndarray:
        batch_norm = hk.BatchNorm(
            create_scale=True,
            create_offset=True,
            decay_rate=0.9,
            data_format="NC...",
        )
        conv = hk.Conv1D(
            output_channels=self._dim if self._dim_out is None else self._dim_out,
            kernel_shape=self._kernel_size,
            padding=(self._kernel_size // 2, self._kernel_size // 2),
            data_format="NCW",
        )

        x = batch_norm(x, is_training=is_training)
        x = gelu_fn(x)
        x = conv(x)
        return x
