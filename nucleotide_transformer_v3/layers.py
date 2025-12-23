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

"""Layer implementations for Nucleotide Transformer v3."""

import logging
from dataclasses import dataclass
from enum import Enum
from functools import partial
from typing import Any, Callable, Sequence

import jax
import jax.nn
import jax.numpy as jnp
import numpy as np
from flax import nnx
from flax.nnx.nn import dtypes
from flax.typing import Dtype
from jax import lax

from nucleotide_transformer_v3.types import AttentionMask, Embedding, Tokens, TransformerOutput
from nucleotide_transformer_v3.utils import get_activation_fn

logger = logging.getLogger(__name__)

# =============================================================================
# Rotary Embeddings
# =============================================================================

# Constant used in Sinusoidal/Rotary Embeddings, reference to this value can be found
# on page 6 of https://arxiv.org/pdf/1706.03762.pdf and page 5 of
# https://arxiv.org/abs/2104.09864
# These rotary positional embeddings are proper to ESM implementation.
# dimensions in key space are rotated 2 by 2. The key difference with
# GPT's one is that in this case each 2 dimensions rotated together are spaced
# by key_size//2
UPPER_FREQ = 10000


@dataclass
class RotaryEmbeddingConfig:
    """
    Parameters to initialize the RotaryEmbedding layer. The rescaling factor allows
    to adapt the rotary embeddings to larger lengths than what was used for training.
    One of this strategy is presented in the Yarn paper: https://arxiv.org/pdf/2309.00071.pdf. # noqa

    Args:
        rescaling_factor: Rescaling factor for the rotary embeddings.
    """

    rescaling_factor: float | None


class RotaryEmbedding(nnx.Module):
    """
    Rotary Positional Embedding inspired by RoFormer:
    https://arxiv.org/abs/2104.09864
    https://github.com/ZhuiyiTechnology/roformer .
    """

    def __init__(
        self,
        key_size: int,
        rotary_embedding_config: RotaryEmbeddingConfig,
    ):
        """
        Args:
            key_size: Dimension of one head.
            rotary_embedding_config: Configuration to specify hyperparameters for
                RotaryEmbedding layer
                (see RoFormer https://arxiv.org/pdf/2104.09864.pdf). It contains
                the hyperparameters specifying the type of rotary embedding applied.
        """
        # Extract argument from the config
        self._rescaling_factor = rotary_embedding_config.rescaling_factor
        self._key_size = key_size

    def get_inv_freq(self) -> jax.Array:
        if self._rescaling_factor is None:
            return 1.0 / (
                UPPER_FREQ ** (jnp.arange(0, self._key_size, 2) / self._key_size)
            )
        else:
            updated_base = UPPER_FREQ * (
                self._rescaling_factor ** (self._key_size / (self._key_size - 2))
            )
            return 1.0 / (
                updated_base ** (jnp.arange(0, self._key_size, 2) / self._key_size)
            )

    def _compute_cos_sin_tables(
        self,
        heads: jax.Array,
    ) -> tuple[jax.Array, jax.Array]:
        """
        Computes the cosine and sine for rotation.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
            key_size).

        Returns:
            Cosine positional embedding of shape (1, seq_len, 1, key_size/2).
            Sine positional embedding of shape (1, seq_len, 1, key_size/2).
        """
        seq_len = heads.shape[1]

        t = jnp.arange(seq_len)
        inv_freq = self.get_inv_freq()
        freqs = jnp.einsum("i,j->ij", t, inv_freq)

        # Compute cos and cast is as (1, seq_len, 1, key_size/2) to be applied to
        # queries of shape (batch_size, seq_len, num_heads, key_size/2)
        cos_cached = jnp.cos(freqs)[None, :, None, :]
        sin_cached = jnp.sin(freqs)[None, :, None, :]

        return cos_cached, sin_cached

    def _apply_rotary_pos_emb(
        self, heads: jax.Array, cos: jax.Array, sin: jax.Array
    ) -> jax.Array:
        """
        Applies the rotary positional embedding to the heads.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
                key_size).
            cos: Cosine values.
            sin: Sine values.

        Returns:
            Embedded heads of shape (batch_size, seq_len, num_heads, key_size).
        """

        # Rotate x
        x_first, x_second = (
            heads[..., : heads.shape[-1] // 2],
            heads[..., heads.shape[-1] // 2 :],
        )
        first_part = x_first * cos - x_second * sin
        second_part = x_second * cos + x_first * sin

        return jnp.concatenate((first_part, second_part), axis=-1, dtype=heads.dtype)

    def __call__(
        self, query_heads: jax.Array, key_heads: jax.Array
    ) -> tuple[jax.Array, jax.Array]:
        """
        Applies rotary embeddings to query_heads and key_heads.

        Args:
            query_heads: Query heads of shape
                (batch_size, seq_len, num_heads, key_size).
            key_heads: Key heads of shape (batch_size, seq_len, num_heads, key_size).

        Returns:
            Embedded query heads.
            Embedded key heads.
        """
        cos, sin = self._compute_cos_sin_tables(query_heads)

        return (
            self._apply_rotary_pos_emb(query_heads, cos, sin),
            self._apply_rotary_pos_emb(key_heads, cos, sin),
        )


# =============================================================================
# Convolution Transpose Layers
# =============================================================================


class Conv1DTranspose(nnx.Module):
    """General n-dimensional transposed convolution (aka. deconvolution)."""

    def __init__(
        self,
        input_channels: int,
        output_channels: int,
        kernel_shape: int,
        stride: int = 1,
        padding: str | Sequence[tuple[int, int]] = "SAME",
        with_bias: bool = True,
        w_init: Callable | None = None,
        b_init: Callable | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        precision: lax.Precision | None = None,
        rngs: nnx.Rngs,
    ):
        self.kernel_shape = (kernel_shape,)
        self.with_bias = with_bias
        self.stride = (stride,)
        self.padding = padding
        self.dtype = dtype
        self.precision = precision

        # Assumes channels_last
        self.dimension_numbers = lax.ConvDimensionNumbers(
            lhs_spec=(0, 2, 1), rhs_spec=(1, 2, 0), out_spec=(0, 2, 1)
        )

        # Create weight parameter
        w_shape = self.kernel_shape + (output_channels, input_channels)
        if w_init is None:
            fan_in_shape = self.kernel_shape + (input_channels,)
            stddev = 1.0 / np.sqrt(np.prod(fan_in_shape))
            w_init = nnx.initializers.truncated_normal(stddev=stddev)

        assert w_init is not None
        self.kernel = nnx.Param(w_init(rngs.params(), w_shape, param_dtype))

        # Create bias parameter
        if self.with_bias:
            bias_shape = (output_channels,)
            if b_init:
                b = nnx.Param(b_init(rngs.params(), bias_shape, param_dtype))
            else:
                b = nnx.Param(jnp.zeros(bias_shape, param_dtype))
            self.bias = b

    def __call__(
        self,
        inputs: jax.Array,
    ) -> jax.Array:
        # Assume inputs of shape (batch_size, seq_len, num_features)
        inputs, kernel, bias = dtypes.promote_dtype(
            (inputs, self.kernel.value, self.bias.value), dtype=self.dtype
        )

        out = lax.conv_transpose(
            lhs=inputs,
            rhs=kernel,
            strides=self.stride,
            padding=self.padding,
            dimension_numbers=self.dimension_numbers,
            precision=self.precision,
        )

        # Add bias if needed
        if self.with_bias:
            out = out + bias

        return out


class UpsampleConv1D(nnx.Conv):
    """Generic upsample convolution that supports different upsampling types."""

    def __init__(
        self,
        input_channels: int,
        output_channels: int,
        kernel_shape: int,
        upsample_factor: int,
        padding: str | Sequence[tuple[int, int]] = "SAME",
        with_bias: bool = True,
        w_init: Callable | None = None,
        b_init: Callable | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        precision: lax.Precision | None = None,
        rngs: nnx.Rngs,
    ):
        self.input_channels = input_channels
        self.output_channels = output_channels
        self.kernel_size_int = kernel_shape
        self.upsample_factor = upsample_factor
        self.padding = padding
        self.with_bias = with_bias
        self.dtype = dtype

        if w_init is None:
            fan_in_shape = (kernel_shape,) + (self.input_channels,)
            stddev = 1.0 / np.sqrt(np.prod(fan_in_shape))
            _w_init = nnx.initializers.truncated_normal(stddev=stddev)
        else:
            _w_init = w_init

        if b_init is None:
            _b_init = nnx.initializers.zeros
        else:
            _b_init = b_init

        super().__init__(
            in_features=self.input_channels,
            out_features=self.output_channels,
            kernel_size=(kernel_shape,),
            strides=(1,),
            padding=self.padding,
            use_bias=self.with_bias,
            kernel_init=_w_init,
            bias_init=_b_init,
            dtype=self.dtype,
            param_dtype=param_dtype,
            precision=precision,
            rngs=rngs,
            feature_group_count=1,
        )

    def __call__(
        self,
        inputs: jax.Array,
    ) -> jax.Array:
        # Assume inputs of shape (batch_size, seq_len, num_features)

        # Upsample along the sequence length dimension (axis=1)
        if self.upsample_factor > 1:
            conv_input = jnp.repeat(inputs, self.upsample_factor, axis=1)
        else:
            conv_input = inputs

        out = super().__call__(conv_input)

        return out


class DeConvUpsampleType(str, Enum):
    CONV_TRANSPOSE = "conv_transpose"
    REPEAT_CONV = "repeat+conv"


def generic_upsample_conv1d(
    input_channels: int,
    output_channels: int,
    kernel_shape: int,
    upsample: DeConvUpsampleType | None = None,
    padding: str | Sequence[tuple[int, int]] = "SAME",
    with_bias: bool = True,
    w_init: Callable | None = None,
    b_init: Callable | None = None,
    *,
    dtype: Dtype | None = None,
    param_dtype: Dtype = jnp.float32,
    precision: lax.Precision | None = None,
    rngs: nnx.Rngs,
) -> Conv1DTranspose | UpsampleConv1D:
    """
    Generic upsample convolution that supports different upsampling types.
    """
    if upsample is None:
        return Conv1DTranspose(
            input_channels=input_channels,
            output_channels=output_channels,
            kernel_shape=kernel_shape,
            stride=1,
            padding=padding,
            with_bias=with_bias,
            w_init=w_init,
            b_init=b_init,
            dtype=dtype,
            param_dtype=param_dtype,
            precision=precision,
            rngs=rngs,
        )
    elif upsample == DeConvUpsampleType.REPEAT_CONV:
        return UpsampleConv1D(
            input_channels=input_channels,
            output_channels=output_channels,
            kernel_shape=kernel_shape,
            upsample_factor=2,
            padding=padding,
            with_bias=with_bias,
            w_init=w_init,
            b_init=b_init,
            dtype=dtype,
            param_dtype=param_dtype,
            precision=precision,
            rngs=rngs,
        )
    elif upsample == DeConvUpsampleType.CONV_TRANSPOSE:
        return Conv1DTranspose(
            input_channels=input_channels,
            output_channels=output_channels,
            kernel_shape=kernel_shape,
            stride=2,
            padding=padding,
            with_bias=with_bias,
            w_init=w_init,
            b_init=b_init,
            dtype=dtype,
            param_dtype=param_dtype,
            precision=precision,
            rngs=rngs,
        )
    else:
        raise ValueError(f"Invalid upsample type: {upsample}")


# =============================================================================
# Multi-Head Attention
# =============================================================================


class LinearProjectionHeInit(nnx.Module):
    def __init__(
        self,
        num_heads: int,
        key_size: int,
        rngs: nnx.Rngs,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
    ):
        initializer = nnx.initializers.variance_scaling(2.0, "fan_in", "uniform")
        self.linear = nnx.Linear(
            in_features=num_heads * key_size,
            out_features=num_heads * key_size,
            rngs=rngs,
            kernel_init=initializer,
            dtype=dtype,
            param_dtype=param_dtype,
        )
        self._num_heads = num_heads
        self._key_size = key_size

    def __call__(self, x: jax.Array) -> jax.Array:
        y = self.linear(x)
        return y.reshape((*x.shape[:-1], self._num_heads, self._key_size))


class MultiHeadAttention(nnx.Module):
    """
    Multi-head attention with masking applied. Modified from the core implementation to
    support biases in keys and values.
    """

    def __init__(
        self,
        num_heads: int,
        key_size: int,
        rotary_embedding_config: RotaryEmbeddingConfig | None = None,
        add_bias_kv: bool = False,
        value_size: int | None = None,
        model_size: int | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            num_heads: Number of independent attention heads.
            key_size: The size of keys and queries used for attention.
            rotary_embedding_config: Configuration to specify hyperparameters for
                RotaryEmbedding layer
                (see RoFormer https://arxiv.org/pdf/2104.09864.pdf). If None,
                rotary embeddings are not used. If specified, it contains the
                hyperparameters specifying the type of rotary embedding applied.
            add_bias_kv: If True, appends biases to key and query heads, used in ESM
                model (https://www.biorxiv.org/content/10.1101/622803v4.full.pdf).
            value_size: Optional size of the value projection. If None, defaults
                to the key size.
            model_size: Optional size of the output embedding. If None, defaults
                to the key size multiplied by the number of heads.
            rngs: Random number generators.
        """
        if not value_size:
            value_size = key_size

        if not model_size:
            model_size = key_size * num_heads

        if add_bias_kv:
            self._bias_k = nnx.Param(
                jnp.zeros([1, 1, num_heads, key_size], dtype=dtype)
            )
            self._bias_v = nnx.Param(
                jnp.zeros([1, 1, num_heads, value_size], dtype=dtype)
            )
        else:
            self._bias_k = None
            self._bias_v = None

        self._rotary_embedding_config = rotary_embedding_config
        self._key_size = key_size

        self.query_head = LinearProjectionHeInit(
            num_heads,
            key_size,
            rngs,
            dtype=dtype,
            param_dtype=param_dtype,
        )
        self.key_head = LinearProjectionHeInit(
            num_heads,
            key_size,
            rngs,
            dtype=dtype,
            param_dtype=param_dtype,
        )
        self.value_head = LinearProjectionHeInit(
            num_heads,
            value_size,
            rngs,
            dtype=dtype,
            param_dtype=param_dtype,
        )

        initializer = nnx.initializers.variance_scaling(2.0, "fan_in", "uniform")
        self.mha_output = nnx.Linear(
            in_features=model_size,
            out_features=model_size,
            kernel_init=initializer,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )

        if self._rotary_embedding_config:
            self.rotary_embedding = RotaryEmbedding(
                key_size,
                rotary_embedding_config,  # type: ignore
            )

    def attention_weights(
        self,
        query: jax.Array,
        key: jax.Array,
        attention_mask: AttentionMask | None = None,
        attention_weight_bias: jax.Array | None = None,
    ) -> jax.Array:
        """
        Computes the attention weights.

        Args:
            query: Embedding sequence to compute queries.
            key: Embedding sequence to compute keys.
            attention_mask: Input attention_mask. Defaults to None.
            attention_weight_bias: Optional bias to add to the attention weights.

        Returns:
            Attention weights.
        """
        query_heads = self.query_head(query)
        key_heads = self.key_head(key)

        if self._bias_k is not None:
            batch_size = key_heads.shape[0]
            attention_bias = jnp.tile(self._bias_k, (batch_size, 1, 1, 1))
            key_heads = jnp.concatenate((key_heads, attention_bias), axis=1)
            if attention_mask is not None:
                attention_mask = jnp.concatenate(
                    (
                        attention_mask,
                        jnp.ones(attention_mask.shape[:-1] + (1,), dtype=jnp.bool_),
                    ),
                    axis=-1,
                )

        if self._rotary_embedding_config:
            query_heads, key_heads = self.rotary_embedding(query_heads, key_heads)

        attention_logits = jnp.einsum("...thd,...Thd->...htT", query_heads, key_heads)
        sqrt_key_size = jnp.sqrt(self._key_size).astype(query.dtype)
        attention_logits = attention_logits / sqrt_key_size

        if attention_mask is not None:
            assert len(attention_mask.shape) == len(attention_logits.shape)
            attention_logits = jnp.where(attention_mask, attention_logits, -1e30)

        if attention_weight_bias is None:
            attention_weights = jax.nn.softmax(attention_logits)
        else:
            attention_weights = jax.nn.softmax(attention_logits + attention_weight_bias)

        return attention_weights

    def compute_embeddings(
        self,
        value: jax.Array,
        attention_weights: jax.Array,
    ) -> jax.Array:
        """
        Computes the output embeddings.

        Args:
            value: Embedding sequence to compute values.
            attention_weights: Attention weights.

        Returns:
            Output embeddings.
        """
        value_heads = self.value_head(value)

        if self._bias_v is not None:
            batch_size = value_heads.shape[0]
            attention_bias = jnp.tile(
                self._bias_v,
                (batch_size, 1, 1, 1),
            )
            value_heads = jnp.concatenate((value_heads, attention_bias), axis=1)

        attention = jnp.einsum("...htT,...Thd->...thd", attention_weights, value_heads)

        # Concatenate attention maps of all heads into a single vector.
        attention_vec = jnp.reshape(attention, (*attention.shape[:-2], -1))
        return self.mha_output(attention_vec)

    def __call__(
        self,
        query: jax.Array,
        key: jax.Array,
        value: jax.Array,
        attention_mask: jax.Array | None = None,
        attention_weight_bias: jax.Array | None = None,
    ) -> TransformerOutput:
        """
        Computes both the embeddings and the attention weights.

        Args:
            query: Embedding sequence to compute queries.
            key: Embedding sequence to compute keys.
            value: Embedding sequence to compute values.
            attention_mask: Mask to be applied during the attention layers.
                Triangular for autoregressive models. Defaults to None.

        Returns:
            Dictionary containing the output embeddings and the attention weights.
        """

        attention_weights = self.attention_weights(
            query,
            key,
            attention_mask=attention_mask,
            attention_weight_bias=attention_weight_bias,
        )
        embeddings = self.compute_embeddings(value, attention_weights)

        return {"embeddings": embeddings, "attention_weights": attention_weights}


# =============================================================================
# Self-Attention Block
# =============================================================================


class SelfAttentionBlock(nnx.Module):
    """
    Attention block made of self-attention.
    """

    def __init__(
        self,
        num_heads: int,
        embed_dim: int,
        ffn_embed_dim: int,
        key_size: int | None = None,
        add_bias_kv: bool = False,
        add_bias_fnn: bool = True,
        ffn_activation_name: str = "gelu-no-approx",
        use_glu_in_ffn: bool = False,
        layer_norm_eps: float = 1e-5,  # this is the default haiku value
        pre_layer_norm: bool = True,
        rotary_embedding_config: RotaryEmbeddingConfig | None = None,
        *,
        ffn_dtype: Dtype | None = None,
        ffn_param_dtype: Dtype = jnp.float32,
        mha_dtype: Dtype | None = None,
        mha_param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            num_heads: Number of attention heads.
            embed_dim: Embedding dimension.
            ffn_embed_dim: Feed forward network embedding dimension.
            key_size: Key size for the attention layer. If None, it is set to
                embed_dim // num_heads.
            add_bias_kv: Whether to add bias to the key and value sequences at
                each attention head.
            add_bias_fnn: Whether to add bias in the feed forward network.
            ffn_activation_name: Activation function name for the feed forward
                network. Default is "gelu-no-approx".
            use_glu_in_ffn: Whether to use GLU in the feed forward network.
            layer_norm_eps: Epsilon value for the layer normalization.
            pre_layer_norm: Whether to apply layer normalization before or after
                the self-attention and feed forward network layers.
            rotary_embedding_config: Configuration for rotary embeddings.
            ffn_dtype: Dtype for the feed forward network.
            ffn_param_dtype: Dtype for the feed forward network parameters.
            mha_dtype: Dtype for the multi-head attention.
            mha_param_dtype: Dtype for the multi-head attention parameters.
            ln_dtype: Dtype for the layer normalization.
            ln_param_dtype: Dtype for the layer normalization parameters.
        """
        # Add checks on dimensions
        if key_size is None:
            if embed_dim % num_heads != 0:
                raise ValueError(
                    f"The embedding dimension should be divisible by the number of "
                    f"heads, however provided embedding dimension is {embed_dim} and "
                    f"the number of heads is {num_heads}."
                )
            else:
                key_size = embed_dim // num_heads

        # Get ffn activation function
        self._pre_layer_norm = pre_layer_norm
        self._ffn_activation_fn = get_activation_fn(activation_name=ffn_activation_name)
        self._use_glu_in_fnn = use_glu_in_ffn

        # Define layers
        self.fc1 = nnx.Linear(
            in_features=embed_dim,
            out_features=2 * ffn_embed_dim if use_glu_in_ffn else ffn_embed_dim,
            use_bias=add_bias_fnn,
            dtype=ffn_dtype,
            param_dtype=ffn_param_dtype,
            rngs=rngs,
        )
        self.fc2 = nnx.Linear(
            in_features=ffn_embed_dim,
            out_features=embed_dim,
            use_bias=add_bias_fnn,
            dtype=ffn_dtype,
            param_dtype=ffn_param_dtype,
            rngs=rngs,
        )
        self.self_attention_layer_norm = nnx.LayerNorm(
            num_features=embed_dim,
            epsilon=layer_norm_eps,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )
        self.final_layer_norm = nnx.LayerNorm(
            num_features=embed_dim,
            epsilon=layer_norm_eps,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )

        self.sa_layer = MultiHeadAttention(
            num_heads=num_heads,
            key_size=key_size,
            add_bias_kv=add_bias_kv,
            model_size=embed_dim,
            rotary_embedding_config=rotary_embedding_config,
            dtype=mha_dtype,
            param_dtype=mha_param_dtype,
            rngs=rngs,
        )

    def self_attention(
        self,
        x: Embedding,
        attention_mask: AttentionMask | None = None,
        attention_weight_bias: jax.Array | None = None,
    ) -> TransformerOutput:
        """
        Applies the self attention mechanism.

        Args:
            x: Input token embeddings of shape (batch_size, seq_len, embed_dim).
            attention_mask: Optional attention mask of
                shape (batch_size, 1, seq_len, seq_len).
            attention_weight_bias: Optional bias to add to the attention weights.

        Returns:
            Dictionary containing the output embeddings and the attention weights.
        """

        return self.sa_layer(
            x,
            x,
            x,
            attention_mask=attention_mask,
            attention_weight_bias=attention_weight_bias,
        )

    def mlp(self, embed: Embedding) -> Embedding:
        """
        Applies one layer-norm, one linear layer, a Gelu activation,
        then a final linear layer.

        Args:
            embed: Embeddings of shape (batch_size, seq_len, key_size * num_heads).

        Returns:
            The transformed sequence embedding.
        """

        if self._pre_layer_norm:
            x = self.final_layer_norm(embed)
        else:
            x = embed

        if self._use_glu_in_fnn:
            x1, x2 = jnp.split(self.fc1(x), indices_or_sections=2, axis=-1)
            x = self._ffn_activation_fn(x1) * x2
        else:
            x = self._ffn_activation_fn(self.fc1(x))

        x = self.fc2(x)

        if not self._pre_layer_norm:
            x = self.final_layer_norm(x + embed)

        return x

    def __call__(
        self,
        x: Tokens,
        attention_mask: AttentionMask | None = None,
        attention_weight_bias: jax.Array | None = None,
    ) -> TransformerOutput:
        """
        Computes the output of the attention layer.

        Args:
            x: Input token embeddings of shape (batch_size,seq_len,embed_dim).
            attention_mask: Optional attention mask of
                shape (batch_size, 1, seq_len, seq_len).
            attention_weight_bias: Optional bias to add to the attention weights.

        Returns:
            A dictionary containing the output embeddings and the attention weights.
        """

        # Self-Attention
        res = x
        if self._pre_layer_norm:
            x = self.self_attention_layer_norm(x)

        output = self.self_attention(
            x=x,
            attention_mask=attention_mask,
            attention_weight_bias=attention_weight_bias,
        )

        if not self._pre_layer_norm:
            output["embeddings"] = self.self_attention_layer_norm(
                output["embeddings"] + res
            )
            x = output["embeddings"]
        else:
            x = output["embeddings"]
            x = res + x

        # MLP
        if not self._pre_layer_norm:
            x = self.mlp(x)
        else:
            x = x + self.mlp(x)

        output["embeddings"] = x
        return output  # type: ignore


# =============================================================================
# Convolution Blocks
# =============================================================================


class ResidualConvBlock(nnx.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        dim_out: int | None = None,
        kernel_size: int = 1,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            rngs: Random number generators.
        """
        self.conv_block = ConvBlock(
            dim=dim,
            dim_out=dim_out,
            kernel_size=kernel_size,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
            **kwargs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        return x + self.conv_block(x)


class ConvBlock(nnx.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: int | None = None,
        kernel_size: int = 1,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            rngs: Random number generators.
        """
        # Extract layer_norm_dim from kwargs with default value of dim
        layer_norm_dim = kwargs.pop("layer_norm_dim", dim)
        self.conv = nnx.Conv(  # assumes NWC format
            in_features=dim,
            out_features=dim if dim_out is None else dim_out,
            kernel_size=kernel_size,
            padding="SAME",
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.layer_norm = nnx.LayerNorm(
            num_features=layer_norm_dim,
            epsilon=1e-5,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.layer_norm(x)
        x = self.conv(x)
        x = jax.nn.gelu(x)
        return x


class ResidualDeConvBlock(nnx.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        dim_out: int | None = None,
        kernel_size: int = 1,
        upsample: DeConvUpsampleType | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            upsample: Upsampling type.
            rngs: Random number generators.
        """
        self.conv_block = DeConvBlock(
            dim=dim,
            dim_out=dim_out,
            kernel_size=kernel_size,
            upsample=upsample,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
            **kwargs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        return x + self.conv_block(x)


class DeConvBlock(nnx.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: int | None = None,
        kernel_size: int = 1,
        upsample: DeConvUpsampleType | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            upsample: Upsampling type.
            rngs: Random number generators.
        """
        # Extract layer_norm_dim from kwargs with default value of dim
        layer_norm_dim = kwargs.pop("layer_norm_dim", dim)
        self.conv = generic_upsample_conv1d(
            input_channels=dim,
            output_channels=dim if dim_out is None else dim_out,
            kernel_shape=kernel_size,
            padding="SAME",
            upsample=upsample,
            dtype=dtype,
            param_dtype=param_dtype,
            rngs=rngs,
        )
        self.layer_norm = nnx.LayerNorm(
            num_features=layer_norm_dim,
            epsilon=1e-5,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.layer_norm(x)
        x = self.conv(x)
        x = jax.nn.gelu(x)
        return x


class Stem(nnx.Module):
    def __init__(
        self,
        token_embed_dim: int,
        conv_init_embed_dim: int,
        use_remat: bool = False,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        self.conv = nnx.Conv(
            in_features=token_embed_dim,
            out_features=conv_init_embed_dim,
            kernel_size=15,
            padding="SAME",
            rngs=rngs,
            param_dtype=param_dtype,
            dtype=dtype,
        )
        self.use_remat = use_remat

    def __call__(self, x: jax.Array) -> jax.Array:
        def apply_conv(conv: nnx.Module, x: jax.Array) -> jax.Array:
            return conv(x)

        if self.use_remat:
            apply_conv = nnx.remat(apply_conv)
        x = apply_conv(self.conv, x)
        x = jax.nn.gelu(x)
        return x


class ConvTowerBlock(nnx.Module):
    def __init__(
        self,
        dim_in: int,
        dim_out: int,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        self.dtype = dtype
        self.conv = ConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
            **kwargs,
        )
        self.res_conv = ResidualConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
            rngs=rngs,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            **kwargs,
        )
        self.avg_pool = partial(
            nnx.avg_pool,
            window_shape=(2,),
            strides=(2,),
            padding="SAME",
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.conv(x)
        x = self.res_conv(x)
        x = self.avg_pool(x)
        return x


class DeconvTowerBlock(nnx.Module):
    def __init__(
        self,
        dim_in: int,
        dim_out: int,
        upsample_type: DeConvUpsampleType,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
        **kwargs: Any,
    ):
        self.conv = DeConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
            upsample=upsample_type,
            rngs=rngs,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            **kwargs,
        )
        self.res_conv = ResidualDeConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
            upsample=None,
            rngs=rngs,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            **kwargs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.conv(x)
        x = self.res_conv(x)
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

    def __call__(self, x: jax.Array) -> jax.Array:
        x = jax.nn.gelu(x)
        for layer in self.hidden_layers:
            x = layer(x)
            x = jax.nn.gelu(x)
        x = self.head(x)
        return x
