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

from typing import Dict, Optional, Tuple

import haiku as hk
import jax
import jax.numpy as jnp
from haiku import initializers

from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)
from nucleotide_transformer.utils import get_activation_fn

# Constant used in Sinusoidal/Rotary Embeddings, reference to this value can be found
# on page 6 of https://arxiv.org/pdf/1706.03762.pdf and page 5 of
# https://arxiv.org/abs/2104.09864
# These rotary positional embeddings are proper to ESM implementation.
# dimensions in key space are rotated 2 by 2. The key difference with
# GPT's one is that in this case each 2 dimensions rotated together are spaced
# by key_size//2
UPPER_FREQ = 10000


class RotaryEmbedding(hk.Module):
    """
    Rotary Positional Embedding inspired by RoFormer:
    https://arxiv.org/abs/2104.09864
    https://github.com/ZhuiyiTechnology/roformer .
    """

    def __init__(
        self,
        key_size: int,
        name: Optional[str] = None,
    ):

        """
        Args:
            key_size: Dimension of one head.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self._inv_freq = 1.0 / (UPPER_FREQ ** (jnp.arange(0, key_size, 2) / key_size))

    def _compute_cos_sin_tables(
        self,
        heads: jnp.ndarray,
    ) -> Tuple[jnp.ndarray, jnp.ndarray]:
        """
        Computes the cosinus and sinus for rotation.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
            key_size).

        Returns:
            Cosinus positional embedding of shape (1, seq_len, 1,
                key_size).
            Sinus positional embedding of shape (1, seq_len, 1,
                key_size).
        """
        seq_len = heads.shape[1]

        self._seq_len_cached = seq_len
        t = jnp.arange(seq_len)
        freqs = jnp.einsum("i,j->ij", t, self._inv_freq)
        emb = jnp.concatenate((freqs, freqs), axis=-1, dtype=heads.dtype)

        # Compute cos and cast is as (1, seq_len, 1, key_size) to be applied to queries
        # of shape (batch_size, seq_len, num_heads, key_size)
        cos_cached = jnp.cos(emb)[None, :, None, :]
        sin_cached = jnp.sin(emb)[None, :, None, :]

        return cos_cached, sin_cached

    def _apply_rotary_pos_emb(
        self, heads: jnp.ndarray, cos: jnp.ndarray, sin: jnp.ndarray
    ) -> jnp.ndarray:
        """
        Applies the rotary positional embedding to the heads.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
                key_size).
            cos: Cosinus values.
            sin: Sinus values.

        Returns:
            Embedded heads of shape (batch_size, seq_len, num_heads,
                key_size).
        """

        # Rotate x
        x1, x2 = heads[..., : heads.shape[-1] // 2], heads[..., heads.shape[-1] // 2 :]
        heads_rotated = jnp.concatenate((-x2, x1), axis=-1)

        embedded_heads = (heads * cos) + (heads_rotated * sin)
        return embedded_heads

    def __call__(
        self, query_heads: jnp.ndarray, key_heads: jnp.ndarray
    ) -> Tuple[jnp.ndarray, jnp.ndarray]:
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


class MultiHeadAttention(hk.MultiHeadAttention):
    """
    Multi-head attention with masking applied. Modified from the core implementation to
    support biases in keys and values.
    """

    def __init__(
        self,
        num_heads: int,
        key_size: int,
        use_rotary_embedding: bool = False,
        add_bias_kv: bool = False,
        value_size: Optional[int] = None,
        model_size: Optional[int] = None,
        name: Optional[str] = None,
    ):
        """
        Args:
            num_heads: Number of independent attention heads.
            key_size: The size of keys and queries used for attention.
            use_rotary_embedding: If true, adds rotary embeddings to the key and query
                heads (see RoFormer https://arxiv.org/pdf/2104.09864.pdf).
            add_bias_kv: If True, appends biases to key and query heads, used in ESM
                model (https://www.biorxiv.org/content/10.1101/622803v4.full.pdf).
            value_size: Optional size of the value projection. If None, defaults
                to the key size.
            model_size: Optional size of the output embedding. If None, defaults
                to the key size multiplied by the number of heads.
            name: Optional name for this module.
        """
        w_init = hk.initializers.VarianceScaling(2.0, "fan_in", "uniform")
        super().__init__(
            num_heads=num_heads,
            key_size=key_size,
            w_init=w_init,
            value_size=value_size,
            model_size=model_size,
            name=name,
        )

        if add_bias_kv:
            self._bias_k = hk.get_parameter(
                "bias_k", [1, 1, self.num_heads, self.key_size], init=jnp.zeros
            )
            self._bias_v = hk.get_parameter(
                "bias_v", [1, 1, self.num_heads, self.value_size], init=jnp.zeros
            )
        else:
            self._bias_k = None
            self._bias_v = None
        self._use_rotary_embedding = use_rotary_embedding

    @hk.transparent
    def attention_weights(
        self,
        query: jnp.ndarray,
        key: jnp.ndarray,
        attention_mask: Optional[AttentionMask] = None,
        attention_weight_bias: Optional[jnp.ndarray] = None,
    ) -> jnp.ndarray:
        """
        Computes the attention weights.

        Args:
            query: Embedding sequence to compute queries.
            key: Embedding sequence to compute keys.
            attention_mask: Input attention_mask. Defaults to None.

        Returns:
            Attention weights.
        """

        query_heads = self._linear_projection_he_init(query, self.key_size, "query")
        key_heads = self._linear_projection_he_init(key, self.key_size, "key")

        # Add bias for key (see ESM architecture)
        jmp_policy = hk.mixed_precision.current_policy()
        if jmp_policy is None:
            # default float32
            compute_dtype = jnp.float32
        else:
            # cast to jmp policy if specified
            compute_dtype = jmp_policy.compute_dtype

        if self._bias_k is not None:
            batch_size = key_heads.shape[0]
            attention_bias = jnp.tile(self._bias_k, (batch_size, 1, 1, 1)).astype(
                dtype=compute_dtype
            )
            key_heads = jnp.concatenate((key_heads, attention_bias), axis=1)
            if attention_mask is not None:
                attention_mask = jnp.concatenate(
                    (
                        attention_mask,
                        jnp.ones(attention_mask.shape[:-1] + (1,), dtype=jnp.bool_),
                    ),
                    axis=-1,
                )

        if self._use_rotary_embedding:
            query_heads, key_heads = RotaryEmbedding(
                self.key_size, name="rotary_embed"
            )(query_heads, key_heads)

        attention_logits = jnp.einsum("...thd,...Thd->...htT", query_heads, key_heads)
        sqrt_key_size = jnp.sqrt(self.key_size).astype(query.dtype)
        attention_logits = attention_logits / sqrt_key_size

        if attention_mask is not None:
            assert len(attention_mask.shape) == len(attention_logits.shape)
            attention_logits = jnp.where(attention_mask, attention_logits, -1e30)

        if attention_weight_bias is None:
            attention_weights = jax.nn.softmax(attention_logits)
        else:
            attention_weights = jax.nn.softmax(attention_logits + attention_weight_bias)

        return attention_weights

    @hk.transparent
    def compute_embeddings(
        self,
        value: jnp.ndarray,
        attention_weights: jnp.ndarray,
    ) -> jnp.ndarray:
        """
        Computes the output embeddings.

        Args:
            value: Embedding sequence to compute values.
            attention_weights: Attention weights.

        Returns:
            Output embeddings.
        """

        # He initialization
        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")

        value_heads = self._linear_projection_he_init(value, self.value_size, "value")

        if self._bias_v is not None:
            batch_size = value_heads.shape[0]
            # Add bias for key (see ESM architecture)
            jmp_policy = hk.mixed_precision.current_policy()
            if jmp_policy is None:
                # default float32
                compute_dtype = jnp.float32
            else:
                # cast to jmp policy if specified
                compute_dtype = jmp_policy.compute_dtype

            attention_bias = jnp.tile(
                self._bias_v,
                (batch_size, 1, 1, 1),
            ).astype(dtype=compute_dtype)
            value_heads = jnp.concatenate((value_heads, attention_bias), axis=1)

        attention = jnp.einsum("...htT,...Thd->...thd", attention_weights, value_heads)

        # Concatenate attention matrix of all heads into a single vector.
        attention_vec = jnp.reshape(attention, (*attention.shape[:-2], -1))
        return hk.Linear(
            self.model_size, w_init=w_init, b_init=b_init, name="mha_output"
        )(attention_vec)

    def __call__(
        self,
        query: jnp.ndarray,
        key: jnp.ndarray,
        value: jnp.ndarray,
        attention_mask: Optional[jnp.ndarray] = None,
        attention_weight_bias: Optional[jnp.ndarray] = None,
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

    @hk.transparent
    def _linear_projection_he_init(
        self, x: jnp.ndarray, head_size: int, name: Optional[str] = None
    ) -> jnp.ndarray:
        """
        Linear layer for multi-head attention mechanism. Initialized with the He method.

        Args:
            x: Input embeddings.
            head_size: Embedding size of each attention head.
            name: Name of the linear layer.

        Returns:
            Multi-head embeddings.
        """

        # He initialization
        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")

        y = hk.Linear(
            self.num_heads * head_size, w_init=w_init, b_init=b_init, name=name
        )(x)
        return y.reshape((*x.shape[:-1], self.num_heads, head_size))


class SelfAttentionBlock(hk.Module):
    """
    Attention block made of self-attention.
    """

    def __init__(
        self,
        num_heads: int,
        embed_dim: int,
        ffn_embed_dim: int,
        key_size: Optional[int] = None,
        use_rotary_embedding: bool = False,
        add_bias_kv: bool = False,
        add_bias_fnn: bool = True,
        ffn_activation_name: str = "gelu-no-approx",
        use_glu_in_ffn: bool = False,
        layer_norm_eps: float = 1e-5,  # this is the default haiku value
        pre_layer_norm: bool = True,
        name: Optional[str] = None,
    ):
        super().__init__(name=name)
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
        if use_glu_in_ffn:
            # user should multiply ffn_embed_dim by 2/3 when using GLU
            # to keep total number of parameters equal
            # see https://arxiv.org/pdf/2002.05202.pdf. for more details
            # we multiply by 2 here as the output will be split in 2 for GLU
            ffn_embed_dim = int(2 * ffn_embed_dim)

        self.fc1 = hk.Linear(ffn_embed_dim, name="fc1", with_bias=add_bias_fnn)
        self.fc2 = hk.Linear(embed_dim, name="fc2", with_bias=add_bias_fnn)

        self.layer_norm_self_attention = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            name="self_attention_layer_norm",
            eps=layer_norm_eps,
        )
        self.layer_norm_mlp = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            name="final_layer_norm",
            eps=layer_norm_eps,
        )
        self.sa_layer = MultiHeadAttention(
            num_heads=num_heads,
            key_size=key_size,
            model_size=embed_dim,
            add_bias_kv=add_bias_kv,
            use_rotary_embedding=use_rotary_embedding,
            name="self_attention",
        )

    @hk.transparent
    def self_attention(
        self,
        x: Embedding,
        attention_mask: Optional[AttentionMask] = None,
        attention_weight_bias: Optional[jnp.ndarray] = None,
    ) -> TransformerOutput:
        """
        Applies the self attention mechanism.

        Args:
            x: Input token embeddings of shape (batch_size, seq_len, embed_dim).
            attention_mask: Attention mask of shape (batch_size, 1, seq_len, seq_len).

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

    @hk.transparent
    def mlp(self, embed: Embedding) -> Embedding:
        """
        Applies one layer-norm, one linear layer, a Gelu activation,
        then a final linear layer.

        Args:
            x: Embeddings of shape (batch_size, seq_len, key_size * num_heads).

        Returns:
            The transformed sequence embedding.
        """

        if self._pre_layer_norm:
            x = self.layer_norm_mlp(embed)
        else:
            x = embed

        if self._use_glu_in_fnn:
            x1, x2 = jnp.split(self.fc1(x), indices_or_sections=2, axis=-1)
            x = self._ffn_activation_fn(x1) * x2
        else:
            x = self._ffn_activation_fn(self.fc1(x))

        x = self.fc2(x)

        if not self._pre_layer_norm:
            x = self.layer_norm_mlp(x + embed)

        return x

    def __call__(
        self,
        x: Tokens,
        attention_mask: Optional[AttentionMask] = None,
        attention_weight_bias: Optional[jnp.ndarray] = None,
    ) -> TransformerOutput:
        """
        Computes the output of the attention layer.

        Args:
            x: Input token embeddings of shape (batch_size,seq_len,embed_dim).
            attention_mask: Attention mask of shape (batch_size, 1,seq_len, seq_len).

        Returns:
            A dictionary containing the output embeddings and the attention weights.
        """

        # Self-Attention
        res = x
        if self._pre_layer_norm:
            x = self.layer_norm_self_attention(x)

        output = self.self_attention(
            x=x,
            attention_mask=attention_mask,
            attention_weight_bias=attention_weight_bias,
        )

        if not self._pre_layer_norm:
            output["embeddings"] = self.layer_norm_self_attention(
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


class RobertaLMHead(hk.Module):
    """
    Roberta Language Model head. Transform final attention layer output into a
    distribution over tokens at each position.
    """

    def __init__(self, embed_dim: int, alphabet_size: int, name: Optional[str] = None):
        """
        Args:
            embed_dim: Embedding dimension.
            alphabet_size: Number of tokens in the alphabet.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self.embed_dim = embed_dim
        self.alphabet_size = alphabet_size

        # Define layers
        self._first_layer_norm = hk.LayerNorm(
            axis=-1, create_scale=True, create_offset=True, name="emb_layer_norm_after"
        )
        self._fc1 = hk.Linear(self.embed_dim, name="lm_head_fc_1")
        self._final_fc = hk.Linear(self.alphabet_size, name="lm_final_fc")
        self._second_layer_norm = hk.LayerNorm(
            axis=-1, create_scale=True, create_offset=True, name="lm_head_layer_norm"
        )

    def __call__(self, x: jnp.ndarray) -> Dict[str, jnp.ndarray]:
        x = self._first_layer_norm(x)
        # Embeddings are computed after the first layer norm to be consistent with ESM
        embeddings = x
        x = self._fc1(x)
        x = jax.nn.gelu(x, approximate=False)
        x = self._second_layer_norm(x)

        # Compute logits
        logits = self._final_fc(x)
        return {"embeddings": embeddings, "logits": logits}


class TokensDropout(hk.Module):
    """
    Tokens dropout layer.
    """

    def __init__(
        self,
        embed_dim: int,
        pad_token_id: int,
        mask_token_id: int,
        masking_ratio: float,
        masking_prob: float,
        name: Optional[str] = None,
    ):
        """
        Args:
            embed_dim: Embedding dimension.
            pad_token_id: ID of the pad token.
            mask_token_id: ID of the pad token.
            masking_ratio: Masking ratio.
            masking_prob: Probability to mask.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self.pad_token_id = pad_token_id
        self.mask_token_id = mask_token_id
        self.masking_ratio = masking_ratio
        self.masking_prob = masking_prob
        self.embed_dim = embed_dim

    def __call__(self, x: jnp.ndarray, tokens: Tokens) -> jnp.ndarray:

        padding_mask_tokens = tokens == self.pad_token_id
        tokens_repeated = jnp.repeat(
            tokens[:, :, None], repeats=self.embed_dim, axis=-1
        )
        x = jnp.where(tokens_repeated == self.mask_token_id, 0.0, x)
        mask_ratio_train = self.masking_ratio * self.masking_prob
        src_lengths = (~padding_mask_tokens).sum(-1)
        mask_ratio_observed = (tokens == self.mask_token_id).sum(-1) / src_lengths
        x = x * (1 - mask_ratio_train) / (1 - mask_ratio_observed)[:, None, None]
        return x


class ESMLearnedPositionalEmbeddings(hk.Module):
    """
    Learned positional embeddings to be added to token embeddings. Specific to ESM as it
    is implemented by shifting the positions by 2 (1 + padding_idx).
    """

    def __init__(
        self,
        vocab_size: int,
        embed_dim: int,
        padding_idx: int,
        name: Optional[str] = None,
    ):
        """
        Args:
            vocab_size: Tokenizer's vocabulary size.
            embed_dim: Embedding size.
            padding_idx: Index attributed to the padding
                token. Defaults to 1.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self.padding_idx = padding_idx
        self._embed_layer = hk.Embed(vocab_size + padding_idx + 1, embed_dim)

    def __call__(self, tokens: jnp.ndarray) -> jnp.ndarray:
        mask = tokens != self.padding_idx
        positions = jnp.cumsum(mask, axis=1) * mask + self.padding_idx

        return self._embed_layer(positions)
