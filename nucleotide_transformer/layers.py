from typing import Dict, Optional

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

        attention_logits = jnp.einsum("...thd,...Thd->...htT", query_heads, key_heads)
        sqrt_key_size = jnp.sqrt(self.key_size).astype(query.dtype)
        attention_logits = attention_logits / sqrt_key_size

        if attention_mask is not None:
            assert len(attention_mask.shape) == len(attention_logits.shape)
            attention_logits = jnp.where(attention_mask, attention_logits, -1e30)

        attention_weights = jax.nn.softmax(attention_logits)

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
            attention_bias = jnp.tile(self._bias_v, (batch_size, 1, 1, 1))
            value_heads = jnp.concatenate((value_heads, attention_bias), axis=1)

        attention = jnp.einsum("...htT,...Thd->...thd", attention_weights, value_heads)

        # Concatenate attention matrix of all heads into a single vector.
        attention_vec = jnp.reshape(attention, (*value.shape[:-1], -1))
        return hk.Linear(
            self.model_size, w_init=w_init, b_init=b_init, name="mha_output"
        )(attention_vec)

    def __call__(
        self,
        query: jnp.ndarray,
        key: jnp.ndarray,
        value: jnp.ndarray,
        attention_mask: Optional[jnp.ndarray] = None,
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

        attention_weights = self.attention_weights(query, key, attention_mask)
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
            key_size = embed_dim // num_heads

        # Define layers
        self.fc1 = hk.Linear(ffn_embed_dim, name="fc1")
        self.fc2 = hk.Linear(embed_dim, name="fc2")

        self.layer_norm_self_attention = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            name="self_attention_layer_norm",
        )
        self.layer_norm_mlp = hk.LayerNorm(
            axis=-1, create_scale=True, create_offset=True, name="final_layer_norm"
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
    ) -> TransformerOutput:
        """
        Applies the self attention mechanism.

        Args:
            x: Input token embeddings of shape (batch_size, seq_len, embed_dim).
            attention_mask: Attention mask of shape (batch_size, 1, seq_len, seq_len).

        Returns:
            Dictionary containing the output embeddings and the attention weights.
        """

        return self.sa_layer(x, x, x, attention_mask=attention_mask)

    @hk.transparent
    def mlp(self, x: Embedding) -> Embedding:
        """
        Applies one layer-norm, one linear layer, a Gelu activation,
        then a final linear layer.

        Args:
            x: Embeddings of shape (batch_size, seq_len, key_size * num_heads).

        Returns:
            The transformed sequence embedding.
        """
        x = self.layer_norm_mlp(x)
        x = jax.nn.gelu(
            self.fc1(x),
            approximate=False,
        )
        x = self.fc2(x)
        return x

    def __call__(
        self,
        x: Tokens,
        attention_mask: Optional[AttentionMask] = None,
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
        x = self.layer_norm_self_attention(x)
        output = self.self_attention(
            x=x,
            attention_mask=attention_mask,
        )
        x = output["embeddings"]
        x = res + x

        # MLP
        x = x + self.mlp(x)

        output["embeddings"] = x
        return output  # type: ignore


class SimpleLMHead(hk.Module):
    """
    Basic Language Model head. Transforms final attention block output
    into a distribution over tokens at each sequence position.
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
        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        self._final_fc = hk.Linear(
            self.alphabet_size, w_init=w_init, b_init=b_init, name="lm_final_fc"
        )

    def __call__(self, x: jnp.ndarray) -> Dict[str, jnp.ndarray]:
        # Compute logits
        logits = self._final_fc(x)
        return {"logits": logits}


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


class ESMSinusoidalPositionalEmbedding(hk.Module):
    """
    Sinusoidal embeddings to be added to token embeddings. Specific to ESM as it
    is implemented by shifting the positions by 2 (1 + padding_idx).
    """

    def __init__(
        self,
        embed_dim: int,
        padding_idx: int,
        name: Optional[str] = None,
    ):

        super().__init__(name=name)
        self.embed_dim = embed_dim
        self.padding_idx = padding_idx

    def __call__(self, tokens: jnp.ndarray) -> jnp.ndarray:
        """
        Create the sinusoidal positional embeddings
        """

        bsz, seq_len = tokens.shape
        max_pos = self.padding_idx + 1 + seq_len
        weights = self._get_embedding(max_pos)
        positions = self._make_positions(tokens)

        return weights[positions.reshape(-1), :].reshape(bsz, seq_len, -1)

    @hk.transparent
    def _make_positions(self, x: jnp.ndarray) -> jnp.ndarray:
        mask = ~jnp.equal(x, self.padding_idx)
        range_buf = (
            jnp.broadcast_to(jnp.arange(x.shape[1]), x.shape) + self.padding_idx + 1
        )
        positions = jnp.broadcast_to(range_buf, x.shape)
        return positions * mask + self.padding_idx * (1 - mask)

    @hk.transparent
    def _get_embedding(self, num_embeddings: int) -> jnp.ndarray:

        half_dim = self.embed_dim // 2
        emb = jnp.log(10000) / (half_dim - 1)
        emb = jnp.exp(jnp.arange(half_dim, dtype=jnp.float32) * -emb)
        emb = jnp.expand_dims(
            jnp.arange(num_embeddings, dtype=jnp.float32), axis=1
        ) * jnp.expand_dims(emb, 0)
        emb = jnp.reshape(
            jnp.concatenate([jnp.sin(emb), jnp.cos(emb)], axis=1), (num_embeddings, -1)
        )

        if self.embed_dim % 2 == 1:
            # zero pad
            emb = jnp.concatenate([emb, jnp.zeros((num_embeddings, 1))], axis=1)

        if self.padding_idx is not None:
            emb = emb.at[self.padding_idx, :].set(0)

        return emb
