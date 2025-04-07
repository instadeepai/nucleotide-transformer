import logging
import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import haiku as hk
import jax
import jax.numpy as jnp
import jmp
from haiku import initializers

from nucleotide_transformer.chatNT.configs import ESMTransformerConfig
from nucleotide_transformer.chatNT.esm_rotary import (
    RotaryEmbedding,
    RotaryEmbeddingConfig,
)
from nucleotide_transformer.chatNT.types import (
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
        rotary_embedding_config: Optional[RotaryEmbeddingConfig] = None,
        add_bias_kv: bool = False,
        value_size: Optional[int] = None,
        model_size: Optional[int] = None,
        name: Optional[str] = None,
    ):
        """
        Args:
            num_heads: Number of independent attention heads.
            key_size: The size of keys and queries used for attention.
            rotary_embedding_config: Configuration to specify hyperparameters for
                RotaryEmbeddig layer
                (see RoFormer https://arxiv.org/pdf/2104.09864.pdf). If None,
                rotary embeddings are not used. If specified, it contains the
                hyperparameters specifying the type of rotary embedding applied.
            add_bias_kv: If True, appends biases to key and query heads, used in ESM
                model (https://www.biorxiv.org/content/10.1101/622803v4.full.pdf).
            value_size: Optional size of the value projection. If None, defaults
                to the key size.
            model_size: Optional size of the output embedding. If None, defaults
                to the key size multiplied by the number of heads.
            name: Optional name for this module.
            rescaling_factor: Scaling factor to use for rotary positional embeddings
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
        self._rotary_embedding_config = rotary_embedding_config

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

        if self._rotary_embedding_config:

            query_heads, key_heads = RotaryEmbedding(
                self.key_size, self._rotary_embedding_config, name="rotary_embed"
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
        add_bias_kv: bool = False,
        add_bias_fnn: bool = True,
        ffn_activation_name: str = "gelu-no-approx",
        use_glu_in_ffn: bool = False,
        layer_norm_eps: float = 1e-5,  # this is the default haiku value
        pre_layer_norm: bool = True,
        name: Optional[str] = None,
        rotary_embedding_config: Optional[RotaryEmbeddingConfig] = None,
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
            add_bias_kv=add_bias_kv,
            model_size=embed_dim,
            name="self_attention",
            rotary_embedding_config=rotary_embedding_config,
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


class SimpleLMHead(hk.Module):
    """
    Basic Language Model head. Transforms final attention block output
    into a distribution over tokens at each sequence position.
    """

    def __init__(
        self,
        embed_dim: int,
        alphabet_size: int,
        add_bias_lm_head: bool = True,
        name: Optional[str] = None,
    ):
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
            self.alphabet_size,
            w_init=w_init,
            b_init=b_init,
            with_bias=add_bias_lm_head,
            name="lm_final_fc",
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


class LearnedPositionalEmbeddings(hk.Module):
    """Position embeddings to be added to token embeddings."""

    def __init__(
        self,
        max_positions: int,
        embed_dim: int,
        name: Optional[str] = None,
    ):
        """
        Args:
            max_positions: Max number of positions in sequence.
            embed_dim: Embedding size.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self._max_positions = max_positions
        self._embed_layer = hk.Embed(max_positions, embed_dim)

    def __call__(
        self, tokens: jnp.ndarray, positions: Optional[jnp.ndarray] = None
    ) -> jnp.ndarray:
        """
        Creates the positional embeddings.

        Args:
            tokens: Input tokens, shape (batch_size, seq_len).

        Returns:
            The positional embeddings.
        """

        batch_size, seq_len = tokens.shape[0], tokens.shape[1]

        if positions is None:
            positions = jnp.tile(jnp.arange(seq_len), (batch_size, 1))

        return self._embed_layer(positions)


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
        # We create this array by summing an array of booleans. This is automatically
        # cast off to jnp.float32. We need to recast it to the desired dtype, otherwise
        # multiplying it with the activations will make the activations jnp.float32
        # again.
        mask_ratio_observed = jnp.asarray(
            (tokens == self.mask_token_id).sum(-1) / src_lengths, dtype=x.dtype
        )
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


class ESMTransformer(hk.Module):
    """
    Jax implementation of ESM models. Covers ESM1, ESM1-b and ESM2 models.
    """

    def __init__(
        self,
        config: ESMTransformerConfig,
        name: Optional[str] = None,
    ):
        """
        Initialize an ESM model.

        Args:
            config: Dataclass containing model hyperparameters.
            name: Name for module (custom will break weight loading).
        """

        self._config = config
        super().__init__(name=name)

        self._embed_layer = hk.Embed(self._config.alphabet_size, self._config.embed_dim)

        if config.bias_word_embedding:
            self.embed_bias = hk.Bias(bias_dims=[-1], name="bias_word_embedding")

        if config.positional_embedding == "learned":
            self._pos_embed_layer = ESMLearnedPositionalEmbeddings(
                config.max_positions, config.embed_dim, config.pad_token_id
            )

        elif config.positional_embedding == "learned_standard":
            self._pos_embed_layer = LearnedPositionalEmbeddings(
                config.max_positions, config.embed_dim
            )

        elif config.positional_embedding == "sinusoidal":
            self._pos_embed_layer = ESMSinusoidalPositionalEmbedding(
                config.embed_dim, config.pad_token_id
            )
        elif self._config.positional_embedding == "alibi_dnabert_2":
            self.build_alibi_tensor_fn = build_alibi_tensor_dnabert2

        if config.lm_head is None:
            self._lm_head = config.lm_head

        elif config.lm_head == "roberta":
            self._lm_head = RobertaLMHead(
                embed_dim=self._config.embed_dim,
                alphabet_size=self._config.alphabet_size,
                name="roberta_lm_head",
            )

        elif config.lm_head == "simple":
            self._lm_head = SimpleLMHead(
                embed_dim=self._config.embed_dim,
                alphabet_size=self._config.alphabet_size,
            )

        if self._config.emb_layer_norm_before:
            self.emb_ln_before = hk.LayerNorm(
                axis=-1,
                create_scale=True,
                create_offset=True,
                name="emb_layer_norm_before",
                eps=self._config.layer_norm_eps,
            )

        if config.use_rotary_embedding:
            self._rotary_embedding_config = RotaryEmbeddingConfig(
                rescaling_factor=config.rescaling_factor
            )
        else:
            self._rotary_embedding_config = None  # type: ignore

        # Process attention maps to save requirement into more suitable format
        attention_maps_to_save = config.attention_maps_to_save
        self._attention_layers_to_save = list({t[0] for t in attention_maps_to_save})
        self._attention_maps_per_layer_to_save = {
            layer: [t[1] for t in attention_maps_to_save if t[0] == layer]
            for layer in self._attention_layers_to_save
        }

        # Checking user request can be executed, raise error otherwise
        max_layer = max(self._attention_layers_to_save + [0])
        if max_layer > config.num_layers:
            raise ValueError(
                f"You are requiring attention maps for layer {max_layer}, "
                f"while the model has {config.num_layers} layers only."
            )

        for layer, maps in self._attention_maps_per_layer_to_save.items():
            max_map = max(maps)
            if max_map > config.attention_heads:
                raise ValueError(
                    f"You are requiring attention maps number {max_map} "
                    f"at layer {layer}, while the model has {config.attention_heads} "
                    f"only."
                )

    @hk.transparent
    def apply_attention_blocks(
        self,
        x: Embedding,
        outs: Dict[str, Embedding],
        attention_mask: Optional[AttentionMask] = None,
    ) -> Tuple[Embedding, Dict[str, Embedding]]:
        """
        Create the blocks of attention layers and applies them.

        Args:
            x: The sequence embedding.
            outs: A dictionary to carry through the attention layers which stores the
                intermediate sequence embedding and attention maps.
            attention_mask: Attention mask of shape (batch_size, 1, seq_len, seq_len).

        Returns:
            The output sequence embedding.
            The optional intermediate results (embeddings of the layer and attention
                weights).
        """

        layers: List[Callable] = [
            self._attention_block(layer_idx)
            for layer_idx in range(self._config.num_layers)
        ]

        if self._config.use_gradient_checkpointing:
            # the remat-ed function cannot take control flow arguments
            layers = [hk.remat(layer) for layer in layers]

        if self._config.positional_embedding and (
            "alibi" in self._config.positional_embedding
        ):
            alibi_bias = self.build_alibi_tensor_fn(
                x.shape[1], self._config.attention_heads
            )
        else:
            alibi_bias = None

        for layer_idx, layer in enumerate(layers):
            output = layer(
                x=x, attention_mask=attention_mask, attention_weight_bias=alibi_bias
            )
            x = output["embeddings"]

            # Save intermediate embeddings if needed
            if (layer_idx + 1) in self._config.embeddings_layers_to_save:
                outs[f"embeddings_{(layer_idx + 1)}"] = output["embeddings"]
            # Save intermediate attention maps if needed
            if (layer_idx + 1) in self._attention_layers_to_save:
                for map_number in self._attention_maps_per_layer_to_save[layer_idx + 1]:
                    dkey = f"attention_map_layer_{layer_idx + 1}_number_{map_number}"
                    outs[dkey] = output["attention_weights"][:, map_number + 1]

        return x, outs

    @hk.transparent
    def _attention_block(self, layer_idx: int) -> SelfAttentionBlock:
        return SelfAttentionBlock(  # type: ignore
            num_heads=self._config.attention_heads,
            embed_dim=self._config.embed_dim,
            key_size=self._config.key_size,
            ffn_embed_dim=self._config.ffn_embed_dim,
            add_bias_kv=self._config.add_bias_kv,
            add_bias_fnn=self._config.add_bias_ffn,
            ffn_activation_name=self._config.ffn_activation_name,
            use_glu_in_ffn=self._config.use_glu_in_ffn,
            rotary_embedding_config=self._rotary_embedding_config,
            layer_norm_eps=self._config.layer_norm_eps,
            pre_layer_norm=self._config.pre_layer_norm,
            name=f"attention_layer_{layer_idx}",
        )

    def __call__(
        self,
        tokens: Tokens,
        attention_mask: Optional[AttentionMask] = None,
    ) -> TransformerOutput:
        """
        Computes the embeddings based on the input tokens.

        Args:
            tokens: Input tokens out of the tokenizer of shape (batch_size, seq_len).
            attention_mask: Attention mask of shape (batch_size, 1, seq_len, seq_len).
                If no mask is provided, a mask by default which equals 1 over all non
                pad tokens and 0 over pad tokens is computed.

        Returns:
            Dictionary containing the final embeddings and logits.
        """

        # Prepare outputs dict
        outs: Dict[str, jnp.ndarray] = {}

        x = self._embed_layer(tokens)

        if self._config.bias_word_embedding:
            x = self.embed_bias(x)

        # Tokens dropout if needed
        if self._config.token_dropout:
            x = TokensDropout(
                embed_dim=self._config.embed_dim,
                mask_token_id=self._config.mask_token_id,
                pad_token_id=self._config.pad_token_id,
                masking_ratio=self._config.masking_ratio,
                masking_prob=self._config.masking_prob,
            )(x, tokens)

        # RoBERTa's mask scaling factor
        x = self._config.embed_scale * x
        if self._config.positional_embedding in [
            "learned",
            "sinusoidal",
            "learned_standard",
        ]:
            if self._config.positional_embedding == "learned":
                # Add check that the sequence fed into the transformer is not longer
                # than the max positions used to instantiate the learned positional
                # embeddings layer
                max_length_authorized = (
                    self._pos_embed_layer._embed_layer.vocab_size
                    - self._pos_embed_layer.padding_idx
                    - 1
                )
                assert tokens.shape[1] <= max_length_authorized, (
                    "Inputs to the learned positional embeddings layer have a length "
                    f"{x.shape[1]} greater than the max positions used to instantiate "
                    f"it: {max_length_authorized}"
                )
            elif self._config.positional_embedding == "learned_standard":
                # Add check that the sequence fed into the transformer is not longer
                # than the max positions used to instantiate the learned positional
                # embeddings layer
                max_length_authorized = self._pos_embed_layer._embed_layer.vocab_size
                assert tokens.shape[1] <= max_length_authorized, (
                    "Inputs to the learned positional embeddings layer have a length "
                    f"{x.shape[1]} greater than the max positions used to instantiate "
                    f"it: {max_length_authorized}"
                )
            x = x + self._pos_embed_layer(tokens)

        if self._config.emb_layer_norm_before:
            x = self.emb_ln_before(x)

        # Attention mask
        if attention_mask is None:
            attention_mask = build_padding_attention_mask(
                tokens=tokens, pad_token_id=self._config.pad_token_id
            )

        # Mask before attention for models ESM1b and ESM2
        if self._config.mask_before_attention:
            x = x - jnp.where(attention_mask == 0, 1, 0)[:, 0, 0][..., None] * x
        # construct a tower of attention layers

        x, outs = self.apply_attention_blocks(
            x=x,
            outs=outs,
            attention_mask=attention_mask,
        )

        # Language Model Head
        if self._lm_head:
            lm_head_outs = self._lm_head(x)
            sequence_mask = attention_mask[:, 0, :, 0][:, :, None]
            outs["logits"] = jnp.where(sequence_mask, lm_head_outs["logits"], 0)

        if self._config.lm_head == "roberta":
            embeddings = lm_head_outs["embeddings"]
        else:
            embeddings = x

        # Save final embeddings if needed
        if self._config.num_layers in self._config.embeddings_layers_to_save:
            outs[f"embeddings_{self._config.num_layers}"] = embeddings

        return outs  # type: ignore


def build_padding_attention_mask(tokens: Tokens, pad_token_id: int) -> AttentionMask:
    """
    Builds a padding mask from a sequence of tokens by masking <pad> in the attention.

    Args:
        tokens: Batch of sequences of shape (batch_size, seq_len).
        pad_token_id: Int corresponding to the <pad> token to mask.

    Returns:
        Batch of attention masks, masking out <pad> tokens.
    """
    padding_mask = tokens != pad_token_id
    padding_mask = padding_mask[:, None, :]
    padding_mask = jnp.einsum("bhT, bht->bhtT", padding_mask, padding_mask)
    return padding_mask


def build_alibi_tensor_dnabert2(sequence_len: int, num_heads: int) -> jnp.ndarray:
    """
    Compute the attention bias tensor for ALIBI-based attention
    ref: https://arxiv.org/pdf/2108.12409.pdf
    BUT, it does in the dnabert2 fashion, which not exactly the way explained in the
    Alibi paper

    Args:
        sequence_len: the sequence_len of the sequence
        num_heads: the number of heads

    Returns:
        The alibi attention bias tensor
    """

    context_position = jnp.arange(
        sequence_len,
    )[:, jnp.newaxis]
    memory_position = jnp.arange(
        sequence_len,
    )[jnp.newaxis, :]
    relative_position = jnp.abs(memory_position - context_position)
    relative_position = jnp.tile(
        jnp.expand_dims(relative_position, 0), (num_heads, 1, 1)
    )
    abstract_num_heads = int(2 ** math.floor(math.log2(num_heads)))
    all_slopes = jnp.power(
        2,
        -(jnp.arange(abstract_num_heads * 2) + 1) * 8 / (abstract_num_heads * 2),
    )
    slopes_a = all_slopes[1::2]
    slopes_b = all_slopes[0::2][: num_heads - abstract_num_heads]

    slopes = jnp.concatenate((slopes_a, slopes_b), axis=0)

    alibi_attention_bias = slopes[:, jnp.newaxis, jnp.newaxis] * -relative_position
    alibi_attention_bias = jnp.expand_dims(alibi_attention_bias, 0)

    return alibi_attention_bias


SUPPORTED_FFN_ACTIVATIONS = ["gelu", "gelu-no-approx", "relu", "swish", "silu", "sin"]


def get_activation_fn(activation_name: str) -> Callable:
    """
    Return activation fn given its name.
    Args:
        activation_name: Activation name.

    Returns:
        activation function.
    """
    if activation_name not in SUPPORTED_FFN_ACTIVATIONS:
        raise NotImplementedError(
            f"Activation {activation_name} not supported yet. "
            f"Supported activations for feed forward "
            f"block are {SUPPORTED_FFN_ACTIVATIONS}"
        )
    if activation_name == "gelu-no-approx":
        activation_fn = lambda x: jax.nn.gelu(x, approximate=False)  # noqa: E731
    elif activation_name == "sin":
        activation_fn = lambda x: jnp.sin(x)  # noqa: E731
    else:
        activation_fn = getattr(jax.nn, activation_name)
    return activation_fn
