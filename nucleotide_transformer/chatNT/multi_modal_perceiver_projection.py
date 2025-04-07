from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
from haiku import initializers
from jax import numpy as jnp

from nucleotide_transformer.chatNT.configs import PerceiverResamplerConfig
from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)

UPPER_FREQ = 10000


@dataclass
class RotaryEmbeddingConfig:
    """
    Parameters to initialize the RotaryEmbedding layer. The rescaling factor allows
    to adapt the rotary embeddings to larger lengths than what was used for training.
    One of this strategy is presented in the Yarn paper: https://arxiv.org/pdf/2309.00071.pdf. # noqa

    Args:

    """

    rescaling_factor: Optional[float]


class RotaryEmbedding(hk.Module):
    """
    Rotary Positional Embedding inspired by RoFormer:
    https://arxiv.org/abs/2104.09864
    https://github.com/ZhuiyiTechnology/roformer .
    """

    def __init__(
        self,
        key_size: int,
        rotary_embedding_config: RotaryEmbeddingConfig,
        name: Optional[str] = None,
    ):
        """
        Args:
            key_size: Dimension of one head.
            rotary_embedding_config: Configuration to specify hyperparameters for
                RotaryEmbeddig layer
                (see RoFormer https://arxiv.org/pdf/2104.09864.pdf). It contains
                the hyperparameters specifying the type of rotary embedding applied.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)

        # Extract argument from the config
        rescaling_factor = rotary_embedding_config.rescaling_factor

        if rescaling_factor is None:
            self._inv_freq = 1.0 / (
                UPPER_FREQ ** (np.arange(0, key_size, 2) / key_size)
            )
        else:
            updated_base = UPPER_FREQ * (
                rescaling_factor ** (key_size / (key_size - 2))
            )
            self._inv_freq = 1.0 / (
                updated_base ** (np.arange(0, key_size, 2) / key_size)
            )

    def _compute_cos_sin_tables(
        self,
        heads: jnp.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the cosinus and sinus for rotation.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
            key_size).

        Returns:
            Cosinus positional embedding of shape (1, seq_len, 1,
                key_size/2).
            Sinus positional embedding of shape (1, seq_len, 1,
                key_size/2).
        """
        seq_len = heads.shape[1]

        self._seq_len_cached = seq_len
        t = np.arange(seq_len)
        freqs = np.einsum("i,j->ij", t, self._inv_freq)

        # Compute cos and cast is as (1, seq_len, 1, key_size/2) to be applied to
        # queries of shape (batch_size, seq_len, num_heads, key_size/2)
        cos_cached = np.cos(freqs)[None, :, None, :]
        sin_cached = np.sin(freqs)[None, :, None, :]

        return cos_cached, sin_cached

    def _apply_rotary_pos_emb(
        self, heads: jnp.ndarray, cos: np.ndarray, sin: np.ndarray
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
        x_first, x_second = (
            heads[..., : heads.shape[-1] // 2],
            heads[..., heads.shape[-1] // 2 :],
        )
        first_part = x_first * cos - x_second * sin
        second_part = x_second * cos + x_first * sin

        return jnp.concatenate((first_part, second_part), axis=-1, dtype=heads.dtype)

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


class MultiModalPerceiverResamplerBlock(hk.Module):
    """
    Adaption of the block of the perceiver resampler architecture to co-encode
    two modalities at the same time (e.g. DNA and English).
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
        name: Optional[str] = None,
    ):
        """
        Documentation generated by Copilot.
        Initializes the MultiModalPerceiverLayer.

        Args:
            num_heads (int): Number of attention heads.
            embed_dim (int): Dimension of the embedding.
            ffn_embed_dim (int): Dimension of the feed-forward network embedding.
            key_size (Optional[int], optional): Size of the key. Defaults to None.
            add_bias_kv (bool, optional): Whether to add bias to key and value. Defaults
                to False.
            add_bias_fnn (bool, optional): Whether to add bias to feed-forward network.
                Defaults to True.
            ffn_activation_name (str, optional): Name of the activation function for
                feed-forward network. Defaults to "gelu-no-approx".
            use_glu_in_ffn (bool, optional): Whether to use Gated Linear Units in
                feed-forward network. Defaults to False.
            name (Optional[str], optional): Name of the layer. Defaults to None.

        Raises:
            ValueError: If the embedding dimension is not divisible by the number of
                heads when key_size is None.
        """
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

        # Hyperparameters internalization
        self._ffn_activation_fn = get_activation_fn(activation_name=ffn_activation_name)
        self._use_glu_in_ffn = use_glu_in_ffn

        # Define layers
        if use_glu_in_ffn:
            # user should multiply ffn_embed_dim by 2/3 when using GLU
            # to keep total number of parameters equal
            # see https://arxiv.org/pdf/2002.05202.pdf. for more details
            # we multiply by 2 here as the output will be split in 2 for GLU
            ffn_embed_dim = int(2 * ffn_embed_dim)

        # Define layers
        self.fc1 = hk.Linear(ffn_embed_dim, name="fc1", with_bias=add_bias_fnn)
        self.fc2 = hk.Linear(embed_dim, name="fc2", with_bias=add_bias_fnn)

        self.layer_norm_cross_attention_1 = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            name="cross_attention_layer_norm_1",
        )
        self.layer_norm_cross_attention_2 = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            name="cross_attention_layer_norm_2",
        )

        self.layer_norm_mlp = hk.LayerNorm(
            axis=-1, create_scale=True, create_offset=True, name="final_layer_norm"
        )

        self.cross_attention_layer_1 = MultiHeadAttention(
            num_heads=num_heads,
            key_size=key_size,
            add_bias_kv=add_bias_kv,
            name="cross_attention_1",
        )

        self.cross_attention_layer_2 = MultiHeadAttention(
            num_heads=num_heads,
            key_size=key_size,
            add_bias_kv=add_bias_kv,
            name="cross_attention_2",
        )

    @hk.transparent
    def mlp(self, x: Embedding) -> Embedding:
        """
        Applies one layer-norm, one linear layer, a Gelu activation,
        then a final linear layer

        Args:
            x: Embeddings of shape (batch_size, seq_len, key_size * num_heads).

        Returns:
            The transformed sequence embedding.
        """
        x = self.layer_norm_mlp(x)
        if self._use_glu_in_ffn:
            x1, x2 = jnp.split(self.fc1(x), indices_or_sections=2, axis=-1)
            x = self._ffn_activation_fn(x1) * x2
        else:
            x = self._ffn_activation_fn(self.fc1(x))
        x = self.fc2(x)
        return x

    def __call__(
        self,
        x: Tokens,
        cross_attention_embeddings_1: Embedding,
        cross_attention_embeddings_2: Embedding,
        attention_mask_1: Optional[AttentionMask] = None,
        attention_mask_2: Optional[AttentionMask] = None,
    ) -> TransformerOutput:
        """
        Computes the output embeddings of the attention block.

        Args:
            x: Input tokens.
            cross_attention_embeddings_1: Embeddings to be used for cross attention
                (in encoder-decoder models, it is the output of the encoder) for the
                first modality.
            cross_attention_embeddings_2: Embeddings to be used for cross attention
                (in encoder-decoder models, it is the output of the encoder) for the
                second modality.
            attention_mask_1: Attention mask of shape (batch_size, 1,seq_len, seq_len)
                for the first modality.
            attention_mask_2: Attention mask of shape (batch_size, 1,seq_len, seq_len)
                for the second modality.

        Returns:
            A dictionary containing the output embeddings and the attention weights.
        """

        # First cross-attention
        res = x
        x = self.layer_norm_cross_attention_1(x)
        output = self.cross_attention_layer_1(
            query=x,
            key=cross_attention_embeddings_1,
            value=cross_attention_embeddings_1,
            attention_mask=attention_mask_1,
        )
        x = output["embeddings"]
        x = res + x

        # Second cross-attention
        res = x
        x = self.layer_norm_cross_attention_2(x)
        output = self.cross_attention_layer_2(
            query=x,
            key=cross_attention_embeddings_2,
            value=cross_attention_embeddings_2,
            attention_mask=attention_mask_2,
        )
        x = output["embeddings"]
        x = res + x

        # MLP
        x = x + self.mlp(x)

        output["embeddings"] = x
        return output  # type: ignore


class MultiModalPerceiverResampler(hk.Module):
    """
    Perceiver Resampler model, made of successive PerceiverResamplerBlocks.
    """

    def __init__(
        self,
        config: PerceiverResamplerConfig,
        name: Optional[str] = None,
    ):
        """
        Initialize a Perceiver Resampler model.

        Args:
            config: Dataclass containing model hyperparameters.
            name: Name for module (custom will break weight loading).
        """

        self._config = config
        super().__init__(name=name)

    @hk.transparent
    def apply_attention_blocks(
        self,
        x: Embedding,
        xf_1: Embedding,
        xf_2: Embedding,
        outs: Dict[str, Embedding],
        attention_mask_1: Optional[AttentionMask] = None,
        attention_mask_2: Optional[AttentionMask] = None,
    ) -> Tuple[Embedding, Dict[str, Embedding]]:
        """
        Create the blocks of attention layers and applies them.

        Args:
            x: Main embedding flowing in the network.
            xf_1: Cross embeddings from first modality to attend too.
            xf_2: Cross embeddings from second modality to attend too.
            outs: A dictionary to carry through the attention layers which stores the
                intermediate sequence embedding and attention maps.
            attention_mask_1: Attention mask of shape (batch_size, 1, seq_len, seq_len)
                for the first modality.
            attention_mask_2: Attention mask of shape (batch_size, 1, seq_len, seq_len)
                for the second modality.

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

        for layer_idx, layer in enumerate(layers):

            concat_input_1 = jnp.concatenate([xf_1, x], axis=1)
            concat_input_2 = jnp.concatenate([xf_2, x], axis=1)

            output = layer(
                x=x,
                cross_attention_embeddings_1=concat_input_1,
                cross_attention_embeddings_2=concat_input_2,
                attention_mask_1=attention_mask_1,
                attention_mask_2=attention_mask_2,
            )
            x = output["embeddings"]

        return x, outs

    @hk.transparent
    def _attention_block(self, layer_idx: int) -> MultiModalPerceiverResamplerBlock:
        return MultiModalPerceiverResamplerBlock(  # type: ignore
            num_heads=self._config.attention_heads,
            embed_dim=self._config.embed_dim,
            key_size=self._config.key_size,
            ffn_embed_dim=self._config.ffn_embed_dim,
            add_bias_kv=self._config.add_bias_kv,
            add_bias_fnn=self._config.add_bias_ffn,
            ffn_activation_name=self._config.ffn_activation_name,
            use_glu_in_ffn=self._config.use_glu_in_ffn,
            name=f"multi_modal_perceiver_layer_{layer_idx}",
        )

    def __call__(
        self,
        input_embeddings_1: Embedding,
        input_embeddings_2: Embedding,
        attention_mask_1: Optional[AttentionMask] = None,
        attention_mask_2: Optional[AttentionMask] = None,
    ) -> TransformerOutput:
        """
        Computes the embeddings based on the input tokens.

        Args:
            input_embeddings_1: Input embeddings from first modality to be resampled.
            input_embeddings_2: Input embeddings from second modality to be resampled.
            attention_mask_1: Attention mask to feed to the attention layers for the
                first modality.
            attention_mask_2: Attention mask to feed to the attention layers for the
                second modality.

        Returns:
            Dictionary containing the final embeddings and logits.
        """
        assert (
            input_embeddings_1.shape[-1] == self._config.embed_dim
        ), "The input embedding dim should match the model embed dim"
        assert (
            input_embeddings_2.shape[-1] == self._config.embed_dim
        ), "The input embedding dim should match the model embed dim"

        batch_size = input_embeddings_1.shape[0]
        # Get latent queries
        w_init = hk.initializers.TruncatedNormal(1.0 / jnp.sqrt(self._config.embed_dim))
        latent_queries = hk.get_parameter(
            "latent_queries",
            shape=[self._config.resampled_length, self._config.embed_dim],
            dtype=input_embeddings_1.dtype,
            init=w_init,
        )
        latent_queries = jnp.expand_dims(latent_queries, axis=0)
        latent_queries = jnp.repeat(latent_queries, repeats=batch_size, axis=0)

        # Prepare outputs dict
        outs: Dict[str, jnp.ndarray] = {}
        x = latent_queries

        # Apply attention blocks
        x, outs = self.apply_attention_blocks(
            x=x,
            xf_1=input_embeddings_1,
            xf_2=input_embeddings_2,
            outs=outs,
            attention_mask_1=attention_mask_1,
            attention_mask_2=attention_mask_2,
        )

        outs["embeddings"] = x

        return outs  # type: ignore


class MultiModalPerceiverResamplerProjection(hk.Module):
    """
    Creates the projection for the DaLlamatian model based on PerceiverResampler.
    (1) DNA embeddings of the sequence will be fed to a linear projection to project
    them to the embedding dimension of the LLM.
    (2) A perceiver resampler is used to reduce the sequence length and create a fixed
    number (resampled_length) of tokens.
    (3) This perceiver resampler also attend to English embeddings to produce a task
    informed embedding of the DNA sequence. Typically, these English embeddings
    would correspond to the question being asked.
    """

    def __init__(
        self,
        config: PerceiverResamplerConfig,
        embed_dim: int,
        bio_pad_token_id: int,
        english_pad_token_id: int,
        english_vocab_size: int,
        name: Optional[str] = None,
    ):
        """
        Initializes the bio to english projection layer.

        Args:
            config: config of the Perceiver Resampler model.
            embed_dim: embedding dimension to project the embeddings to.
            bio_pad_token_id: Index of the bio pad token
            english_pad_token_id: Index of the english pad token.
            name: haiku module name
        """
        super().__init__(name=name)
        self._config = config
        self._embed_dim = embed_dim
        self._bio_pad_token_id = bio_pad_token_id
        self._english_pad_token_id = english_pad_token_id
        self._english_vocab_size = english_vocab_size

    def __call__(
        self,
        bio_token_ids: jnp.ndarray,
        bio_embeddings: jnp.ndarray,
        english_tokens_ids: jnp.ndarray,
    ) -> jnp.ndarray:
        # project bio embeddings to language embed_dim
        projected_bio_embeddings = hk.Linear(
            output_size=self._config.embed_dim,
            name="bio_to_english_projection",
        )(
            bio_embeddings
        )  # (batch_size, dna_seq_len, embed_dim)

        perceiver_resampler = MultiModalPerceiverResampler(config=self._config)
        bio_attention_mask = build_perceiver_padding_attention_mask(
            tokens=bio_token_ids,
            resampled_length=self._config.resampled_length,
            pad_token_id=self._bio_pad_token_id,
        )
        english_attention_mask = build_perceiver_padding_attention_mask(
            tokens=english_tokens_ids,
            resampled_length=self._config.resampled_length,
            pad_token_id=self._english_pad_token_id,
        )

        english_embeddings = hk.Embed(
            vocab_size=self._english_vocab_size,
            embed_dim=self._embed_dim,
            name="token_embed",
        )(english_tokens_ids)
        projected_embeddings = perceiver_resampler(
            input_embeddings_1=projected_bio_embeddings,
            attention_mask_1=bio_attention_mask,
            input_embeddings_2=english_embeddings,
            attention_mask_2=english_attention_mask,
        )[
            "embeddings"
        ]  # (batch_size, resampled_length, embed_dim)

        return projected_embeddings


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


def build_perceiver_padding_attention_mask(
    tokens: Tokens, resampled_length: int, pad_token_id: int
) -> AttentionMask:
    """
    Builds a padding mask from a sequence of tokens by masking <pad> in the attention,
    specific to perceiver resampler.

    Args:
        tokens: Batch of sequences of shape (batch_size, seq_len).
        resampled_length: New length after a PerceiverResampler module. Will be used to
            define the shape of the mask.
        pad_token_id: Int corresponding to the <pad> token to mask.

    Returns:
        Batch of attention masks, masking out <pad> tokens, of shape (batch_size, 1,
            resampled_length, seq_len+resampled_length)
    """
    batch_size, _ = tokens.shape
    padding_mask = tokens != pad_token_id  # (batch_size, seq_len)

    # add 1 for resampled_length
    # (in perceiver we concat [token_embeddings, latent query_embeddings])
    padding_mask = jnp.concatenate(
        (padding_mask, jnp.ones((batch_size, resampled_length))), axis=1
    )  # (batch_size, seq_len + resampled_length)

    padding_mask = padding_mask[
        :, None, None, :
    ]  # (batch_size, 1, seq_len + resampled_length)
    # repeat the mask resampled_length times for latent_query attention
    padding_mask = jnp.tile(padding_mask, (1, 1, resampled_length, 1))
    return padding_mask
