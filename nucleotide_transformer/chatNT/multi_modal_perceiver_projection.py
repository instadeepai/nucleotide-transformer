from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import haiku as hk
import jax.numpy as jnp
import numpy as np
from jax import numpy as jnp

from nucleotide_transformer.layers import MultiHeadAttention
from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)
from nucleotide_transformer.utils import get_activation_fn

UPPER_FREQ = 10000


@dataclass
class PerceiverResamplerConfig:
    """
    Parameters to initialize an PerceiverResampler model.

    Args:
        emb_layer_norm_before: Whether to use layer norm before the first attention
            layer.
        attention_heads: Number of attention heads.
        key_size: The dimension of the query, key, and values within each attention
            head, if not specified, it is set to attention_heads//embed_dim.
            It can be useful to set a custom key size if we want to impose the size of
            the query, key and value tensor ( for example, tensors shaped with
            power of 2 are more efficiently handled on TPUs ).
            Note: Parametrizing the model with a custom key size has been done in :
            Brown, Tom, et al. "Language models are few-shot learners."
            Advances in neural information processing systems 33 (2020): 1877-1901.
        embed_dim: Embedding dimension.
        ffn_embed_dim: Feed forward embedding dimension.
        num_layers: Number of attention blocks.
        ffn_activation_name: Activation function to be used in FFN block. Supported
            names are "gelu", "relu", "swish".
        use_glu_in_ffn: Whether to use Gated Linear Unit (GLU) in Feed
            Forward Network (FFN) block. To do a swiGLU (gated-swish) put this arg
            to True and use swish as ffn_activation_name.
            Same principle for a gated-relu. To keep the same number of parameters in
            the FFN block, one should multiply by 2/3 the ffn_embed_dim when using GLU.
            See https://arxiv.org/pdf/2002.05202.pdf for more details.
        resampled_length: length of the resampled output of the module
        use_gradient_checkpointing: Whether to use gradient checkpointing (checkpoint
            gradients in the forward pass to reduce the computation in the backward).
    """

    # architecture
    emb_layer_norm_before: bool = False
    attention_heads: int = 20
    key_size: Optional[int] = None
    embed_dim: int = 1280
    ffn_embed_dim: int = 5120
    num_layers: int = 24
    add_bias_kv: bool = False
    add_bias_ffn: bool = True
    ffn_activation_name: str = "gelu-no-approx"
    use_glu_in_ffn: bool = False
    resampled_length: int = 64

    # performance
    use_gradient_checkpointing: bool = False

    def __post_init__(self) -> None:
        """
        Checks that the given values are compatible.
        """

        if self.key_size is None:
            if not self.embed_dim % self.attention_heads == 0:
                raise ValueError(
                    f"When no key size is provided, the embedding dimension should be "
                    f"divisible by the number of heads, however provided embedding "
                    f"dimension is {self.embed_dim} and the number of heads is "
                    f"{self.attention_heads}."
                )
            self.key_size = self.embed_dim // self.attention_heads


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
