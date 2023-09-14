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

"""Implementation of the Nucleotide Transformer model in Jax."""
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import haiku as hk
import jax.numpy as jnp
import jmp

from nucleotide_transformer.layers import (
    ESMLearnedPositionalEmbeddings,
    RobertaLMHead,
    SelfAttentionBlock,
    TokensDropout,
)
from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)


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


@dataclass
class NucleotideTransformerConfig:
    """
    Parameters to initialize a Nucleotide Transformer model.

    Args:
        alphabet_size: Token vocabulary.
        pad_token_id: ID of pad token.
        mask_token_id: ID of mask token.
        max_positions: Maximum sequence length.
        embed_scale: Correction ratio applied to the embeddings to make up for the
            norm difference between the input during training and inference.
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
        token_dropout: Token dropout.
        masking_ratio: Masking ratio (used if token dropout is enabled).
        masking_prob: Masking probability (used if token dropout is enabled).
        use_gradient_checkpointing: Whether to use gradient checkpointing (checkpoint
            gradients in the forward pass to reduce the computation in the backward).
    """

    alphabet_size: int
    pad_token_id: int
    mask_token_id: int

    max_positions: int = 1000
    embed_scale: float = 1.0

    # architecture
    emb_layer_norm_before: bool = False
    attention_heads: int = 20
    key_size: Optional[int] = None
    embed_dim: int = 1280
    ffn_embed_dim: int = 5120
    num_layers: int = 24
    positional_embedding: Optional[str] = "learned"
    add_bias_kv: bool = False
    add_bias_ffn: bool = True
    use_rotary_embedding: bool = False
    ffn_activation_name: str = "gelu-no-approx"
    use_glu_in_ffn: bool = False
    layer_norm_eps: float = 1e-5
    pre_layer_norm: bool = True

    # dropout
    token_dropout: bool = False
    masking_ratio: float = 0.1
    masking_prob: float = 0.8

    # logging
    use_gradient_checkpointing: bool = False

    # return
    embeddings_layers_to_save: Tuple[int, ...] = ()
    attention_maps_to_save: List[Tuple[int, int]] = field(default_factory=list)

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


class NucleotideTransformer(hk.Module):
    """
    Jax implementation of Nucleotide Transformer models.
    """

    def __init__(
        self,
        config: NucleotideTransformerConfig,
        name: Optional[str] = None,
    ):
        """
        Initialize a Nucleotide Transformer model.

        Args:
            config: Dataclass containing model hyperparameters.
            name: Name for module (custom will break weight loading).
        """

        self._config = config
        super().__init__(name=name)

        self._embed_layer = hk.Embed(self._config.alphabet_size, self._config.embed_dim)

        if config.positional_embedding == "learned":
            self._pos_embed_layer = ESMLearnedPositionalEmbeddings(
                config.max_positions, config.embed_dim, config.pad_token_id
            )

        self._lm_head = RobertaLMHead(
            embed_dim=self._config.embed_dim,
            alphabet_size=self._config.alphabet_size,
            name="roberta_lm_head",
        )

        if self._config.emb_layer_norm_before:
            self.emb_ln_before = hk.LayerNorm(
                axis=-1,
                create_scale=True,
                create_offset=True,
                name="emb_layer_norm_before",
            )

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
        for layer_idx, layer in enumerate(layers):
            output = layer(
                x=x, attention_mask=attention_mask, attention_weight_bias=None
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
            use_rotary_embedding=self._config.use_rotary_embedding,
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

        # Compute embeddings
        x = self._embed_layer(tokens)
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
            x = x + self._pos_embed_layer(tokens)

        if self._config.emb_layer_norm_before:
            x = self.emb_ln_before(x)

        # Attention mask
        if attention_mask is None:
            attention_mask = build_padding_attention_mask(
                tokens=tokens, pad_token_id=self._config.pad_token_id
            )

        # construct a tower of attention layers
        x, outs = self.apply_attention_blocks(
            x=x,
            outs=outs,
            attention_mask=attention_mask,
        )

        # Language Model Head
        lm_head_outs = self._lm_head(x)
        sequence_mask = attention_mask[:, 0, :, 0][:, :, None]
        outs["logits"] = jnp.where(sequence_mask, lm_head_outs["logits"], 0)

        embeddings = lm_head_outs["embeddings"]
        # Save final embeddings if needed
        if self._config.num_layers in self._config.embeddings_layers_to_save:
            outs[f"embeddings_{self._config.num_layers}"] = embeddings

        return outs  # type: ignore


def build_nucleotide_transformer_fn(
    model_config: NucleotideTransformerConfig,
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    model_name: Optional[str] = None,
) -> Callable:
    """
    Creates the model's forward pass.

    Args:
        model_config: Model hyperparameters.
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too. During
            training, it is preferable to keep parameters in float32 for better
            numerical stability.
        output_dtype: the output type of the model. it determines the float precioson
            of the gradient when training the model.
        model_name: Model's name.

    Returns:
        Nucleotide Transformer model forward function.
    """
    policy = jmp.Policy(
        compute_dtype=compute_dtype, param_dtype=param_dtype, output_dtype=output_dtype
    )
    hk.mixed_precision.set_policy(NucleotideTransformer, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32, param_dtype=param_dtype, output_dtype=compute_dtype
    )
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)

    def nucleotide_transformer_fn(
        tokens: Tokens, attention_mask: Optional[AttentionMask] = None
    ) -> TransformerOutput:
        """Forward pass."""
        # Run the encoder over the inputs.
        encoder = NucleotideTransformer(config=model_config, name=model_name)
        outs = encoder(
            tokens=tokens,
            attention_mask=attention_mask,
        )
        return outs

    return nucleotide_transformer_fn
