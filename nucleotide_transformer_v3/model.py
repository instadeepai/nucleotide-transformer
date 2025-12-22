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

"""Model implementations for Nucleotide Transformer v3."""

import copy
import dataclasses
from dataclasses import asdict, dataclass, field
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx
from flax.typing import Dtype

from nucleotide_transformer_v3.layers import (
    ConvTowerBlock,
    DeconvTowerBlock,
    DeConvUpsampleType,
    LMHead,
    RotaryEmbeddingConfig,
    SelfAttentionBlock,
    Stem,
)
from nucleotide_transformer_v3.adaptive_layers import (
    AdaptiveSelfAttentionBlock,
    ConditionedConvTowerBlock,
    ConditionedDeConvTowerBlock,
)
from nucleotide_transformer_v3.heads import ClassificationHead, LinearHead
from nucleotide_transformer_v3.types import AttentionMask, Embedding, Tokens


# =============================================================================
# Configuration Classes
# =============================================================================


@dataclass(kw_only=True)
class NTv3PreTrainedConfig:
    """
    This architecture used a convolution tower to downsample the sequence length,
    followed by a Transformer torso and a deconvolution tower to upsample the sequence
    length back to its input size.

    Args:
        alphabet_size: number of possible tokens.
        pad_token_id: id of pad token.
        mask_token_id: id of mask token.
        num_downsamples: number of times the sequences length is divided by two
            through convolutions before flowing in the Transformer torso. The sequence
            length seen by the Transformer will be
            initial_seq_length / 2**num_downsamples. E.g. for a sequence length
            of 1M tokens and 8 downsamples, the Transformer will process
            roughly 4k tokens.
        attention_heads: number of heads in the Transformer torso.
        key_size: key size in the Transformer torso.
        token_embed_dim: token embedding dimension.
        conv_init_embed_dim: Embedding dimension of first conv layer.
        embed_dim: Embedding dimension in the Transformer torso.
        ffn_embed_dim: feed forward dimension in the Transformer torso.
        num_layers: number of Transformer layers.
        layer_norm_eps: epsilon for layer norm.
        num_hidden_layers_head: number of hidden layers in head.
        use_remat_in_transformer: whether to use gradient checkpointing in Transformer.
        use_remat_in_convs: whether to use gradient checkpointing in conv towers.
        embeddings_layers_to_save: indices of Transformer layers to save embeddings for.
        attention_maps_to_save: indices of Transformer layers to save attention map for.
        deconv_layers_to_save: indices of de-convolution layers to save embeddings for.
        embedding_compute_dtype: dtype for the embedding layer computation.
        embedding_param_dtype: dtype for the embedding layer parameters.
        stem_compute_dtype: dtype for the stem layer computation.
        stem_param_dtype: dtype for the stem layer parameters.
        up_convolution_compute_dtype: dtype for the up convolution layer computation.
        up_convolution_param_dtype: dtype for the up convolution layer parameters.
        down_convolution_compute_dtype: dtype for the down convolution layer computation.
        down_convolution_param_dtype: dtype for the down convolution layer parameters.
        layernorm_compute_dtype: dtype for the layer normalization computation.
        layernorm_param_dtype: dtype for the layer normalization parameters.
        transformer_qkvo_compute_dtype: dtype for the multi-head attention computation.
        transformer_qkvo_param_dtype: dtype for the multi-head attention parameters.
        transformer_ffn_compute_dtype: dtype for the feed forward network computation.
        transformer_ffn_param_dtype: dtype for the feed forward network parameters.
        lmhead_compute_dtype: dtype for the output layer computation.
        lmhead_param_dtype: dtype for the output layer parameters.
    """

    alphabet_size: int
    pad_token_id: int
    mask_token_id: int

    num_downsamples: int = 8

    # architecture
    attention_heads: int = 16
    key_size: int | None = None
    token_embed_dim: int = 16
    conv_init_embed_dim: int = 512
    embed_dim: int = 512
    ffn_embed_dim: int = 2048
    num_layers: int = 12
    layer_norm_eps: float = 1e-5
    num_hidden_layers_head: int = 0
    use_skip_connection: bool = True

    # optimization
    use_remat_in_transformer: bool = False
    use_remat_in_convs: bool = False

    # return
    embeddings_layers_to_save: tuple[int, ...] = ()
    attention_maps_to_save: list = field(default_factory=list)
    deconv_layers_to_save: tuple[int, ...] = ()

    deconv_upsample_type: DeConvUpsampleType = DeConvUpsampleType.CONV_TRANSPOSE

    # dtypes
    embedding_param_dtype: Dtype = "float32"
    embedding_compute_dtype: Dtype = "float32"
    stem_param_dtype: Dtype = "float32"
    stem_compute_dtype: Dtype = "float32"
    down_convolution_param_dtype: Dtype = "float32"
    down_convolution_compute_dtype: Dtype = "float32"
    up_convolution_param_dtype: Dtype = "float32"
    up_convolution_compute_dtype: Dtype = "float32"
    layernorm_param_dtype: Dtype = "float32"
    layernorm_compute_dtype: Dtype = "float32"
    transformer_qkvo_param_dtype: Dtype = "float32"
    transformer_qkvo_compute_dtype: Dtype = "float32"
    transformer_ffn_param_dtype: Dtype = "float32"
    transformer_ffn_compute_dtype: Dtype = "float32"
    lmhead_param_dtype: Dtype = "float32"
    lmhead_compute_dtype: Dtype = "float32"
    modulation_param_dtype: Dtype = "float32"
    modulation_compute_dtype: Dtype = "float32"

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

        filter_list = list(
            np.linspace(
                self.conv_init_embed_dim, self.embed_dim, self.num_downsamples + 1
            ).astype(int)
        )
        self.filter_list = filter_list


@dataclass(kw_only=True)
class DiscreteConditionedNTv3PreTrainedConfig(NTv3PreTrainedConfig):
    conditions_vocab_size: list[int]
    conditions_names: None | list[str] = None

    @classmethod
    def from_base_config(
        cls,
        config: NTv3PreTrainedConfig,
        conditions_vocab_size: list[int],
        conditions_names: None | list[str] = None,
    ) -> "DiscreteConditionedNTv3PreTrainedConfig":
        """
        Instantiate a DiscreteConditionedNTv3PreTrainedConfig from a NTv3PreTrainedConfig.
        """
        return cls(
            **dataclasses.asdict(config),
            conditions_vocab_size=conditions_vocab_size,
            conditions_names=conditions_names,
        )


@dataclass(kw_only=True)
class NTv3PostTrainedConfig(DiscreteConditionedNTv3PreTrainedConfig):
    """
    Configuration for the NTv3PostTrained model.

    Attributes:
        bigwigs_per_species: A dictionary mapping species names to a list of
            bigwig track names. Used to store the tracks the model was trained on.
        bed_elements_names: A list of BED element names to predict.
        species_to_token_id: A dictionary mapping species names to their token IDs.
        keep_target_center_fraction: The fraction of the sequence to predict the tracks
            for. Set to 1.0 to predict the tracks for the entire sequence.
        num_species_special_tokens: Number of special tokens reserved for species vocabulary.
    """

    bigwigs_per_species: dict[str, list[str]]
    bed_elements_names: list[str]
    species_to_token_id: dict[str, int]
    keep_target_center_fraction: float = 1.0
    num_species_special_tokens: int = 6

    def __post_init__(self) -> None:
        super().__post_init__()

    @classmethod
    def from_base_config(
        cls,
        config: DiscreteConditionedNTv3PreTrainedConfig,
        bigwigs_per_species: dict[str, list[str]],
        bed_elements_names: list[str],
        species_to_token_id: dict[str, int],
        keep_target_center_fraction: float = 1.0,
        num_species_special_tokens: int = 6,
    ) -> "NTv3PostTrainedConfig":
        return cls(
            **asdict(config),
            bigwigs_per_species=bigwigs_per_species,
            bed_elements_names=bed_elements_names[:],
            species_to_token_id=species_to_token_id,
            keep_target_center_fraction=keep_target_center_fraction,
            num_species_special_tokens=num_species_special_tokens,
        )

class NTv3PreTrained(nnx.Module):
    def __init__(
        self,
        config: NTv3PreTrainedConfig,
        *,
        rngs: nnx.Rngs,
    ):
        self.config = config

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

        # --- Define layers
        # Embedding layers
        self.embed_layer = nnx.Embed(
            self.config.alphabet_size,
            self.config.token_embed_dim,
            rngs=rngs,
            dtype=self.config.embedding_compute_dtype,
            param_dtype=self.config.embedding_param_dtype,
        )
        self._rotary_embedding_config = RotaryEmbeddingConfig(rescaling_factor=None)

        # Stem layer
        self.stem = Stem(
            token_embed_dim=self.config.token_embed_dim,
            conv_init_embed_dim=self.config.conv_init_embed_dim,
            use_remat=self.config.use_remat_in_convs,
            dtype=self.config.stem_compute_dtype,
            param_dtype=self.config.stem_param_dtype,
            rngs=rngs,
        )

        # Convolution tower
        filter_list = copy.deepcopy(self.config.filter_list)
        self.conv_tower_blocks = []

        for dim_in, dim_out in zip(filter_list[:-1], filter_list[1:]):
            tower_block = ConvTowerBlock(
                dim_in=dim_in,
                dim_out=dim_out,
                dtype=self.config.down_convolution_compute_dtype,
                param_dtype=self.config.down_convolution_param_dtype,
                ln_dtype=self.config.layernorm_compute_dtype,
                ln_param_dtype=self.config.layernorm_param_dtype,
                rngs=rngs,
            )
            self.conv_tower_blocks.append(tower_block)

        # Transformer tower
        self.transformer_blocks = []
        for _ in range(self.config.num_layers):
            attention_block = self._attention_block(rngs)
            self.transformer_blocks.append(attention_block)

        # Deconvolution tower
        filter_list = copy.deepcopy(self.config.filter_list)
        filter_list.reverse()
        self.deconv_tower_blocks = []

        for dim_in, dim_out in zip(filter_list[:-1], filter_list[1:]):
            tower_block = DeconvTowerBlock(
                dim_in=dim_in,
                dim_out=dim_out,
                upsample_type=self.config.deconv_upsample_type,
                dtype=self.config.up_convolution_compute_dtype,
                param_dtype=self.config.up_convolution_param_dtype,
                ln_dtype=self.config.layernorm_compute_dtype,
                ln_param_dtype=self.config.layernorm_param_dtype,
                rngs=rngs,
            )
            self.deconv_tower_blocks.append(tower_block)

        # LM head
        self.lm_head = LMHead(
            embed_dim=self.config.conv_init_embed_dim,
            num_hidden_layers=self.config.num_hidden_layers_head,
            alphabet_size=self.config.alphabet_size,
            dtype=self.config.lmhead_compute_dtype,
            param_dtype=self.config.lmhead_param_dtype,
            rngs=rngs,
        )

    def transformer_tower(
        self,
        x: Embedding,
        outs: dict[str, Embedding],
        attention_mask: AttentionMask | None = None,
    ) -> tuple[Embedding, dict[str, Embedding]]:
        def apply_attention_block(
            block: SelfAttentionBlock,
            x: Embedding,
            attention_mask: AttentionMask | None,
            attention_weight_bias: jax.Array | None,
        ) -> dict[str, Embedding]:
            return block(  # type: ignore
                x=x,
                attention_mask=attention_mask,
                attention_weight_bias=attention_weight_bias,
            )

        if self.config.use_remat_in_transformer:
            apply_attention_block = nnx.remat(apply_attention_block)

        for layer_idx, layer in enumerate(self.transformer_blocks):
            output = apply_attention_block(
                layer, x=x, attention_mask=attention_mask, attention_weight_bias=None
            )
            x = output["embeddings"]

            # Save intermediate embeddings if needed
            if (layer_idx + 1) in self.config.embeddings_layers_to_save:
                outs[f"embeddings_{(layer_idx + 1)}"] = output["embeddings"]
            # Save intermediate attention maps if needed
            if (layer_idx + 1) in self._attention_layers_to_save:
                for map_number in self._attention_maps_per_layer_to_save[layer_idx + 1]:
                    dkey = f"attention_map_layer_{layer_idx + 1}_number_{map_number}"
                    outs[dkey] = output["attention_weights"][:, map_number + 1]

        return x, outs

    def conv_tower(self, x: jax.Array) -> tuple[jax.Array, list[jax.Array]]:
        residuals = []

        def apply_conv_block(block: nnx.Module, x: jax.Array) -> jax.Array:
            return block(x)

        if self.config.use_remat_in_convs:
            apply_conv_block = nnx.remat(apply_conv_block)

        for tower_block in self.conv_tower_blocks:
            residuals.append(x)
            x = apply_conv_block(tower_block, x)

        return x, residuals

    def deconv_tower(
        self,
        x: jax.Array,
        residuals: list[jax.Array],
        outs: dict[str, Embedding],
    ) -> tuple[jax.Array, dict[str, jax.Array]]:
        residuals_generator = reversed(residuals)

        def apply_deconv_block(block: DeconvTowerBlock, x: jax.Array) -> jax.Array:
            return block(x)

        if self.config.use_remat_in_convs:
            apply_deconv_block = nnx.remat(apply_deconv_block)

        for i, tower_block in enumerate(self.deconv_tower_blocks):
            x = apply_deconv_block(tower_block, x)

            if self.config.use_skip_connection:
                residuals = next(residuals_generator)
                x = x + residuals

            # Save intermediate embeddings if needed
            if (i + 1) in self.config.deconv_layers_to_save:
                outs[f"embeddings_deconv_{(i + 1)}"] = x

        return x, outs

    def _attention_block(self, rngs: nnx.Rngs) -> SelfAttentionBlock:
        return SelfAttentionBlock(  # type: ignore
            num_heads=self.config.attention_heads,
            embed_dim=self.config.embed_dim,
            key_size=self.config.key_size,
            ffn_embed_dim=self.config.ffn_embed_dim,
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            rotary_embedding_config=self._rotary_embedding_config,
            layer_norm_eps=self.config.layer_norm_eps,
            pre_layer_norm=True,
            ffn_dtype=self.config.transformer_ffn_compute_dtype,
            ffn_param_dtype=self.config.transformer_ffn_param_dtype,
            mha_dtype=self.config.transformer_qkvo_compute_dtype,
            mha_param_dtype=self.config.transformer_qkvo_param_dtype,
            ln_dtype=self.config.layernorm_compute_dtype,
            ln_param_dtype=self.config.layernorm_param_dtype,
            rngs=rngs,
        )

    def _attention_mask(self, tokens: Tokens) -> AttentionMask:
        """
        Computes the attention mask for the transformer.

        Args:
            tokens: tokens of shape (batch_size, sequence_length,).

        Returns:
            Attention mask of shape
                (batch_size, 1, block_size, block_size),
                where block_size is 2**num_downsamples.
        """
        block_size = 2**self.config.num_downsamples
        batch_size, seq_len = tokens.shape
        assert (
            seq_len % block_size == 0
        ), f"Sequence length ({seq_len}) must be divisible by block_size {block_size}"

        seq_len_reduced = seq_len // block_size

        # Identify non-pad tokens
        is_not_pad = (tokens != self.config.pad_token_id).astype(jnp.int32)

        # Group into blocks of size block_size
        is_not_pad_blocks = is_not_pad.reshape(batch_size, seq_len_reduced, block_size)

        # A block is valid if it contains any non-pad token
        block_valid = (is_not_pad_blocks.sum(axis=-1) > 0).astype(
            jnp.int32
        )  # (batch, seq_len_reduced)

        # Broadcast into attention mask
        # Shape: (batch, 1, seq_len_reduced, seq_len_reduced)
        attention_mask = block_valid[:, None, None, :]
        attention_mask = jnp.broadcast_to(
            attention_mask, (batch_size, 1, seq_len_reduced, seq_len_reduced)
        )

        return attention_mask

    def __call__(
        self, tokens: Tokens, **kwargs: dict[str, Any]
    ) -> dict[str, Embedding]:
        """
        Forward pass in NTv3PreTrained.

        Args:
            tokens: tokens of shape (batch_size, sequence_length,)

        Returns:
            outs: dictionary of outputs:
                - "embedding": the final embedding of shape
                    (batch_size, sequence_length, embed_dim)
                - "after_transformer_embedding": the embedding after transformer
                    of shape (batch_size, sequence_length, embed_dim)
                - "logits": the logits of
                    shape (batch_size, sequence_length, alphabet_size)
        """
        _, seq_len = tokens.shape
        # Assert that the sequence length is divisible by 2**num_downsamples
        assert seq_len % (2**self.config.num_downsamples) == 0, (
            f"Sequence length ({seq_len}) "
            "should be divisible by 2 to the power of the number of downsample "
            f" layers ({self.config.num_downsamples})"
        )
        # Prepare outputs dict
        outs: dict[str, Embedding] = {}

        # Compute embeddings
        x = self.embed_layer(tokens)

        # Main body
        x = self.stem(x)

        x, residuals = self.conv_tower(x)

        x, outs = self.transformer_tower(
            x, outs=outs, attention_mask=self._attention_mask(tokens)
        )
        outs["after_transformer_embedding"] = x

        x, outs = self.deconv_tower(x, residuals=residuals, outs=outs)
        outs["embedding"] = x

        x = self.lm_head(x)

        outs["logits"] = x

        return outs  # type: ignore


class ConditionedNTv3PreTrained(NTv3PreTrained):
    """
    NTv3PreTrained with conditions support. Assume that the conditions are a list of
    arrays of shape (batch_size, condition_dim). The model uses adaptive layer norm
    to condition the layers on the conditions. Conditioning is done through an additive
    mechanism that ensures that the initial forward pass of the model is not modified
    by adding one or several conditions. This is done by rescaling and shifting layers
    with rescaling factors initialized to ones and shifts initialized to zeros.
    """

    def __init__(
        self,
        config: NTv3PreTrainedConfig,
        conditions_dims: list[int],
        *,
        rngs: nnx.Rngs,
    ):
        self.conditions_dims = conditions_dims
        super().__init__(config=config, rngs=rngs)

        # --- Redefine layers
        # Convolution tower
        filter_list = copy.deepcopy(self.config.filter_list)
        self.conv_tower_blocks = []

        for dim_in, dim_out in zip(filter_list[:-1], filter_list[1:]):
            tower_block = ConditionedConvTowerBlock(
                dim_in=dim_in,
                dim_out=dim_out,
                conditions_dims=conditions_dims,
                dtype=self.config.down_convolution_compute_dtype,
                param_dtype=self.config.down_convolution_param_dtype,
                ln_dtype=self.config.layernorm_compute_dtype,
                ln_param_dtype=self.config.layernorm_param_dtype,
                modulation_dtype=self.config.modulation_compute_dtype,
                modulation_param_dtype=self.config.modulation_param_dtype,
                rngs=rngs,
            )
            self.conv_tower_blocks.append(tower_block)

        # Transformer tower
        self.transformer_blocks = []
        for _ in range(self.config.num_layers):
            attention_block = self._attention_block(rngs)
            self.transformer_blocks.append(attention_block)

        # Deconvolution tower
        filter_list = copy.deepcopy(self.config.filter_list)
        filter_list.reverse()
        self.deconv_tower_blocks = []

        for dim_in, dim_out in zip(filter_list[:-1], filter_list[1:]):
            tower_block = ConditionedDeConvTowerBlock(
                dim_in=dim_in,
                dim_out=dim_out,
                conditions_dims=conditions_dims,
                upsample_type=self.config.deconv_upsample_type,
                dtype=self.config.up_convolution_compute_dtype,
                param_dtype=self.config.up_convolution_param_dtype,
                ln_dtype=self.config.layernorm_compute_dtype,
                ln_param_dtype=self.config.layernorm_param_dtype,
                modulation_dtype=self.config.modulation_compute_dtype,
                modulation_param_dtype=self.config.modulation_param_dtype,
                rngs=rngs,
            )
            self.deconv_tower_blocks.append(tower_block)

    def conv_tower(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> tuple[jax.Array, list[jax.Array]]:
        residuals = []

        # Build layers
        def apply_block(
            block: ConditionedConvTowerBlock,
            x: jax.Array,
            conditions: list[jax.Array],
            conditions_masks: list[jax.Array] | None = None,
        ) -> jax.Array:
            return block(x, conditions, conditions_masks)

        if self.config.use_remat_in_convs:
            apply_block = nnx.remat(apply_block)

        for tower_block in self.conv_tower_blocks:
            residuals.append(x)
            x = apply_block(tower_block, x, conditions, conditions_masks)

        return x, residuals

    def deconv_tower(  # type: ignore
        self,
        x: jax.Array,
        residuals: list[jax.Array],
        outs: dict[str, Embedding],
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> tuple[jax.Array, dict[str, Embedding]]:
        filter_list = copy.deepcopy(self.config.filter_list)
        filter_list.reverse()
        residuals_generator = reversed(residuals)

        def apply_conv_block(
            block: ConditionedDeConvTowerBlock,
            x: jax.Array,
            conditions: list[jax.Array],
            conditions_masks: list[jax.Array] | None = None,
        ) -> jax.Array:
            return block(x, conditions, conditions_masks)

        if self.config.use_remat_in_convs:
            apply_conv_block = nnx.remat(apply_conv_block)

        # Build layers
        for i, tower_block in enumerate(self.deconv_tower_blocks):
            x = apply_conv_block(tower_block, x, conditions, conditions_masks)

            if self.config.use_skip_connection:
                current_residual = next(residuals_generator)
                x = x + current_residual
            # Save intermediate embeddings if needed
            if (i + 1) in self.config.deconv_layers_to_save:
                outs[f"embeddings_deconv_{(i + 1)}"] = x

        return x, outs

    def transformer_tower(  # type: ignore
        self,
        x: Embedding,
        outs: dict[str, Embedding],
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        attention_mask: AttentionMask | None = None,
    ) -> tuple[Embedding, dict[str, Embedding]]:
        def apply_attention_block(
            block: AdaptiveSelfAttentionBlock,
            x: Embedding,
            conditions: list[jax.Array],
            conditions_masks: list[jax.Array] | None = None,
            attention_mask: AttentionMask | None = None,
        ) -> Embedding:
            return block(
                x=x,
                conditions=conditions,
                conditions_masks=conditions_masks,
                attention_mask=attention_mask,
            )

        if self.config.use_remat_in_transformer:
            apply_attention_block = nnx.remat(apply_attention_block)
        for layer_idx, layer in enumerate(self.transformer_blocks):
            output = apply_attention_block(
                layer,
                x=x,
                conditions=conditions,
                attention_mask=attention_mask,
                conditions_masks=conditions_masks,
            )
            x = output["embeddings"]

            # Save intermediate embeddings if needed
            if (layer_idx + 1) in self.config.embeddings_layers_to_save:
                outs[f"embeddings_{(layer_idx + 1)}"] = output["embeddings"]
            # Save intermediate attention maps if needed
            if (layer_idx + 1) in self._attention_layers_to_save:
                for map_number in self._attention_maps_per_layer_to_save[layer_idx + 1]:
                    dkey = f"attention_map_layer_{layer_idx + 1}_number_{map_number}"
                    outs[dkey] = output["attention_weights"][:, map_number + 1]

        return x, outs

    def _attention_block(self, rngs: nnx.Rngs) -> AdaptiveSelfAttentionBlock:
        return AdaptiveSelfAttentionBlock(  # type: ignore
            num_heads=self.config.attention_heads,
            embed_dim=self.config.embed_dim,
            key_size=self.config.key_size,
            ffn_embed_dim=self.config.ffn_embed_dim,
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            rotary_embedding_config=self._rotary_embedding_config,
            layer_norm_eps=self.config.layer_norm_eps,
            pre_layer_norm=True,
            conditions_dims=self.conditions_dims,
            ffn_dtype=self.config.transformer_ffn_compute_dtype,
            ffn_param_dtype=self.config.transformer_ffn_param_dtype,
            mha_dtype=self.config.transformer_qkvo_compute_dtype,
            mha_param_dtype=self.config.transformer_qkvo_param_dtype,
            ln_dtype=self.config.layernorm_compute_dtype,
            ln_param_dtype=self.config.layernorm_param_dtype,
            modulation_dtype=self.config.modulation_compute_dtype,
            modulation_param_dtype=self.config.modulation_param_dtype,
            rngs=rngs,
        )

    def __call__(  # type: ignore
        self,
        tokens: Tokens,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        **kwargs: Any,
    ) -> dict[str, Embedding]:
        """
        Args:
            tokens: tokens of shape (batch_size, sequence_length,).
            conditions: list of conditions of shape (batch_size, condition_dim).
            conditions_masks: list of masks of shape (batch_size,). If None,
                no mask is applied. Mask is used to dropout conditions tokens where
                here dropout means that the condition is not applied in the adaptive
                layer norm.
        """
        seq_len = tokens.shape[1]
        # Assert that the sequence length is divisible by 2**num_downsamples
        assert seq_len % (2**self.config.num_downsamples) == 0, (
            f"Sequence length ({seq_len}) "
            "should be divisible by 2 to the power of the number of downsample "
            f" layers ({self.config.num_downsamples})"
        )
        # Prepare outputs dict
        outs: dict[str, Embedding] = {}

        # Compute embeddings
        x = self.embed_layer(tokens)

        # Main body
        with jax.named_scope("stem"):
            x = self.stem(x)

        with jax.named_scope("conv_tower"):
            x, residuals = self.conv_tower(
                x=x, conditions=conditions, conditions_masks=conditions_masks
            )

        outs["before_transformer_embedding"] = x
        with jax.named_scope("transformer_tower"):
            x, outs = self.transformer_tower(
                x,
                conditions=conditions,
                conditions_masks=conditions_masks,
                outs=outs,
                attention_mask=None,
            )
        outs["after_transformer_embedding"] = x

        with jax.named_scope("deconv_tower"):
            x, outs = self.deconv_tower(
                x=x,
                conditions=conditions,
                conditions_masks=conditions_masks,
                residuals=residuals,
                outs=outs,
            )
        outs["embedding"] = x

        x = self.lm_head(x=x)

        outs["logits"] = x

        return outs  # type: ignore


class DiscreteConditionedNTv3PreTrained(ConditionedNTv3PreTrained):
    """
    NTv3PreTrained with conditioning on one or several discrete tokens. These tokens can
    for instance represent a species, a tissue, an assay, an activity level, etc.
    It leverages the conditioned NTv3PreTrained model and calls it after embedding the
    discrete token(s). This model returns logits over the sequence as well as
    logits over each discrete token.
    """

    def __init__(
        self,
        config: DiscreteConditionedNTv3PreTrainedConfig,
        *,
        rngs: nnx.Rngs,
    ):
        super().__init__(
            config=config,
            conditions_dims=[
                config.token_embed_dim for _ in config.conditions_vocab_size
            ],
            rngs=rngs,
        )

        conditions_names = config.conditions_names
        conditions_vocab_size = config.conditions_vocab_size

        if conditions_names is None:
            conditions_names = [
                f"condition_{i}" for i in range(len(conditions_vocab_size))
            ]

        if not len(conditions_names) == len(conditions_vocab_size):
            raise ValueError(
                "conditions_names and conditions_vocab_size must have the same length"
            )

        self._conditions_vocab_size = conditions_vocab_size
        self._conditions_names = conditions_names

        # Define condition embed layers
        self.conditions_embed_layers = [
            nnx.Embed(
                num_embeddings=vocab_size,
                features=self.config.token_embed_dim,
                rngs=rngs,
            )
            for vocab_size in conditions_vocab_size
        ]

        # Define conditions heads
        self.conditions_heads = [
            nnx.Linear(
                in_features=config.embed_dim,
                out_features=vocab_size,
                rngs=rngs,
            )
            for vocab_size in self._conditions_vocab_size
        ]

    def __call__(
        self,
        tokens: Tokens,
        conditions_tokens: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        **kwargs: Any,
    ) -> dict[str, Embedding]:
        """
        Args:
            tokens: tokens of shape (batch_size, sequence_length,)
            conditions_tokens: list of condition tokens of shape (batch_size,)
            conditions_masks: list of boolean masks of shape (batch_size,). If None,
                no mask is applied. Mask is used to dropout conditions tokens where
                here dropout means that the condition is not applied in the adaptive
                layer norm.
        """
        # Embed each condition token, concatenate embeddings and do forward pass
        conditions = [
            embed_layer(condition_tokens)
            for (condition_tokens, embed_layer) in zip(
                conditions_tokens, self.conditions_embed_layers
            )
        ]
        # Pass in adaptive NTv3PreTrained model
        outs = super().__call__(tokens, conditions, conditions_masks)  # type: ignore

        # Define additional heads to generate probabilities over conditions
        # Can be used for joint modelling of tokens and conditions
        # These heads are plugged on top of the transformer output
        transformer_outs = outs["after_transformer_embedding"]
        transformer_outs = jnp.mean(transformer_outs, axis=-2)
        transformer_outs = jax.nn.gelu(transformer_outs)
        heads_outs = {
            f"{name}_logits": cond_head(transformer_outs)
            for name, cond_head in zip(self._conditions_names, self.conditions_heads)
        }
        outs.update(heads_outs)

        return outs  # type: ignore


class MultiSpeciesHead(nnx.Module):
    """
    A multi-species prediction head that predicts a set of tracks for each species.

    The head has a linear head for each species. The number of tracks is padded to a
    constant value, that is the maximum number of tracks in all species. The head
    returns the logits for the tracks of the selected species.
    """

    def __init__(
        self,
        num_tracks_per_head: list[int],
        embed_dim: int,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            num_tracks_per_head: A list of the number of tracks that will be predicted
                from each head.
            embed_dim: The dimension of the embeddings.
            rngs: The random number generators.

        The head predicts logits with shape (..., max(num_tracks_per_head)). The first,
        num_tracks_per_head[i] tracks are the actual predictions for the i-th head, and
        the remaining tracks are padded with 0.
        """
        self.num_output_tracks = max(num_tracks_per_head)
        self.species_heads = [
            LinearHead(
                embed_dim=embed_dim,
                num_labels=num_tracks_per_head_i,
                rngs=rngs,
            )
            if num_tracks_per_head_i > 0
            else None
            for num_tracks_per_head_i in num_tracks_per_head
        ]

    def __call__(self, x: jax.Array, head_idx: jax.Array) -> jax.Array:
        logits = jnp.zeros((*x.shape[:-1], self.num_output_tracks))
        for i, head in enumerate(self.species_heads):
            if head is None:
                continue
            head_logits = head(x)
            head_padding = self.num_output_tracks - head_logits.shape[-1]
            padding_tuple = ((0, 0), (0, 0), (0, head_padding))
            padded_logits = jnp.pad(
                head_logits, padding_tuple, mode="constant", constant_values=0
            )
            mask = head_idx == i
            mask = mask[:, None, None]  # add dimensions for sequence length and tracks
            logits = jnp.where(mask, padded_logits, logits)

        return logits


def crop_sequence_center(
    sequences_array: jax.Array, keep_target_center_fraction: float
) -> jax.Array:
    """
    Crops an array to keep only the center portion specified by keep_center_fraction.

    Args:
        sequences_array: Input array of shape (batch_size, sequence_length, ...).
        keep_target_center_fraction: Fraction of the sequence to keep in the center
            (0.0 to 1.0). For example, 0.375 means keeping the middle 37.5% of the seq.

    Returns:
        Cropped array with the same number of dimensions but reduced sequence length.
    """
    assert (
        0.0 < keep_target_center_fraction <= 1.0
    ), f"{keep_target_center_fraction=} must be in (0, 1]"
    sequence_length = sequences_array.shape[1]
    crop_length = int((1.0 - keep_target_center_fraction) * sequence_length)
    start_idx = crop_length // 2
    end_idx = sequence_length - crop_length // 2

    return sequences_array[:, start_idx:end_idx, ...]


class NTv3PostTrained(DiscreteConditionedNTv3PreTrained):
    config: NTv3PostTrainedConfig  # type: ignore

    def __init__(
        self,
        config: NTv3PostTrainedConfig,
        rngs: nnx.Rngs,
    ):
        super().__init__(config, rngs=rngs)

        # Bigwig head (functional track predictions)
        if config.bigwigs_per_species:
            # We take the convention that the order of the species heads is the order 
            # of the species token IDs
            sorted_species = sorted(
                config.bigwigs_per_species.keys(),
                key=lambda s: config.species_to_token_id[s],
            )
            sorted_num_tracks_per_head = [
                len(config.bigwigs_per_species[s]) for s in sorted_species
            ]
            self.bigwig_head = MultiSpeciesHead(
                num_tracks_per_head=sorted_num_tracks_per_head,
                embed_dim=config.embed_dim,
                rngs=rngs,
            )

        # BED head (genome annotation predictions)
        if len(config.bed_elements_names) > 0:
            self.bed_head = ClassificationHead(
                embed_dim=config.embed_dim,
                num_labels=len(config.bed_elements_names),
                num_classes=2,
                rngs=rngs,
            )

    def __call__(
        self,
        tokens: Tokens,
        species_tokens: jax.Array,
    ) -> dict[str, jax.Array]:
        """
        Forward pass for the model.

        Args:
            tokens: The tokens to encode.
            species_tokens: The tokens for the species, shape (batch_size, ).
        """
        # We do not support special tokens for species; this is used to index
        # the species heads
        assert (
            species_tokens.min() >= self.config.num_species_special_tokens 
        ), (
            f"No special tokens are supported for the species, species_tokens must be "
            f"greater than or equal to {self.config.num_species_special_tokens}"
        )
        assert (
            species_tokens.shape == (tokens.shape[0],)
        ), (
            f"species_tokens must be shape (batch_size,), got {species_tokens.shape}"
        )
        output_dict: dict[str, jax.Array] = super().__call__(
            tokens, [species_tokens], None
        )

        # Crop the sequence center
        embedding = crop_sequence_center(
            output_dict["embedding"],  # type: ignore
            self.config.keep_target_center_fraction,  # type: ignore
        )

        # Bigwig head (functional track predictions)
        if hasattr(self, "bigwig_head"):
            # Only call the forward pass on the selected head
            head_indices = species_tokens - self.config.num_species_special_tokens
            bigwig_logits = self.bigwig_head(embedding, head_indices)
            output_dict["bigwig_tracks_logits"] = bigwig_logits

        # BED head (genome annotation predictions)
        if hasattr(self, "bed_head"):
            output_dict["bed_tracks_logits"] = self.bed_head(embedding)
        return output_dict  # type: ignore
    @property
    def supported_species(self) -> list[str]:
        """List of supported species names (excludes special tokens)."""
        return sorted([k for k in self.config.species_to_token_id.keys() if not k.startswith("<")])

    def encode_species(self, species: str | list[str]) -> Tokens:
        """
        Encode species name(s) to token IDs.

        Args:
            species: Species name(s) (e.g., "human" or ["human", "mouse"]).
                Use `model.supported_species` to see all valid options.

        Returns:
            Tokens of shape (len(species),) containing token IDs.

        Raises:
            ValueError: If any species name is not supported.

        Example:
            >>> print(model.supported_species)  # See valid species
            >>> species_ids = model.encode_species(["human", "mouse"])
            >>> out = model(input_ids=tokens, species_ids=species_ids)
        """
        if isinstance(species, str):
            species = [species]
        token_ids = []
        for s in species:
            if s not in self.supported_species:
                raise ValueError(
                    f"Unknown species '{s}'. Supported species: {self.supported_species}"
                )
            token_ids.append(self.config.species_to_token_id[s])
        return jnp.array(token_ids, dtype=jnp.int32)
