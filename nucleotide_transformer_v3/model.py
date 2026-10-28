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

"""NTv3 base model: U-Net encoder/decoder with transformer tower."""

from typing import Any

import jax
import jax.numpy as jnp
from flax import nnx
from flax.linen.spmd import with_logical_constraint
from jax.sharding import PartitionSpec

from nucleotide_transformer_v3.heads import LMHead
from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.layers import (
    DownResBlock,
    GQABlock,
    MLPBlock,
    TokenDnaEmbedder,
    UpResBlock,
    pool,
)
from nucleotide_transformer_v3.types import AttentionMask, Embedding, Tokens


class NTv3(nnx.Module):
    """U-Net encoder/decoder with a transformer tower at the bottleneck.

    Pipeline::

        TokenEmbedding -> TokenDnaEmbedder -> MaxPool 2x
        -> DownResBlock x num_downsamples (encoder, 2x downsample each)
        -> [GQABlock + MLPBlock] x num_layers (transformer tower)
        -> RMSNorm (final trunk norm)
        -> UpResBlock x (num_downsamples + 1) (decoder, 2x upsample each, U-Net skips)
        -> LMHead

    Sequence length must be divisible by ``2 ** (num_downsamples + 1)`` (the
    initial 2x pool plus one halving per encoder block). Decoder symmetry
    restores the input length at the LM head.
    """

    def __init__(
        self,
        config: NTv3Config,
        *,
        rngs: nnx.Rngs,
    ) -> None:
        self.config = config

        channels = config.encoder_channels
        n_down = config.num_downsamples

        self._encoder_bin_sizes = tuple(2 ** (i + 1) for i in range(n_down))
        self._decoder_bin_sizes = tuple(2**i for i in range(n_down, -1, -1))

        # --- Token embedding ---
        self.embed_layer = nnx.Embed(
            num_embeddings=config.alphabet_size,
            features=config.token_embed_dim,
            rngs=rngs,
        )

        # --- Stem ---
        self.stem = TokenDnaEmbedder(
            token_embed_dim=config.token_embed_dim,
            conv_init_embed_dim=config.conv_init_embed_dim,
            dtype=config.stem_compute_dtype,
            param_dtype=config.stem_param_dtype,
            rngs=rngs,
        )

        # --- Encoder ---
        self.encoder_blocks = [
            DownResBlock(
                in_channels=channels[i],
                out_channels=channels[i + 1],
                dtype=config.down_convolution_compute_dtype,
                param_dtype=config.down_convolution_param_dtype,
                rngs=rngs,
            )
            for i in range(n_down)
        ]

        # --- Transformer tower ---
        # FiLM conditioning (empty conditions_dims => unconditioned).
        conditions_dims = config.conditions_dims
        self.mha_blocks = [
            GQABlock(
                embed_dim=config.embed_dim,
                attention_heads=config.attention_heads,
                num_kv_heads=config.num_kv_heads,
                key_size=config.key_size,
                value_size=config.value_size,
                attention_implementation=config.attention_implementation,
                conditions_dims=conditions_dims,
                dtype=config.transformer_qkvo_compute_dtype,
                param_dtype=config.transformer_qkvo_param_dtype,
                rngs=rngs,
            )
            for _ in range(config.num_layers)
        ]
        self.mlp_blocks = [
            MLPBlock(
                embed_dim=config.embed_dim,
                conditions_dims=conditions_dims,
                dtype=config.transformer_ffn_compute_dtype,
                param_dtype=config.transformer_ffn_param_dtype,
                rngs=rngs,
            )
            for _ in range(config.num_layers)
        ]

        # --- Final trunk norm ---
        # Bounds the residual stream before the decoder / LM head consume it.
        # Identity-init: gamma=1. Always fp32 (dtype=None default) regardless
        # of the transformer/decoder compute dtypes, like every other norm.
        self.final_norm = nnx.RMSNorm(num_features=config.embed_dim, rngs=rngs)

        # --- Decoder ---
        self.decoder_blocks = [
            UpResBlock(
                x_channels=x_ch,
                unet_skip_channels=skip_ch,
                conditions_dims=conditions_dims,
                dtype=config.up_convolution_compute_dtype,
                param_dtype=config.up_convolution_param_dtype,
                rngs=rngs,
            )
            for x_ch, skip_ch in zip(
                config.decoder_x_channels, config.decoder_skip_channels, strict=True
            )
        ]

        # --- LM head ---
        self.lm_head = LMHead(
            embed_dim=config.conv_init_embed_dim,
            num_hidden_layers=0,
            alphabet_size=config.alphabet_size,
            rngs=rngs,
        )

    def _block_attention_mask(self, tokens: Tokens) -> AttentionMask:
        """Build a block-level pad-aware attention mask.

        Each transformer position covers ``total_downsampling`` input tokens.
        A position is valid iff its block contains at least one non-pad
        token; the mask is broadcast to ``(batch, 1, T, T)`` over the K axis.
        ``True`` = attend, ``False`` = mask out (the convention
        ``jax.nn.dot_product_attention`` expects).
        """
        block_size = self.config.total_downsampling
        batch_size, seq_len = tokens.shape
        seq_len_reduced = seq_len // block_size

        is_not_pad = tokens != self.config.pad_token_id
        is_not_pad_blocks = is_not_pad.reshape(batch_size, seq_len_reduced, block_size)
        block_valid = is_not_pad_blocks.any(axis=-1)

        return jnp.broadcast_to(
            block_valid[:, None, None, :],
            (batch_size, 1, seq_len_reduced, seq_len_reduced),
        )

    def __call__(
        self,
        tokens: Tokens,
        conditions: list[jax.Array] | None = None,
        conditions_masks: list[jax.Array] | None = None,
        stop_after_transformer: bool = False,
        skip_lm_head: bool = False,
        **kwargs: Any,
    ) -> dict[str, Embedding]:
        """Forward pass.

        Args:
            tokens: integer token ids of shape ``(batch_size, sequence_length)``.
            conditions: optional list of per-example condition embeddings, each
                of shape ``(batch_size, cond_dim_i)``. Its length must match
                ``config.conditions_dims``. ``None`` (the default) runs the
                model unconditioned (required when ``conditions_dims`` is empty).
            conditions_masks: optional list of per-example masks, each of shape
                ``(batch_size,)``. Where ``False``, the corresponding condition
                contributes zero FiLM modulation (condition dropout).
            stop_after_transformer: Return as soon as
                ``"after_transformer_embedding"`` is available, skipping the
                decoder tower and the LM head. For heads that only consume the
                trunk output (e.g. a pooled gene-level regression head) the
                decoder is dead compute, and at long context lengths its
                1bp-resolution activations dominate peak memory. The returned
                ``"after_transformer_embedding"`` is identical either way, so
                this changes cost only, never values. ``"embedding"`` and
                ``"logits"`` are absent from the output when True.
            skip_lm_head: Omit ``"logits"`` and skip the LM head projection,
                while still running the decoder tower. The head applies a GELU
                over the full 1 bp embedding, so it allocates another
                ``(batch_size, sequence_length, conv_init_embed_dim)`` tensor
                alongside the one the track heads still need - worth skipping
                whenever the logits are unused, as in tracks inference. Unlike
                ``stop_after_transformer`` this keeps the decoder outputs, which
                the track heads consume.


        Returns:
            outs: dictionary of outputs:
                - ``"embedding"``: the final 1bp-resolution embedding of shape
                    ``(batch_size, sequence_length, conv_init_embed_dim)``
                - ``"before_transformer_embedding"``: the trunk embedding before
                    the transformer of shape
                    ``(batch_size, sequence_length / total_downsampling, embed_dim)``
                - ``"after_transformer_embedding"``: same shape as above, after
                    the transformer and final norm
                - ``"logits"``: the LM head logits of shape
                    ``(batch_size, sequence_length, alphabet_size)``
                - additional intermediate outputs as configured by the
                    ``conv_layers_to_save`` / ``embeddings_layers_to_save`` /
                    ``deconv_layers_to_save`` config fields.

            When ``stop_after_transformer`` is True, ``"embedding"``,
            ``"logits"`` and the ``embeddings_deconv_*`` entries are omitted.
            When ``skip_lm_head`` is True, only ``"logits"`` is omitted.
        """
        del kwargs  # accepted for forward compatibility, currently unused
        ds = self.config.total_downsampling
        seq_len = tokens.shape[1]
        if seq_len % ds != 0:
            raise ValueError(
                f"Sequence length ({seq_len}) must be divisible by "
                f"total_downsampling={ds} (2^(num_downsamples+1))"
            )

        outs: dict[str, Embedding] = {}

        # 1. Token embedding -> stem
        x = self.embed_layer(tokens)
        x = with_logical_constraint(x, PartitionSpec("data", "sequence"))
        x = self.stem(x)
        x = with_logical_constraint(x, PartitionSpec("data", "sequence"))

        # 2. Encoder: initial max-pool, then N blocks each followed by max-pool
        intermediates: dict[str, jax.Array] = {}
        intermediates["bin_size_1"] = x
        x = pool(x, window=2, reduce="max")

        def apply_encoder_block(block: nnx.Module, x: jax.Array) -> jax.Array:
            return block(x)

        if self.config.use_remat_in_convs:
            apply_encoder_block = nnx.remat(apply_encoder_block)

        for i, block in enumerate(self.encoder_blocks):
            x = apply_encoder_block(block, x)
            x = with_logical_constraint(x, PartitionSpec("data", "sequence"))
            intermediates[f"bin_size_{self._encoder_bin_sizes[i]}"] = x
            if (i + 1) in self.config.conv_layers_to_save:
                outs[f"embeddings_conv_{i + 1}"] = x
            x = pool(x, window=2, reduce="max")

        outs["before_transformer_embedding"] = x

        # 3. Transformer tower
        attention_mask = self._block_attention_mask(tokens)

        def apply_transformer_layer(
            mha_block: GQABlock,
            mlp_block: MLPBlock,
            x: jax.Array,
            attention_mask: jax.Array,
            conditions: list[jax.Array] | None,
            conditions_masks: list[jax.Array] | None,
        ) -> jax.Array:
            x = x + mha_block(
                x,
                attention_mask=attention_mask,
                conditions=conditions,
                conditions_masks=conditions_masks,
            )
            x = x + mlp_block(x, conditions, conditions_masks)
            return x

        if self.config.use_remat_in_transformer:
            apply_transformer_layer = nnx.remat(apply_transformer_layer)

        for i in range(self.config.num_layers):
            x = apply_transformer_layer(
                self.mha_blocks[i],
                self.mlp_blocks[i],
                x,
                attention_mask,
                conditions,
                conditions_masks,
            )
            x = with_logical_constraint(x, PartitionSpec("data", "sequence"))

            if (i + 1) in self.config.embeddings_layers_to_save:
                outs[f"embeddings_{i + 1}"] = x

        x = self.final_norm(x)
        outs["after_transformer_embedding"] = x

        if stop_after_transformer:
            return outs

        # 4. Decoder
        def apply_decoder_block(
            block: UpResBlock,
            x: jax.Array,
            skip: jax.Array,
            conditions: list[jax.Array] | None,
            conditions_masks: list[jax.Array] | None,
        ) -> jax.Array:
            return block(x, skip, conditions, conditions_masks)

        if self.config.use_remat_in_convs:
            apply_decoder_block = nnx.remat(apply_decoder_block)

        for i, block in enumerate(self.decoder_blocks):
            bin_size = self._decoder_bin_sizes[i]
            skip = intermediates[f"bin_size_{bin_size}"]
            x = apply_decoder_block(block, x, skip, conditions, conditions_masks)
            x = with_logical_constraint(x, PartitionSpec("data", "sequence"))
            if (i + 1) in self.config.deconv_layers_to_save:
                outs[f"embeddings_deconv_{i + 1}"] = x

        # 5. LM head + outputs
        if not skip_lm_head:
            outs["logits"] = self.lm_head(x)
        outs["embedding"] = x
        return outs
