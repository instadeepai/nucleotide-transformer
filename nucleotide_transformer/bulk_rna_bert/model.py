from typing import Callable, Optional

import haiku as hk
import jax.numpy as jnp
import jmp

from nucleotide_transformer.bulk_rna_bert.config import BulkRNABertConfig
from nucleotide_transformer.bulk_rna_bert.layers import SimpleLMHead
from nucleotide_transformer.layers import SelfAttentionBlock
from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)


class BulkRNABert(hk.Module):
    """
    Jax implementation of BulkRNABert model.
    """

    def __init__(
        self,
        config: BulkRNABertConfig,
        name: Optional[str] = None,
    ):
        super().__init__(name=name)
        self._config = config

        self._gene_embedding_layer = hk.Embed(
            self._config.n_genes,
            self._config.init_gene_embed_dim,
            name="gene_embedding",
        )
        self._fc_gene_embedding = hk.Linear(self._config.embed_dim)
        self._expression_embedding_layer = hk.Embed(
            self._config.n_expressions_bins,
            self._config.embed_dim,
            name="expression_embedding",
        )
        self._lm_head = SimpleLMHead(
            embed_dim=self._config.embed_dim,
            alphabet_size=self._config.n_expressions_bins,
        )

    @hk.transparent
    def apply_attention_blocks(
        self,
        x: Embedding,
        outs: dict[str, Embedding],
        attention_mask: Optional[AttentionMask] = None,
    ) -> tuple[Embedding, dict[str, Embedding]]:
        """
        Creates the blocks of attention layers and applies them.

        Args:
            x: Sequence embedding of shape (batch,seq_len,embed_dim).
            outs: A dictionary to carry through the attention layers which stores the
                intermediate sequence embedding and attention maps.
            attention_mask: attention mask of shape (batch_size, 1, seq_len, seq_len).

        Returns:
            Output sequence embedding.
            Dictionary of optional intermediate results (embeddings of the layer and
            attention weights).
        """

        layers: list[Callable] = [
            self._self_attention(layer_idx)
            for layer_idx in range(self._config.num_layers)
        ]

        if self._config.use_gradient_checkpointing:
            # the remat-ed function cannot take control flow arguments
            layers = [hk.remat(layer) for layer in layers]

        for layer_idx, layer in enumerate(layers):
            output = layer(
                x=x,
                attention_mask=attention_mask,
            )
            x = output["embeddings"]

            if (layer_idx + 1) in self._config.embeddings_layers_to_save:
                outs[f"embeddings_{(layer_idx + 1)}"] = output["embeddings"]
            if (layer_idx + 1) in self._config.attention_layers_to_save:
                outs[f"attention_map_layer_{layer_idx + 1}"] = output[
                    "attention_weights"
                ]

        return x, outs

    @hk.transparent
    def _self_attention(self, layer_idx: int) -> SelfAttentionBlock:
        return SelfAttentionBlock(  # type: ignore
            num_heads=self._config.num_attention_heads,
            embed_dim=self._config.embed_dim,
            key_size=self._config.key_size,
            ffn_embed_dim=self._config.ffn_embed_dim,
            name=f"self_attention_block_{layer_idx}",
        )

    def __call__(
        self,
        tokens: Tokens,
        attention_mask: AttentionMask = None,
    ) -> TransformerOutput:
        batch_size, seq_len = tokens.shape
        outs: dict[str, jnp.ndarray] = {}

        x = self._expression_embedding_layer(tokens)
        if self._config.use_gene_embedding:
            gene_embedding = self._gene_embedding_layer(
                jnp.arange(self._config.n_genes)
            )
            if self._config.project_gene_embedding:
                gene_embedding = self._fc_gene_embedding(gene_embedding)
            x = x + gene_embedding

        if attention_mask is None:
            attention_mask = jnp.ones((batch_size, 1, seq_len, seq_len))

        x, outs = self.apply_attention_blocks(
            x=x,
            outs=outs,
            attention_mask=attention_mask,
        )
        lm_head_outs = self._lm_head(x)
        outs["logits"] = lm_head_outs["logits"]
        return outs


def build_bulk_rna_bert_forward_fn(
    model_config: BulkRNABertConfig,
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    model_name: Optional[str] = None,
) -> Callable:
    assert {compute_dtype, param_dtype, output_dtype}.issubset(
        {
            jnp.bfloat16,
            jnp.float32,
            jnp.float16,
        }
    ), f"Please provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    policy = jmp.Policy(
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
    )
    hk.mixed_precision.set_policy(BulkRNABert, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32,
        param_dtype=param_dtype,
        output_dtype=compute_dtype,
    )
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.RMSNorm, norm_policy)

    def forward_fn(
        tokens: jnp.ndarray,
        attention_mask: Optional[jnp.ndarray] = None,
    ):
        """Forward pass"""
        model = BulkRNABert(
            config=model_config,
            name=model_name,
        )
        outs = model(
            tokens=tokens,
            attention_mask=attention_mask,
        )
        return outs

    return forward_fn
