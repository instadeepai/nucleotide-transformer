import copy
from typing import Callable, Optional

import haiku as hk
import jax
import jax.numpy as jnp
import jmp

from nucleotide_transformer.layers import RotaryEmbeddingConfig, SelfAttentionBlock
from nucleotide_transformer.mojo.config import MOJOConfig
from nucleotide_transformer.mojo.layers import (
    ConvBlock,
    DeConvBlock,
    ResidualConvBlock,
    ResidualDeConvBlock,
)
from nucleotide_transformer.types import AttentionMask, Embedding, Tokens


class MOJO(hk.Module):
    def __init__(
        self,
        config: MOJOConfig,
        name: Optional[str] = None,
    ):
        super().__init__(name=name)
        self._config = config
        self._embedding_layers = {}
        self._embedding_layers = {
            omic: hk.Embed(
                self._config.alphabet_size[omic],
                self._config.token_embed_dim,
                name=f"{omic}_embedding",
            )
            for omic in self._config.alphabet_size
        }

        self._gene_embedding_layer = hk.Embed(
            self._config.fixed_sequence_length,
            self._config.init_gene_embed_dim,
            name="gene_embedding",
        )
        self._fc_gene_embedding = hk.Linear(self._config.token_embed_dim)

        self._rotary_embedding_config = RotaryEmbeddingConfig(rescaling_factor=None)

    @hk.transparent
    def stem(self, x: jnp.ndarray) -> jnp.ndarray:
        with hk.experimental.name_scope("stem"):
            conv = hk.Conv1D(
                output_channels=self._config.conv_init_embed_dim,
                kernel_shape=self._config.stem_kernel_shape,
                padding="SAME",
                data_format="NWC",
            )
            if self._config.use_remat_in_convs:
                conv = hk.remat(conv)

        x = conv(x)
        x = jax.nn.gelu(x)

        return x

    @hk.transparent
    def conv_tower(self, x: jnp.ndarray) -> tuple[jnp.ndarray, list[jnp.ndarray]]:
        filter_list = copy.deepcopy(self._config.filter_list)
        residuals = []

        for i, (dim_in, dim_out) in enumerate(zip(filter_list[:-1], filter_list[1:])):
            with hk.experimental.name_scope(f"conv_block_{i}"):
                conv, res_conv = self._conv_block(dim_in, dim_out)
                avg_pool = hk.AvgPool(window_shape=2, strides=2, padding="SAME")
                if self._config.use_remat_in_convs:
                    conv = hk.remat(conv)
                    res_conv = hk.remat(res_conv)
                    avg_pool = hk.remat(avg_pool)

            residuals.append(x)
            x = conv(x)
            x = res_conv(x)
            x = avg_pool(x)

        return x, residuals

    @hk.transparent
    def deconv_tower(self, x: jnp.ndarray, residuals: list[jnp.ndarray]) -> jnp.ndarray:
        filter_list = copy.deepcopy(self._config.filter_list)
        filter_list.reverse()
        residuals_generator = reversed(residuals)

        for i, (dim_in, dim_out) in enumerate(zip(filter_list[:-1], filter_list[1:])):
            with hk.experimental.name_scope(f"deconv_block_{i}"):
                conv, res_conv = self._deconv_block(dim_in, dim_out)
                if self._config.use_remat_in_convs:
                    conv = hk.remat(conv)
                    res_conv = hk.remat(res_conv)

            x = conv(x)
            x = res_conv(x)

            if self._config.use_skip_connection:
                residuals = next(residuals_generator)
                x = x + residuals

        return x

    @hk.transparent
    def transformer_tower(
        self,
        x: Embedding,
        outs: dict[str, Embedding],
        attention_mask: Optional[AttentionMask] = None,
    ) -> tuple[Embedding, dict[str, Embedding]]:

        layers: list[Callable] = [
            self._attention_block(layer_idx)
            for layer_idx in range(self._config.num_layers)
        ]

        if self._config.use_remat_in_transformer:
            layers = [hk.remat(layer) for layer in layers]

        for layer_idx, layer in enumerate(layers):
            output = layer(
                x=x, attention_mask=attention_mask, attention_weight_bias=None
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
    def omic_lm_heads(self, x: jnp.ndarray, omic: str) -> jnp.ndarray:
        x = jax.nn.gelu(x)
        for i in range(self._config.num_hidden_layers_head):
            name = f"{omic}_head_linear_{i}"
            x = hk.Linear(self._config.embed_dim, name=name)(x)
            x = jax.nn.gelu(x)
        name = f"{omic}_head_linear_final"
        head = hk.Linear(self._config.alphabet_size[omic], name=name)
        return head(x)

    def get_embeddings(
        self,
        tokens: dict[str, Tokens],
        attention_masks: Optional[dict[str, Optional[AttentionMask]]],
    ) -> dict[str, Embedding]:
        omic_embeddings = {}
        for omic, omic_tokens in tokens.items():
            omic_embeddings[omic] = self._embedding_layers[omic](omic_tokens)

        return omic_embeddings

    @hk.transparent
    def _conv_block(self, dim_in: int, dim_out: int) -> tuple[hk.Module, hk.Module]:
        conv = ConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
        )
        res_conv = ResidualConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
        )
        return conv, res_conv

    @hk.transparent
    def _deconv_block(self, dim_in: int, dim_out: int) -> tuple[hk.Module, hk.Module]:
        conv = DeConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
            stride=2,
        )
        res_conv = ResidualDeConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
        )
        return conv, res_conv

    @hk.transparent
    def _attention_block(self, layer_idx: int) -> SelfAttentionBlock:
        return SelfAttentionBlock(  # type: ignore
            num_heads=self._config.num_attention_heads,
            embed_dim=self._config.embed_dim,
            key_size=self._config.key_size,
            ffn_embed_dim=self._config.ffn_embed_dim,
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            rotary_embedding_config=self._rotary_embedding_config,
            layer_norm_eps=self._config.layer_norm_eps,
            pre_layer_norm=True,
            name=f"attention_layer_{layer_idx}",
        )

    def __call__(
        self,
        tokens_dict: dict[str, Tokens],
        attention_masks: Optional[dict[str, Optional[AttentionMask]]] = None,
    ) -> dict[str, Embedding]:
        outs: dict[str, dict[str, jnp.ndarray] | jnp.ndarray | list[jnp.ndarray]] = {}

        embeddings = self.get_embeddings(tokens_dict, attention_masks)
        outs["omic_embeddings"] = embeddings

        x = jnp.sum(jnp.array(list(embeddings.values())), axis=0)
        outs["embeddings"] = x

        if self._config.use_gene_embedding:
            gene_embedding = self._gene_embedding_layer(
                jnp.arange(self._config.fixed_sequence_length)
            )
            if self._config.project_gene_embedding:
                gene_embedding = self._fc_gene_embedding(gene_embedding)
            x = x + gene_embedding
            outs["embeddings_with_gene_embedding"] = x

        x = self.stem(x)
        outs["stem"] = x

        x, residuals = self.conv_tower(x)
        outs["conv_tower"] = x
        outs["conv_tower_residuals"] = residuals

        x, outs = self.transformer_tower(x, outs=outs, attention_mask=None)
        outs["after_transformer_embedding"] = x

        x = self.deconv_tower(x, residuals)
        outs["deconv_tower"] = x

        outs["logits"] = {
            omic: self.omic_lm_heads(x, omic) for omic in self._config.alphabet_size
        }

        return outs


def build_mojo_fn(
    model_config: MOJOConfig,
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
    ), f"provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    policy = jmp.Policy(
        compute_dtype=compute_dtype, param_dtype=param_dtype, output_dtype=output_dtype
    )
    hk.mixed_precision.set_policy(MOJO, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32, param_dtype=param_dtype, output_dtype=compute_dtype
    )
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)

    def multiomics_lm_fn(
        tokens_dict: dict[str, Tokens],
        attention_masks: Optional[dict[str, Optional[AttentionMask]]] = None,
    ) -> dict[str, jnp.ndarray]:
        model = MOJO(model_config, name=model_name)
        return model(tokens_dict, attention_masks)

    return multiomics_lm_fn
