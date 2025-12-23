import jax
import jax.nn
import jax.numpy as jnp
from flax import nnx
from flax.typing import Dtype

from nucleotide_transformer_v3.layers import (
    ConvBlock,
    ConvTowerBlock,
    DeConvBlock,
    DeconvTowerBlock,
    DeConvUpsampleType,
    RotaryEmbeddingConfig,
    SelfAttentionBlock,
)
from nucleotide_transformer_v3.types import AttentionMask, Embedding, Tokens, TransformerOutput


class AdaptiveLayerNorm(nnx.LayerNorm):
    """
    Adaptive layer norm. Rescale a layer norm with one or several conditions.
    """

    def __init__(
        self,
        num_features: int,
        conditions_dims: list[int],
        epsilon: float = 1e-5,
        *,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            num_features: input dimension.
            conditions_dims: conditions dimensions.
            epsilon: epsilon value.
            rngs: random number generator.
        """
        super().__init__(
            num_features=num_features,
            epsilon=epsilon,
            dtype=ln_dtype,
            param_dtype=ln_param_dtype,
            rngs=rngs,
        )
        self.modulation_layers = [
            nnx.Linear(
                in_features=cond_dim,
                out_features=2 * num_features,
                kernel_init=nnx.initializers.constant(0.0),
                bias_init=nnx.initializers.constant(0.0),
                dtype=modulation_dtype,
                param_dtype=modulation_param_dtype,
                rngs=rngs,
            )
            for cond_dim in conditions_dims
        ]
        self._num_conditions = len(conditions_dims)
        self._dim = num_features

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = super().__call__(x=x)

        assert len(conditions) == self._num_conditions, "num conditions mismatch"

        if conditions_masks is None:
            conditions_masks = [
                jnp.ones(shape=(x.shape[0]), dtype=x.dtype) for _ in conditions
            ]

        assert (
            len(conditions_masks) == self._num_conditions
        ), "num conditions masks mismatch"

        scale = 1.0
        shift = 0.0
        for i, x_c_i in enumerate(conditions):
            tmp = self.modulation_layers[i](x_c_i)[:, None, :]
            shift_i, scale_i = jnp.split(tmp, indices_or_sections=2, axis=-1)
            condition_mask = conditions_masks[i]
            shift_i = jnp.where(condition_mask[:, None, None], shift_i, 0.0)
            scale_i = jnp.where(condition_mask[:, None, None], scale_i, 0.0)
            scale = scale * (1.0 + scale_i)
            shift = shift + shift_i
        x = x * scale + shift
        return x


class AdaptiveConvBlock(ConvBlock):
    """
    Adaptive version of the Conv Block.
    """

    def __init__(
        self,
        dim: int,
        conditions_dims: list[int],
        dim_out: int | None = None,
        kernel_size: int = 1,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            dim: input dimension.
            conditions_dims: conditions dimensions.
            dim_out: output dimension.
            kernel_size: kernel's size.
            rngs: random number generator.
        """
        super().__init__(
            dim=dim,
            dim_out=dim_out,
            kernel_size=kernel_size,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
        )
        self.layer_norm = AdaptiveLayerNorm(
            num_features=dim,
            conditions_dims=conditions_dims,
            epsilon=1e-5,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = self.layer_norm(x, conditions, conditions_masks)
        x = self.conv(x)
        x = jax.nn.gelu(x)
        return x


class AdaptiveResidualConvBlock(nnx.Module):
    """
    Adaptive version of the Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        conditions_dims: list[int],
        dim_out: int | None = None,
        kernel_size: int = 1,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            dim: input dimension.
            conditions_dims: conditions dimensions.
            dim_out: output dimension.
            kernel_size: kernel's size.
            rngs: random number generator.
        """
        self.conv_block = AdaptiveConvBlock(
            dim=dim,
            conditions_dims=conditions_dims,
            dim_out=dim_out,
            kernel_size=kernel_size,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.modulation_layers = [
            nnx.Linear(
                in_features=cond_dim,
                out_features=dim,
                kernel_init=nnx.initializers.constant(0.0),
                bias_init=nnx.initializers.constant(0.0),
                dtype=modulation_dtype,
                param_dtype=modulation_param_dtype,
                rngs=rngs,
            )
            for cond_dim in conditions_dims
        ]

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        assert len(conditions) == len(self.modulation_layers), "num conditions mismatch"

        if conditions_masks is None:
            conditions_masks = [
                jnp.ones(shape=(x.shape[0]), dtype=x.dtype) for _ in conditions
            ]
        assert len(conditions_masks) == len(conditions), "num conditions masks mismatch"

        gate = 1.0
        for i, condition_i in enumerate(conditions):
            gate_i = self.modulation_layers[i](condition_i)[:, None, :]
            gate_i = jnp.where(conditions_masks[i][:, None, None], gate_i, 0.0)
            gate = gate * (1.0 + gate_i)
        return x + gate * self.conv_block(x, conditions, conditions_masks)


class AdaptiveDeConvBlock(DeConvBlock):
    """
    Adaptive version of the Conv Block.
    """

    def __init__(
        self,
        dim: int,
        conditions_dims: list[int],
        dim_out: int | None = None,
        kernel_size: int = 1,
        upsample: DeConvUpsampleType | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            dim: input dimension.
            conditions_dims: conditions dimensions.
            dim_out: output dimension.
            kernel_size: kernel's size.
            upsample: Upsampling type.
            rngs: random number generator.
        """
        super().__init__(
            dim=dim,
            dim_out=dim_out,
            kernel_size=kernel_size,
            upsample=upsample,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
        )
        self.layer_norm = AdaptiveLayerNorm(
            num_features=dim,
            epsilon=1e-5,
            conditions_dims=conditions_dims,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = self.layer_norm(x, conditions, conditions_masks)
        x = self.conv(x)
        x = jax.nn.gelu(x)
        return x


class AdaptiveResidualDeConvBlock(nnx.Module):
    """
    Adaptive version of the Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        conditions_dims: list[int],
        dim_out: int | None = None,
        kernel_size: int = 1,
        upsample: DeConvUpsampleType | None = None,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        """
        Args:
            dim: input dimension.
            conditions_dims: conditions dimensions.
            dim_out: output dimension.
            kernel_size: kernel's size.
            upsample: Upsampling type.
            rngs: random number generator.
        """
        self.conv_block = AdaptiveDeConvBlock(
            dim=dim,
            conditions_dims=conditions_dims,
            dim_out=dim_out,
            kernel_size=kernel_size,
            upsample=upsample,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.modulation_layers = [
            nnx.Linear(
                in_features=cond_dim,
                out_features=dim,
                kernel_init=nnx.initializers.constant(0.0),
                bias_init=nnx.initializers.constant(0.0),
                dtype=modulation_dtype,
                param_dtype=modulation_param_dtype,
                rngs=rngs,
            )
            for cond_dim in conditions_dims
        ]

    def __call__(
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        assert len(conditions) == len(self.modulation_layers), "num conditions mismatch"
        if conditions_masks is None:
            conditions_masks = [
                jnp.ones(shape=(x.shape[0]), dtype=x.dtype) for _ in conditions
            ]
        assert len(conditions_masks) == len(conditions), "num conditions masks mismatch"
        gate = 1.0
        for i, condition_i in enumerate(conditions):
            gate_i = self.modulation_layers[i](condition_i)[:, None, :]
            gate_i = jnp.where(conditions_masks[i][:, None, None], gate_i, 0.0)
            gate = gate * (1.0 + gate_i)
        return x + gate * self.conv_block(x, conditions, conditions_masks)


class AdaptiveSelfAttentionBlock(SelfAttentionBlock):
    """
    Adaptive version of the Attention block made of self-attention.
    """

    def __init__(
        self,
        num_heads: int,
        embed_dim: int,
        ffn_embed_dim: int,
        conditions_dims: list[int],
        key_size: int | None = None,
        add_bias_kv: bool = False,
        add_bias_fnn: bool = True,
        ffn_activation_name: str = "gelu-no-approx",
        use_glu_in_ffn: bool = False,
        layer_norm_eps: float = 1e-5,  # this is the default haiku value
        pre_layer_norm: bool = True,
        rotary_embedding_config: RotaryEmbeddingConfig | None = None,
        *,
        ffn_dtype: Dtype | None = None,
        ffn_param_dtype: Dtype = jnp.float32,
        mha_dtype: Dtype | None = None,
        mha_param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        super().__init__(
            num_heads=num_heads,
            embed_dim=embed_dim,
            ffn_embed_dim=ffn_embed_dim,
            key_size=key_size,
            add_bias_kv=add_bias_kv,
            add_bias_fnn=add_bias_fnn,
            ffn_activation_name=ffn_activation_name,
            use_glu_in_ffn=use_glu_in_ffn,
            layer_norm_eps=layer_norm_eps,
            pre_layer_norm=pre_layer_norm,
            rotary_embedding_config=rotary_embedding_config,
            ffn_dtype=ffn_dtype,
            ffn_param_dtype=ffn_param_dtype,
            mha_dtype=mha_dtype,
            mha_param_dtype=mha_param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            rngs=rngs,
        )

        self.self_attention_layer_norm = AdaptiveLayerNorm(
            num_features=embed_dim,
            conditions_dims=conditions_dims,
            epsilon=layer_norm_eps,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.final_layer_norm = AdaptiveLayerNorm(
            num_features=embed_dim,
            conditions_dims=conditions_dims,
            epsilon=layer_norm_eps,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.modulation_layers = [
            nnx.Linear(
                in_features=cond_dim,
                out_features=embed_dim,
                kernel_init=nnx.initializers.constant(0.0),
                bias_init=nnx.initializers.constant(0.0),
                dtype=modulation_dtype,
                param_dtype=modulation_param_dtype,
                rngs=rngs,
            )
            for cond_dim in conditions_dims
        ]

    def mlp(
        self,
        embed: Embedding,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> Embedding:
        if self._pre_layer_norm:
            x = self.final_layer_norm(embed, conditions, conditions_masks)
        else:
            x = embed

        if self._use_glu_in_fnn:
            x1, x2 = jnp.split(self.fc1(x), indices_or_sections=2, axis=-1)
            x = self._ffn_activation_fn(x1) * x2
        else:
            x = self._ffn_activation_fn(self.fc1(x))

        x = self.fc2(x)

        if not self._pre_layer_norm:
            x = self.final_layer_norm(x + embed, conditions, conditions_masks)

        return x

    def __call__(  # type: ignore
        self,
        x: Tokens,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        attention_mask: AttentionMask | None = None,
        attention_weight_bias: jax.Array | None = None,
    ) -> TransformerOutput:
        """
        Computes the output of the attention layer.

        Args:
            x: Input token embeddings of shape (batch_size,seq_len,embed_dim).
            conditions: List of condition embeddings of shape
                (batch_size, cond_embed_dim).
            conditions_masks: List of condition masks of shape (batch_size,).
            attention_mask: Attention mask of shape (batch_size, 1,seq_len, seq_len).

        Returns:
            A dictionary containing the output embeddings and the attention weights.
        """
        if conditions_masks is None:
            conditions_masks = [
                jnp.ones(shape=(x.shape[0]), dtype=x.dtype) for _ in conditions
            ]
        assert len(conditions_masks) == len(conditions), "num conditions masks mismatch"

        # Self-Attention
        res = x
        if self._pre_layer_norm:
            x = self.self_attention_layer_norm(x, conditions, conditions_masks)

        output = self.self_attention(
            x=x,
            attention_mask=attention_mask,
            attention_weight_bias=attention_weight_bias,
        )

        if not self._pre_layer_norm:
            output["embeddings"] = self.self_attention_layer_norm(
                output["embeddings"] + res, conditions, conditions_masks
            )
            x = output["embeddings"]
        else:
            x = output["embeddings"]
            x = res + x

        # MLP
        if not self._pre_layer_norm:
            x = self.mlp(x, conditions, conditions_masks)
        else:
            assert len(conditions) == len(
                self.modulation_layers
            ), "num conditions mismatch"
            gate = 1.0
            for i, x_c_i in enumerate(conditions):
                gate_i = self.modulation_layers[i](x_c_i)[:, None, :]
                gate_i = jnp.where(conditions_masks[i][:, None, None], gate_i, 0.0)
                gate = gate * (1.0 + gate_i)
            x = x + gate * self.mlp(x, conditions, conditions_masks)

        output["embeddings"] = x
        return output  # type: ignore


class ConditionedConvTowerBlock(ConvTowerBlock):
    def __init__(
        self,
        dim_in: int,
        dim_out: int,
        conditions_dims: list[int],
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        super().__init__(dim_in=dim_in, dim_out=dim_out, rngs=rngs)
        self.conv = AdaptiveConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.res_conv = AdaptiveResidualConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = self.conv(x, conditions, conditions_masks)
        x = self.res_conv(x, conditions, conditions_masks)
        x = self.avg_pool(x)
        return x


class ConditionedDeConvTowerBlock(DeconvTowerBlock):
    def __init__(
        self,
        dim_in: int,
        dim_out: int,
        conditions_dims: list[int],
        upsample_type: DeConvUpsampleType,
        *,
        dtype: Dtype | None = None,
        param_dtype: Dtype = jnp.float32,
        ln_dtype: Dtype | None = None,
        ln_param_dtype: Dtype = jnp.float32,
        modulation_dtype: Dtype | None = None,
        modulation_param_dtype: Dtype = jnp.float32,
        rngs: nnx.Rngs,
    ):
        super().__init__(
            dim_in=dim_in, dim_out=dim_out, upsample_type=upsample_type, rngs=rngs
        )
        self.conv = AdaptiveDeConvBlock(
            dim=dim_in,
            dim_out=dim_out,
            kernel_size=5,
            upsample=upsample_type,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )
        self.res_conv = AdaptiveResidualDeConvBlock(
            dim=dim_out,
            dim_out=dim_out,
            kernel_size=1,
            upsample=None,
            conditions_dims=conditions_dims,
            dtype=dtype,
            param_dtype=param_dtype,
            ln_dtype=ln_dtype,
            ln_param_dtype=ln_param_dtype,
            modulation_dtype=modulation_dtype,
            modulation_param_dtype=modulation_param_dtype,
            rngs=rngs,
        )

    def __call__(  # type: ignore
        self,
        x: jax.Array,
        conditions: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
    ) -> jax.Array:
        x = self.conv(x, conditions, conditions_masks)
        x = self.res_conv(x, conditions, conditions_masks)
        return x

