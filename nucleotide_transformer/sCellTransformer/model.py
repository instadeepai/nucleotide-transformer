import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import haiku as hk
import jax
import jax.nn
import jax.numpy as jnp
import jmp
import numpy as np
from einops import rearrange

from nucleotide_transformer.layers import RotaryEmbeddingConfig, SelfAttentionBlock
from nucleotide_transformer.types import (
    AttentionMask,
    Embedding,
    SequenceMask,
    Tokens,
    TransformerOutput,
)

logger = logging.getLogger(__name__)


class ResidualConvBlock(hk.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim_out: int,
        kernel_size: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim_out = dim_out
        self._kernel_size = kernel_size

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv_block = ConvBlock(dim_out=self._dim_out, kernel_size=self._kernel_size)
        y = conv_block(x)

        return x.reshape(y.shape) + y


class ConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim_out: int,
        kernel_size: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)

        self._dim_out = dim_out
        self._kernel_size = kernel_size

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv = hk.Conv1D(
            output_channels=self._dim_out,
            kernel_shape=self._kernel_size,
            padding="SAME",
            data_format="NCW",
        )
        layer_norm = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            eps=1e-5,
            param_axis=-1,
        )
        x = layer_norm(x)
        x = x.reshape((x.shape[0], x.shape[1], -1))
        x = conv(x)
        x = jax.nn.gelu(x)
        return x


class ResidualDeConvBlock(hk.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim_out: int,
        kernel_size: int = 1,
        stride: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim_out: output dimension.
            kernel_size: kernel's size.
            stride: kernel's stride.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim_out = dim_out
        self._kernel_size = kernel_size
        self._stride = stride

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv_block = DeConvBlock(
            dim_out=self._dim_out,
            kernel_size=self._kernel_size,
            stride=self._stride,
        )
        y = conv_block(x)
        return x.reshape(y.shape) + y


class DeConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim_out: int,
        kernel_size: int = 1,
        stride: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim_out: output dimension.
            kernel_size: kernel's size.
            stride: kernel's stride.
            name: model's name.
        """
        super().__init__(name=name)

        self._dim_out = dim_out
        self._kernel_size = kernel_size
        self._stride = stride

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv = hk.Conv1DTranspose(
            output_channels=self._dim_out,
            kernel_shape=self._kernel_size,
            padding="SAME",
            data_format="NCW",
            stride=self._stride,
        )
        layer_norm = hk.LayerNorm(
            axis=-1,
            create_scale=True,
            create_offset=True,
            eps=1e-5,
            param_axis=-1,
        )
        x = layer_norm(x)
        x = x.reshape((x.shape[0], x.shape[1], -1))

        x = conv(x)
        x = jax.nn.gelu(x)
        return x


class SpatialEncoding(hk.Module):
    """
    Spatial coordinates encoding module
    """

    def __init__(
        self,
        embed_dim: int,
        num_scales: int = 10,
        sigma_min: float = 1.0,
        sigma_max: float = 10.0,
        name: str | None = None,
    ):
        """
        Spatial Encoding Layer built on positional
        sinusoidal embeddings. Implemented from
        https://proceedings.mlr.press/v206/klemmer23a/klemmer23a.pdf # noqa: E501
        Args:
            embed_dim: spatial embedding dimensions.
                Should be the same as token embedding
                dimensions.
            num_scales: number of scales for scaling the position
                embedding. helps learn scal-invariant transforms.
            sigma_min: minimum scalling parameter for PE
            sigma_max: maximum scaling parameter for PE
            name: Embedding layer name
        Returns:
            None
        """
        super().__init__(name=name)
        self._num_scales = num_scales
        self._sigma_min = sigma_min
        self._sigma_max = sigma_max
        self._g = self._sigma_max / self._sigma_min
        self._scales = jnp.linspace(self._sigma_min, self._sigma_max, self._num_scales)

        self.fc_layer = hk.Linear(embed_dim, name="spatial_fully_connected")

    def scale_specific_encoder(
        self, coordinates: jnp.ndarray, scale: float
    ) -> jnp.ndarray:
        x, y = coordinates[..., 0], coordinates[..., 1]
        constant = self._sigma_min * self._g ** (scale / (self._num_scales - 1))
        x_transform = jnp.cos(x / constant)
        y_transform = jnp.sin(y / constant)
        transformed_coordinates = jnp.stack([x_transform, y_transform], axis=-1)
        return transformed_coordinates

    def __call__(self, coordinates: jnp.ndarray) -> jnp.ndarray:

        transformed_coordinates = [
            self.scale_specific_encoder(coordinates, scale) for scale in self._scales
        ]
        transformed_coordinates = jnp.concatenate(transformed_coordinates, axis=-1)
        return self.fc_layer(transformed_coordinates)


def upsample_x2(x: jnp.ndarray, method: str) -> jnp.ndarray:
    """
    Upsample the sequence length by a factor of 2.

    Args:
        x: input tensor of shape (batch, dim, seq_len)
        method: method used to upsample the sequence length.
    """

    def upsample_x2_1d(x: jnp.ndarray) -> jnp.ndarray:
        if method == "linear":
            x_int = x + jnp.diff(x, append=x[-1]) / 2
            y = jnp.stack([x, x_int]).reshape(-1, order="F")
        elif method == "nearest":
            y = jnp.repeat(x, 2)
        else:
            raise NotImplementedError(f"Upsampling method {method} not implemented.")

        return y

    return jax.vmap(jax.vmap(upsample_x2_1d))(x)


@dataclass
class sCTConfig:
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
        use_gradient_checkpointing: whether to use gradient checkpointing.
        embeddings_layers_to_save: indices of Transformer layers to save embeddings for.
        attention_maps_to_save: indices of Transformer layers to save attention map for.
        use_spatial_information: Adds a spatial encoding to the input tokens for spatial
            tx data
        num_scales: number of scales for scaling the position embedding. (see spatial
            encoding layer for more details.).
        sigma_min: minimum scalling parameter for PE
        sigma_max: maximum scaling parameter for PE


    """

    alphabet_size: int
    pad_token_id: int
    mask_token_id: int

    num_downsamples: int = 8

    # architecture
    attention_heads: int = 16
    key_size: Optional[int] = None
    token_embed_dim: int = 16
    embed_dim: int = 1024
    ffn_embed_dim: int = 2048
    num_layers: int = 4
    layer_norm_eps: float = 1e-5
    interpolation_method: str = "nearest"

    # bad hack to satisfy cellnt_celltype_annotation.py:312
    max_positions: int = 20480
    num_cells: int = 50
    num_hidden_layers_head: int = 1

    use_skip_connection: bool = True

    # logging
    use_gradient_checkpointing: bool = False

    # return
    embeddings_layers_to_save: Tuple[int, ...] = ()
    attention_maps_to_save: List[Tuple[int, int]] = field(default_factory=list)

    # Spatial info configuration
    use_spatial_information: bool = False
    num_scales: int = 10
    sigma_min: float = 1.0
    sigma_max: float = 10.0

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


class sCT(hk.Module):
    def __init__(
        self,
        config: sCTConfig,
        name: Optional[str] = None,
    ):
        super().__init__(name=name)
        self._config = config
        logger.info(
            f"gradient checkpointing: {self._config.use_gradient_checkpointing}"
        )

        # Default embedding layer used when spatial information is not used. If
        # `use_spatial_information` is set to True, the spatial embedding layer
        # will be used instead.
        self._embed_layer = hk.Embed(
            self._config.alphabet_size, self._config.token_embed_dim
        )
        self._rotary_embedding_config = RotaryEmbeddingConfig(rescaling_factor=None)

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

        filter_list = np.linspace(
            self._config.token_embed_dim,
            self._config.embed_dim,
            self._config.num_downsamples + 1,
        )

        filter_list = np.ceil(filter_list / 32) * 32
        filter_list = filter_list.astype(int).tolist()

        self._filter_list = filter_list

    @hk.transparent
    def stem(self, x: jnp.ndarray) -> jnp.ndarray:
        with hk.experimental.name_scope("stem"):
            conv = hk.Conv1D(
                output_channels=self._config.token_embed_dim,
                kernel_shape=15,
                padding="SAME",
                data_format="NCW",
            )
        x = conv(x)
        x = jax.nn.gelu(x)

        return x

    @hk.transparent
    def conv_tower(self, x: jnp.ndarray) -> Tuple[jnp.ndarray, List[jnp.ndarray]]:

        filter_list = self._filter_list[1:]
        residuals = []
        for i, dim_out in enumerate(filter_list):
            with hk.experimental.name_scope(f"conv_block_{i}"):
                conv = ConvBlock(dim_out=dim_out, kernel_size=5)
                res_conv = ResidualConvBlock(dim_out=dim_out, kernel_size=1)
                avg_pool = hk.AvgPool(window_shape=2, strides=2, padding="SAME")
            residuals.append(x)
            x = x.reshape((x.shape[0], x.shape[1], self._config.num_cells, -1))
            x = conv(x)
            x = x.reshape((x.shape[0], x.shape[1], self._config.num_cells, -1))
            x = res_conv(x)
            x = x.transpose(0, 2, 1)
            x = avg_pool(x)
            x = x.transpose(0, 2, 1)

        return x, residuals

    @hk.transparent
    def deconv_tower(self, x: jnp.ndarray, residuals: List[jnp.ndarray]) -> jnp.ndarray:

        residuals = residuals[::-1]
        filter_list = self._filter_list[::-1]
        filter_list = filter_list[1:]

        for i, (dim_out, res) in enumerate(zip(filter_list, residuals)):
            with hk.experimental.name_scope(f"deconv_block_{i}"):
                conv = DeConvBlock(dim_out=dim_out, kernel_size=5, stride=2)
                res_conv = ResidualDeConvBlock(dim_out=dim_out, kernel_size=1)

            x = x.reshape((x.shape[0], x.shape[1], self._config.num_cells, -1))
            x = conv(x)
            x = x.reshape((x.shape[0], x.shape[1], self._config.num_cells, -1))
            x = res_conv(x)

            if self._config.use_skip_connection:

                x = x + res

        return x

    @hk.transparent
    def lm_head(self, x: jnp.ndarray) -> jnp.ndarray:
        x = jax.nn.gelu(x)
        for _ in range(self._config.num_hidden_layers_head):
            x = hk.Linear(self._config.embed_dim)(x)
            x = jax.nn.gelu(x)
        head = hk.Linear(self._config.alphabet_size)
        return head(x)

    @hk.transparent
    def transformer_tower(
        self,
        x: Embedding,
        outs: Dict[str, Embedding],
        attention_mask: Optional[AttentionMask] = None,
    ) -> Tuple[Embedding, Dict[str, Embedding]]:

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
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            rotary_embedding_config=self._rotary_embedding_config,
            layer_norm_eps=self._config.layer_norm_eps,
            pre_layer_norm=True,
            name=f"attention_layer_{layer_idx}",
        )

    @hk.transparent
    def _spatial_embedding_layer(self, tokens: tuple[Tokens]) -> Embedding:

        tokens, coordinates = tokens  # type: ignore
        batch, seq_len = tokens.shape  # type: ignore
        x = self._embed_layer(tokens)  # type: ignore
        coordinates_embed_layer = SpatialEncoding(
            embed_dim=self._config.token_embed_dim,
            num_scales=self._config.num_scales,
            sigma_min=self._config.sigma_min,
            sigma_max=self._config.sigma_max,
        )
        x_coordinates = coordinates_embed_layer(coordinates)  # type: ignore
        x = rearrange(
            x,
            "batch (cells genes) dim -> batch cells genes dim",
            cells=self._config.num_cells,
            genes=seq_len // self._config.num_cells,
        )
        x = x + x_coordinates[:, :, jnp.newaxis, :]
        x = rearrange(x, "batch cells genes dim -> batch (cells genes) dim")
        return x

    def __call__(self, tokens: Tokens | tuple[Tokens]) -> Dict[str, Embedding]:

        # Prepare outputs dict
        outs: Dict[str, Embedding] = {}

        # Compute embeddings

        if self._config.use_spatial_information:
            x = self._spatial_embedding_layer(tokens)
        else:
            x = self._embed_layer(tokens)

        # Main body

        x = jnp.transpose(x, (0, 2, 1))
        x = self.stem(x)
        x, residuals = self.conv_tower(x)
        x = jnp.transpose(x, (0, 2, 1))
        outs["residuals"] = residuals
        outs["conv_out"] = x

        x, outs = self.transformer_tower(x, outs=outs, attention_mask=None)
        outs["transformer_out"] = x

        x = jnp.transpose(x, (0, 2, 1))
        x = self.deconv_tower(x, residuals)
        outs["deconv_out"] = x
        x = jnp.transpose(x, (0, 2, 1))
        outs["embedding"] = x.reshape(
            (x.shape[0], self._config.num_cells, -1, x.shape[-1])
        )

        x = self.lm_head(x)
        outs["logits"] = x

        return outs  # type: ignore


def build_sct_fn(
    config: sCTConfig,
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    name: Optional[str] = None,
) -> Callable:
    """
    Create the model's forward pass.

    Args:
        config: Configuration data class containing the hyperparameters for the GPT
            forward function.
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too
        output_dtype: the output type of the model. it determines the float precision
            of the gradient when training the model.
            NOTE: when training, the gradient is often accumulated in fp32, therefore
            output_dtype need to be in fp32.
        name: the name of the model. example: gpt_j_decoder.


        # NOTE: in inference, the model could be in fp16 without too much degradation
        # NOTE: on NVIDIA accelerator, XLA inter-device operation ( psum, all_gather,
        etc ... ) are not always implemented for bf16. but on TPU hardware yes

    Returns:
        Enformer model forward function.
    """

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
    hk.mixed_precision.set_policy(sCT, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32, param_dtype=param_dtype, output_dtype=compute_dtype
    )
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)

    def sct_fn(
        tokens: jnp.ndarray, attention_mask: Optional[AttentionMask] = None
    ) -> Dict[str, jnp.ndarray]:
        model = sCT(config, name=name)
        return model(tokens)

    return sct_fn


def build_sct_with_head_fn(
    model_config: sCTConfig,
    head_fn: Callable[
        [], Callable[[jnp.ndarray, SequenceMask], Dict[str, jnp.ndarray]]
    ],
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    model_name: Optional[str] = None,
) -> Callable:
    """
    Creates a forward pass for that sCT and adds the input head.

    Args:
        model_config: Model hyperparameters.
        head_fn: Wrapper initializing a Classification/Regression head. The head cannot
            be passed directly as haiku modules cannot be initialized outside
            hk.transform.
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
        model_name: Optional name of the model.

    Example of the function being used with a classification head:
        The classification head is wrapped inside head_fn because
        haiku modules cannot be instantiated outside hk.transform.
        def head_fn():
            return SimpleClassificationHead(num_classes=num_classes)
        finetune_forward_fn = build_esm_ia3_rescaling_with_head_fn(
            model_config=config, head_fn=head_fn, model_name=model_name,
        )
        finetune_forward_fn = hk.transform(finetune_forward_fn)

    Returns:
        sCT model forward function with indicated head.
    """
    policy = jmp.Policy(
        compute_dtype=compute_dtype, param_dtype=param_dtype, output_dtype=output_dtype
    )
    hk.mixed_precision.set_policy(sCT, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=compute_dtype, param_dtype=param_dtype, output_dtype=output_dtype
    )
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)

    def sct_fn(
        tokens: jnp.ndarray,
        attention_mask: Optional[AttentionMask] = None,
        sequence_mask: Optional[SequenceMask] = None,
    ) -> TransformerOutput:
        """Forward pass."""
        # Run the encoder over the inputs.
        encoder = sCT(config=model_config, name=model_name)
        outs: TransformerOutput = encoder(
            tokens=tokens,
            # attention_mask=attention_mask,
        )
        embeddings = outs["embedding"]

        # Define head.
        head = head_fn()

        if sequence_mask is None:
            sequence_mask = jnp.ones_like(tokens)

        head_outs = head(  # type: ignore[call-arg]
            x=embeddings, sequence_mask=sequence_mask
        )
        outs.update(head_outs)
        return outs

    return sct_fn
