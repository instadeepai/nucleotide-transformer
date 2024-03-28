from typing import Dict, Optional, Tuple

import haiku as hk
import jax
import jax.numpy as jnp
from haiku import initializers

from nucleotide_transformer.types import SequenceMask
from nucleotide_transformer.utils import get_activation_fn


class DownSample1D(hk.Module):
    """
    1D-UNET downsampling block.
    """

    def __init__(
        self,
        output_channels: int,
        activation_fn: str = "swish",
        num_layers: int = 2,
        name: Optional[str] = None,
    ):
        """
        Args:
            output_channels: number of output channels.
            activation_fn: name of the activation function to use.
                Should be one of "gelu",
                "gelu-no-approx", "relu", "swish", "silu", "sin".
            num_layers: number of convolution layers.
            name: module name.
        """

        super().__init__(name=name)

        self._conv_layers = [
            hk.Conv1D(
                output_channels=output_channels,
                kernel_shape=3,
                stride=1,
                rate=1,
                padding="SAME",
                data_format="NWC",
                name=f"conv{i}",
            )
            for i in range(num_layers)
        ]

        self._avg_pool = hk.AvgPool(
            window_shape=2,
            strides=2,
            padding="SAME",
            channel_axis=-1,
        )

        self._activation_fn = get_activation_fn(activation_name=activation_fn)

    def __call__(self, x: jnp.ndarray) -> Tuple[jnp.ndarray, jnp.ndarray]:
        for _, conv_layer in enumerate(self._conv_layers):
            x = self._activation_fn(conv_layer(x))
        hidden = x
        x = self._avg_pool(hidden)
        return x, hidden


class UpSample1D(hk.Module):
    """
    1D-UNET upsampling block.
    """

    def __init__(
        self,
        output_channels: int,
        activation_fn: str = "swish",
        num_layers: int = 2,
        interpolation_method: str = "nearest",
        name: Optional[str] = None,
    ):
        """
        Args:
            output_channels: number of output channels.
            activation_fn: name of the activation function to use.
                Should be one of "gelu",
                "gelu-no-approx", "relu", "swish", "silu", "sin".
            interpolation_method: Method to be used for upsampling interpolation.
                Should be one of "nearest", "linear", "cubic", "lanczos3", "lanczos5".
            num_layers: number of convolution layers.
            name: module name.
        """
        super().__init__(name=name)

        self._conv_layers = [
            hk.Conv1DTranspose(
                output_channels=output_channels,
                kernel_shape=3,
                stride=1,
                padding="SAME",
                data_format="NWC",
                name=f"conv_transpose{i}",
            )
            for i in range(num_layers)
        ]

        self._interpolation_method = interpolation_method
        self._activation_fn = get_activation_fn(activation_name=activation_fn)

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        for _, conv_layer in enumerate(self._conv_layers):
            x = self._activation_fn(conv_layer(x))
        x = jax.image.resize(
            x,
            shape=(x.shape[0], 2 * x.shape[1], x.shape[2]),
            method=self._interpolation_method,
        )
        return x


class FinalConv1D(hk.Module):
    """
    Final output block of the 1D-UNET.
    """

    def __init__(
        self,
        output_channels: int,
        activation_fn: str = "swish",
        num_layers: int = 2,
        name: Optional[str] = None,
    ):
        """
        Args:
            output_channels: number of output channels.
            activation_fn: name of the activation function to use.
                Should be one of "gelu",
                "gelu-no-approx", "relu", "swish", "silu", "sin".
            num_layers: number of convolution layers.
            name: module name.
        """
        super().__init__(name=name)

        self._conv_layers = [
            hk.Conv1D(
                output_channels=output_channels,
                kernel_shape=3,
                stride=1,
                rate=1,
                padding="SAME",
                data_format="NWC",
                name=f"conv{i}",
            )
            for i in range(num_layers)
        ]

        self._activation_fn = get_activation_fn(activation_name=activation_fn)

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        for i, conv_layer in enumerate(self._conv_layers):
            x = conv_layer(x)
            if i < len(self._conv_layers) - 1:
                x = self._activation_fn(x)
        return x


class UNET1DSegmentationHead(hk.Module):
    """
    1D-UNET based head to be plugged on top of a pretrained model to perform
    semantic segmentation.
    """

    def __init__(
        self,
        num_classes: int,
        output_channels_list: Tuple[int, ...] = (64, 128, 256),
        activation_fn: str = "swish",
        num_conv_layers_per_block: int = 2,
        upsampling_interpolation_method: str = "nearest",
        name: Optional[str] = None,
    ):
        """
        Args:
            num_classes: number of classes to segment
            output_channels_list: list of the number of output channel at each level of
                the UNET
            activation_fn: name of the activation function to use.
                Should be one of "gelu",
                "gelu-no-approx", "relu", "swish", "silu", "sin".
            num_conv_layers_per_block: number of convolution layers per block.
            upsampling_interpolation_method: Method to be used for
                interpolation in upsampling blocks. Should be one of "nearest",
                "linear", "cubic", "lanczos3", "lanczos5".
            name: module name.
        """
        super().__init__(name=name)
        self._num_pooling_layers = len(output_channels_list)
        self._downsample_blocks = [
            DownSample1D(
                output_channels=output_channels,
                activation_fn=activation_fn,
                num_layers=num_conv_layers_per_block,
                name=f"downsample_block_{i}",
            )
            for i, output_channels in enumerate(output_channels_list)
        ]

        self._upsample_blocks = [
            UpSample1D(
                output_channels=output_channels,
                activation_fn=activation_fn,
                num_layers=num_conv_layers_per_block,
                interpolation_method=upsampling_interpolation_method,
                name=f"upsample_block_{i}",
            )
            for i, output_channels in enumerate(reversed(output_channels_list))
        ]

        self._final_block = FinalConv1D(
            activation_fn=activation_fn,
            output_channels=num_classes * 2,
            num_layers=num_conv_layers_per_block,
        )

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:

        if x.shape[1] % 2**self._num_pooling_layers:
            raise ValueError(
                "Input length must be divisible by the 2 to the power of"
                " number of poolign layers."
            )

        hiddens = []
        for downsample_block in self._downsample_blocks:
            x, hidden = downsample_block(x)
            hiddens.append(hidden)

        for upsample_block, hidden in zip(self._upsample_blocks, reversed(hiddens)):
            x = upsample_block(x) + hidden

        x = self._final_block(x)
        return x


class UNetHead(hk.Module):
    """
    Returns a probability between 0 and 1 over a target feature presence
    for each nucleotide in the input sequence. Assumes the sequence has been tokenized
    with non-overlapping 6-mers.
    """

    def __init__(
        self,
        num_features: int,
        embed_dimension: int = 1024,
        num_layers: int = 2,
        name: Optional[str] = None,
    ):
        """
        Args:
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self._num_features = num_features

        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        unet = UNET1DSegmentationHead(
            num_classes=embed_dimension // 2,
            output_channels_list=tuple(
                embed_dimension * (2**i) for i in range(num_layers)
            ),
        )
        fc = hk.Linear(
            6 * 2 * self._num_features, w_init=w_init, b_init=b_init, name="fc"
        )
        self._fc = hk.Sequential([unet, jax.nn.swish, fc])

    def __call__(
        self, x: jnp.ndarray, sequence_mask: SequenceMask
    ) -> Dict[str, jnp.ndarray]:
        """
        Input shape: (batch_size, sequence_length + 1, embed_dim)
        Output_shape: (batch_size, 6 * sequence_length, 2)
        """
        batch_size, seq_len = x.shape[0], x.shape[1] - 1
        logits = self._fc(x[:, 1:])  # remove CLS token
        logits = jnp.reshape(logits, (batch_size, seq_len * 6, self._num_features, 2))
        return {"logits": logits}
