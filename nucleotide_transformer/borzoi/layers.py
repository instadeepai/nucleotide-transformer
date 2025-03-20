from collections.abc import Sequence
from typing import Optional, Tuple, Union

import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


def _prepend_dims(x: np.ndarray, num_dims: int) -> np.ndarray:
    return jnp.reshape(x, tuple([1] * num_dims) + x.shape)


def get_positional_features_central_mask_borzoi(
    positions: jnp.ndarray, feature_size: int, seq_length: int
) -> np.ndarray:
    """Positional features using a central mask (allow only central features)."""
    pow_rate = jnp.exp(jnp.log(seq_length + 1) / feature_size).astype("float32")
    center_widths = jnp.power(pow_rate, jnp.arange(1, feature_size + 1, 1))
    center_widths = center_widths - 1
    center_widths = _prepend_dims(center_widths, positions.ndim)
    outputs = jnp.asarray(
        center_widths > jnp.abs(positions)[..., jnp.newaxis], jnp.float32
    )

    return outputs


def get_positional_embed_borzoi(seq_len: int, feature_size: int) -> jnp.ndarray:
    """
    Compute positional embedding for Borzoi. Note that it is different than the one
    used in Enformer
    """
    distances = jnp.arange(-seq_len + 1, seq_len)

    num_components = 2

    if (feature_size % num_components) != 0:
        raise ValueError(
            f"feature size is not divisible by number of components ({num_components})"
        )

    num_basis_per_class = feature_size // num_components

    embeddings = []

    embeddings.append(
        get_positional_features_central_mask_borzoi(
            distances, num_basis_per_class, seq_len
        )
    )

    embeddings = jnp.concatenate(embeddings, axis=-1)
    embeddings = jnp.concatenate(
        (embeddings, jnp.sign(distances)[..., None] * embeddings), axis=-1
    )
    return embeddings


class SeparableDepthwiseConv1D(hk.Module):
    """Separable 2-D Depthwise Convolution Module."""

    def __init__(
        self,
        channel_multiplier: int,
        kernel_shape: int,
        stride: int = 1,
        padding: Union[str, Sequence[Tuple[int, int]]] = "SAME",
        with_bias: bool = True,
        w_init: Optional[hk.initializers.Initializer] = None,
        b_init: Optional[hk.initializers.Initializer] = None,
        data_format: str = "NWC",
        name: Optional[str] = None,
    ):
        """Construct a Separable 2D Depthwise Convolution module.

        Args:
          channel_multiplier: Multiplicity of output channels. To keep the number of
            output channels the same as the number of input channels, set 1.
          kernel_shape: The shape of the kernel. Either an integer or a sequence of
            length ``num_spatial_dims``.
          stride: Optional stride for the kernel. Either an integer or a sequence of
            length ``num_spatial_dims``. Defaults to 1.
          padding: Optional padding algorithm. Either ``VALID``, ``SAME`` or a
            sequence of ``before, after`` pairs. Defaults to ``SAME``. See:
            https://www.tensorflow.org/xla/operation_semantics#conv_convolution.
          with_bias: Whether to add a bias. By default, true.
          w_init: Optional weight initialization. By default, truncated normal.
          b_init: Optional bias initialization. By default, zeros.
          data_format: The data format of the input.  Can be either
            ``channels_first``, ``channels_last``, ``N...C`` or ``NC...``. By
            default, ``channels_last``.
          name: The name of the module.
        """
        super().__init__(name=name)
        self._conv1 = hk.DepthwiseConv1D(
            channel_multiplier=channel_multiplier,
            kernel_shape=kernel_shape,
            stride=stride,
            padding=padding,
            with_bias=False,
            w_init=w_init,
            b_init=b_init,
            data_format=data_format,
        )

        self._conv2 = hk.Conv1D(
            output_channels=1536,
            kernel_shape=1,
            stride=1,
            padding=padding,
            with_bias=with_bias,
            w_init=w_init,
            b_init=b_init,
            data_format=data_format,
        )

    def __call__(self, inputs: jax.Array) -> jax.Array:

        x = self._conv1(inputs)
        x = self._conv2(x)

        return x


class ConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._dim_out = dim_out
        self._kernel_size = kernel_size

    def __call__(self, x: jnp.ndarray, is_training: bool = False) -> jnp.ndarray:
        batch_norm = hk.BatchNorm(
            create_scale=True,
            create_offset=True,
            decay_rate=0.99,
            eps=0.001,
            data_format="NC...",
        )
        conv = hk.Conv1D(
            output_channels=self._dim if self._dim_out is None else self._dim_out,
            kernel_shape=self._kernel_size,
            padding=(self._kernel_size // 2, self._kernel_size // 2),
            data_format="NCW",
        )

        x = batch_norm(x, is_training=is_training)
        x = jax.nn.gelu(x)
        x = conv(x)
        return x


class UNetConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim

    def __call__(
        self, x: jnp.ndarray, unet_repr: jnp.ndarray, is_training: bool = False
    ) -> jnp.ndarray:
        batch_norm = hk.BatchNorm(
            create_scale=True,
            create_offset=True,
            decay_rate=0.99,
            eps=0.001,
            data_format="NC...",
            name="batch_norm",
        )
        batch_norm_unet_repr = hk.BatchNorm(
            create_scale=True,
            create_offset=True,
            decay_rate=0.99,
            eps=0.001,
            data_format="NC...",
            name="batch_norm_unet_repr",
        )

        linear = hk.Linear(output_size=self._dim, with_bias=True, name="linear")
        linear_unet = hk.Linear(
            output_size=self._dim, with_bias=True, name="linear_unet"
        )
        separable_conv = SeparableDepthwiseConv1D(
            channel_multiplier=1, kernel_shape=3, with_bias=False, data_format="NCW"
        )

        x = batch_norm(x, is_training=is_training)
        unet_repr = batch_norm_unet_repr(unet_repr, is_training=is_training)

        x = jax.nn.gelu(x)
        unet_repr = jax.nn.gelu(unet_repr)

        x = jnp.transpose(x, axes=(0, 2, 1))
        unet_repr = jnp.transpose(unet_repr, axes=(0, 2, 1))

        x = linear(x)
        unet_repr = linear_unet(unet_repr)

        x = jnp.transpose(x, axes=(0, 2, 1))
        unet_repr = jnp.transpose(unet_repr, axes=(0, 2, 1))

        # Upsample x
        x = jax.image.resize(
            x,
            shape=(x.shape[0], x.shape[1], x.shape[2] * 2),
            method="nearest",
        )

        x = x + unet_repr

        x = separable_conv(x)

        return x
