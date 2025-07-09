from typing import Optional

import haiku as hk
import jax
import jax.numpy as jnp


class ConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        layer_norm_axis: int = -1,
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
        self._layer_norm_axis = layer_norm_axis

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv = hk.Conv1D(
            output_channels=self._dim if self._dim_out is None else self._dim_out,
            kernel_shape=self._kernel_size,
            padding="SAME",
            data_format="NWC",
        )

        layer_norm = hk.LayerNorm(
            axis=self._layer_norm_axis,
            create_scale=True,
            create_offset=True,
            eps=1e-5,
            param_axis=self._layer_norm_axis,
        )

        x = layer_norm(x)
        x = conv(x)
        x = jax.nn.gelu(x)
        return x


class ResidualConvBlock(hk.Module):
    """
    Conv Block with Residual connection.
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

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv_block = ConvBlock(
            dim=self._dim,
            dim_out=self._dim_out,
            kernel_size=self._kernel_size,
        )
        return x + conv_block(x)


class ResidualDeConvBlock(hk.Module):
    """
    Conv Block with Residual connection.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        stride: int = 1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            stride: kernel's stride.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._dim_out = dim_out
        self._kernel_size = kernel_size
        self._stride = stride

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv_block = DeConvBlock(
            dim=self._dim,
            dim_out=self._dim_out,
            kernel_size=self._kernel_size,
            stride=self._stride,
        )
        return x + conv_block(x)


class DeConvBlock(hk.Module):
    """
    Conv Block.
    """

    def __init__(
        self,
        dim: int,
        dim_out: Optional[int] = None,
        kernel_size: int = 1,
        stride: int = 1,
        layer_norm_axis: int = -1,
        name: Optional[str] = None,
    ):
        """
        Args:
            dim: input dimension.
            dim_out: output dimension.
            kernel_size: kernel's size.
            stride: kernel's stride.
            name: model's name.
        """
        super().__init__(name=name)
        self._dim = dim
        self._dim_out = dim_out
        self._kernel_size = kernel_size
        self._stride = stride
        self._layer_norm_axis = layer_norm_axis

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        conv = hk.Conv1DTranspose(
            output_channels=self._dim if self._dim_out is None else self._dim_out,
            kernel_shape=self._kernel_size,
            padding="SAME",
            data_format="NWC",
            stride=self._stride,
        )

        layer_norm = hk.LayerNorm(
            axis=self._layer_norm_axis,
            create_scale=True,
            create_offset=True,
            eps=1e-5,
            param_axis=self._layer_norm_axis,
        )

        x = layer_norm(x)
        x = conv(x)
        x = jax.nn.gelu(x)
        return x
