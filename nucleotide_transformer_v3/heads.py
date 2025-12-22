"""Head modules for Nucleotide Transformer v3."""

import jax
import jax.numpy as jnp
from flax import nnx


class LinearHead(nnx.Module):
    """A linear head that predicts one scalar value per track."""

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        *,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        self.layer_norm = nnx.LayerNorm(embed_dim, rngs=rngs)
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels,
            rngs=rngs,
        )

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.layer_norm(x)
        x = self.head(x)
        x = jax.nn.softplus(x)
        return x


class ClassificationHead(nnx.Module):
    """A linear head that predicts num_classes values per label."""

    def __init__(
        self,
        embed_dim: int,
        num_labels: int,
        num_classes: int,
        *,
        rngs: nnx.Rngs,
    ):
        super().__init__()
        self.layer_norm = nnx.LayerNorm(embed_dim, rngs=rngs)
        self.head = nnx.Linear(
            in_features=embed_dim,
            out_features=num_labels * num_classes,
            rngs=rngs,
        )
        self.num_labels = num_labels
        self.num_classes = num_classes

    def __call__(self, x: jax.Array) -> jax.Array:
        x = self.layer_norm(x)
        x = self.head(x)
        x = jnp.reshape(x, x.shape[:-1] + (self.num_labels, self.num_classes))
        return x

