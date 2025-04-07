from dataclasses import dataclass
from typing import Optional, Tuple

import haiku as hk
import jax.numpy as jnp
import numpy as np

UPPER_FREQ = 10000


@dataclass
class RotaryEmbeddingConfig:
    """
    Parameters to initialize the RotaryEmbedding layer. The rescaling factor allows
    to adapt the rotary embeddings to larger lengths than what was used for training.
    One of this strategy is presented in the Yarn paper: https://arxiv.org/pdf/2309.00071.pdf. # noqa

    Args:

    """

    rescaling_factor: Optional[float]


class RotaryEmbedding(hk.Module):
    """
    Rotary Positional Embedding inspired by RoFormer:
    https://arxiv.org/abs/2104.09864
    https://github.com/ZhuiyiTechnology/roformer .
    """

    def __init__(
        self,
        key_size: int,
        rotary_embedding_config: RotaryEmbeddingConfig,
        name: Optional[str] = None,
    ):
        """
        Args:
            key_size: Dimension of one head.
            rotary_embedding_config: Configuration to specify hyperparameters for
                RotaryEmbeddig layer
                (see RoFormer https://arxiv.org/pdf/2104.09864.pdf). It contains
                the hyperparameters specifying the type of rotary embedding applied.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)

        # Extract argument from the config
        rescaling_factor = rotary_embedding_config.rescaling_factor

        if rescaling_factor is None:
            self._inv_freq = 1.0 / (
                UPPER_FREQ ** (np.arange(0, key_size, 2) / key_size)
            )
        else:
            updated_base = UPPER_FREQ * (
                rescaling_factor ** (key_size / (key_size - 2))
            )
            self._inv_freq = 1.0 / (
                updated_base ** (np.arange(0, key_size, 2) / key_size)
            )

    def _compute_cos_sin_tables(
        self,
        heads: jnp.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the cosinus and sinus for rotation.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
            key_size).

        Returns:
            Cosinus positional embedding of shape (1, seq_len, 1,
                key_size/2).
            Sinus positional embedding of shape (1, seq_len, 1,
                key_size/2).
        """
        seq_len = heads.shape[1]

        self._seq_len_cached = seq_len
        t = np.arange(seq_len)
        freqs = np.einsum("i,j->ij", t, self._inv_freq)

        # Compute cos and cast is as (1, seq_len, 1, key_size/2) to be applied to
        # queries of shape (batch_size, seq_len, num_heads, key_size/2)
        cos_cached = np.cos(freqs)[None, :, None, :]
        sin_cached = np.sin(freqs)[None, :, None, :]

        return cos_cached, sin_cached

    def _apply_rotary_pos_emb(
        self, heads: jnp.ndarray, cos: np.ndarray, sin: np.ndarray
    ) -> jnp.ndarray:
        """
        Applies the rotary positional embedding to the heads.

        Args:
            heads: Query or key heads of shape (batch_size, seq_len, num_heads,
                key_size).
            cos: Cosinus values.
            sin: Sinus values.

        Returns:
            Embedded heads of shape (batch_size, seq_len, num_heads,
                key_size).
        """

        # Rotate x
        x_first, x_second = (
            heads[..., : heads.shape[-1] // 2],
            heads[..., heads.shape[-1] // 2 :],
        )
        first_part = x_first * cos - x_second * sin
        second_part = x_second * cos + x_first * sin

        return jnp.concatenate((first_part, second_part), axis=-1, dtype=heads.dtype)

    def __call__(
        self, query_heads: jnp.ndarray, key_heads: jnp.ndarray
    ) -> Tuple[jnp.ndarray, jnp.ndarray]:
        """
        Applies rotary embeddings to query_heads and key_heads.

        Args:
            query_heads: Query heads of shape
                (batch_size, seq_len, num_heads, key_size).
            key_heads: Key heads of shape (batch_size, seq_len, num_heads, key_size).

        Returns:
            Embedded query heads.
            Embedded key heads.
        """
        cos, sin = self._compute_cos_sin_tables(query_heads)

        return (
            self._apply_rotary_pos_emb(query_heads, cos, sin),
            self._apply_rotary_pos_emb(key_heads, cos, sin),
        )
