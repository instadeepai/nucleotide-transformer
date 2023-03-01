from typing import Dict

import jax.numpy as jnp
from typing_extensions import TypeAlias

RNGKey: TypeAlias = jnp.ndarray
Embedding: TypeAlias = jnp.ndarray
Tokens: TypeAlias = jnp.ndarray
Labels: TypeAlias = jnp.ndarray
Images: TypeAlias = jnp.ndarray
AttentionMask: TypeAlias = jnp.ndarray
SequenceMask: TypeAlias = jnp.ndarray
AttentionWeights: TypeAlias = jnp.ndarray
TransformerOutput: TypeAlias = Dict[str, jnp.ndarray]  # type: ignore
Metrics: TypeAlias = Dict[str, jnp.ndarray]
