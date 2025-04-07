from typing import Dict, Tuple

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
TransformerOutput: TypeAlias = Dict[str, jnp.ndarray]
MultiOmicsTokens: TypeAlias = Tuple[jnp.ndarray, jnp.ndarray]
