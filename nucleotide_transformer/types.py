from typing import Dict

import jax.numpy as jnp
from typing_extensions import TypeAlias

Embedding: TypeAlias = jnp.ndarray
Tokens: TypeAlias = jnp.ndarray
AttentionMask: TypeAlias = jnp.ndarray
TransformerOutput: TypeAlias = Dict[str, jnp.ndarray]  # type: ignore
