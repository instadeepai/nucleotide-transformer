# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Dict

import jax.numpy as jnp
from typing_extensions import TypeAlias

Embedding: TypeAlias = jnp.ndarray
Tokens: TypeAlias = jnp.ndarray
AttentionMask: TypeAlias = jnp.ndarray
TransformerOutput: TypeAlias = Dict[str, jnp.ndarray]  # type: ignore
