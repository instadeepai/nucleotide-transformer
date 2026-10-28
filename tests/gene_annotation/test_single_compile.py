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

"""The annotation forward is traced and compiled once per process.

Two regressions guarded here. ``model.eval()`` writes each ``RMSBatchNorm`` mode
flag as an uncommitted host array, and ``nnx.jit`` hands it back committed, so
without ``commit_model_state`` the second call recompiled the whole forward. And
the single-strand path inlined the network once per strand, compiling it twice.
"""

import contextlib
from collections.abc import Callable, Iterator

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx

from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
    build_annotation_inference_fn,
    commit_model_state,
)
from nucleotide_transformer_v3.layers import RMSBatchNorm

NUM_CLASSES = 15
VOCAB = 16
BACKEND_COMPILE_EVENT = "/jax/core/compile/backend_compile_duration"
_TRACES: list[int] = []


class _NormedStub(nnx.Module):
    """Single-strand head behind a real RMSBatchNorm, counting its traces."""

    def __init__(self) -> None:
        self.norm = RMSBatchNorm(NUM_CLASSES, ndim=3)

    def __call__(
        self,
        tokens: jax.Array,
        conditions_tokens: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        file_assembly_indices: jax.Array | None = None,
    ) -> dict[str, jax.Array]:
        _TRACES.append(1)
        logits = jax.nn.one_hot(tokens % NUM_CLASSES, NUM_CLASSES) * 10.0
        return {"bed_tracks_logits": self.norm(logits)[:, :, None, :]}


@contextlib.contextmanager
def _count_backend_compiles() -> Iterator[list[str]]:
    events: list[str] = []

    def listener(name: str, *_args: object, **_kwargs: object) -> None:
        if name == BACKEND_COMPILE_EVENT:
            events.append(name)

    jax.monitoring.register_event_duration_secs_listener(listener)
    try:
        yield events
    finally:
        # Private JAX API, stable under the jax<0.6 pin in requirements-jax.txt.
        jax._src.monitoring._unregister_event_duration_listener_by_callback(listener)


def _eval_model_on_mesh() -> _NormedStub:
    """An eval-mode stub whose weights and running stats sit on a mesh, as loaded."""
    model = _NormedStub()
    mesh = jax.sharding.Mesh(np.asarray(jax.devices()[:1]), ("dp",))
    replicated = jax.sharding.NamedSharding(mesh, jax.sharding.PartitionSpec())
    weights = nnx.state(model, nnx.Any(nnx.Param, nnx.BatchStat))
    nnx.update(model, jax.device_put(weights, replicated))
    model.eval()
    return model


def _inference(
    model: nnx.Module,
) -> tuple[Callable[..., dict[str, jax.Array]], jax.Array]:
    inference_fn = build_annotation_inference_fn(
        model,
        decoders=["argmax"],
        center_offset=0,
        rc_token_lut=jnp.arange(VOCAB),
    )
    tokens = jnp.asarray([[1, 2, 3, 4, 5, 6, 7, 8]], dtype=jnp.int32)
    return inference_fn, tokens


def _state_arrays(model: nnx.Module) -> list[jax.Array]:
    return [x for x in jax.tree.leaves(nnx.state(model)) if isinstance(x, jax.Array)]


def test_eval_model_compiles_once_across_calls() -> None:
    inference_fn, tokens = _inference(_eval_model_on_mesh())
    with _count_backend_compiles() as compiles:
        first = np.asarray(inference_fn(tokens, [], None)["argmax"])
        second = np.asarray(inference_fn(tokens, [], None)["argmax"])
    n_compiles = len(compiles)
    np.testing.assert_array_equal(first, second)
    assert n_compiles == 1, f"expected one compile, got {n_compiles}"


def test_both_strands_trace_the_network_once() -> None:
    inference_fn, tokens = _inference(_eval_model_on_mesh())
    _TRACES.clear()
    inference_fn(tokens, [], None)
    assert len(_TRACES) == 1, f"network traced {len(_TRACES)} times"


def test_commit_model_state_joins_the_parameters_mesh() -> None:
    model = _eval_model_on_mesh()
    replicated = next(x.sharding for x in _state_arrays(model) if x.committed)
    commit_model_state(model)
    arrays = _state_arrays(model)
    assert all(x.committed for x in arrays)
    assert all(x.sharding == replicated for x in arrays)


def test_commit_model_state_with_nothing_placed_moves_nothing() -> None:
    model = _NormedStub()
    model.eval()
    before = [x.devices() for x in _state_arrays(model)]
    commit_model_state(model)
    arrays = _state_arrays(model)
    assert all(x.committed for x in arrays)
    assert [x.devices() for x in arrays] == before
