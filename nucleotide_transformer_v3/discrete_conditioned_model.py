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

"""NTv3 conditioned on discrete tokens (e.g. species) via FiLM trunk norms."""

import dataclasses
from typing import Any

import jax
import jax.numpy as jnp
from flax import nnx

from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.model import NTv3
from nucleotide_transformer_v3.types import Embedding, Tokens


@dataclasses.dataclass(kw_only=True)
class DiscreteConditionedNTv3Config(NTv3Config):
    """NTv3 config extended with discrete-condition (e.g. species) metadata.

    A discrete condition is a single integer token per example (a species id, a
    tissue id, ...). Each condition gets its own embedding table of dimension
    ``species_embed_dim``; the resulting embeddings are the FiLM conditions
    consumed by the trunk's adaptive norms (via ``NTv3Config.conditions_dims``).

    Args:
        conditions_vocab_size: vocabulary size of each discrete condition.
        conditions_names: optional human-readable names, one per condition.
        species_embed_dim: embedding dimension of every condition table (also
            the FiLM condition dim).
        enable_trunk_conditioning: gates trunk FiLM ONLY. When ``True`` and
            conditions exist, ``conditions_dims`` is populated so the NTv3
            trunk's adaptive norms become FiLM-conditioned. When ``False``
            ``conditions_dims`` stays ``()`` and the trunk is shared/unconditioned.
            This is independent of the tracks heads: a head in
            :class:`DiscreteConditionedNTv3WithTracksHeads` is FiLM-conditioned
            when conditions exist (``conditions_vocab_size`` non-empty) and its
            own ``condition_bigwig_head`` / ``condition_annotation_head`` flag
            is set, which is what enables head-only conditioning.
    """

    conditions_vocab_size: list[int]
    conditions_names: list[str] | None = None
    species_embed_dim: int = 16
    enable_trunk_conditioning: bool = True

    def __post_init__(self) -> None:
        super().__post_init__()
        # conditions_dims controls TRUNK conditioning only, gated by the flag,
        # which must be authoritative for trunk conditioning.
        if self.enable_trunk_conditioning:
            if not self.conditions_dims:
                self.conditions_dims = tuple(
                    self.species_embed_dim for _ in self.conditions_vocab_size
                )
        elif self.conditions_dims:
            raise ValueError(
                "enable_trunk_conditioning=False but conditions_dims is non-empty; "
                "the flag must be authoritative for trunk conditioning."
            )

    @classmethod
    def from_base_config(
        cls,
        config: NTv3Config,
        conditions_vocab_size: list[int],
        conditions_names: list[str] | None = None,
        species_embed_dim: int = 16,
        enable_trunk_conditioning: bool = True,
    ) -> "DiscreteConditionedNTv3Config":
        """Build a :class:`DiscreteConditionedNTv3Config` from a base config."""
        return cls(
            **dataclasses.asdict(config),
            conditions_vocab_size=conditions_vocab_size[:],
            conditions_names=None if conditions_names is None else conditions_names[:],
            species_embed_dim=species_embed_dim,
            enable_trunk_conditioning=enable_trunk_conditioning,
        )


class DiscreteConditionedNTv3(NTv3):
    """NTv3 conditioned on one or more discrete tokens.

    Each discrete condition token is embedded via its own ``nnx.Embed`` table;
    the embeddings are passed as FiLM ``conditions`` into the NTv3 trunk. The
    trunk sees this conditioning iff ``conditions_dims`` is non-empty, which is
    gated by ``enable_trunk_conditioning`` (trunk FiLM only). With trunk
    conditioning off (or no conditions) the base ``NTv3.__call__`` runs with a
    shared/unconditioned trunk; the tracks heads in
    :class:`DiscreteConditionedNTv3WithTracksHeads` are conditioned separately.

    A linear head per condition also predicts logits over that condition's
    vocabulary from the mean-pooled trunk output (``{name}_logits``), enabling
    joint modelling of tokens and conditions (e.g. masked-species prediction),
    mirroring the LongRangeNT conditioned model.
    """

    config: DiscreteConditionedNTv3Config  # type: ignore[assignment]

    def __init__(
        self,
        config: DiscreteConditionedNTv3Config,
        *,
        rngs: nnx.Rngs,
    ) -> None:
        super().__init__(config=config, rngs=rngs)

        conditions_names = config.conditions_names
        conditions_vocab_size = config.conditions_vocab_size
        if conditions_names is None:
            conditions_names = [
                f"condition_{i}" for i in range(len(conditions_vocab_size))
            ]
        if len(conditions_names) != len(conditions_vocab_size):
            raise ValueError(
                "conditions_names and conditions_vocab_size must have equal length"
            )
        self._conditions_names = conditions_names

        # One embedding table per discrete condition.
        self.conditions_embed_layers = [
            nnx.Embed(
                num_embeddings=vocab_size,
                features=config.species_embed_dim,
                rngs=rngs,
            )
            for vocab_size in conditions_vocab_size
        ]

        # One linear head per condition, predicting logits over its vocabulary
        # from the mean-pooled trunk output (joint token/condition modelling).
        self.conditions_heads = [
            nnx.Linear(
                in_features=config.embed_dim,
                out_features=vocab_size,
                rngs=rngs,
            )
            for vocab_size in conditions_vocab_size
        ]

    def _condition_logits(
        self, after_transformer_embedding: jax.Array
    ) -> dict[str, jax.Array]:
        """Per-condition vocabulary logits from the mean-pooled trunk output.

        Args:
            after_transformer_embedding: trunk output of shape
                ``(batch_size, reduced_length, embed_dim)``.

        Returns:
            A ``{f"{name}_logits": (batch_size, condition_vocab_size)}`` dict,
            one entry per discrete condition.
        """
        pooled = jnp.mean(after_transformer_embedding, axis=-2)
        pooled = jax.nn.gelu(pooled)
        return {
            f"{name}_logits": head(pooled)
            for name, head in zip(
                self._conditions_names, self.conditions_heads, strict=True
            )
        }

    def _embed_conditions(self, conditions_tokens: list[jax.Array]) -> list[jax.Array]:
        """Look up condition embeddings: ``[(B,)] -> [(B, species_embed_dim)]``."""
        return [
            embed_layer(cond_tokens)
            for cond_tokens, embed_layer in zip(
                conditions_tokens, self.conditions_embed_layers, strict=True
            )
        ]

    def _trunk_conditions(self, conditions: list[jax.Array]) -> list[jax.Array] | None:
        """Conditions to pass to the trunk (``None`` when trunk conditioning off)."""
        return conditions if self.config.conditions_dims else None

    def __call__(  # type: ignore[override]
        self,
        tokens: Tokens,
        conditions_tokens: list[jax.Array],
        conditions_masks: list[jax.Array] | None = None,
        stop_after_transformer: bool = False,
        skip_lm_head: bool = False,
        **kwargs: Any,
    ) -> dict[str, Embedding]:
        """Forward pass.

        Args:
            tokens: integer token ids of shape ``(batch_size, sequence_length)``.
            conditions_tokens: list of discrete condition token ids, each of
                shape ``(batch_size,)``.
            conditions_masks: optional list of per-example masks, each of shape
                ``(batch_size,)``. Where ``False``, the condition contributes
                zero FiLM modulation (condition dropout).
            skip_lm_head: Omit ``"logits"`` and skip the LM head projection;
                see :meth:`NTv3.__call__`. Threaded explicitly rather than left
                to ``**kwargs``, which is deleted below - a silently dropped
                flag here would leave the LM head running with no error and no
                symptom beyond the memory it was meant to save.
            stop_after_transformer: Skip the decoder tower and LM head; see
                :meth:`NTv3.__call__`. The per-condition logits are unaffected,
                since they are derived from the trunk output.
        """
        del kwargs  # accepted for forward compatibility
        conditions = self._embed_conditions(conditions_tokens)
        outs = super().__call__(
            tokens,
            conditions=self._trunk_conditions(conditions),
            conditions_masks=conditions_masks,
            stop_after_transformer=stop_after_transformer,
            skip_lm_head=skip_lm_head,
        )
        outs.update(self._condition_logits(outs["after_transformer_embedding"]))
        return outs
