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

"""Species-conditioned NTv3 with BigWig / BED (genome-annotation) track heads."""

from dataclasses import asdict, dataclass

import jax
from flax import nnx
from flax.typing import Dtype

from nucleotide_transformer_v3.heads import (
    AdaptiveClassificationHead,
    AdaptiveLinearHead,
)
from nucleotide_transformer_v3.discrete_conditioned_model import (
    DiscreteConditionedNTv3,
    DiscreteConditionedNTv3Config,
)
from nucleotide_transformer_v3.model import NTv3
from nucleotide_transformer_v3.layers import crop_sequence_center
from nucleotide_transformer_v3.types import AttentionMask, Tokens


@dataclass(kw_only=True)
class DiscreteConditionedNTv3WithTracksHeadsConfig(DiscreteConditionedNTv3Config):
    """Config for :class:`DiscreteConditionedNTv3WithTracksHeads`.

    Args:
        bigwigs_per_file_assembly: mapping from file assembly name to the list
            of bigwig file ids it predicts. Fixes the per-species head sizes.
        number_bigwig_tracks: maximum number of bigwig tracks across species
            (the output width; shorter species heads are zero-padded).
        bed_elements_names: BED / genome-annotation element class names.
        bed_num_classes: number of classes the BED head predicts per element.
        keep_target_center_fraction: fraction of the sequence center kept for
            target prediction (1.0 keeps the whole sequence).
        condition_bigwig_head: if ``True`` (and conditions exist) the bigwig head
            is FiLM-conditioned on the species embeddings; if ``False`` it is a
            plain (non-adaptive) head with no FiLM/bias params. Gates the bigwig
            head independently of the trunk (``enable_trunk_conditioning``) and of
            ``condition_annotation_head``.
        condition_annotation_head: same as ``condition_bigwig_head`` but for the
            BED / genome-annotation head. Defaults to ``False`` (annotation-head
            conditioning is not yet validated), so by default the bigwig head is
            conditioned while the annotation head is left unconditioned; set to
            ``True`` to also FiLM-condition the annotation head.
        bigwig_head_compute_dtype: compute dtype for the bigwig head's output
            projection. The head norm and output-bias layers stay float32, and
            the logits are cast back to float32, so no output dtype changes.
            Worth lowering: at ``number_bigwig_tracks=7362`` this projection is
            ~0.93x the FLOPs of the whole trunk. This is the one
            ``build_bf16_compute_model``'s ``cast_head`` flips.
        bigwig_head_param_dtype: storage dtype for the bigwig head's projection.
        annotation_head_compute_dtype: same for the BED head. Left out of
            ``cast_head`` (~0.19% of trunk FLOPs, and argmax reads its logits
            directly), so only ever lowered by setting this field.
        annotation_head_param_dtype: storage dtype for the BED head's projection.
    """

    bigwigs_per_file_assembly: dict[str, list[str]]
    number_bigwig_tracks: int
    bed_elements_names: list[str]
    bed_num_classes: int = 2
    keep_target_center_fraction: float = 1.0
    condition_bigwig_head: bool = True
    condition_annotation_head: bool = False

    bigwig_head_param_dtype: Dtype = "float32"
    bigwig_head_compute_dtype: Dtype = "float32"
    annotation_head_param_dtype: Dtype = "float32"
    annotation_head_compute_dtype: Dtype = "float32"

    @classmethod
    def from_base_config(  # type: ignore[override]
        cls,
        config: DiscreteConditionedNTv3Config,
        bigwigs_per_file_assembly: dict[str, list[str]],
        bed_elements_names: list[str],
        bed_num_classes: int = 2,
        keep_target_center_fraction: float = 1.0,
        condition_bigwig_head: bool = True,
        condition_annotation_head: bool = False,
    ) -> "DiscreteConditionedNTv3WithTracksHeadsConfig":
        number_bw_tracks = (
            max(len(bws) for bws in bigwigs_per_file_assembly.values())
            if bigwigs_per_file_assembly
            else 0
        )
        return cls(
            **asdict(config),
            bigwigs_per_file_assembly=bigwigs_per_file_assembly,
            number_bigwig_tracks=number_bw_tracks,
            bed_elements_names=bed_elements_names[:],
            bed_num_classes=bed_num_classes,
            keep_target_center_fraction=keep_target_center_fraction,
            condition_bigwig_head=condition_bigwig_head,
            condition_annotation_head=condition_annotation_head,
        )


class DiscreteConditionedNTv3WithTracksHeads(DiscreteConditionedNTv3):
    """Species-conditioned NTv3 backbone plus BigWig and BED track heads.

    Trunk and head conditioning are decoupled, giving three modes:

    * Full conditioning (``enable_trunk_conditioning=True`` + conditions): the
      NTv3 trunk is FiLM-conditioned AND the tracks heads are FiLM-conditioned.
    * Head-only conditioning (``enable_trunk_conditioning=False`` + conditions):
      the trunk is shared/unconditioned (``conditions_dims == ()``) while the
      bigwig/BED heads are still FiLM-conditioned on the species embeddings.
    * No conditioning (no conditions, ``conditions_vocab_size`` empty): the
      heads are plain (non-adaptive) and the trunk is unconditioned.

    The heads consume the 1bp decoder embedding (``conv_init_embed_dim``
    channels). Each head is FiLM-conditioned iff conditions exist AND its own
    per-head flag (``condition_bigwig_head`` / ``condition_annotation_head``) is
    ``True``, independent of the trunk flag and of the other head; a head with
    conditioning disabled is built with an empty ``conditions_dims`` and carries
    no FiLM/bias params. This allows conditioning one head but not the other
    (e.g. bigwig-only conditioning). All FiLM (trunk and heads) is zero-init
    identity, so at step 0 the outputs match the unconditioned heads regardless
    of the mode.
    """

    config: DiscreteConditionedNTv3WithTracksHeadsConfig  # type: ignore[assignment]

    def __init__(
        self,
        config: DiscreteConditionedNTv3WithTracksHeadsConfig,
        rngs: nnx.Rngs,
    ) -> None:
        super().__init__(config, rngs=rngs)

        head_dim = config.conv_init_embed_dim
        # Head FiLM dims are populated whenever conditions exist, independent of
        # the trunk (config.conditions_dims, which is gated by
        # enable_trunk_conditioning). This enables head-only conditioning. Each
        # head is gated by its own flag: a disabled head gets empty conditions_dims
        # (no FiLM/bias params), so conditioning can target one head only.
        head_cond_dims = [
            config.species_embed_dim for _ in config.conditions_vocab_size
        ]
        bigwig_cond_dims = head_cond_dims if config.condition_bigwig_head else []
        annotation_cond_dims = (
            head_cond_dims if config.condition_annotation_head else []
        )

        # --- BigWig head ---
        if config.number_bigwig_tracks > 0:
            self.bigwig_head = AdaptiveLinearHead(
                embed_dim=head_dim,
                num_labels=config.number_bigwig_tracks,
                conditions_dims=bigwig_cond_dims,
                dtype=config.bigwig_head_compute_dtype,
                param_dtype=config.bigwig_head_param_dtype,
                rngs=rngs,
            )

        # --- BED / genome-annotation head ---
        if len(config.bed_elements_names) > 0:
            self.bed_head = AdaptiveClassificationHead(
                embed_dim=head_dim,
                num_labels=len(config.bed_elements_names),
                num_classes=config.bed_num_classes,
                conditions_dims=annotation_cond_dims,
                use_output_bias=True,
                dtype=config.annotation_head_compute_dtype,
                param_dtype=config.annotation_head_param_dtype,
                rngs=rngs,
            )

    def __call__(  # type: ignore[override]
        self,
        tokens: Tokens,
        conditions_tokens: list[Tokens],
        file_assembly_indices: jax.Array | None = None,
        conditions_masks: list[AttentionMask] | None = None,
        *,
        single_head_idx: int | None = None,
        track_indices: jax.Array | None = None,
        skip_bed_head: bool = False,
        skip_lm_head: bool = False,
    ) -> dict[str, jax.Array]:
        """Forward pass.

        Args:
            tokens: integer token ids of shape ``(batch_size, sequence_length)``.
            conditions_tokens: list of discrete condition token ids, each of
                shape ``(batch_size,)`` (e.g. species ids).
            file_assembly_indices: unused; kept for interface compatibility with
                the multi-species trainer step functions.
            conditions_masks: optional per-example condition masks
                ``[(batch_size,)]``.
            single_head_idx: unused; kept for interface compatibility - this
                architecture has a single flat bigwig head rather than one head
                per assembly, so track indices are absolute and no head needs
                selecting.
            track_indices: when given, compute only these bigwig tracks. The
                head is ``number_bigwig_tracks`` wide, so leaving this unset
                materialises the entire panel at 1 bp resolution even when a
                caller wants a handful of tracks. Indices are positions in
                ``bigwigs_per_file_assembly[<assembly>]``, which are also head
                columns because the single flat head is shared positionally
                across assemblies - there is no per-assembly offset to apply,
                and adding one would be wrong rather than merely redundant,
                since no canonical cross-species assay schema backs those
                columns.
            skip_bed_head: if ``True``, skip the BED head.
            skip_lm_head: if ``True``, skip the LM head (no ``"logits"`` output).
        """
        del file_assembly_indices, single_head_idx
        # Embed conditions once and reuse for both the trunk and the heads.
        conditions = self._embed_conditions(conditions_tokens)
        # Call NTv3.__call__ directly with the already-embedded conditions to
        # avoid re-embedding in DiscreteConditionedNTv3.__call__.
        output_dict: dict[str, jax.Array] = NTv3.__call__(
            self,
            tokens,
            conditions=self._trunk_conditions(conditions),
            conditions_masks=conditions_masks,
            skip_lm_head=skip_lm_head,
        )
        # Per-condition vocabulary logits (reuses the trunk output already
        # computed above; no extra trunk pass).
        output_dict.update(
            self._condition_logits(output_dict["after_transformer_embedding"])
        )

        embedding = crop_sequence_center(
            output_dict["embedding"],
            self.config.keep_target_center_fraction,
        )

        # Each head is conditioned iff conditions exist AND its own per-head flag
        # is set, independent of the trunk's enable_trunk_conditioning (which
        # drives config.conditions_dims) and of the other head. A head that is
        # not conditioned gets an empty conditions list (no FiLM modulation),
        # matching its empty conditions_dims from __init__.
        conditions_exist = bool(self.config.conditions_vocab_size)

        if hasattr(self, "bigwig_head"):
            bigwig_conditioned = conditions_exist and self.config.condition_bigwig_head
            output_dict["bigwig_tracks_logits"] = self.bigwig_head(
                embedding,
                conditions if bigwig_conditioned else [],
                conditions_masks if bigwig_conditioned else None,
                track_indices=track_indices,
            )

        if hasattr(self, "bed_head") and not skip_bed_head:
            bed_conditioned = conditions_exist and self.config.condition_annotation_head
            output_dict["bed_tracks_logits"] = self.bed_head(
                embedding,
                conditions if bed_conditioned else [],
                conditions_masks if bed_conditioned else None,
            )

        return output_dict
