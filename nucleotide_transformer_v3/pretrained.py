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

"""Load a released NTv3 model from a directory of Flax artifacts.

A model directory holds::

    hyperparameters.json          the model config
    flax_model.msgpack            the NNX variable state
    tokenizer_config.json         the sequence tokenizer, optional
    conditions_tokenizers.json    one tokenizer per discrete condition, required
                                  iff the config declares conditions

Nothing here reads Orbax, so the package needs neither ``orbax-checkpoint`` nor
``tensorstore``. Every loader takes either such a directory or the name of a
released model (``RELEASED_MODELS``), which is fetched from the Hugging Face Hub
into the local cache. ``huggingface_hub`` is imported only then, so loading
from a local directory needs nothing beyond the package's own dependencies.
"""

import dataclasses
import json
from pathlib import Path
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx, serialization, traverse_util
from flax.nnx import filterlib

from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.dtypes import build_bf16_compute_model
from nucleotide_transformer_v3.logging import get_logger
from nucleotide_transformer_v3.model import NTv3
from nucleotide_transformer_v3.tokenizers import (
    BaseTokenizer,
    NucleotideTokenizer,
    get_ntv3_tokenizer,
)
from nucleotide_transformer_v3.with_tracks_heads import (
    DiscreteConditionedNTv3WithTracksHeads,
    DiscreteConditionedNTv3WithTracksHeadsConfig,
)

logger = get_logger()

CONFIG_FILE = "hyperparameters.json"
WEIGHTS_FILE = "flax_model.msgpack"
TOKENIZER_FILE = "tokenizer_config.json"
CONDITIONS_TOKENIZERS_FILE = "conditions_tokenizers.json"

# Released model name -> Hugging Face repository. The Flax artifacts live in the
# repository's ARTIFACTS_SUBFOLDER, beside the PyTorch weights.
RELEASED_MODELS: dict[str, str] = {
    "NTv3_pretrained": "InstaDeepAI/NTv3_pretrained",
    "NTv3_tracks": "InstaDeepAI/NTv3_tracks",
    "NTv3_tracks_1mb": "InstaDeepAI/NTv3_tracks_1mb",
    "NTv3_annot_mammals": "InstaDeepAI/NTv3_annot_mammals",
    "NTv3_annot_plants": "InstaDeepAI/NTv3_annot_plants",
    "NTv3_annot_insects": "InstaDeepAI/NTv3_annot_insects",
}
ARTIFACTS_SUBFOLDER = "jax_model"

# Params alone would drop RMSBatchNorm.var_ema, which eval-mode normalization
# divides by. The model would then run against a ones-init and be quietly wrong.
VARIABLES = filterlib.Any(nnx.Param, nnx.BatchStat)


def parse_hyperparameters(
    values: dict[str, Any],
    config_class: type[NTv3Config] = DiscreteConditionedNTv3WithTracksHeadsConfig,
) -> Any:
    """Build a config dataclass from a hyperparameters dict.

    Keys the dataclass does not declare are dropped, so an artifact carrying
    provenance or derived fields still loads. They are named in a warning,
    because a misspelled field would otherwise fall back to its default with no
    sign: a wrong crop, for one, changes no tensor shape and so passes every
    check on the weights.
    """
    accepted = {f.name for f in dataclasses.fields(config_class) if f.init}
    if ignored := sorted(values.keys() - accepted - {"provenance"}):
        logger.warning(
            f"{CONFIG_FILE} keys ignored by {config_class.__name__}: {ignored}"
        )
    return config_class(**{k: v for k, v in values.items() if k in accepted})


def _leaf_shapes(tree: dict[str, Any]) -> dict[str, tuple[int, ...]]:
    """Map each leaf's ``/``-joined path to its shape."""
    flat = traverse_util.flatten_dict(tree)
    return {"/".join(map(str, path)): np.shape(leaf) for path, leaf in flat.items()}


def load_variables(model: nnx.Module, weights_path: Path) -> None:
    """Restore params *and* batch statistics from a msgpack blob, in place.

    The blob must match the model's variable tree exactly. ``from_state_dict``
    only enforces that no expected key is missing: it discards tensors the model
    has no slot for and takes a mis-shaped tensor as is, so an 11-token
    embedding would load into a 16-token model and every lowercase base would
    gather out of range, silently.

    Raises:
        ValueError: If the blob is missing a variable, carries one the model
            does not have, or stores one at a different shape.
    """
    target = nnx.state(model, VARIABLES).to_pure_dict()
    state = serialization.msgpack_restore(weights_path.read_bytes())
    try:
        restored = serialization.from_state_dict(target, state)
    except ValueError as error:
        raise ValueError(
            f"{weights_path} does not match this model's variable tree. The "
            "usual cause is an export written with a nnx.Param-only filter, "
            "which omits the RMSBatchNorm var_ema buffers; re-export with "
            f"filterlib.Any(nnx.Param, nnx.BatchStat). Underlying error: {error}"
        ) from error

    expected, stored = _leaf_shapes(target), _leaf_shapes(state)
    if unexpected := sorted(stored.keys() - expected.keys()):
        raise ValueError(
            f"{weights_path} carries {len(unexpected)} variables this model has "
            f"no place for, starting with {unexpected[:3]}. The weights belong to "
            "a different model than the one configured."
        )
    if reshaped := {
        path: (stored[path], shape)
        for path, shape in expected.items()
        if stored[path] != shape
    }:
        path, (found, wanted) = next(iter(reshaped.items()))
        raise ValueError(
            f"{weights_path} stores {len(reshaped)} variables at a shape this "
            f"model does not expect, e.g. {path}: {found} instead of {wanted}. "
            f"{CONFIG_FILE} does not describe these weights."
        )
    nnx.update(model, restored)


def _warn_on_ones_batch_stats(model: nnx.Module) -> None:
    """Warn when every var_ema came back at its ones-init.

    ``RMSBatchNorm`` divides by ``var_ema`` in eval mode, so a model that lost
    those buffers runs and is quietly wrong: measured on a tiny model, a 0.83
    shift in the annotation logits.

    A warning rather than an error, because an all-ones state is also what a
    freshly initialised model legitimately has, and test fixtures are built
    that way. The case this is really guarding -- a msgpack written with a
    ``nnx.Param``-only filter -- is already caught hard by ``load_variables``,
    where the key sets cannot match.
    """
    stats = jax.tree_util.tree_leaves(nnx.state(model, nnx.BatchStat).to_pure_dict())
    if stats and all(bool(jnp.all(s == 1.0)) for s in stats):
        logger.warning(
            f"All {len(stats)} batch-statistic buffers are exactly 1.0. That is "
            "expected for an untrained model and wrong for a released one, "
            "where it would mean the export dropped them."
        )


def load_tokenizer(model_dir: Path, alphabet_size: int) -> NucleotideTokenizer:
    """Rebuild the tokenizer, falling back to the vocabulary the size implies."""
    config_path = model_dir / TOKENIZER_FILE
    if config_path.exists():
        tokenizer = BaseTokenizer.from_config(json.loads(config_path.read_text()))
    else:
        tokenizer = get_ntv3_tokenizer(include_softmasked_tokens=alphabet_size == 16)
        logger.warning(f"No {TOKENIZER_FILE}; built a {alphabet_size}-token default.")
    if tokenizer.vocabulary_size != alphabet_size:
        raise ValueError(
            f"Tokenizer has {tokenizer.vocabulary_size} tokens but the model "
            f"expects {alphabet_size}. Soft-masked checkpoints need the "
            "16-token vocabulary (get_ntv3_tokenizer(include_softmasked_tokens=True))."
        )
    return tokenizer  # type: ignore[return-value]


def load_conditions_tokenizers(model_dir: Path, config: Any) -> list[BaseTokenizer]:
    """Rebuild one tokenizer per discrete condition, in ``conditions_names`` order.

    The species ids a tracks model was trained on exist nowhere but in this
    file, and a wrong id changes every track through the FiLM-conditioned head
    without raising. So a model that declares conditions refuses to load
    without it, rather than loading half-usable.

    Raises:
        FileNotFoundError: If the config declares conditions and the file is
            missing.
        ValueError: If a condition has no tokenizer, or one whose vocabulary
            size differs from the config's.
    """
    vocab_sizes = list(getattr(config, "conditions_vocab_size", None) or [])
    if not vocab_sizes:
        return []
    path = model_dir / CONDITIONS_TOKENIZERS_FILE
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. The config declares {len(vocab_sizes)} discrete "
            "condition(s), and without their vocabulary no condition id can be "
            "trusted."
        )
    configs = json.loads(path.read_text())["conditions_tokenizers"]
    names = config.conditions_names or [
        f"condition_{i}" for i in range(len(vocab_sizes))
    ]
    tokenizers = []
    for name, vocab_size in zip(names, vocab_sizes, strict=True):
        if name not in configs:
            raise ValueError(f"{path} has no tokenizer for condition {name!r}.")
        tokenizer = BaseTokenizer.from_config(configs[name])
        if tokenizer.vocabulary_size != vocab_size:
            raise ValueError(
                f"The {name!r} tokenizer has {tokenizer.vocabulary_size} tokens "
                f"but the model's embedding table has {vocab_size}."
            )
        tokenizers.append(tokenizer)
    return tokenizers


def resolve_model_dir(model: str | Path, *, revision: str | None = None) -> Path:
    """Turn a released model name or a local directory into a local directory.

    A released name downloads the Flax artifacts of ``RELEASED_MODELS[model]``
    into the Hugging Face cache, or reuses them if already there.

    Args:
        model: A local directory, or a key of ``RELEASED_MODELS``.
        revision: Branch, tag or commit of the Hub repository. Ignored for a
            local directory.

    Raises:
        FileNotFoundError: If ``model`` is neither a directory nor a released
            name.
    """
    if Path(model).is_dir():
        return Path(model)
    repo_id = RELEASED_MODELS.get(str(model))
    if repo_id is None:
        raise FileNotFoundError(
            f"{model} is not a directory, nor one of the released models "
            f"{sorted(RELEASED_MODELS)}."
        )
    # Imported here so that loading from a directory never needs it.
    from huggingface_hub import snapshot_download

    snapshot = snapshot_download(
        repo_id,
        revision=revision,
        # The .npz is the export's round-trip fixture, which check_release
        # replays: 164 KB for the pretrained model, 40 MB for tracks.
        allow_patterns=[
            f"{ARTIFACTS_SUBFOLDER}/*.json",
            f"{ARTIFACTS_SUBFOLDER}/*.msgpack",
            f"{ARTIFACTS_SUBFOLDER}/*.npz",
        ],
    )
    return Path(snapshot) / ARTIFACTS_SUBFOLDER


def _reject_head_fields(values: dict[str, Any], config_class: type[NTv3Config]) -> None:
    """Refuse a config that declares heads or conditioning the class lacks.

    ``parse_hyperparameters`` drops keys the dataclass does not declare, so a
    tracks config parsed as ``NTv3Config`` would build a bare backbone.
    """
    accepted = {f.name for f in dataclasses.fields(config_class) if f.init}
    tracks_only = {
        f.name
        for f in dataclasses.fields(DiscreteConditionedNTv3WithTracksHeadsConfig)
        if f.init
    }
    if foreign := sorted((tracks_only - accepted) & values.keys()):
        raise ValueError(
            f"{CONFIG_FILE} declares {foreign[:3]}, which {config_class.__name__} "
            "cannot build. These are tracks or annotation weights; load them "
            "with load_tracks_model or load_annotation_model."
        )


def _load(
    model: str | Path,
    *,
    model_class: type[NTv3],
    config_class: type[NTv3Config],
    use_bf16_compute: bool,
    revision: str | None,
) -> tuple[Any, NucleotideTokenizer, list[BaseTokenizer]]:
    """Load a model of ``model_class`` and its tokenizers."""
    model_dir = resolve_model_dir(model, revision=revision)
    config_path, weights_path = model_dir / CONFIG_FILE, model_dir / WEIGHTS_FILE
    for path in (config_path, weights_path):
        if not path.exists():
            raise FileNotFoundError(f"{path} not found in {model_dir}")

    values = json.loads(config_path.read_text())
    _reject_head_fields(values, config_class)
    config = parse_hyperparameters(values, config_class=config_class)
    tokenizer = load_tokenizer(model_dir, config.alphabet_size)
    conditions_tokenizers = load_conditions_tokenizers(model_dir, config)

    loaded = model_class(config=config, rngs=nnx.Rngs(0))
    loaded.eval()
    load_variables(loaded, weights_path)
    _warn_on_ones_batch_stats(loaded)

    if use_bf16_compute:
        loaded = build_bf16_compute_model(loaded)
        loaded.eval()
    return loaded, tokenizer, conditions_tokenizers


def load_pretrained_model(
    model_dir: str | Path,
    *,
    use_bf16_compute: bool = True,
    revision: str | None = None,
) -> tuple[NTv3, NucleotideTokenizer]:
    """Load a pretrained NTv3 backbone and its tokenizer.

    The backbone always carries its LM head, so ``model(tokens)["logits"]`` is
    available; pass ``skip_lm_head=True`` when only embeddings are wanted.

    Its 11-token vocabulary has no lowercase bases: soft-masked sequence
    tokenizes to ``<unk>``, so uppercase it first.

    Args:
        model_dir: Directory holding the artifacts listed in the module
            docstring, or a released name such as ``"NTv3_pretrained"``.
        use_bf16_compute: Rebuild with bfloat16 compute dtypes. Parameter
            storage stays as exported.
        revision: Hub revision, when ``model_dir`` is a released name.

    Returns:
        The model in eval mode, and its sequence tokenizer. Its ``provenance``
        block, which the config dataclass does not carry, is available from
        :func:`read_provenance`.

    Raises:
        FileNotFoundError: If the config or weights file is missing.
        ValueError: If the artifacts carry tracks heads or conditioning, if the
            tokenizer disagrees with the model, or if the weights do not match
            the model's variable tree.
    """
    model, tokenizer, _ = _load(
        model_dir,
        model_class=NTv3,
        config_class=NTv3Config,
        use_bf16_compute=use_bf16_compute,
        revision=revision,
    )
    logger.info(
        f"Loaded {model_dir}: pretrained backbone, "
        f"alphabet {model.config.alphabet_size}, bf16={use_bf16_compute}."
    )
    return model, tokenizer


def load_tracks_model(
    model_dir: str | Path,
    *,
    use_bf16_compute: bool = True,
    revision: str | None = None,
) -> tuple[
    DiscreteConditionedNTv3WithTracksHeads, NucleotideTokenizer, list[BaseTokenizer]
]:
    """Load an NTv3 functional-tracks model and its tokenizers.

    The bigwig head is ``number_bigwig_tracks`` wide at 1 bp resolution, so a
    forward pass without ``track_indices`` materialises the whole panel. At the
    released 7362 tracks and a 1 Mb window that is the difference between a
    working call and an out-of-memory error.

    The model takes one condition token per example and condition, in
    ``config.conditions_names`` order, and the returned condition tokenizers
    come in that same order. For the released model that is a single
    ``species`` condition: ``conditions_tokenizers[0].token_to_id("human")``.

    Args:
        model_dir: Directory holding the artifacts listed in the module
            docstring, or a released name such as ``"NTv3_tracks"``.
        use_bf16_compute: Rebuild with bfloat16 compute dtypes. Parameter
            storage stays as exported, and the bigwig head keeps float32
            compute.
        revision: Hub revision, when ``model_dir`` is a released name.

    Returns:
        The model in eval mode, its sequence tokenizer, and one tokenizer per
        discrete condition.

    Raises:
        FileNotFoundError: If the config, weights or condition tokenizers file
            is missing.
        ValueError: If the artifacts carry no bigwig head, declare a conditioned
            bigwig head with no conditioning to feed it, disagree with the
            tokenizer, or do not match the model's variable tree.
    """
    model, tokenizer, conditions_tokenizers = _load(
        model_dir,
        model_class=DiscreteConditionedNTv3WithTracksHeads,
        config_class=DiscreteConditionedNTv3WithTracksHeadsConfig,
        use_bf16_compute=use_bf16_compute,
        revision=revision,
    )
    config = model.config
    if not config.number_bigwig_tracks:
        raise ValueError(
            "These artifacts declare no bigwig tracks. A genome-annotation "
            "checkpoint loads with load_annotation_model."
        )
    # The species embedding reaches the bigwig head through FiLM. A config that
    # claims a conditioned head but has no condition tables is an export that
    # dropped them in both places, which loads cleanly and predicts wrong.
    if config.condition_bigwig_head and not config.conditions_vocab_size:
        raise ValueError(
            "condition_bigwig_head is set but conditions_vocab_size is empty: "
            "the export dropped the species conditioning the bigwig head consumes."
        )
    logger.info(
        f"Loaded {model_dir}: {config.number_bigwig_tracks} bigwig "
        f"tracks, alphabet {config.alphabet_size}, "
        f"crop {config.keep_target_center_fraction}, bf16={use_bf16_compute}."
    )
    return model, tokenizer, conditions_tokenizers


def read_provenance(
    model_dir: str | Path, *, revision: str | None = None
) -> dict[str, Any]:
    """The ``provenance`` block an export records in its hyperparameters.

    Kept out of the model config on purpose: ``parse_hyperparameters`` keeps
    only the fields the config dataclass declares, and provenance is not one
    of them. It carries two things nothing else does -- which released model
    these weights are an export of, and the crop the head was *trained* at,
    which is not the crop it is *served* at.

    Args:
        model_dir: Directory of artifacts, or a released name.
        revision: Hub revision, when ``model_dir`` is a released name.

    Returns:
        The provenance block, or an empty dict for an export without one.
    """
    config_path = resolve_model_dir(model_dir, revision=revision) / CONFIG_FILE
    values: dict[str, Any] = json.loads(config_path.read_text())
    provenance = values.get("provenance", {})
    return provenance if isinstance(provenance, dict) else {}


def load_annotation_model(
    model_dir: str | Path,
    *,
    use_bf16_compute: bool = True,
    revision: str | None = None,
) -> tuple[DiscreteConditionedNTv3WithTracksHeads, NucleotideTokenizer]:
    """Load a genome-annotation model and its tokenizer.

    Args:
        model_dir: Directory holding the artifacts listed in the module
            docstring, or a released name such as ``"NTv3_annot_mammals"``.
        use_bf16_compute: Rebuild with bfloat16 compute dtypes. Parameter
            storage stays as exported.
        revision: Hub revision, when ``model_dir`` is a released name.

    Returns:
        The model in eval mode, and its sequence tokenizer.

    Raises:
        FileNotFoundError: If the config or weights file is missing.
        ValueError: If the artifacts carry no genome-annotation (BED) head, if
            the tokenizer disagrees with the model, or if the weights do not
            match the model's variable tree. Batch statistics that are all at
            their ones-init only warn, since that is also a fresh model's state.
    """
    model, tokenizer, _ = _load(
        model_dir,
        model_class=DiscreteConditionedNTv3WithTracksHeads,
        config_class=DiscreteConditionedNTv3WithTracksHeadsConfig,
        use_bf16_compute=use_bf16_compute,
        revision=revision,
    )
    config = model.config
    if not config.bed_elements_names:
        raise ValueError(
            "These artifacts declare no genome-annotation elements. A tracks "
            "checkpoint loads with load_tracks_model."
        )
    logger.info(
        f"Loaded {model_dir}: {config.bed_num_classes}-class annotation head, "
        f"alphabet {config.alphabet_size}, "
        f"crop {config.keep_target_center_fraction}, bf16={use_bf16_compute}."
    )
    return model, tokenizer


def count_parameters(model: nnx.Module) -> int:
    """Total number of scalar parameters, for reporting."""
    leaves = jax.tree_util.tree_leaves(nnx.state(model, nnx.Param).to_pure_dict())
    return int(sum(np.prod(np.shape(leaf)) for leaf in leaves))
