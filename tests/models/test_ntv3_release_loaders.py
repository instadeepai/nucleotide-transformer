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

"""The pretrained and tracks loaders rebuild exactly the model that was exported.

Every model here is tiny and randomly initialised, so the suite needs no weights
and no network. What it pins is the export contract and the ways a released
artifact can load cleanly and still be wrong: batch statistics left behind, the
species conditioning stripped, or one family's weights read as another's.
"""

import dataclasses
import json
import logging
import sys
import types
from pathlib import Path
from typing import Any, TypeVar

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from flax import nnx, serialization

from nucleotide_transformer_v3.check_release import FIXTURE_FILE, check
from nucleotide_transformer_v3.config import NTv3Config
from nucleotide_transformer_v3.model import NTv3
from nucleotide_transformer_v3.pretrained import (
    ARTIFACTS_SUBFOLDER,
    CONDITIONS_TOKENIZERS_FILE,
    CONFIG_FILE,
    RELEASED_MODELS,
    TOKENIZER_FILE,
    VARIABLES,
    WEIGHTS_FILE,
    load_annotation_model,
    load_pretrained_model,
    load_tracks_model,
    load_variables,
)
from nucleotide_transformer_v3.tokenizers import (
    StandardTokenizer,
    get_ntv3_tokenizer,
)
from nucleotide_transformer_v3.with_tracks_heads import (
    DiscreteConditionedNTv3WithTracksHeads,
    DiscreteConditionedNTv3WithTracksHeadsConfig,
)

# Divisible by the total downsampling, 2 ** (num_downsamples + 1) = 8.
SEQUENCE_LENGTH = 32
BATCH_SIZE = 2
# Ids 0-5 are the tokenizer's special tokens, so these two are 6 and 7.
SPECIES = ["human", "mouse"]
SPECIES_VOCAB_SIZE = 6 + len(SPECIES)

ModuleT = TypeVar("ModuleT", bound=nnx.Module)

TINY_BACKBONE: dict[str, Any] = {
    "token_embed_dim": 8,
    "conv_init_embed_dim": 16,
    "embed_dim": 32,
    "num_downsamples": 2,
    "num_layers": 1,
    "attention_heads": 2,
    "num_kv_heads": 1,
    "key_size": 8,
    "value_size": 8,
}


def _pretrained_config(alphabet_size: int = 11) -> NTv3Config:
    return NTv3Config(alphabet_size=alphabet_size, **TINY_BACKBONE)


def _tracks_config(
    condition_bigwig_head: bool = True,
    bed_elements_names: tuple[str, ...] = (),
) -> DiscreteConditionedNTv3WithTracksHeadsConfig:
    """The released tracks topology, shrunk: head-only species conditioning."""
    return DiscreteConditionedNTv3WithTracksHeadsConfig(
        alphabet_size=16,
        **TINY_BACKBONE,
        conditions_vocab_size=[SPECIES_VOCAB_SIZE],
        conditions_names=["species"],
        species_embed_dim=4,
        enable_trunk_conditioning=False,
        bigwigs_per_file_assembly={"hg38": ["a", "b", "c"], "mm10": ["a", "b"]},
        number_bigwig_tracks=3,
        bed_elements_names=list(bed_elements_names),
        bed_num_classes=3,
        keep_target_center_fraction=0.375,
        condition_bigwig_head=condition_bigwig_head,
    )


def _randomized(model: ModuleT, seed: int = 0) -> ModuleT:
    """Overwrite every param and batch statistic with random values.

    At init the FiLM projections are zero and ``var_ema`` is all ones, so a
    dropped conditioning table or batch statistic would round-trip unnoticed.
    """
    key = jax.random.key(seed)
    for kind, offset in ((nnx.Param, 0.0), (nnx.BatchStat, 1.0)):
        leaves, treedef = jax.tree_util.tree_flatten(nnx.state(model, kind))
        key, *keys = jax.random.split(key, len(leaves) + 1)
        # Batch statistics are variances, so they stay positive.
        noisy = [
            offset + 0.2 * jax.random.normal(k, leaf.shape, leaf.dtype)
            for k, leaf in zip(keys, leaves, strict=True)
        ]
        nnx.update(model, jax.tree_util.tree_unflatten(treedef, noisy))
    return model


def _state_dict(model: nnx.Module) -> dict[str, Any]:
    return jax.tree_util.tree_map(
        np.asarray, nnx.state(model, VARIABLES).to_pure_dict()
    )


def _write_artifacts(
    out_dir: Path,
    config: NTv3Config,
    variables: dict[str, Any],
    *,
    alphabet_size: int | None = None,
) -> Path:
    """Write the files the loaders read, as the release export does."""
    out_dir.mkdir(parents=True, exist_ok=True)
    values = dataclasses.asdict(config)
    (out_dir / CONFIG_FILE).write_text(json.dumps(values))
    (out_dir / WEIGHTS_FILE).write_bytes(serialization.to_bytes(variables))
    tokenizer = get_ntv3_tokenizer(
        include_softmasked_tokens=(alphabet_size or config.alphabet_size) == 16
    )
    (out_dir / TOKENIZER_FILE).write_text(json.dumps(tokenizer.to_config()))
    if values.get("conditions_vocab_size"):
        species = StandardTokenizer(standard_tokens=SPECIES).to_config()
        conditions = {"conditions_tokenizers": {"species": species}}
        (out_dir / CONDITIONS_TOKENIZERS_FILE).write_text(json.dumps(conditions))
    return out_dir


def _tokens(alphabet_size: int, seed: int = 1) -> jax.Array:
    # Ids 6 and up are nucleotides in both vocabularies; 0-5 are specials.
    return jax.random.randint(
        jax.random.key(seed), (BATCH_SIZE, SEQUENCE_LENGTH), 6, alphabet_size
    )


def _species(*ids: int) -> list[jax.Array]:
    return [jnp.asarray(ids, dtype=jnp.int32)]


def _assert_same_outputs(
    reference: dict[str, jax.Array], reloaded: dict[str, jax.Array]
) -> None:
    assert reference.keys() == reloaded.keys()
    for key in reference:
        np.testing.assert_array_equal(
            np.asarray(reference[key]), np.asarray(reloaded[key]), err_msg=key
        )


@pytest.fixture
def pretrained_model() -> NTv3:
    model = _randomized(NTv3(_pretrained_config(), rngs=nnx.Rngs(0)))
    model.eval()
    return model


def _tracks_model(**overrides: Any) -> DiscreteConditionedNTv3WithTracksHeads:
    model = DiscreteConditionedNTv3WithTracksHeads(
        _tracks_config(**overrides), rngs=nnx.Rngs(0)
    )
    model = _randomized(model)
    model.eval()
    return model


def test_pretrained_round_trip(tmp_path: Path, pretrained_model: NTv3) -> None:
    """A pretrained export reloads to bitwise-identical outputs."""
    out_dir = _write_artifacts(
        tmp_path, pretrained_model.config, _state_dict(pretrained_model)
    )
    reloaded, tokenizer = load_pretrained_model(out_dir, use_bf16_compute=False)

    assert isinstance(reloaded, NTv3)
    assert tokenizer.vocabulary_size == 11
    tokens = _tokens(11)
    _assert_same_outputs(pretrained_model(tokens), reloaded(tokens))


def test_pretrained_bf16_compute_stays_close(
    tmp_path: Path, pretrained_model: NTv3
) -> None:
    """The default bf16 rebuild keeps the loaded weights, so logits stay close."""
    out_dir = _write_artifacts(
        tmp_path, pretrained_model.config, _state_dict(pretrained_model)
    )
    reloaded, _ = load_pretrained_model(out_dir)

    tokens = _tokens(11)
    reference, served = pretrained_model(tokens), reloaded(tokens)
    assert served["logits"].dtype == jnp.float32
    np.testing.assert_allclose(served["logits"], reference["logits"], atol=1e-2)


# (True, ()) is the released checkpoint's corner, and the one nothing else in
# the repo exercises: a conditioned bigwig head with no BED head.
@pytest.mark.parametrize(
    ("condition_bigwig_head", "bed_elements_names"),
    [(True, ()), (False, ()), (True, ("exon", "intron")), (False, ("exon",))],
)
def test_tracks_round_trip(
    tmp_path: Path,
    condition_bigwig_head: bool,
    bed_elements_names: tuple[str, ...],
) -> None:
    """A tracks export reloads to bitwise-identical outputs on every head."""
    model = _tracks_model(
        condition_bigwig_head=condition_bigwig_head,
        bed_elements_names=bed_elements_names,
    )
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    reloaded, tokenizer, (species,) = load_tracks_model(out_dir, use_bf16_compute=False)

    assert tokenizer.vocabulary_size == 16
    assert species.token_to_id("mouse") == 7
    assert hasattr(reloaded, "bed_head") == bool(bed_elements_names)
    tokens, species_ids = _tokens(16), _species(6, 7)
    reference = model(tokens, species_ids)
    assert "species_logits" in reference
    _assert_same_outputs(reference, reloaded(tokens, species_ids))


def test_tracks_species_reaches_the_bigwig_head(tmp_path: Path) -> None:
    """The reloaded bigwig head still depends on the species, through FiLM."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    reloaded, _, _ = load_tracks_model(out_dir, use_bf16_compute=False)

    tokens = jnp.tile(_tokens(16)[:1], (2, 1))
    tracks = reloaded(tokens, _species(6, 7))["bigwig_tracks_logits"]
    assert not np.allclose(tracks[0], tracks[1])


def test_tracks_track_indices_select_head_columns(tmp_path: Path) -> None:
    """``track_indices`` computes a subset of the full panel, column for column."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    reloaded, _, _ = load_tracks_model(out_dir, use_bf16_compute=False)

    tokens, species, indices = _tokens(16), _species(6, 7), jnp.asarray([2, 0])
    full = reloaded(tokens, species)["bigwig_tracks_logits"]
    subset = reloaded(tokens, species, track_indices=indices)["bigwig_tracks_logits"]
    kept = int(SEQUENCE_LENGTH * 0.375)
    assert full.shape == (BATCH_SIZE, kept, 3)
    np.testing.assert_allclose(subset, full[..., indices], rtol=1e-6, atol=1e-6)


def test_param_only_export_is_refused(tmp_path: Path, pretrained_model: NTv3) -> None:
    """Trap: a ``nnx.Param`` filter drops ``var_ema``, which eval mode divides by."""
    params = jax.tree_util.tree_map(
        np.asarray, nnx.state(pretrained_model, nnx.Param).to_pure_dict()
    )
    out_dir = _write_artifacts(tmp_path, pretrained_model.config, params)

    with pytest.raises(ValueError, match=r"nnx\.Param-only filter"):
        load_pretrained_model(out_dir, use_bf16_compute=False)


def test_tracks_with_conditioning_tables_stripped_is_refused(tmp_path: Path) -> None:
    """Trap: dropping the species tables, as the annotation export does, fails."""
    model = _tracks_model()
    variables = _state_dict(model)
    for subtree in ("conditions_embed_layers", "conditions_heads"):
        variables.pop(subtree)
    out_dir = _write_artifacts(tmp_path, model.config, variables)

    with pytest.raises(ValueError, match="does not match this model's variable tree"):
        load_tracks_model(out_dir, use_bf16_compute=False)


def test_tracks_with_conditioning_dropped_everywhere_is_refused(
    tmp_path: Path,
) -> None:
    """Trap: config and tables blanked together still loads, and must not."""
    model = _tracks_model()
    variables = _state_dict(model)
    for subtree in ("conditions_embed_layers", "conditions_heads"):
        variables.pop(subtree)
    # The FiLM params sized by the condition count go with them.
    unconditioned = DiscreteConditionedNTv3WithTracksHeads(
        dataclasses.replace(
            model.config, conditions_vocab_size=[], conditions_names=[]
        ),
        rngs=nnx.Rngs(0),
    )
    variables["bigwig_head"] = _state_dict(unconditioned)["bigwig_head"]
    out_dir = _write_artifacts(tmp_path, unconditioned.config, variables)

    with pytest.raises(ValueError, match="dropped the species conditioning"):
        load_tracks_model(out_dir, use_bf16_compute=False)


def test_tracks_artifacts_refused_by_the_pretrained_loader(tmp_path: Path) -> None:
    """A tracks config never builds a bare backbone that discards its heads."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))

    with pytest.raises(ValueError, match="load_tracks_model"):
        load_pretrained_model(out_dir, use_bf16_compute=False)


def test_unexpected_variables_are_refused(
    tmp_path: Path, pretrained_model: NTv3
) -> None:
    """Tensors the model has no slot for fail loudly instead of being dropped."""
    variables = _state_dict(pretrained_model)
    variables["bigwig_head"] = {"kernel": np.zeros((2, 2), np.float32)}
    weights = tmp_path / WEIGHTS_FILE
    weights.write_bytes(serialization.to_bytes(variables))

    with pytest.raises(ValueError, match="no place for"):
        load_variables(pretrained_model, weights)


def test_misshapen_variables_are_refused(tmp_path: Path) -> None:
    """A 16-token embedding never lands in an 11-token model."""
    wide = NTv3(_pretrained_config(alphabet_size=16), rngs=nnx.Rngs(0))
    narrow = NTv3(_pretrained_config(alphabet_size=11), rngs=nnx.Rngs(0))
    weights = tmp_path / WEIGHTS_FILE
    weights.write_bytes(serialization.to_bytes(_state_dict(wide)))

    with pytest.raises(ValueError, match="at a shape this model does not expect"):
        load_variables(narrow, weights)


def test_no_bigwig_head_is_refused_by_the_tracks_loader(tmp_path: Path) -> None:
    """An annotation-only checkpoint is pointed at its own loader."""
    config = dataclasses.replace(
        _tracks_config(bed_elements_names=("exon",)),
        bigwigs_per_file_assembly={},
        number_bigwig_tracks=0,
    )
    model = DiscreteConditionedNTv3WithTracksHeads(config, rngs=nnx.Rngs(0))
    out_dir = _write_artifacts(tmp_path, config, _state_dict(model))

    with pytest.raises(ValueError, match="load_annotation_model"):
        load_tracks_model(out_dir, use_bf16_compute=False)


@pytest.mark.parametrize("alphabet_size", [11, 16])
def test_tokenizer_defaults_to_the_vocabulary_the_config_implies(
    tmp_path: Path, alphabet_size: int
) -> None:
    """With no tokenizer file, 11 tokens is plain and 16 is soft-masked.

    The plain vocabulary has no lowercase bases at all: soft-masked sequence
    tokenizes to ``<unk>``, so callers of the 11-token model must uppercase.
    """
    model = NTv3(_pretrained_config(alphabet_size), rngs=nnx.Rngs(0))
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    (out_dir / TOKENIZER_FILE).unlink()

    _, tokenizer = load_pretrained_model(out_dir, use_bf16_compute=False)

    assert tokenizer.vocabulary_size == alphabet_size
    upper, lower = tokenizer.tokenize("ACGT")[1], tokenizer.tokenize("acgt")[1]
    assert upper == [6, 8, 9, 7]
    assert lower == ([11, 13, 14, 12] if alphabet_size == 16 else [0, 0, 0, 0])


def test_mismatched_tokenizer_is_refused(
    tmp_path: Path, pretrained_model: NTv3
) -> None:
    """A 16-token tokenizer shipped with an 11-token model fails loudly."""
    out_dir = _write_artifacts(
        tmp_path,
        pretrained_model.config,
        _state_dict(pretrained_model),
        alphabet_size=16,
    )

    with pytest.raises(ValueError, match="Tokenizer has 16 tokens"):
        load_pretrained_model(out_dir, use_bf16_compute=False)


def test_tracks_without_conditions_tokenizers_is_refused(tmp_path: Path) -> None:
    """Without the species vocabulary no condition id can be trusted."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    (out_dir / CONDITIONS_TOKENIZERS_FILE).unlink()

    # Matched on the guard's own message: the OS error a bare read would raise
    # also names the file.
    with pytest.raises(FileNotFoundError, match=r"declares 1 discrete condition"):
        load_tracks_model(out_dir, use_bf16_compute=False)


@pytest.mark.parametrize(
    ("conditions", "message"),
    [
        ({"tissue": StandardTokenizer(SPECIES)}, "no tokenizer for condition"),
        ({"species": StandardTokenizer(["human"])}, "has 7 tokens"),
    ],
)
def test_conditions_tokenizers_must_match_the_config(
    tmp_path: Path, conditions: dict[str, StandardTokenizer], message: str
) -> None:
    """The file names each condition, at the embedding table's size."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    configs = {name: tokenizer.to_config() for name, tokenizer in conditions.items()}
    (out_dir / CONDITIONS_TOKENIZERS_FILE).write_text(
        json.dumps({"conditions_tokenizers": configs})
    )

    with pytest.raises(ValueError, match=message):
        load_tracks_model(out_dir, use_bf16_compute=False)


def test_local_directory_needs_no_huggingface_hub(
    tmp_path: Path, pretrained_model: NTv3, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Loading from disk never imports the Hub client, installed or not."""
    monkeypatch.setitem(sys.modules, "huggingface_hub", None)
    out_dir = _write_artifacts(
        tmp_path, pretrained_model.config, _state_dict(pretrained_model)
    )

    load_pretrained_model(out_dir, use_bf16_compute=False)


def test_unknown_name_is_refused(tmp_path: Path) -> None:
    """A name that is neither a directory nor a release fails before any download."""
    with pytest.raises(FileNotFoundError, match="NTv3_tracks"):
        load_tracks_model(tmp_path / "missing", use_bf16_compute=False)


def test_released_name_downloads_the_flax_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A released name fetches its Hub repo's Flax subfolder, and loads from it."""
    model = _tracks_model()
    _write_artifacts(tmp_path / ARTIFACTS_SUBFOLDER, model.config, _state_dict(model))
    calls: list[dict[str, Any]] = []

    def snapshot_download(repo_id: str, **kwargs: Any) -> str:
        calls.append({"repo_id": repo_id, **kwargs})
        return str(tmp_path)

    # A stand-in module, so this runs where huggingface_hub is not installed.
    hub = types.ModuleType("huggingface_hub")
    setattr(hub, "snapshot_download", snapshot_download)
    monkeypatch.setitem(sys.modules, "huggingface_hub", hub)

    reloaded, _, _ = load_tracks_model(
        "NTv3_tracks", use_bf16_compute=False, revision="abc123"
    )

    (call,) = calls
    assert call["repo_id"] == RELEASED_MODELS["NTv3_tracks"]
    assert call["revision"] == "abc123"
    # The Flax subfolder only, fixture included, so check_release can run.
    patterns = call["allow_patterns"]
    assert all(pattern.startswith(ARTIFACTS_SUBFOLDER) for pattern in patterns)
    assert any(pattern.endswith(".npz") for pattern in patterns)
    tokens, species = _tokens(16), _species(6, 7)
    _assert_same_outputs(model(tokens, species), reloaded(tokens, species))


def _annotation_model() -> DiscreteConditionedNTv3WithTracksHeads:
    """The released annotation topology, shrunk: a BED head and no conditions."""
    config = dataclasses.replace(
        _tracks_config(bed_elements_names=("plus",)),
        conditions_vocab_size=[],
        conditions_names=[],
        bigwigs_per_file_assembly={},
        number_bigwig_tracks=0,
        condition_bigwig_head=False,
        keep_target_center_fraction=0.75,
    )
    model = _randomized(
        DiscreteConditionedNTv3WithTracksHeads(config, rngs=nnx.Rngs(0))
    )
    model.eval()
    return model


def test_annotation_round_trip(tmp_path: Path) -> None:
    """An annotation export reloads bitwise, with no conditions tokenizers file."""
    model = _annotation_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    assert not (out_dir / CONDITIONS_TOKENIZERS_FILE).exists()

    reloaded, tokenizer = load_annotation_model(out_dir, use_bf16_compute=False)

    assert tokenizer.vocabulary_size == 16
    tokens = _tokens(16)
    _assert_same_outputs(model(tokens, []), reloaded(tokens, []))


def test_tracks_artifacts_refused_by_the_annotation_loader(tmp_path: Path) -> None:
    """A checkpoint with no BED head is pointed at its own loader."""
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))

    with pytest.raises(ValueError, match="load_tracks_model"):
        load_annotation_model(out_dir, use_bf16_compute=False)


def test_ignored_config_keys_are_named(
    tmp_path: Path, pretrained_model: NTv3, caplog: pytest.LogCaptureFixture
) -> None:
    """A misspelled field falls back to its default, so it must not go unseen."""
    out_dir = _write_artifacts(
        tmp_path, pretrained_model.config, _state_dict(pretrained_model)
    )
    values = json.loads((out_dir / CONFIG_FILE).read_text())
    values["num_layer"] = values.pop("num_layers")
    values["provenance"] = {"source_model_name": "tiny"}
    (out_dir / CONFIG_FILE).write_text(json.dumps(values))

    with caplog.at_level(logging.WARNING):
        # The default num_layers builds a deeper model than the weights hold.
        with pytest.raises(ValueError):
            load_pretrained_model(out_dir, use_bf16_compute=False)

    ignored = [
        r.getMessage() for r in caplog.records if "keys ignored" in r.getMessage()
    ]
    assert ignored == ["hyperparameters.json keys ignored by NTv3Config: ['num_layer']"]


def test_ones_batch_stats_warn_only_when_untrained(
    tmp_path: Path, pretrained_model: NTv3, caplog: pytest.LogCaptureFixture
) -> None:
    """An all-ones var_ema, what a dropped BatchStat leaves, is flagged."""
    fresh = NTv3(_pretrained_config(), rngs=nnx.Rngs(0))
    for name, model in (("fresh", fresh), ("trained", pretrained_model)):
        out_dir = _write_artifacts(tmp_path / name, model.config, _state_dict(model))
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            load_pretrained_model(out_dir, use_bf16_compute=False)
        warned = any(
            "batch-statistic buffers" in r.getMessage() for r in caplog.records
        )
        assert warned == (name == "fresh"), name


def _write_tracks_fixture(
    out_dir: Path, model: DiscreteConditionedNTv3WithTracksHeads
) -> None:
    """The export's fixture: condition id 0 with the condition masked out."""
    tokens = _tokens(16)
    rows = tokens.shape[0]
    outputs = model(
        tokens,
        [jnp.zeros((rows,), jnp.int32)],
        conditions_masks=[jnp.zeros((rows,), bool)],
    )
    np.savez(
        out_dir / FIXTURE_FILE,
        tokens=np.asarray(tokens),
        bigwig_tracks_logits=np.asarray(outputs["bigwig_tracks_logits"]),
    )


def test_check_release_passes_a_conditioned_tracks_model(tmp_path: Path) -> None:
    model = _tracks_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    _write_tracks_fixture(out_dir, model)

    result = check("tracks", str(out_dir), tolerance=1e-5)

    assert result["species_effect"] > 0.01


def test_check_release_refuses_dead_species_conditioning(tmp_path: Path) -> None:
    """The fixture masks the species out, so this check is the only one that sees it."""
    model = _tracks_model(condition_bigwig_head=False)
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    _write_tracks_fixture(out_dir, model)

    with pytest.raises(ValueError, match="does not reach the bigwig head"):
        check("tracks", str(out_dir), tolerance=1e-5)


def test_check_release_replays_an_annotation_fixture(tmp_path: Path) -> None:
    model = _annotation_model()
    out_dir = _write_artifacts(tmp_path, model.config, _state_dict(model))
    tokens = _tokens(16)
    np.savez(
        out_dir / FIXTURE_FILE,
        tokens=np.asarray(tokens),
        bed_tracks_logits=np.asarray(model(tokens, [])["bed_tracks_logits"]),
    )

    result = check("annotation", str(out_dir), tolerance=1e-5)

    assert result["argmax_agreement"] == 1.0
