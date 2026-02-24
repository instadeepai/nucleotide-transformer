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

import dataclasses
import json
from typing import Any

from flax import nnx, serialization
from huggingface_hub import hf_hub_download

from nucleotide_transformer_v3.model import (
    NTv3PreTrained,
    NTv3PreTrainedConfig,
    NTv3PostTrained,
    NTv3PostTrainedConfig,
)
from nucleotide_transformer_v3.layers import DeConvUpsampleType
from nucleotide_transformer_v3.tokenizers import get_ntv3_tokenizer, NucleotideTokenizer


SUPPORTED_PRETRAINED_MODEL_LIST = [
    # Production models
    "NTv3_8M_pre",
    "NTv3_100M_pre",
    "NTv3_650M_pre",
    # 8kb context intermediate checkpoints
    "NTv3_8M_pre_8kb",
    "NTv3_100M_pre_8kb",
    "NTv3_650M_pre_8kb",
    # 5 downsample ablations
    "NTv3_5downsample_pre",
    "NTv3_5downsample_pre_8kb",
]

SUPPORTED_POSTTRAINED_MODEL_LIST = [
    # Production models
    "NTv3_100M_post",
    "NTv3_650M_post",
    # 131kb context checkpoints
    "NTv3_100M_post_131kb",
    "NTv3_650M_post_131kb",
    # 5 downsample ablations
    "NTv3_5downsample_post",
    "NTv3_5downsample_post_131kb",
]

def set_bfloat16_dtypes(cfg: NTv3PreTrainedConfig | NTv3PostTrainedConfig) -> None:
    """Mutate Flax config in-place: all *compute* and *param* dtypes -> bfloat16."""
    cfg.embedding_param_dtype = "bfloat16"
    cfg.embedding_compute_dtype = "bfloat16"
    cfg.stem_param_dtype = "bfloat16"
    cfg.stem_compute_dtype = "bfloat16"
    cfg.down_convolution_param_dtype = "bfloat16"
    cfg.down_convolution_compute_dtype = "bfloat16"
    cfg.up_convolution_param_dtype = "bfloat16"
    cfg.up_convolution_compute_dtype = "bfloat16"
    cfg.layernorm_param_dtype = "bfloat16"
    cfg.layernorm_compute_dtype = "bfloat16"
    cfg.transformer_qkvo_param_dtype = "bfloat16"
    cfg.transformer_qkvo_compute_dtype = "bfloat16"
    cfg.transformer_ffn_param_dtype = "bfloat16"
    cfg.transformer_ffn_compute_dtype = "bfloat16"
    cfg.lmhead_param_dtype = "bfloat16"
    cfg.lmhead_compute_dtype = "bfloat16"
    cfg.modulation_param_dtype = "bfloat16"
    cfg.modulation_compute_dtype = "bfloat16"


def _parse_jax_hyperparams(
    cfg_dict: dict[str, Any],
    config_cls: type[NTv3PreTrainedConfig] = NTv3PreTrainedConfig,
) -> NTv3PreTrainedConfig | NTv3PostTrainedConfig:
    """Parse a hyperparameters.json dict into a config dataclass.

    Handles type conversions that JSON cannot represent natively (enums, tuples)
    and filters out any keys not declared on *config_cls*.
    """
    d = dict(cfg_dict)
    if "deconv_upsample_type" in d:
        d["deconv_upsample_type"] = DeConvUpsampleType(d["deconv_upsample_type"])
    if "embeddings_layers_to_save" in d:
        d["embeddings_layers_to_save"] = tuple(d["embeddings_layers_to_save"])
    if "deconv_layers_to_save" in d:
        d["deconv_layers_to_save"] = tuple(d["deconv_layers_to_save"])
    valid_fields = {f.name for f in dataclasses.fields(config_cls)}
    d = {k: v for k, v in d.items() if k in valid_fields}
    return config_cls(**d)


def _load_jax_weights(model: NTv3PreTrained | NTv3PostTrained, weights_path: str) -> None:
    """Load Flax NNX weights from a msgpack file into *model*."""
    with open(weights_path, "rb") as f:
        weights_bytes = f.read()
    target = nnx.state(model, nnx.Param).to_pure_dict()
    restored = serialization.from_bytes(target, weights_bytes)
    nnx.update(model, restored)


def get_pretrained_ntv3_model(
    model_name: str,
    embeddings_layers_to_save: tuple[int, ...] = (),
    attention_maps_to_save: tuple[tuple[int, int], ...] = (),
    deconv_layers_to_save: tuple[int, ...] = (),
    use_bfloat16: bool = False,
) -> tuple[NTv3PreTrained, NucleotideTokenizer, NTv3PreTrainedConfig]:
    """
    Download and load a pre-trained NTv3 model from Hugging Face.
    
    Pre-trained models are trained on DNA sequences for masked language modeling.

    Args:
        model_name: Model name (e.g., "NTv3_8M_pre").
        embeddings_layers_to_save: Transformer layer indices to save embeddings for.
        attention_maps_to_save: Attention maps to save as (layer, head) tuples.
        deconv_layers_to_save: Deconv layer indices to save embeddings for.
        use_bfloat16: Whether to use bfloat16 for the model.

    Returns:
        Tuple of (model, tokenizer, config).

    Example:
        >>> model, tokenizer, config = get_pretrained_ntv3_model(
        ...     model_name="NTv3_8M_pre",
        ...     embeddings_layers_to_save=(5, 10),
        ... )
    """
    assert model_name in SUPPORTED_PRETRAINED_MODEL_LIST, f"Model {model_name} not supported"

    repo_id = f"InstaDeepAI/{model_name}"
    cfg_path = hf_hub_download(repo_id, "jax_model/hyperparameters.json")
    wt_path = hf_hub_download(repo_id, "jax_model/flax_model.msgpack")

    with open(cfg_path) as f:
        config = _parse_jax_hyperparams(json.load(f))

    config.embeddings_layers_to_save = embeddings_layers_to_save
    config.attention_maps_to_save = list(attention_maps_to_save)
    config.deconv_layers_to_save = deconv_layers_to_save

    if use_bfloat16:
        set_bfloat16_dtypes(config)

    tokenizer = get_ntv3_tokenizer()

    model = NTv3PreTrained(config=config, rngs=nnx.Rngs(0))
    _load_jax_weights(model, wt_path)

    return model, tokenizer, config


def get_posttrained_ntv3_model(
    model_name: str,
    embeddings_layers_to_save: tuple[int, ...] = (),
    attention_maps_to_save: tuple[tuple[int, int], ...] = (),
    deconv_layers_to_save: tuple[int, ...] = (),
    use_bfloat16: bool = False,
) -> tuple[NTv3PostTrained, NucleotideTokenizer, NTv3PostTrainedConfig]:
    """
    Download and load a post-trained NTv3 model from Hugging Face.

    Post-trained models extend pre-trained NTv3 with:
    - Species/assembly conditioning via adaptive layer normalization
    - Bigwig prediction heads (histone modifications, chromatin accessibility, etc.)
    - BED element classification head (regulatory element prediction)

    Args:
        model_name: Model name (e.g., "NTv3_100M_post").
        embeddings_layers_to_save: Transformer layer indices to save embeddings for.
        attention_maps_to_save: Attention maps to save as (layer, head) tuples.
        deconv_layers_to_save: Deconv layer indices to save embeddings for.
        use_bfloat16: Whether to use bfloat16 for the model.

    Returns:
        Tuple of (model, tokenizer, config).

    Example:
        >>> model, tokenizer, config = get_posttrained_ntv3_model(
        ...     model_name="NTv3_100M_post",
        ... )
        >>> species_tokens = model.encode_species("human")
    """
    assert model_name in SUPPORTED_POSTTRAINED_MODEL_LIST, f"Model {model_name} not supported"

    repo_id = f"InstaDeepAI/{model_name}"
    cfg_path = hf_hub_download(repo_id, "jax_model/hyperparameters.json")
    wt_path = hf_hub_download(repo_id, "jax_model/flax_model.msgpack")

    with open(cfg_path) as f:
        config = _parse_jax_hyperparams(json.load(f), config_cls=NTv3PostTrainedConfig)

    config.embeddings_layers_to_save = embeddings_layers_to_save
    config.attention_maps_to_save = list(attention_maps_to_save)
    config.deconv_layers_to_save = deconv_layers_to_save

    if use_bfloat16:
        set_bfloat16_dtypes(config)

    tokenizer = get_ntv3_tokenizer()

    model = NTv3PostTrained(config=config, rngs=nnx.Rngs(0))
    _load_jax_weights(model, wt_path)

    return model, tokenizer, config