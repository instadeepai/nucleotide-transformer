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

from typing import Any

import jax.numpy as jnp
import numpy as np
from flax import nnx
from transformers import AutoModel, AutoModelForMaskedLM

from nucleotide_transformer_v3.model import (
    NTv3PreTrained,
    NTv3PreTrainedConfig,
    DiscreteConditionedNTv3PreTrainedConfig,
    NTv3PostTrained,
    NTv3PostTrainedConfig,
)
from nucleotide_transformer_v3.layers import DeConvUpsampleType
from nucleotide_transformer_v3.tokenizers import StandardTokenizer, get_ntv3_tokenizer, NucleotideTokenizer


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

def parse_hf_config_to_ntv3_config(
    config_dict: dict[str, Any],
    tokenizer: NucleotideTokenizer,
) -> NTv3PreTrainedConfig:
    """
    Parse Hugging Face config.json into NTv3PreTrainedConfig.

    Args:
        config_dict: Dictionary loaded from config.json.
        tokenizer: Tokenizer instance to get pad_token_id and mask_token_id.

    Returns:
        NTv3PreTrainedConfig instance.
    """
    # Handle deconv_upsample_type conversion
    deconv_type_str = config_dict.get("deconv_upsample_type", "conv_transpose")
    deconv_upsample_type = (
        DeConvUpsampleType.REPEAT_CONV
        if deconv_type_str == "repeat+conv"
        else DeConvUpsampleType.CONV_TRANSPOSE
    )

    return NTv3PreTrainedConfig(
        alphabet_size=config_dict.get(
            "alphabet_size",
            config_dict.get("vocab_size", tokenizer.vocabulary_size),
        ),
        pad_token_id=config_dict.get("pad_token_id", tokenizer.pad_token_id),
        mask_token_id=config_dict.get("mask_token_id", tokenizer.mask_token_id),
        num_downsamples=config_dict.get("num_downsamples", 8),
        attention_heads=config_dict.get(
            "attention_heads",
            config_dict.get("num_attention_heads", 16),
        ),
        key_size=config_dict.get("key_size", None),
        token_embed_dim=config_dict.get("token_embed_dim", 16),
        conv_init_embed_dim=config_dict.get("conv_init_embed_dim", 512),
        embed_dim=config_dict.get("embed_dim", config_dict.get("hidden_size", 512)),
        ffn_embed_dim=config_dict.get(
            "ffn_embed_dim",
            config_dict.get("intermediate_size", 2048),
        ),
        num_layers=config_dict.get(
            "num_layers",
            config_dict.get("num_hidden_layers", 12),
        ),
        layer_norm_eps=config_dict.get("layer_norm_eps", 1e-5),
        num_hidden_layers_head=config_dict.get("num_hidden_layers_head", 0),
        use_skip_connection=config_dict.get("use_skip_connection", True),
        deconv_upsample_type=deconv_upsample_type,
    )


def parse_hf_config_to_posttrained_config(
    config_dict: dict[str, Any],
    tokenizer: NucleotideTokenizer,
) -> NTv3PostTrainedConfig:
    """
    Parse Hugging Face config.json into NTv3PostTrainedConfig.

    Post-trained models add species/assembly conditioning and prediction heads.

    Args:
        config_dict: Dictionary loaded from config.json.
        tokenizer: Tokenizer instance to get pad_token_id and mask_token_id.

    Returns:
        NTv3PostTrainedConfig instance.
    """
    # Parse base config
    base_config = parse_hf_config_to_ntv3_config(config_dict, tokenizer)
    # Extract species/assembly information
    # Post-trained models have ONE condition: species/assembly with each assembly as a token
    conditions_names = config_dict.get("conditions_names", None)
    conditions_vocab_size = config_dict.get("conditions_vocab_size", None)

    # Create conditioned config with species condition
    conditioned_config = DiscreteConditionedNTv3PreTrainedConfig.from_base_config(
        config=base_config,
        conditions_vocab_size=conditions_vocab_size,
        conditions_names=conditions_names,
    )

    # Extract post-trained specific fields
    # Note: These fields are specific to post-trained models and may need
    # refactoring if the config structure changes
    bigwigs_per_species: dict[str, list[str]] = config_dict.get(
        "bigwigs_per_species", {}
    )
    bed_elements_names: list[str] = config_dict.get("bed_elements_names", [])
    keep_target_center_fraction = config_dict.get("keep_target_center_fraction", 1.0)
    species_to_token_id = config_dict.get("species_to_token_id", {})
    num_species_special_tokens = config_dict.get("num_species_special_tokens", 6)

    return NTv3PostTrainedConfig.from_base_config(
        config=conditioned_config,
        bigwigs_per_species=bigwigs_per_species,
        bed_elements_names=bed_elements_names,
        species_to_token_id=species_to_token_id,
        keep_target_center_fraction=keep_target_center_fraction,
        num_species_special_tokens=num_species_special_tokens,
    )


def _convert_conv_weight(weight: np.ndarray) -> np.ndarray:
    """Convert PyTorch Conv1d/ConvTranspose1d weight to JAX format."""
    # PyTorch: [out_channels, in_channels, kernel_size]
    # JAX: [kernel_size, in_channels, out_channels]
    return weight.transpose(2, 1, 0)


def _load_layer_norm_with_modulation(
    model_ln: Any,
    state_ln: Any,
    tensors: dict[str, np.ndarray],
    key_prefix: str,
) -> None:
    """
    Load layer norm weights, with optional modulation layers for post-trained models.
    
    Pre-trained models: Only load scale and bias.
    Post-trained models: Also load modulation layers that condition on species/assembly.
    """
    # Load basic layer norm parameters
    state_ln.scale.value = jnp.array(tensors[f"{key_prefix}.weight"])
    state_ln.bias.value = jnp.array(tensors[f"{key_prefix}.bias"])
    
    # Load modulation layers if they exist (post-trained models only)
    if hasattr(model_ln, "modulation_layers"):
        for j in range(len(model_ln.modulation_layers)):
            k = tensors[f"{key_prefix}.modulation_layers.{j}.weight"]
            state_ln.modulation_layers[j].kernel.value = jnp.array(k.T)
            state_ln.modulation_layers[j].bias.value = jnp.array(
                tensors[f"{key_prefix}.modulation_layers.{j}.bias"]
            )


def _load_conv_block(
    model_block: Any,
    state_block: Any,
    tensors: dict[str, np.ndarray],
    key_prefix: str,
) -> None:
    """
    Load convolutional block: main conv + layer norm + residual conv + layer norm.
    
    Handles both pre-trained (no modulation) and post-trained (with modulation) models.
    """
    # Main conv
    k = tensors[f"{key_prefix}.conv.conv.weight"]
    state_block.conv.conv.kernel.value = jnp.array(_convert_conv_weight(k))
    state_block.conv.conv.bias.value = jnp.array(tensors[f"{key_prefix}.conv.conv.bias"])
    
    # Layer norm in main conv (with optional modulation)
    _load_layer_norm_with_modulation(
        model_block.conv.layer_norm,
        state_block.conv.layer_norm,
        tensors,
        f"{key_prefix}.conv.layer_norm",
    )
    
    # Residual conv
    k = tensors[f"{key_prefix}.res_conv.conv_block.conv.weight"]
    state_block.res_conv.conv_block.conv.kernel.value = jnp.array(_convert_conv_weight(k))
    state_block.res_conv.conv_block.conv.bias.value = jnp.array(
        tensors[f"{key_prefix}.res_conv.conv_block.conv.bias"]
    )
    
    # Layer norm in residual conv (with optional modulation)
    _load_layer_norm_with_modulation(
        model_block.res_conv.conv_block.layer_norm,
        state_block.res_conv.conv_block.layer_norm,
        tensors,
        f"{key_prefix}.res_conv.conv_block.layer_norm",
    )
    
    # Residual block modulation layers (post-trained only)
    if hasattr(model_block.res_conv, "modulation_layers"):
        for j in range(len(model_block.res_conv.modulation_layers)):
            k = tensors[f"{key_prefix}.res_conv.modulation_layers.{j}.weight"]
            state_block.res_conv.modulation_layers[j].kernel.value = jnp.array(k.T)
            state_block.res_conv.modulation_layers[j].bias.value = jnp.array(
                tensors[f"{key_prefix}.res_conv.modulation_layers.{j}.bias"]
            )


def _load_transformer_block(
    model_block: Any,
    state_block: Any,
    tensors: dict[str, np.ndarray],
    key_prefix: str,
) -> None:
    """
    Load transformer block: layer norms + attention + FFN.
    
    Handles both pre-trained (no modulation) and post-trained (with modulation) models.
    """
    # Self-attention layer norm (with optional modulation)
    _load_layer_norm_with_modulation(
        model_block.self_attention_layer_norm,
        state_block.self_attention_layer_norm,
        tensors,
        f"{key_prefix}.self_attention_layer_norm",
    )
    
    # Final layer norm (with optional modulation)
    _load_layer_norm_with_modulation(
        model_block.final_layer_norm,
        state_block.final_layer_norm,
        tensors,
        f"{key_prefix}.final_layer_norm",
    )
    
    # Block-level modulation layers (post-trained only)
    if hasattr(model_block, "modulation_layers"):
        for j in range(len(model_block.modulation_layers)):
            k = tensors[f"{key_prefix}.modulation_layers.{j}.weight"]
            state_block.modulation_layers[j].kernel.value = jnp.array(k.T)
            state_block.modulation_layers[j].bias.value = jnp.array(
                tensors[f"{key_prefix}.modulation_layers.{j}.bias"]
            )
    
    # Attention heads (query, key, value)
    for head_name in ["query_head", "key_head", "value_head"]:
        k = tensors[f"{key_prefix}.sa_layer.{head_name}.linear.weight"]
        getattr(state_block.sa_layer, head_name).linear.kernel.value = jnp.array(k.T)
        getattr(state_block.sa_layer, head_name).linear.bias.value = jnp.array(
            tensors[f"{key_prefix}.sa_layer.{head_name}.linear.bias"]
        )
    
    # Multi-head attention output
    k = tensors[f"{key_prefix}.sa_layer.mha_output.weight"]
    state_block.sa_layer.mha_output.kernel.value = jnp.array(k.T)
    state_block.sa_layer.mha_output.bias.value = jnp.array(
        tensors[f"{key_prefix}.sa_layer.mha_output.bias"]
    )
    
    # Feed-forward network
    state_block.fc1.kernel.value = jnp.array(tensors[f"{key_prefix}.fc1.weight"].T)
    state_block.fc2.kernel.value = jnp.array(tensors[f"{key_prefix}.fc2.weight"].T)

def convert_pytorch_to_jax_params(
    pytorch_model: Any,
    model: NTv3PreTrained | NTv3PostTrained,
) -> None:
    """
    Convert PyTorch model state_dict to JAX params and load into model.
    
    Handles both pre-trained and post-trained models:
    - Pre-trained: Base architecture only
    - Post-trained: Base + species conditioning + bigwig/bed prediction heads
    
    Args:
        pytorch_model: The PyTorch model loaded from Hugging Face.
        model: The JAX model instance (NTv3PreTrained or NTv3PostTrained).
    """
    # Extract and convert state dict
    state_dict = pytorch_model.state_dict()
    tensors = {k: v.cpu().numpy() for k, v in state_dict.items()}
    
    # All models now use "core." prefix on HuggingFace
    prefix = "core."
    
    # Get model's parameter state
    state = nnx.state(model, nnx.Param)
    
    # Embedding layer
    state.embed_layer.embedding.value = jnp.array(tensors[f"{prefix}embed_layer.weight"])
    
    # Stem (initial convolution)
    k = tensors[f"{prefix}stem.conv.weight"]
    state.stem.conv.kernel.value = jnp.array(_convert_conv_weight(k))
    state.stem.conv.bias.value = jnp.array(tensors[f"{prefix}stem.conv.bias"])
    
    # Convolutional tower (downsampling)
    for i, block in enumerate(model.conv_tower_blocks):
        _load_conv_block(block, state.conv_tower_blocks[i], tensors, f"{prefix}conv_tower_blocks.{i}")
    
    # Deconvolutional tower (upsampling)
    for i, block in enumerate(model.deconv_tower_blocks):
        _load_conv_block(block, state.deconv_tower_blocks[i], tensors, f"{prefix}deconv_tower_blocks.{i}")
    
    # Transformer blocks
    for i, block in enumerate(model.transformer_blocks):
        _load_transformer_block(block, state.transformer_blocks[i], tensors, f"{prefix}transformer_blocks.{i}")
    
    # Language modeling head
    k = tensors[f"{prefix}lm_head.head.weight"]
    state.lm_head.head.kernel.value = jnp.array(k.T)
    state.lm_head.head.bias.value = jnp.array(tensors[f"{prefix}lm_head.head.bias"])
    
    print("  ✓ Loaded base model weights")
    

    if isinstance(model, NTv3PostTrained):
        # Species/assembly embeddings
        num_conditions = len(model.conditions_embed_layers)
        for i in range(num_conditions):
            state.conditions_embed_layers[i].embedding.value = jnp.array(
                tensors[f"{prefix}cond_tables.{i}.weight"]
            )
        
        # Species condition heads
        for i in range(num_conditions):
            k = tensors[f"{prefix}conditions_heads.{i}.weight"]
            state.conditions_heads[i].kernel.value = jnp.array(k.T)
            state.conditions_heads[i].bias.value = jnp.array(
                tensors[f"{prefix}conditions_heads.{i}.bias"]
            )
        
        print(f"  ✓ Loaded {num_conditions} species condition embedding(s)")
        
        # Bigwig head (per-species prediction heads for histone marks, etc.)
        if hasattr(model, "bigwig_head"):
            num_species_with_bigwigs = 0
            for i, head in enumerate(model.bigwig_head.species_heads):
                if head is not None:
                    head_key = f"{prefix}bigwig_head.species_heads.{i}.layer_norm.weight"
                    if head_key in tensors:
                        # Layer norm
                        state.bigwig_head.species_heads[i].layer_norm.scale.value = jnp.array(
                            tensors[f"{prefix}bigwig_head.species_heads.{i}.layer_norm.weight"]
                        )
                        state.bigwig_head.species_heads[i].layer_norm.bias.value = jnp.array(
                            tensors[f"{prefix}bigwig_head.species_heads.{i}.layer_norm.bias"]
                        )
                        # Linear head
                        k = tensors[f"{prefix}bigwig_head.species_heads.{i}.head.weight"]
                        state.bigwig_head.species_heads[i].head.kernel.value = jnp.array(k.T)
                        state.bigwig_head.species_heads[i].head.bias.value = jnp.array(
                            tensors[f"{prefix}bigwig_head.species_heads.{i}.head.bias"]
                        )
                        num_species_with_bigwigs += 1
            
            print(f"  ✓ Loaded bigwig prediction heads for {num_species_with_bigwigs} species")
        
        # BED head (genomic element classification)
        if hasattr(model, "bed_head"):
            # Layer norm
            state.bed_head.layer_norm.scale.value = jnp.array(
                tensors[f"{prefix}bed_head.layer_norm.weight"]
            )
            state.bed_head.layer_norm.bias.value = jnp.array(
                tensors[f"{prefix}bed_head.layer_norm.bias"]
            )
            # Linear head
            k = tensors[f"{prefix}bed_head.head.weight"]
            state.bed_head.head.kernel.value = jnp.array(k.T)
            state.bed_head.head.bias.value = jnp.array(tensors[f"{prefix}bed_head.head.bias"])
            
            num_bed_elements = len(model.config.bed_elements_names)
            print(f"  ✓ Loaded BED element classification head ({num_bed_elements} elements)")
    
    # Update the model with loaded weights
    nnx.update(model, state)
    print(f"✓ Model loaded successfully\n")


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
    
    # Load PyTorch model and config
    pytorch_model = AutoModelForMaskedLM.from_pretrained(
        f"InstaDeepAI/{model_name}", trust_remote_code=True
    )
    config_dict = pytorch_model.config.to_dict()

    # Get tokenizer
    tokenizer = get_ntv3_tokenizer()

    # Parse config
    config = parse_hf_config_to_ntv3_config(config_dict, tokenizer)

    # Update config with layer save options
    config.embeddings_layers_to_save = embeddings_layers_to_save
    config.attention_maps_to_save = list(attention_maps_to_save)
    config.deconv_layers_to_save = deconv_layers_to_save

    if use_bfloat16:
        set_bfloat16_dtypes(config)

    # Create and load model
    rngs = nnx.Rngs(0)
    model = NTv3PreTrained(config=config, rngs=rngs)
    convert_pytorch_to_jax_params(pytorch_model, model)

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

    # Load PyTorch model and config
    pytorch_model = AutoModel.from_pretrained(
        f"InstaDeepAI/{model_name}", trust_remote_code=True
    )
    config_dict = pytorch_model.config.to_dict()
    
    tokenizer = get_ntv3_tokenizer()
    
    # Parse config
    config = parse_hf_config_to_posttrained_config(config_dict, tokenizer)

    # Update config with layer save options
    config.embeddings_layers_to_save = embeddings_layers_to_save
    config.attention_maps_to_save = list(attention_maps_to_save)
    config.deconv_layers_to_save = deconv_layers_to_save

    if use_bfloat16:
        set_bfloat16_dtypes(config)

    # Create and load model
    rngs = nnx.Rngs(0)
    model = NTv3PostTrained(config=config, rngs=rngs)
    convert_pytorch_to_jax_params(pytorch_model, model)

    return model, tokenizer, config

