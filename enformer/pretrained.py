"""Implementation of utilities to load a pretrained Enformer model in Trix."""

import gc
from collections import defaultdict
from typing import Callable, Dict, Tuple

import haiku as hk
import jax.numpy as jnp
import numpy as np
import torch
from enformer_pytorch import Enformer

from enformer.features import FEATURES
from enformer.heads import UNetHead
from enformer.model import (
    EnformerConfig,
    build_enformer_fn,
    build_enformer_with_head_fn,
)
from enformer.params import download_ckpt
from enformer.tokenizer import NucleotidesKmersTokenizer

__all__ = [
    "build_enformer_fn",
    "get_pretrained_segment_enformer_model",
    "ENFORMER_MODEL_NAME",
]


ENFORMER_MODEL_NAME = "Enformer"


def translate_weights(
    torch_params: torch.ParameterDict,
    num_conv_tower_layers: int,
    num_transformer_tower_layers: int,
) -> hk.Params:

    translate_dict: Dict[str, Tuple[str, str]] = {}

    def _translate_stem_weights(
        translate_dict: Dict[str, Tuple[str, str]]
    ) -> Dict[str, Tuple[str, str]]:
        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_stem/"
        translate_dict["stem.0.weight"] = (
            trix_prefix + "conv1_d",
            "w",
        )
        translate_dict["stem.0.bias"] = (
            trix_prefix + "conv1_d",
            "b",
        )
        translate_dict["stem.1.fn.0.weight"] = (
            trix_prefix + "residual_conv_block/conv_block/batch_norm",
            "scale",
        )
        translate_dict["stem.1.fn.0.bias"] = (
            trix_prefix + "residual_conv_block/conv_block/batch_norm",
            "offset",
        )
        translate_dict["stem.1.fn.2.weight"] = (
            trix_prefix + "residual_conv_block/conv_block/conv1_d",
            "w",
        )
        translate_dict["stem.1.fn.2.bias"] = (
            trix_prefix + "residual_conv_block/conv_block/conv1_d",
            "b",
        )
        translate_dict["stem.2.to_attn_logits.weight"] = (
            trix_prefix + "attention_pool/~/conv2_d",
            "w",
        )
        return translate_dict

    def _translate_conv_tower_weights(
        translate_dict: Dict[str, Tuple[str, str]],
        layer_num: int,
    ) -> Dict[str, Tuple[str, str]]:

        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_conv_tower/layer_{layer_num}/"
        torch_prefix = f"conv_tower.{layer_num}."

        translate_dict[torch_prefix + "0.0.weight"] = (
            trix_prefix + "conv_block/batch_norm",
            "scale",
        )
        translate_dict[torch_prefix + "0.0.bias"] = (
            trix_prefix + "conv_block/batch_norm",
            "offset",
        )

        translate_dict[torch_prefix + "0.2.weight"] = (
            trix_prefix + "conv_block/conv1_d",
            "w",
        )
        translate_dict[torch_prefix + "0.2.bias"] = (
            trix_prefix + "conv_block/conv1_d",
            "b",
        )

        translate_dict[torch_prefix + "1.fn.0.weight"] = (
            trix_prefix + "residual_conv_block/conv_block/batch_norm",
            "scale",
        )
        translate_dict[torch_prefix + "1.fn.0.bias"] = (
            trix_prefix + "residual_conv_block/conv_block/batch_norm",
            "offset",
        )

        translate_dict[torch_prefix + "1.fn.2.weight"] = (
            trix_prefix + "residual_conv_block/conv_block/conv1_d",
            "w",
        )
        translate_dict[torch_prefix + "1.fn.2.bias"] = (
            trix_prefix + "residual_conv_block/conv_block/conv1_d",
            "b",
        )
        translate_dict[torch_prefix + "2.to_attn_logits.weight"] = (
            trix_prefix + "attention_pool/~/conv2_d",
            "w",
        )
        return translate_dict

    def _translate_transformer_tower_weights(
        translate_dict: Dict[str, Tuple[str, str]],
        layer_num: int,
    ) -> Dict[str, Tuple[str, str]]:

        trix_prefix = (
            f"{ENFORMER_MODEL_NAME}/~_transformer_attention/layer_{layer_num}/"
        )
        torch_prefix = f"transformer.{layer_num}."

        translate_dict[torch_prefix + "0.fn.1.rel_content_bias"] = (
            trix_prefix + "attention",
            "rel_content_bias",
        )
        translate_dict[torch_prefix + "0.fn.1.rel_pos_bias"] = (
            trix_prefix + "attention",
            "rel_pos_bias",
        )
        translate_dict[torch_prefix + "0.fn.1.to_k.weight"] = (
            trix_prefix + "attention/~/attn_k",
            "w",
        )
        translate_dict[torch_prefix + "0.fn.1.to_q.weight"] = (
            trix_prefix + "attention/~/attn_q",
            "w",
        )
        translate_dict[torch_prefix + "0.fn.1.to_v.weight"] = (
            trix_prefix + "attention/~/attn_v",
            "w",
        )
        translate_dict[torch_prefix + "0.fn.1.to_rel_k.weight"] = (
            trix_prefix + "attention/~/attn_to_rel_k",
            "w",
        )
        translate_dict[torch_prefix + "0.fn.1.to_out.weight"] = (
            trix_prefix + "attention/~/attn_o",
            "w",
        )
        translate_dict[torch_prefix + "0.fn.1.to_out.bias"] = (
            trix_prefix + "attention/~/attn_o",
            "b",
        )
        translate_dict[torch_prefix + "0.fn.0.weight"] = (
            trix_prefix + "attn_layer_norm",
            "scale",
        )
        translate_dict[torch_prefix + "0.fn.0.bias"] = (
            trix_prefix + "attn_layer_norm",
            "offset",
        )

        trix_prefix = (
            f"{ENFORMER_MODEL_NAME}/~_transformer_ffn_block/layer_{layer_num}/"
        )
        translate_dict[torch_prefix + "1.fn.0.weight"] = (
            trix_prefix + "ffn_layer_norm",
            "scale",
        )
        translate_dict[torch_prefix + "1.fn.0.bias"] = (
            trix_prefix + "ffn_layer_norm",
            "offset",
        )
        translate_dict[torch_prefix + "1.fn.1.weight"] = (
            trix_prefix + "ffn_1",
            "w",
        )
        translate_dict[torch_prefix + "1.fn.1.bias"] = (
            trix_prefix + "ffn_1",
            "b",
        )
        translate_dict[torch_prefix + "1.fn.4.weight"] = (
            trix_prefix + "ffn_2",
            "w",
        )
        translate_dict[torch_prefix + "1.fn.4.bias"] = (
            trix_prefix + "ffn_2",
            "b",
        )

        return translate_dict

    def _translate_final_pointwise_weights(
        translate_dict: Dict[str, Tuple[str, str]]
    ) -> Dict[str, Tuple[str, str]]:

        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_final_pointwise/"
        torch_prefix = "final_pointwise.1."

        translate_dict[torch_prefix + "0.weight"] = (
            trix_prefix + "conv_block/batch_norm",
            "scale",
        )
        translate_dict[torch_prefix + "0.bias"] = (
            trix_prefix + "conv_block/batch_norm",
            "offset",
        )

        translate_dict[torch_prefix + "2.weight"] = (
            trix_prefix + "conv_block/conv1_d",
            "w",
        )
        translate_dict[torch_prefix + "2.bias"] = (
            trix_prefix + "conv_block/conv1_d",
            "b",
        )

        return translate_dict

    def _translate_heads_weights(
        translate_dict: Dict[str, Tuple[str, str]]
    ) -> Dict[str, Tuple[str, str]]:

        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_heads/"
        torch_prefix = "_heads."

        translate_dict[torch_prefix + "human.0.weight"] = (
            trix_prefix + "human_head",
            "w",
        )
        translate_dict[torch_prefix + "human.0.bias"] = (
            trix_prefix + "human_head",
            "b",
        )

        translate_dict[torch_prefix + "mouse.0.weight"] = (
            trix_prefix + "mouse_head",
            "w",
        )
        translate_dict[torch_prefix + "mouse.0.bias"] = (
            trix_prefix + "mouse_head",
            "b",
        )
        return translate_dict

    translate_dict = _translate_stem_weights(translate_dict)

    for i in range(num_conv_tower_layers):
        translate_dict = _translate_conv_tower_weights(translate_dict, layer_num=i)

    for i in range(num_transformer_tower_layers):
        translate_dict = _translate_transformer_tower_weights(
            translate_dict, layer_num=i
        )

    translate_dict = _translate_final_pointwise_weights(translate_dict)

    translate_dict = _translate_heads_weights(translate_dict)

    params: Dict[str, Dict[str, np.ndarray]] = defaultdict(dict)
    for torch_key, (trix_key, weight_key) in translate_dict.items():
        if "weight" in torch_key and not ("embed" in torch_key):
            # in pytorch, the weights of dense matrices indexation is transposes
            # compared to haiku, except for word-token-embedding
            params[trix_key][weight_key] = np.array(
                torch_params[torch_key], dtype=np.float32
            ).transpose()
        else:
            params[trix_key][weight_key] = np.array(
                torch_params[torch_key], dtype=np.float32
            )

        if "conv1_d" in trix_key and weight_key == "b":
            params[trix_key][weight_key] = np.expand_dims(
                params[trix_key][weight_key], axis=-1
            )

        if "batch_norm" in trix_key:
            params[trix_key][weight_key] = np.expand_dims(
                np.expand_dims(params[trix_key][weight_key], axis=-1), axis=0
            )

    return params  # type: ignore


def translate_state(
    torch_params: torch.ParameterDict,
    num_conv_tower_layers: int,
) -> hk.State:

    state: Dict[str, Dict[str, np.ndarray]] = {}

    def _get_one_state(
        state: Dict[str, Dict[str, np.ndarray]], trix_prefix: str, torch_prefix: str
    ) -> Dict[str, Dict[str, np.ndarray]]:
        for stat_type in ["mean", "var"]:
            state[trix_prefix + f"{stat_type}_ema"] = {
                "counter": np.array(
                    torch_params[torch_prefix + "num_batches_tracked"], dtype=np.int32
                ),
                "hidden": np.expand_dims(
                    np.expand_dims(
                        np.array(
                            torch_params[torch_prefix + f"running_{stat_type}"],
                            dtype=np.float32,
                        ),
                        axis=-1,
                    ),
                    axis=0,
                ),
                "average": np.expand_dims(
                    np.expand_dims(
                        np.array(
                            torch_params[torch_prefix + f"running_{stat_type}"],
                            dtype=np.float32,
                        ),
                        axis=-1,
                    ),
                    axis=0,
                ),
            }
        return state

    def _get_stem_state(
        state: Dict[str, Dict[str, np.ndarray]]
    ) -> Dict[str, Dict[str, np.ndarray]]:

        trix_prefix = (
            f"{ENFORMER_MODEL_NAME}/~_stem/residual_conv_block/conv_block/batch_norm/~/"
        )
        torch_prefix = "stem.1.fn.0."

        state = _get_one_state(state, trix_prefix, torch_prefix)
        return state

    def _get_pointwise_state(
        state: Dict[str, Dict[str, np.ndarray]]
    ) -> Dict[str, Dict[str, np.ndarray]]:

        trix_prefix = (
            f"{ENFORMER_MODEL_NAME}/~_final_pointwise/conv_block/batch_norm/~/"
        )
        torch_prefix = "final_pointwise.1.0."

        state = _get_one_state(state, trix_prefix, torch_prefix)
        return state

    def _get_conv_tower_state(
        state: Dict[str, Dict[str, np.ndarray]], layer_num: int
    ) -> Dict[str, Dict[str, np.ndarray]]:

        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_conv_tower/layer_{layer_num}/"
        trix_prefix += "conv_block/batch_norm/~/"
        torch_prefix = f"conv_tower.{layer_num}.0.0."
        state = _get_one_state(state, trix_prefix, torch_prefix)

        trix_prefix = f"{ENFORMER_MODEL_NAME}/~_conv_tower/layer_{layer_num}/"
        trix_prefix += "residual_conv_block/conv_block/batch_norm/~/"
        torch_prefix = f"conv_tower.{layer_num}.1.fn.0."

        state = _get_one_state(state, trix_prefix, torch_prefix)
        return state

    # update state
    state = _get_stem_state(state)
    for i in range(num_conv_tower_layers):
        state = _get_conv_tower_state(state, layer_num=i)
    state = _get_pointwise_state(state)

    return state  # type: ignore


def get_pretrained_enformer_model(
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
) -> Tuple[hk.Params, hk.State, Callable, NucleotidesKmersTokenizer, EnformerConfig]:
    """
    Create a Haiku Enformer model by downloading the pytorch weights hosted by
    HuggingFace and translating them. Note that the Trix Enformer code has been
    directly adapted from the enformer pytorch code and not from the tensorflow code.


    Args:
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too
        output_dtype: the output type of the model. It determines the float precioson
            of the gradient when training the model.
            NOTE: when training, the gradient is often accumulated in fp32, therefore
            output_dtype need to be in fp32.

    Returns:
        Model parameters.
        Model state.
        Haiku function to call the model.
        Tokenizer.
        Model config.

    Usage:
        >>> params, state, hk_fn, tokenizer, config = get_pretrained_enformer_model()
        >>> enformer_fn = hk.transform_with_state(hk_fn)
        >>> apply_fn = enformer_fn.apply
        >>> rng_key = jax.random.PRNGKey(0)
        >>> seq_length = 196_608
        >>> seq = ''.join(np.random.choice(list('ATCGN'), size=(seq_length,)))
        >>> tokens_ids = [b[1] for b in tokenizer.batch_tokenize([seq])]
        >>> tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)
        >>> outs, state = apply_fn(params, state, rng_key, tokens, is_training=False)
    """

    torch_enformer = Enformer.from_pretrained("EleutherAI/enformer-official-rough")
    torch_params = torch_enformer.state_dict()

    # Default params correspond to params from HF and official ones from the paper
    config = EnformerConfig()

    # Get parameters
    hk_params = translate_weights(
        torch_params=torch_params,
        num_conv_tower_layers=config.num_downsamples - 1,
        num_transformer_tower_layers=config.num_transformer_layers,
    )

    # Get state
    hk_state = translate_state(
        torch_params=torch_params, num_conv_tower_layers=config.num_downsamples - 1
    )

    # Remove the torch parameters from the RAM.
    del torch_params
    gc.collect()
    torch.cuda.empty_cache()

    enformer_fn = build_enformer_fn(
        config=config,
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
        name=ENFORMER_MODEL_NAME,
    )

    # NOTE: the tokenizer should have tokens for ONLY A,T,C,G,N nucleotides.
    # The token IDs of A,T,C,G,N should be the ones by default
    # of NucleotidesKmersTokenizer: A:10 / C:12 / G:13 / T:11 / N:14.
    # If the token IDs are different, the one-hot encoded vectors from Enformer
    # will not match the nucleotides and it will fail.

    tokenizer = NucleotidesKmersTokenizer(
        k_mers=1,
        prepend_bos_token=False,
        prepend_cls_token=False,
        append_eos_token=False,
        tokens_to_ids=None,
    )

    # check if token IDs of nucleotides are correct
    nucleotide_to_token_ids = {"A": 10, "C": 12, "G": 13, "T": 11, "N": 14}
    assert all(
        [tokenizer.token_to_id(n) == i for n, i in nucleotide_to_token_ids.items()]
    ), (
        "Token IDs of nucleotides don't match, they should be: "
        "A:10 / C:12 / G:13 / T:11 / N:14"
    )

    return hk_params, hk_state, enformer_fn, tokenizer, config


def get_pretrained_segment_enformer_model(
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
) -> Tuple[hk.Params, hk.State, Callable, NucleotidesKmersTokenizer, EnformerConfig]:
    _, _, _, tokenizer, config = get_pretrained_enformer_model(
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
    )

    def head_fn() -> hk.Module:
        return UNetHead(
            features=FEATURES,
            embed_dimension=config.embed_dim,
            nucl_per_token=config.dim_divisible_by,
            remove_cls_token=False,
        )

    forward_fn = build_enformer_with_head_fn(
        config=config,
        head_fn=head_fn,
        embedding_name="embedding_transformer_tower",
        name="Enformer",
        compute_dtype=jnp.float32,
    )

    parameters, state = download_ckpt("segment_enformer")

    return parameters, state, forward_fn, tokenizer, config
