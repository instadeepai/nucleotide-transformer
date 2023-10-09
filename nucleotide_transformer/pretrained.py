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

import json
import os
from typing import Any, Callable, Dict, Optional, Tuple

import boto3
import haiku as hk
import jax.numpy as jnp
import joblib
import tqdm
from botocore import UNSIGNED
from botocore.config import Config

from nucleotide_transformer.model import (
    NucleotideTransformerConfig,
    build_nucleotide_transformer_fn,
)
from nucleotide_transformer.tokenizers import (
    FixedSizeNucleotidesKmersTokenizer,
    compute_tokens_to_ids_v2,
)

ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"


def _get_dir() -> str:
    """
    Get directory to save files on user machine.
    """
    return os.path.expanduser(
        os.path.join(
            os.getenv(ENV_XDG_CACHE_HOME, DEFAULT_CACHE_DIR), "nucleotide_transformer"
        )
    )


def download_from_s3_bucket(
    s3_client: boto3.session.Session, bucket: str, key: str, filename: str
) -> None:
    """
    Download data from the s3 bucket and display downloading progression bar.

    Args:
        s3_client: Boto3 s3 client
        bucket: Bucket name.
        key: Path towards file in the bucket.
        filename: Path to save file locally.
    """
    kwargs = {
        "Bucket": bucket,
        "Key": key,
    }
    object_size = s3_client.head_object(**kwargs)["ContentLength"]
    with tqdm.tqdm(total=object_size, unit="B", unit_scale=True, desc=filename) as pbar:
        with open(filename, "wb") as f:
            s3_client.download_fileobj(
                Bucket=bucket,
                Key=key,
                ExtraArgs=None,
                Fileobj=f,
                Callback=lambda bytes_transferred: pbar.update(bytes_transferred),
            )


def download_ckpt_and_hyperparams(model_name: str) -> Tuple[hk.Params, Dict[str, Any]]:
    """
    Download checkpoint and hyperparams on kao datacenter.

    Args:
        model_name: Name of the model.

    Returns:
        Model parameters.
        Model hyperparameters' dict.

    """
    # Get directories
    save_dir = os.path.join(_get_dir(), model_name)

    params_save_dir = os.path.join(save_dir, "ckpt.joblib")
    hyperparams_save_dir = os.path.join(save_dir, "hyperparams.json")

    if os.path.exists(hyperparams_save_dir) and os.path.exists(params_save_dir):
        # Load locally
        with open(hyperparams_save_dir, "rb") as f:
            hyperparams = json.load(f)

        with open(params_save_dir, "rb") as f:
            params = joblib.load(f)

        return params, hyperparams

    else:
        os.makedirs(save_dir, exist_ok=True)

        s3_endpoint = "https://s3.kao-prod.instadeep.io"

        session = boto3.Session()
        s3_client = session.client(
            service_name="s3",
            endpoint_url=s3_endpoint,
            config=Config(signature_version=UNSIGNED),
        )

        # Download params and hyperparams
        bucket = "nucleotide-transformer"
        print(f"checkpoints/{model_name}/hyperparams.json")
        download_from_s3_bucket(
            s3_client=s3_client,
            bucket=bucket,
            key=f"checkpoints/{model_name}/hyperparams.json",
            filename=hyperparams_save_dir,
        )

        download_from_s3_bucket(
            s3_client=s3_client,
            bucket=bucket,
            key=f"checkpoints/{model_name}/ckpt.joblib",
            filename=params_save_dir,
        )

        # Load locally
        with open(hyperparams_save_dir, "rb") as f:
            hyperparams = json.load(f)

        with open(params_save_dir, "rb") as f:
            params = joblib.load(f)

        return params, hyperparams


def get_pretrained_model(
    model_name: str,
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    embeddings_layers_to_save: Tuple[int, ...] = (),
    attention_maps_to_save: Optional[Tuple[Tuple[int, int], ...]] = None,
    max_positions: int = 1024,
) -> Tuple[
    hk.Params, Callable, FixedSizeNucleotidesKmersTokenizer, NucleotideTransformerConfig
]:
    """
    Create a Haiku Nucleotide Transformer
    model by downloading pre-trained weights and hyperparameters.
    Nucleotide Transformer Models have ESM-like architectures.

    Args:
        model_name: Name of the model.
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too. During
            training, it is preferable to keep parameters in float32 for better
            numerical stability.
        output_dtype: the output type of the model. it determines the float precioson
            of the gradient when training the model.
        embeddings_layers_to_save: Intermediate embeddings to return in the output.
        attention_maps_to_save: Intermediate attention maps to return in the output.
        max_positions: Maximum length of a token (for padding).

    Returns:
        Model parameters.
        Haiku function to call the model.
        Tokenizer.
        Model config (hyperparameters).

    Example:
        parameters, forward_fn, tokenizer, config = get_pretrained_model(
            model_name="500M_1000G",
            mixed_precision=False,
            # Get embedding at layers 5 and 20
            embeddings_layers_to_save=(5, 20,),
            # Get attention map number 4 at layer 1 and attention map number 14
            # at layer 12
            attention_maps_to_save=((1,4), (12, 14)),
            max_positions=128,
        )
    """
    if attention_maps_to_save is None:
        attention_maps_to_save = ()

    supported_models = [
        "500M_human_ref",
        "500M_1000G",
        "2B5_1000G",
        "2B5_multi_species",
        "50M_multi_species_v2",
        "100M_multi_species_v2",
        "250M_multi_species_v2",
        "500M_multi_species_v2",
    ]

    if not (model_name in supported_models):
        raise NotImplementedError(
            f"Unknown {model_name} model. " f"Supported models are {supported_models}"
        )

    # Download weights and hyperparams
    parameters, hyperparams = download_ckpt_and_hyperparams(model_name)

    if "v2" in model_name:
        tokens_to_ids, _ = compute_tokens_to_ids_v2(k_mers=hyperparams["k_for_kmers"])
        tokenizer = FixedSizeNucleotidesKmersTokenizer(
            k_mers=hyperparams["k_for_kmers"],
            fixed_length=max_positions,
            prepend_cls_token=True,
            tokens_to_ids=tokens_to_ids,
        )

        # Adapat hyperparameters from config
        alphabet_size = len(tokenizer.vocabulary)
    else:
        tokenizer = FixedSizeNucleotidesKmersTokenizer(
            k_mers=hyperparams["k_for_kmers"],
            fixed_length=max_positions,
            prepend_cls_token=True,
        )
        alphabet_size = len(tokenizer.vocabulary) - 2

    if "add_bias_ffn" in hyperparams.keys():
        add_bias_ffn = hyperparams["add_bias_ffn"]
    else:
        add_bias_ffn = True

    if "ffn_activation_name" in hyperparams.keys():
        ffn_activation_name = hyperparams["ffn_activation_name"]
    else:
        ffn_activation_name = "gelu-no-approx"

    if "use_glu_in_ffn" in hyperparams.keys():
        use_glu_in_ffn = hyperparams["use_glu_in_ffn"]
    else:
        use_glu_in_ffn = False

    if "add_bias_kv" in hyperparams.keys():
        add_bias_kv = hyperparams["add_bias_kv"]
    else:
        add_bias_kv = False

    if "use_rotary_embedding" in hyperparams.keys():
        use_rotary_embedding = hyperparams["use_rotary_embedding"]
    else:
        use_rotary_embedding = False

    if "positional_embedding" in hyperparams.keys():
        positional_embedding = hyperparams["positional_embedding"]
    else:
        positional_embedding = "learned"

    # Get config
    config = NucleotideTransformerConfig(
        alphabet_size=alphabet_size,
        pad_token_id=tokenizer.pad_token_id,
        mask_token_id=tokenizer.mask_token_id,
        max_positions=hyperparams["max_positions"],
        embed_scale=hyperparams["embed_scale"],
        # architecture
        emb_layer_norm_before=hyperparams["emb_layer_norm_before"],
        key_size=hyperparams["key_dim"] if "key_dim" in hyperparams.keys() else None,
        attention_heads=hyperparams["attention_heads"],
        embed_dim=hyperparams["embed_dim"],
        ffn_embed_dim=hyperparams["ffn_embed_dim"],
        num_layers=hyperparams["num_layers"],
        positional_embedding=positional_embedding,
        add_bias_kv=add_bias_kv,
        add_bias_ffn=add_bias_ffn,
        use_glu_in_ffn=use_glu_in_ffn,
        ffn_activation_name=ffn_activation_name,
        use_rotary_embedding=use_rotary_embedding,
        # bert
        token_dropout=hyperparams["token_dropout"],
        masking_ratio=hyperparams["masking_ratio"],
        masking_prob=hyperparams["masking_prob"],
        # embeddings to save
        embeddings_layers_to_save=embeddings_layers_to_save,
        attention_maps_to_save=attention_maps_to_save,
    )

    # NOTE: module names are changed here, to validate !
    full_model_name = "nucleotide_transformer" + model_name
    parameters = rename_modules(parameters, full_model_name)

    forward_fn = build_nucleotide_transformer_fn(
        model_config=config,
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
        model_name=full_model_name,
    )

    return parameters, forward_fn, tokenizer, config


def rename_modules(parameters: hk.Params, model_name: str) -> hk.Params:
    """
    Adjusts the names of the modules from checkpoints to NucleotideTransformer.

    Args:
        parameters: Parameters loaded from .joblib archive.
        model_name: Name of the loaded model.

    Returns:
        Parameters with updated names.
    """
    for layer_name in list(parameters.keys()):
        new_name = layer_name.replace("esm_transformer", model_name)

        if "attention_layer" in new_name:
            if new_name.split("/")[3] == "mha":
                new_name = "/".join(
                    new_name.split("/")[:3]
                    + ["self_attention"]
                    + new_name.split("/")[4:]
                )
        if "mha_layer_norm" in new_name:
            new_name = new_name.replace("mha_layer_norm", "self_attention_layer_norm")
        if "esm_roberta_lm_head" in new_name:
            new_name = new_name.replace("esm_roberta_lm_head", "roberta_lm_head")

        parameters[new_name] = parameters.pop(layer_name)

    return parameters
