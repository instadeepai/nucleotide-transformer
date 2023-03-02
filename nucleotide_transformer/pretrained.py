import json
import os
from typing import Any, Callable, Dict, Optional, Tuple

import boto3
import haiku as hk
import joblib
import tqdm

from nucleotide_transformer.model import (
    NucleotideTransformerConfig,
    build_nucleotide_transformer_fn,
)
from nucleotide_transformer.tokenizers import FixedSizeNucleotidesKmersTokenizer

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

        # This part assumes the user has already exported required env variables
        s3_endpoint = "https://s3.kao-prod.instadeep.io"

        session = boto3.Session()
        s3_client = session.client(
            service_name="s3",
            aws_access_key_id="BG39IQWKYHL08M02WQ55",
            aws_secret_access_key="0iA69KYqGMwQf7t8F8HChXODLbZhgWmmcDL0bATw",
            endpoint_url=s3_endpoint,
        )

        # Download params and hyperparams
        bucket = "nucleotide-transformer"

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
    mixed_precision: bool = False,
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
        mixed_precision: Whether to use mixed precision.
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
    ]

    if not (model_name in supported_models):
        raise NotImplementedError(
            f"Unknown {model_name} model. " f"Supported models are {supported_models}"
        )

    # Download weights and hyperparams
    parameters, hyperparams = download_ckpt_and_hyperparams(model_name)

    tokenizer = FixedSizeNucleotidesKmersTokenizer(
        k_mers=hyperparams["k_for_kmers"],
        fixed_length=max_positions,
        prepend_cls_token=True,
    )

    # Get config
    config = NucleotideTransformerConfig(
        alphabet_size=len(tokenizer.vocabulary) - 2,
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
        model_config=config, mixed_precision=mixed_precision, model_name=full_model_name
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
