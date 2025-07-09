import json
import os
from typing import Callable

import haiku as hk
import jax.numpy as jnp
import joblib
from huggingface_hub import hf_hub_download

from nucleotide_transformer.bulk_rna_bert.config import BulkRNABertConfig
from nucleotide_transformer.bulk_rna_bert.model import build_bulk_rna_bert_forward_fn
from nucleotide_transformer.bulk_rna_bert.tokenizer import BinnedOmicTokenizer

ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"


def _get_dir(model_name: str) -> str:
    """
    Get directory to save files on user machine.
    """
    return os.path.expanduser(
        os.path.join(os.getenv(ENV_XDG_CACHE_HOME, DEFAULT_CACHE_DIR), model_name)
    )


def download_bulkrnabert_ckpt() -> tuple[hk.Params, BulkRNABertConfig]:
    """
    Download BulkRNABert checkpoint from Hugging Face.


    Returns:
        Model parameters.
        Model configuration
    """

    save_dir = os.path.join(_get_dir("bulkrnabert"), "bulkrnabert")

    repo_id = "InstaDeepAI/BulkRNABert"

    # Download parameters
    print("Downloading model's weights...")
    params = joblib.load(
        hf_hub_download(
            repo_id=repo_id,
            filename="jax_params/params.joblib",
            cache_dir=save_dir,
        )
    )

    config_path = hf_hub_download(
        repo_id=repo_id,
        filename="jax_params/config.json",
    )
    with open(config_path) as f:
        config_dict = json.load(f)
    config = BulkRNABertConfig(**config_dict)

    return params, config


def get_pretrained_bulkrnabert_model(
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    embeddings_layers_to_save: tuple[int, ...] = (),
) -> tuple[hk.Params, Callable, BinnedOmicTokenizer, BulkRNABertConfig]:
    """
    Create a Haiku BulkRNABert model.
    model by downloading pre-trained weights and hyperparameters.

    Args:
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too. During
            training, it is preferable to keep parameters in float32 for better
            numerical stability.
        output_dtype: the output type of the model. it determines the float precision
            of the gradient when training the model.
        embeddings_layers_to_save: Intermediate embeddings to return in the output.

    Returns:
        Model parameters.
        Haiku function to call the model.
        Tokenizer.
        Model config (hyperparameters).

    """
    parameters, config = download_bulkrnabert_ckpt()
    tokenizer = BinnedOmicTokenizer(
        n_expressions_bins=config.n_expressions_bins,
        use_max_normalization=config.use_max_normalization,
        normalization_factor=config.normalization_factor,  # type: ignore
        prepend_cls_token=False,
    )

    config.embeddings_layers_to_save = embeddings_layers_to_save

    forward_fn = build_bulk_rna_bert_forward_fn(
        model_config=config,
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
        model_name="bulk_bert",
    )

    return parameters, forward_fn, tokenizer, config
