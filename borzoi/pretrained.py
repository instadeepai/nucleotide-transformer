"""Implementation of utilities to load a pretrained Borzoi model in Trix."""

from typing import Callable, Tuple

import haiku as hk
import jax.numpy as jnp
from trix.tokenizers.language_models.bio import NucleotidesKmersTokenizer

from borzoi.model import BorzoiConfig, build_borzoi_fn

BORZOI_MODEL_NAME = "Borzoi"


def get_borzoi_config_and_tokenizer() -> Tuple[
    BorzoiConfig,
    NucleotidesKmersTokenizer,
]:
    """
    Downloads the weights from the Borzoi models in keras.

    Args:

    Returns:
        Model parameters.
    """

    config = BorzoiConfig()

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

    return config, tokenizer


def get_pretrained_borzoi_model(
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
) -> Tuple[hk.Params, hk.State, Callable, NucleotidesKmersTokenizer, BorzoiConfig]:
    """
    Create a Haiku Borzoi model

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
        Haiku function to call the model.
        Tokenizer.
        Model config.

    """

    config, tokenizer = get_borzoi_config_and_tokenizer()
    borzoi_fn = build_borzoi_fn(
        config=config,
        compute_dtype=compute_dtype,
        param_dtype=param_dtype,
        output_dtype=output_dtype,
        name=BORZOI_MODEL_NAME,
    )

    return borzoi_fn, tokenizer, config
