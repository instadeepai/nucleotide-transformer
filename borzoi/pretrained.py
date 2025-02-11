"""Implementation of utilities to load a pretrained Borzoi model in Trix."""

from typing import Callable, Tuple

import haiku as hk
import jax.numpy as jnp

from borzoi.model import BorzoiConfig, build_borzoi_fn, build_borzoi_fn_with_head_fn
from enformer.features import FEATURES
from enformer.heads import UNetHead
from enformer.params import download_ckpt
from enformer.tokenizer import NucleotidesKmersTokenizer

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


def get_pretrained_segment_borzoi_model():
    config, tokenizer = get_borzoi_config_and_tokenizer()

    def head_fn() -> hk.Module:
        return UNetHead(
            features=FEATURES,
            embed_dimension=config.embed_dim,
            nucl_per_token=config.dim_divisible_by,
            remove_cls_token=False,
        )

    forward_fn = build_borzoi_fn_with_head_fn(
        config=config,
        head_fn=head_fn,
        embedding_name="embedding",
        name="Borzoi",
        compute_dtype=jnp.float32,
    )

    parameters, state = download_ckpt("segment_borzoi")

    return parameters, state, forward_fn, tokenizer, config
