"""Implementation of utilities to load a pretrained Enformer model in Trix."""

from typing import Callable, Tuple

import haiku as hk
import jax.numpy as jnp

from nucleotide_transformer.enformer.features import FEATURES
from nucleotide_transformer.enformer.heads import UNetHead
from nucleotide_transformer.enformer.model import (
    EnformerConfig,
    build_enformer_with_head_fn,
)
from nucleotide_transformer.enformer.params import download_ckpt
from nucleotide_transformer.enformer.tokenizer import NucleotidesKmersTokenizer


def get_pretrained_segment_enformer_model() -> (
    Tuple[hk.Params, hk.State, Callable, NucleotidesKmersTokenizer, EnformerConfig]
):
    """
    Loads the pretrained SegmentEnformer model.

    Returns:
        hk.Params: Model parameters
        hk.State: Model state
        Callable: Haiku forward function
        NucleotidesKmersTokenizer: Tokenizer
        EnformerConfig: Configuration of the Enformer model

    Example:
        >>> import jax
        >>> import haiku as hk
        >>> parameters, state, forward_fn, tokenizer, config = get_pretrained_segment_enformer_model()
        >>> apply_fn = hk.transform_with_state(forward_fn).apply
        >>> random_key = jax.random.PRNGKey(seed=0)
        >>> sequences = ["A" * 196608]
        >>> tokens = jax.numpy.asarray([b[1] for b in tokenizer.batch_tokenize(sequences)])
        >>> outs, _ = apply_fn(parameters, state, random_key, tokens)
    """
    config = EnformerConfig()
    tokenizer = NucleotidesKmersTokenizer(
        k_mers=1,
        prepend_bos_token=False,
        prepend_cls_token=False,
        append_eos_token=False,
        tokens_to_ids=None,
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
