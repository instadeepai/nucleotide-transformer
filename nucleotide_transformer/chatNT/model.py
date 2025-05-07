from dataclasses import replace
from typing import Callable, Optional, Tuple

import haiku as hk
import jax
import jax.numpy as jnp
import jmp

from nucleotide_transformer.chatNT.gpt_decoder import GptConfig, GptDecoder
from nucleotide_transformer.chatNT.multi_modal_perceiver_projection import (
    MultiModalPerceiverResamplerProjection,
    PerceiverResamplerConfig,
)
from nucleotide_transformer.chatNT.types import MultiOmicsTokens, TransformerOutput
from nucleotide_transformer.model import (
    NucleotideTransformer,
    NucleotideTransformerConfig,
)


class ChatNTDecoder(hk.Module):
    """
    Creates the decoder for the ChatNT model.
    Mainly the GPT decoder but wrapped inside a module to properly
    control mixed precision.
    """

    def __init__(
        self,
        gpt_config: GptConfig,
        seq_token_id: int,
        gpt_name: Optional[str] = None,
        name: Optional[str] = None,
    ):
        """
        Initializes the bio to english projection layer.

        Args:
            gpt_config: config of the GPT model
            gpt_name: name for the GPT model
            seq_token_id: Index of the SEQ token
            name: haiku module name
        """
        super().__init__(name=name)
        self._gpt_config = gpt_config
        self._gpt_name = gpt_name
        self._seq_token_id = seq_token_id

    def __call__(
        self, english_token_ids: jnp.ndarray, projected_bio_embeddings: jnp.ndarray
    ) -> jnp.ndarray:

        # get GPT model
        outs = {}
        gpt_model = GptDecoder(self._gpt_config, name=self._gpt_name)

        # compute english tokens embeddings
        tokens_embeddings = gpt_model.token_embed(english_token_ids)

        if projected_bio_embeddings is not None:

            # bio token ids assumed to be stacked (after padding) with shape
            # (batch_size, num_bio_sequences, bio_sequence_length, bio_embed_dim)
            num_bio_sequences = projected_bio_embeddings.shape[1]

            # replace english tok embedding by bio embeddings at seq token positions
            def insert_embeddings(
                tokens: jnp.ndarray,
                input_embeddings: jnp.ndarray,
                resampled_embeddings: jnp.ndarray,
                bio_seq_num: int,
            ) -> Tuple[jnp.ndarray, jnp.ndarray]:
                """
                Insert resampled embeddings in input_embeddings, starting at
                seq_token_id position.

                NOTE: this function only accepts one DNA sequence in the input text, we
                could modify this by adding a max_num_dna_sequences and using
                size=max_num_dna_sequences (if there are less, then idx will be 0 for
                those values).
                """

                def _insert(tokens, input_embeddings, resampled_embeddings):  # type: ignore # noqa
                    idx = jnp.argwhere(tokens == self._seq_token_id, size=1)[0][0]
                    x = jnp.insert(
                        input_embeddings,
                        idx + bio_seq_num * resampled_embeddings.shape[-2],
                        resampled_embeddings,
                        axis=0,
                    )[: tokens.shape[0] + 1]
                    # remove embedding of seq_token_id from the tokens
                    x = jnp.roll(jnp.roll(x, idx)[:-1], -idx)
                    tokens = tokens.at[idx].set(-1)
                    return x, tokens

                return jax.vmap(_insert)(  # type: ignore
                    tokens, input_embeddings, resampled_embeddings
                )

            resampled_length = projected_bio_embeddings.shape[-2]

            def cleanup_logits(
                tokens: jnp.ndarray, logits: jnp.ndarray
            ) -> Tuple[jnp.ndarray, jnp.ndarray]:
                """
                Removes the logits corresponding to the unused embeddings.

                Args:
                    tokens: Input english tokens.
                    logits: Input logits.

                Returns:
                    Cleaned logits, last values will be equal to 0.
                """

                def _cleanup(
                    tokens: jnp.ndarray, logit: jnp.ndarray
                ) -> Tuple[jnp.ndarray, jnp.ndarray]:
                    idx = jnp.argwhere(tokens == self._seq_token_id, size=1)[0][0]
                    mask_idx = jnp.arange(logit.shape[0] - resampled_length) > idx
                    mask_idx = jnp.expand_dims(mask_idx, axis=1)
                    # remove values corresponding to bio tokens
                    logit = (
                        logit[:-resampled_length] * (1 - mask_idx)
                        + logit[resampled_length:] * mask_idx
                    )
                    # append zeros at the end
                    logit = jnp.concatenate(
                        (
                            logit,
                            jnp.zeros(
                                shape=(resampled_length, logit.shape[1]),
                                dtype=logit.dtype,
                            ),
                        )
                    )
                    # replace seq_token_id by -1 to apply this function sequentially to
                    # clean one by one the bio embeddings logits
                    tokens = tokens.at[idx].set(-1)
                    return logit, tokens

                return jax.vmap(_cleanup)(tokens, logits)  # type: ignore

            processed_tokens_ids = english_token_ids
            for bio_seq_num in range(num_bio_sequences):
                tokens_embeddings, processed_tokens_ids = insert_embeddings(
                    processed_tokens_ids,
                    tokens_embeddings,
                    resampled_embeddings=projected_bio_embeddings[:, bio_seq_num],
                    bio_seq_num=bio_seq_num,
                )

        # regular pass through GPT
        embeddings = gpt_model.apply_transformer_layers(tokens_embeddings)
        embeddings = gpt_model.final_norm(embeddings)

        # compute logits
        logits = gpt_model.lm_head(embeddings)["logits"]

        if projected_bio_embeddings is not None:
            # clean logits sequentially
            processed_tokens_ids = english_token_ids
            for _ in range(num_bio_sequences):
                logits, processed_tokens_ids = cleanup_logits(
                    tokens=processed_tokens_ids, logits=logits
                )

        outs["logits"] = logits

        return outs  # type: ignore


class ChatNTEncoder(hk.Module):
    """
    Creates the bio encoder for the ChatNT model.
    Mainly the Nucleotide Transformer but wrapped inside a module to properly
    control mixed precision.
    """

    def __init__(
        self,
        nt_config: GptConfig,
        nt_name: Optional[str] = None,
        name: Optional[str] = None,
    ):
        """
        Args:
            nt_config: config of the NT model
            nt_name: name for the NT model
            name: haiku module name
        """
        super().__init__(name=name)
        self._config = nt_config
        self._nt_config = nt_config
        self._nt_name = nt_name

    def __call__(
        self,
        bio_token_ids: jnp.ndarray,
    ) -> jnp.ndarray:
        # get NT model
        nt_model = NucleotideTransformer(self._nt_config, name=self._nt_name)

        # compute bio seq embeddings
        nt_outs = nt_model(tokens=bio_token_ids)
        num_nt_layers = self._config.num_layers
        bio_embeddings = nt_outs[f"embeddings_{num_nt_layers}"]

        return bio_embeddings


def build_chat_nt_fn(
    nt_config: NucleotideTransformerConfig,
    gpt_config: GptConfig,
    seq_token_id: int,
    bio_pad_token_id: int,
    english_pad_token_id: int,
    perceiver_resampler_config: PerceiverResamplerConfig,
    nt_compute_dtype: jnp.dtype = jnp.float32,
    nt_param_dtype: jnp.dtype = jnp.float32,
    nt_output_dtype: jnp.dtype = jnp.float32,
    gpt_compute_dtype: jnp.dtype = jnp.float16,
    gpt_param_dtype: jnp.dtype = jnp.float16,
    gpt_output_dtype: jnp.dtype = jnp.float16,
    nt_name: Optional[str] = None,
    gpt_name: Optional[str] = None,
) -> Callable:
    """
    Creates ChatNT forward pass.
    """

    assert {nt_compute_dtype, nt_param_dtype, nt_output_dtype}.issubset(
        {
            jnp.bfloat16,
            jnp.float32,
            jnp.float16,
        }
    ), f"provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    assert {gpt_compute_dtype, gpt_param_dtype, gpt_output_dtype}.issubset(
        {
            jnp.bfloat16,
            jnp.float32,
            jnp.float16,
        }
    ), f"provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    if perceiver_resampler_config.embed_dim != gpt_config.embed_dim:
        raise ValueError(
            "when using perceiver projection function, its embed_dim must match "
            "the one of gpt_config"
        )

    # set policy for GPTDecoder
    gpt_policy = jmp.Policy(
        compute_dtype=gpt_compute_dtype,
        param_dtype=gpt_param_dtype,
        output_dtype=gpt_output_dtype,
    )
    hk.mixed_precision.set_policy(ChatNTDecoder, gpt_policy)

    # set policy for NT
    nt_policy = jmp.Policy(
        compute_dtype=nt_compute_dtype,
        param_dtype=nt_param_dtype,
        output_dtype=nt_output_dtype,
    )
    hk.mixed_precision.set_policy(ChatNTEncoder, nt_policy)

    # set policy for projection
    # output dtype should match GPT compute dtype
    projection_policy = jmp.Policy(
        compute_dtype=nt_compute_dtype,
        param_dtype=nt_param_dtype,
        output_dtype=gpt_compute_dtype,
    )
    hk.mixed_precision.set_policy(
        MultiModalPerceiverResamplerProjection, projection_policy
    )

    # remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32,  # TO BE CHECKED
        param_dtype=gpt_param_dtype,
        output_dtype=gpt_compute_dtype,
    )
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.RMSNorm, norm_policy)

    # make sure to get final embeddings with NT
    num_nt_layers = nt_config.num_layers
    emb_layers_to_save = nt_config.embeddings_layers_to_save + (num_nt_layers,)
    nt_config = replace(nt_config, embeddings_layers_to_save=emb_layers_to_save)

    def chat_nt_fn(
        multi_omics_tokens_ids: MultiOmicsTokens,
        projection_english_tokens_ids: jnp.ndarray,
        projected_bio_embeddings: Optional[jnp.ndarray] = None,
    ) -> TransformerOutput:

        # unpack tokens
        english_token_ids, bio_token_ids = multi_omics_tokens_ids

        # create decoder
        decoder = ChatNTDecoder(
            gpt_config=gpt_config, gpt_name=gpt_name, seq_token_id=seq_token_id
        )

        if bio_token_ids is None:
            projected_bio_embeddings = None
        else:
            # bio token ids assumed to be stacked (after padding) with shape
            # (batch_size, num_bio_sequences, bio_sequence_length)
            num_bio_sequences = bio_token_ids.shape[1]

            # therefore, we might expect this function to be compiled several times
            # for the different possible number of bio sequences

            # create models
            bio_encoder = ChatNTEncoder(nt_config=nt_config, nt_name=nt_name)

            # project bio-embeddings to same shape than English token embeddings
            projection_model = MultiModalPerceiverResamplerProjection(
                config=perceiver_resampler_config,
                embed_dim=gpt_config.embed_dim,
                bio_pad_token_id=bio_pad_token_id,
                english_pad_token_id=english_pad_token_id,
                english_vocab_size=gpt_config.vocab_size,
            )

            if projected_bio_embeddings is None:
                # compute bio seq embeddings
                # for loop over different bio seqs on purpose - will slow down
                # compilation but worth it for inference and training afterward
                bio_embeddings_list = [
                    bio_encoder(bio_token_ids=bio_token_ids[:, bio_seq_num])
                    for bio_seq_num in range(num_bio_sequences)
                ]

                projected_bio_embeddings = [
                    projection_model(
                        bio_token_ids=bio_token_ids[:, bio_seq_num],
                        bio_embeddings=bio_embeddings,
                        english_tokens_ids=projection_english_tokens_ids,
                    )
                    for bio_seq_num, bio_embeddings in enumerate(bio_embeddings_list)
                ]
                projected_bio_embeddings = jnp.stack(projected_bio_embeddings, axis=1)

        # decode
        outs = decoder(
            english_token_ids=english_token_ids,
            projected_bio_embeddings=projected_bio_embeddings,
        )

        # Insert projected_bio_embeddings in the outputs
        outs["projected_bio_embeddings"] = projected_bio_embeddings
        return outs  # type: ignore

    return chat_nt_fn
