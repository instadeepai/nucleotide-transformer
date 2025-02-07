"""
This file gathers layers used to finetune pre-trained Borzoi models.

NOTE: need to implement also IA3 scaling weights!!
"""

from typing import Callable, Dict, Optional

import haiku as hk
import jax.numpy as jnp
import jmp
from typing_extensions import TypeAlias

from borzoi.model import Borzoi, BorzoiConfig

SequenceMask: TypeAlias = jnp.ndarray


def build_borzoi_fn_with_head_fn(
    config: BorzoiConfig,
    head_fn: Callable[
        [], Callable[[jnp.ndarray, SequenceMask], Dict[str, jnp.ndarray]]
    ],
    embedding_name: str = "embedding",
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    name: Optional[str] = None,
) -> Callable:
    """
    Create the model's forward pass for that Borzoi and adds the input head.

    Args:
        config: Configuration data class containing the hyperparameters for the Borzoi
            forward function.
        head_fn: Wrapper initializing a Classification/Regression head. The head cannot
            be passed directly as haiku modules cannot be initialized outside
            hk.transform.
        embedding_name: Which embeddings to use from the Borzoi pre-trained model as
            input to the head. It should be the name of the key inside model
            predictions ( outs ). Default is "embedding".
        compute_dtype: the type of the activations. fp16 runs faster and is lighter in
            memory. bf16 handles better large int, and is hence more stable ( it avoids
            float overflows ).
        param_dtype: if compute_dtype is fp16, the model weights will be cast to fp16
            during the forward pass anyway. So in inference mode ( not training mode ),
            it is better to use params in fp16 if compute_dtype is fp16 too
        output_dtype: the output type of the model. it determines the float precision
            of the gradient when training the model.
            NOTE: when training, the gradient is often accumulated in fp32, therefore
            output_dtype need to be in fp32.
        name: the name of the model. example: borzoi.

        Example of the function being used with a classification head:
        The classification head is wrapped inside head_fn because
        haiku modules cannot be instantiated outside hk.transform.
        def head_fn():
            return SimpleClassificationHead(num_classes=num_classes)
        finetune_forward_fn = build_esm_ia3_rescaling_with_head_fn(
            model_config=config, head_fn=head_fn, model_name=model_name,
        )
        finetune_forward_fn = hk.transform(finetune_forward_fn)

        # NOTE: the input tokens for borzoi of shape (batch_size, seq_len)
            are expected to be vectors with token IDs corresponding
            to ONLY A,T,C,G,N nucleotides.
            The token IDs of A,T,C,G,N should be the ones by default
            of NucleotidesKmersTokenizer: A:10 / C:12 / G:13 / T:11 / N:14.
            If the token IDs are different, the one-hot encoded vectors
            will not match the nucleotides and it will fail.
            Sequences cannot have the pad token - to pad the sequences
            you can add the nucleotide N.

    Returns:
        Borzoi model forward function with indicated head.
    """

    assert {compute_dtype, param_dtype, output_dtype}.issubset(
        {
            jnp.bfloat16,
            jnp.float32,
            jnp.float16,
        }
    ), f"provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32, param_dtype=param_dtype, output_dtype=compute_dtype
    )
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)

    def borzoi_fn(
        tokens: jnp.ndarray,
        is_training: bool = False,
        sequence_mask: Optional[SequenceMask] = None,
    ) -> Dict[str, jnp.ndarray]:
        """Forward pass."""
        # Run the encoder over the inputs.
        encoder = Borzoi(config, name)
        outs = encoder(tokens, is_training)
        # NOTE: For now we don't remove the Borzoi human/mouse prediction heads
        # since they will not be finetuned. But we could remove them to save space

        # Get embeddings to use as input for head
        embeddings = outs[embedding_name]

        # Define head
        head = head_fn()

        if sequence_mask is None:
            # I should not get the last (embedding) dimension,
            # because the mask should have the dimensions of the input sequence
            sequence_mask = jnp.ones_like(embeddings[:, :, 0])

        head_outs = head(  # type: ignore[call-arg]
            x=embeddings, sequence_mask=sequence_mask
        )
        outs.update(head_outs)
        return outs

    return borzoi_fn
