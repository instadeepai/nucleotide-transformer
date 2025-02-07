"""
This file gathers layers used to finetune pre-trained Enformer models.
It adds IA3 rescailing weights and input head.
"""

from typing import Callable, Dict, Optional

import haiku as hk
import jax
import jax.numpy as jnp
import jmp
from einops import rearrange
from typing_extensions import TypeAlias

from enformer.layers import Attention, relative_shift
from enformer.model import Enformer, EnformerConfig

SequenceMask: TypeAlias = jnp.ndarray


class AttentionIA3Rescaling(Attention):
    """
    This class adds rescaling weights following the IA³ methodology to the
    Enformer Attention Layer. This aims to be used for fine-tuning
    pre-trained Enformer models.
    """

    def __init__(
        self,
        dim: int,
        *,
        num_rel_pos_features: int,
        heads: int = 8,
        dim_key: int = 64,
        dim_value: int = 64,
        name: Optional[str] = None,
    ):
        super().__init__(
            dim=dim,
            num_rel_pos_features=num_rel_pos_features,
            heads=heads,
            dim_key=dim_key,
            dim_value=dim_value,
            name=name,
        )

    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        n = x.shape[-2]

        q = self.to_q(x)
        k = self.to_k(x)
        v = self.to_v(x)

        def _rearrange(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "b n (h d) -> b h n d", h=self.heads)

        q = _rearrange(q)
        k = _rearrange(k)
        v = _rearrange(v)

        # add IA³ rescaling after rearranging
        key_ia3_rescaling = hk.get_parameter(
            "key_ia3_rescaling",
            shape=[k.shape[-3], 1, k.shape[-1]],
            dtype=k.dtype,
            init=hk.initializers.Constant(1.0),
        )
        k = k * key_ia3_rescaling

        # add IA³ rescaling after rearranging
        value_ia3_rescaling = hk.get_parameter(
            "value_ia3_rescaling",
            shape=[v.shape[-3], 1, v.shape[-1]],
            dtype=v.dtype,
            init=hk.initializers.Constant(1.0),
        )
        v = v * value_ia3_rescaling

        # query self-scaling
        q = q * self.scale

        content_logits = jnp.einsum(
            "b h i d, b h j d -> b h i j", q + self.rel_content_bias, k
        )

        positions = self.get_positional_embed(n, self.num_rel_pos_features)
        rel_k = self.to_rel_k(positions)

        def _rearrange_1(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "n (h d) -> h n d", h=self.heads)

        rel_k = _rearrange_1(rel_k)
        rel_logits = jnp.einsum(
            "b h i d, h j d -> b h i j", q + self.rel_pos_bias, rel_k
        )
        rel_logits = relative_shift(rel_logits)

        logits = content_logits + rel_logits
        attn = jax.nn.softmax(logits, axis=-1)

        out = jnp.einsum("b h i j, b h j d -> b h i d", attn, v)

        def _rearrange_2(x: jnp.ndarray) -> jnp.ndarray:
            return rearrange(x, "b h n d -> b n (h d)")

        out = _rearrange_2(out)

        return self.to_out(out)


class EnformerIA3Rescaling(Enformer):
    """
    This class adds rescaling weights following the IA³ methodology to the
    transformer model of Enformer, only in the attention layers.
    This aims to be used for fine-tuning pre-trained Enformer models.
    """

    def __init__(
        self,
        config: EnformerConfig,
        name: Optional[str] = None,
    ):
        super().__init__(config=config, name=name)

    @hk.transparent
    def _attention_layer(self) -> Attention:
        return AttentionIA3Rescaling(  # type: ignore
            dim=self._config.embed_dim,
            heads=self._config.num_attention_heads,
            dim_key=self._config.attention_dim_key,
            dim_value=self._config.embed_dim // self._config.num_attention_heads,
            num_rel_pos_features=self._config.embed_dim
            // self._config.num_attention_heads,
            name="attention",
        )

    def _transformer_ffn_block(self, x: jnp.ndarray, layer_num: int) -> jnp.ndarray:

        with hk.experimental.name_scope(f"layer_{layer_num}"):
            layer_norm = hk.LayerNorm(
                axis=-1, create_scale=True, create_offset=True, name="ffn_layer_norm"
            )
            linear_1 = hk.Linear(output_size=self._config.embed_dim * 2, name="ffn_1")
            linear_2 = hk.Linear(output_size=self._config.embed_dim, name="ffn_2")

            # moved inside hk.experimental.name_scope(f"layer_{layer_num}")
            # for ffn_ia3_rescaling inherit the layer name
            y = layer_norm(x)
            y = linear_1(y)
            y = jax.nn.relu(y)

            # IA3 rescaling
            ffn_ia3_rescaling = hk.get_parameter(
                "ffn_ia3_rescaling",
                shape=[y.shape[-1]],
                dtype=y.dtype,
                init=hk.initializers.Constant(1.0),
            )
            y = y * ffn_ia3_rescaling

            y = linear_2(y)

        return x + y

    def __call__(
        self, tokens: jnp.ndarray, is_training: bool = False
    ) -> Dict[str, jnp.ndarray]:
        """
        Calls the Enformer model.

        Args:
            tokens: Input tokens out of the NucleotidesKmersTokenizer
                of shape (batch_size, seq_len). Token IDs should ONLY
                contain tokens for A,T,C,G,N nucleotides.
            is_training: Pass to True during training, will change the
                behaviour of Convolutional layers.

            # NOTE: the input tokens are expected to be vectors with
                token IDs corresponding to ONLY A,T,C,G,N nucleotides.
                The token IDs of A,T,C,G,N should be the ones by default
                of NucleotidesKmersTokenizer: A:10 / C:12 / G:13 / T:11 / N:14.
                If the token IDs are different, the one-hot encoded vectors
                will not match the nucleotides and it will fail.
                Sequences cannot have the pad token - to pad the sequences
                you can add the nucleotide N.

        Returns:
            Dictionary containing the final embeddings and logits
                for the human and mouse head.

        """

        # one-hot encoding
        x = self._batch_one_hot_encode(tokens)

        # Enformer model
        x = jnp.transpose(x, (0, 2, 1))
        x = self._stem(x, is_training=is_training)
        x = self._conv_tower(x, is_training=is_training)
        x = jnp.transpose(x, (0, 2, 1))
        x = self._transformer_tower(x)
        # return also embeddings before _target_length_crop and _final_pointwise
        outs = {}
        outs["embedding_transformer_tower"] = x

        x = self._target_length_crop(x)
        x = self._final_pointwise(x, is_training=is_training)

        outs["embedding"] = x

        # human and mouse heads
        heads_outs = self._heads(x)
        outs.update(heads_outs)

        return outs


def build_enformer_fn_ia3_rescaling_with_head_fn(
    config: EnformerConfig,
    head_fn: Callable[
        [], Callable[[jnp.ndarray, SequenceMask], Dict[str, jnp.ndarray]]
    ],
    embedding_name: str = "embedding_transformer_tower",
    compute_dtype: jnp.dtype = jnp.float32,
    param_dtype: jnp.dtype = jnp.float32,
    output_dtype: jnp.dtype = jnp.float32,
    name: Optional[str] = None,
) -> Callable:
    """
    Create the model's forward pass for that Enformer,
        adds rescaling weights in the attention layers and the input head.

    Args:
        config: Configuration data class containing the hyperparameters for the Enformer
            forward function.
        head_fn: Wrapper initializing a Classification/Regression head. The head cannot
            be passed directly as haiku modules cannot be initialized outside
            hk.transform.
        embedding_name: Which embeddings to use from the enformer pre-trained model as
            input to the head. It should be the name of the key inside model
            predictions ( outs ). Default is "embedding_transformer_tower".
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
        name: the name of the model. example: enformer.

        Example of the function being used with a classification head:
        The classification head is wrapped inside head_fn because
        haiku modules cannot be instantiated outside hk.transform.
        def head_fn():
            return SimpleClassificationHead(num_classes=num_classes)
        finetune_forward_fn = build_esm_ia3_rescaling_with_head_fn(
            model_config=config, head_fn=head_fn, model_name=model_name,
        )
        finetune_forward_fn = hk.transform(finetune_forward_fn)

        # NOTE: the input tokens for enformer of shape (batch_size, seq_len)
            are expected to be vectors with token IDs corresponding
            to ONLY A,T,C,G,N nucleotides.
            The token IDs of A,T,C,G,N should be the ones by default
            of NucleotidesKmersTokenizer: A:10 / C:12 / G:13 / T:11 / N:14.
            If the token IDs are different, the one-hot encoded vectors
            will not match the nucleotides and it will fail.
            Sequences cannot have the pad token - to pad the sequences
            you can add the nucleotide N.

    Returns:
        Enformer model forward function with IA³ rescaling and indicated head.
    """

    assert {compute_dtype, param_dtype, output_dtype}.issubset(
        {
            jnp.bfloat16,
            jnp.float32,
            jnp.float16,
        }
    ), f"provide a dtype in {jnp.bfloat16, jnp.float32, jnp.float16}"

    policy = jmp.Policy(
        compute_dtype=compute_dtype, param_dtype=param_dtype, output_dtype=output_dtype
    )
    hk.mixed_precision.set_policy(EnformerIA3Rescaling, policy)

    # Remove it in batch norm to avoid instabilities
    norm_policy = jmp.Policy(
        compute_dtype=jnp.float32, param_dtype=param_dtype, output_dtype=compute_dtype
    )
    hk.mixed_precision.set_policy(hk.LayerNorm, norm_policy)
    hk.mixed_precision.set_policy(hk.BatchNorm, norm_policy)

    def enformer_fn(
        tokens: jnp.ndarray,
        is_training: bool = False,
        sequence_mask: Optional[SequenceMask] = None,
    ) -> Dict[str, jnp.ndarray]:
        """Forward pass."""
        # Run the encoder over the inputs.
        encoder = EnformerIA3Rescaling(config, name)
        outs = encoder(tokens, is_training)
        # NOTE: For now we don't remove the Enformer human/mouse prediction heads
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

    return enformer_fn
