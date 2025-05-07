from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import haiku as hk
import jax
import jax.numpy as jnp
from haiku import initializers

from nucleotide_transformer.chatNT.gpt_rotary import (
    RotaryEmbedding,
    RotaryEmbeddingConfig,
)
from nucleotide_transformer.chatNT.types import (
    AttentionMask,
    Embedding,
    Tokens,
    TransformerOutput,
)


@dataclass
class GptConfig:
    """
    Parameters to initialize a Gpt model.

    NOTE: the pad token is not defined

    Args:
        vocab_size: Token vocabulary.
        eos_token_id: used to stop sentence generation
        embed_dim: Embedding dimension.
        ffn_embed_dim: Feed forward embedding dimension.
        num_heads: Number of attention heads.
        num_kv_heads: Number of key and value heads to support Grouped-Query and
            Multi-Query Attention. If None, the number of key and value heads is
            equal to the number of attention heads.
        num_layers: Number of Decoder layer_stack
        rope_config: The configuration for the rotary positional embeddings
        add_bias_ffn: Add bias in feed forward network block.
        ffn_activation_name: Activation function to be used in FFN block. Supported
            names are "gelu", "gelu-no-approx", "relu", "swish".
        use_glu_in_ffn: whether to use Gated Linear Unit (GLU) in Feed
            Forward Network (FFN) block.
            example: To do a swiGLU (gated-swish) put this arg
            to True and use swish as ffn_activation_name.
            Same principle for a gated-relu.
        add_bias_lm_head: whether to use bias in the final LM layer
        norm_type: The type of norm used ( pre normalization scheme ) used. can be
            one of ["layer_norm", "RMS_norm"]
        parallel_attention_ff: Whether to do the attention and the MLP in parallel,
            and then sum up the results as it is done in Gpt-NeoX :
            Black, Sid, et al. "Gpt-neox-20b: An open-source autoregressive
            language model." arXiv preprint arXiv:2204.06745 (2022).
            It is said to improve the training time of 15% when compiling with JAX
        use_gradient_checkpointing: Whether to use gradient checkpointing (checkpoint
            gradients in the forward pass to reduce the computation in the backward).
        add_bias_attn: Add bias to the attention mechanism (key, query, value, and
            output projections).
    """

    # vocabulary
    vocab_size: int
    eos_token_id: int

    # architecture
    embed_dim: int = 16
    ffn_embed_dim: int = 64
    num_heads: int = 2
    num_kv_heads: Optional[int] = None
    num_layers: int = 2
    rope_config: RotaryEmbeddingConfig = field(
        default_factory=lambda: RotaryEmbeddingConfig(
            max_seq_len=512, dim=16 // 2, theta=10000.0
        )
    )
    add_bias_ffn: bool = False
    ffn_activation_name: str = "swish"
    use_glu_in_ffn: bool = True
    add_bias_lm_head: bool = False
    norm_type: str = "RMS_norm"
    rms_norm_eps: float = 1e-6
    parallel_attention_ff: bool = True

    # inference / backward behavior
    use_gradient_checkpointing: bool = False

    # architecture params with default values
    add_bias_attn: bool = False

    def __post_init__(self) -> None:
        """
        Checks that the given values are compatible.
        """
        if not self.embed_dim % self.num_heads == 0:
            raise ValueError(
                f"The embedding dimension should be "
                f"divisible by the number of heads, however provided embedding "
                f"dimension is {self.embed_dim} and the number of heads is "
                f"{self.num_heads}."
            )

        if not self.embed_dim // self.num_heads > 1:
            raise ValueError(
                "embed_dim / num_heads must be higher than 2 to apply rotary embeddings"
            )

        if not self.embed_dim // self.num_heads >= self.rope_config.dim:
            raise ValueError(
                "embed_dim // num_heads must be higher than rope_config.dim "
                "to apply rotary embeddings"
            )


class SimpleLMHead(hk.Module):
    """
    Basic Language Model head. Transforms final attention block output
    into a distribution over tokens at each sequence position.
    """

    def __init__(
        self,
        embed_dim: int,
        alphabet_size: int,
        add_bias_lm_head: bool = True,
        name: Optional[str] = None,
    ):
        """
        Args:
            embed_dim: Embedding dimension.
            alphabet_size: Number of tokens in the alphabet.
            name: Name of the layer. Defaults to None.
        """
        super().__init__(name=name)
        self.embed_dim = embed_dim
        self.alphabet_size = alphabet_size

        # Define layers
        w_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        b_init = initializers.VarianceScaling(2.0, "fan_in", "uniform")
        self._final_fc = hk.Linear(
            self.alphabet_size,
            w_init=w_init,
            b_init=b_init,
            with_bias=add_bias_lm_head,
            name="lm_final_fc",
        )

    def __call__(self, x: jnp.ndarray) -> Dict[str, jnp.ndarray]:
        # Compute logits
        logits = self._final_fc(x)
        return {"logits": logits}


class GptGroupedQueryAttention(hk.Module):
    """
    Generalised attention module supporting Multi-Head Attention, Multi-Query Attention
    and Grouped-Query Attention with masking applied. Computes the keys, queries, and
    values from the input embeddings. Future versions could compute and store these
    to prevent redundant calculations during autoregressive inference.
    """

    def __init__(
        self,
        embed_dim: int,
        num_heads: int,
        rope_config: RotaryEmbeddingConfig,
        num_kv_heads: Optional[int] = None,
        head_dim: Optional[int] = None,
        name: Optional[str] = "attention",
        add_bias_attn: bool = False,
    ):
        """
        Initializes the attention layer.

        Args:
            embed_dim: Length of the token embedding at each position in the sequence.
            num_heads: Number of independent attention heads.
            rope_config: The configuration for the rotary positional embeddings
            num_kv_heads: Number of independent key-value heads.
            head_dim: The dimension of the heads in the attention mechanism.
            name: Optional name for this module.
            add_bias_attn: Add bias to the attention mechanism (key, query, value, and
                output projections).
        """
        super().__init__(name=name)
        self.num_heads = num_heads
        self.num_kv_heads = num_kv_heads or num_heads
        self.embed_dim = embed_dim
        self.head_dim = head_dim or self.embed_dim // self.num_heads
        self.rope_config = rope_config
        self.add_bias_attn = add_bias_attn

        self.query_linear = hk.Linear(
            output_size=self.head_dim * self.num_heads,
            with_bias=self.add_bias_attn,
            name="query_linear",
        )

        self.key_linear = hk.Linear(
            output_size=self.head_dim * self.num_kv_heads,
            with_bias=self.add_bias_attn,
            name="key_linear",
        )

        self.value_linear = hk.Linear(
            output_size=self.head_dim * self.num_kv_heads,
            with_bias=self.add_bias_attn,
            name="value_linear",
        )

        self.out_linear = hk.Linear(
            output_size=embed_dim,
            with_bias=self.add_bias_attn,
            name="out_linear",
        )

        self.rotary_embedding = RotaryEmbedding(self.rope_config)

    def __call__(
        self,
        query_inputs: jnp.ndarray,
        key_inputs: jnp.ndarray,
        value_inputs: jnp.ndarray,
        attention_mask: Optional[jnp.ndarray],
    ) -> jnp.ndarray:
        """
        Computes the result of multiheaded dot-product attention, using
        pre-computed projections for the queries, keys, and values.

        Args:
            query_inputs: Embeddings that will be projected to become the queries.
            key_inputs: Embeddings that will be projected to become the keys.
            value_inputs: Embeddings that will be projected to become the values.
            attention_mask: Mask to be applied in the attention layers.
                Triangular for autoregressive models.
                shape : (1, 1, seq_len, seq_len)

        Returns:
            The standard output of multi-headed attention
        """

        # Compute queries, keys and values:
        queries = self.query_linear(
            query_inputs
        )  # (batch_size, seq_len, num_heads * embed_dim)
        keys = self.key_linear(
            key_inputs
        )  # (batch_size, seq_len, num_kv_heads * embed_dim)
        values = self.value_linear(
            value_inputs
        )  # (batch_size, seq_len, num_kv_heads * embed_dim)

        queries = queries.reshape(
            queries.shape[0], queries.shape[1], self.num_heads, self.head_dim
        )
        keys = keys.reshape(
            keys.shape[0], keys.shape[1], self.num_kv_heads, self.head_dim
        )
        values = values.reshape(
            values.shape[0], values.shape[1], self.num_kv_heads, self.head_dim
        )

        # Compute RoPE for queries and keys
        keys, queries = self.rotary_embedding(keys, queries)

        # Repeat keys and values to match the number of query-heads
        n_rep = self.num_heads // self.num_kv_heads
        keys = jnp.repeat(keys, n_rep, axis=2)
        values = jnp.repeat(values, n_rep, axis=2)

        # Compute attention
        attention_logits = jnp.einsum("...thd,...Thd->...htT", queries, keys)
        sqrt_key_size = jnp.sqrt(self.head_dim).astype(queries.dtype)
        attention_logits = attention_logits / sqrt_key_size
        attention_logits = jnp.where(attention_mask, attention_logits, -1e30)
        attention_weights = jax.nn.softmax(attention_logits, axis=-1)
        values = jnp.einsum("...htT,...Thd->...thd", attention_weights, values)
        values = jnp.reshape(values, (values.shape[0], values.shape[1], -1))

        return self.out_linear(values)


class GptDecoderLayer(hk.Module):
    """
    Single layer in the encoder, including self-attention and feed-forward operations.
    The feed-forward network uses a ReLU activation and has no biases.
    """

    def __init__(
        self,
        embed_dim: int,
        ffn_embed_dim: int,
        num_heads: int,
        rope_config: RotaryEmbeddingConfig,
        norm_type: str,
        parallel_attention_ff: bool,
        add_bias_ffn: bool,
        ffn_activation_name: str,
        use_glu_in_ffn: bool,
        num_kv_heads: Optional[int],
        rms_norm_eps: float = 1e-6,
        add_bias_attn: bool = False,
        name: Optional[str] = None,
    ):
        """
        Initializes the encoder layer, including the projections needed for
        self-attention and the linear layers applied in the fully connected portion

        Args:
            embed_dim: Dimension of the embeddings
            ffn_embed_dim: Dimension of the hidden layer in the MLP
            num_heads: Number of independent attention heads.
            rope_config: The configuration for the rotary positional embeddings
            norm_type: The type of norm used ( pre normalization scheme ) used. can be
                one of ["layer_norm", "RMS_norm"]
            parallel_attention_ff: Whether to do the attention and the MLP in parallel,
                and then sum up the results as it is done in Gpt-NeoX :
                Black, Sid, et al. "Gpt-neox-20b: An open-source autoregressive
                language model." arXiv preprint arXiv:2204.06745 (2022).
                It is said to improve the training time of 15% when compiling with JAX
            add_bias_ffn: Add bias in feed forward network block.
            ffn_activation_name: Activation function to be used in FFN block. Supported
                names are "gelu", "gelu-no-approx", "relu", "swish", and "silu"
            use_glu_in_ffn: Whether to use Gated Linear Unit (GLU) in Feed
                Forward Network (FFN) block. To do a swiGLU (gated-swish) put this arg
                to True and use swish as ffn_activation_name.
                Same principle for a gated-relu.
            num_kv_heads: Number of independent key-value heads.
            name: Optional name for this module.
            add_bias_attn: Add bias to the attention mechanism (key, query, value, and
                output projections).
        """
        super().__init__(name=name)

        self.num_heads = num_heads
        self.parallel_attention_ff = parallel_attention_ff

        self.sa_layer = GptGroupedQueryAttention(
            embed_dim=embed_dim,
            num_heads=num_heads,
            num_kv_heads=num_kv_heads,
            name="self_attn",
            rope_config=rope_config,
            add_bias_attn=add_bias_attn,
        )

        if norm_type == "layer_norm":
            self.attn_norm = hk.LayerNorm(
                axis=-1, create_scale=True, create_offset=True, name="attn_layer_norm"
            )
            if not (self.parallel_attention_ff):
                self.ffn_norm = hk.LayerNorm(
                    axis=-1,
                    create_scale=True,
                    create_offset=True,
                    name="ffn_layer_norm",
                )
        elif norm_type == "RMS_norm":
            self.attn_norm = hk.RMSNorm(
                axis=-1, create_scale=True, name="attn_RMS_norm", eps=rms_norm_eps
            )
            if not (self.parallel_attention_ff):
                self.ffn_norm = hk.RMSNorm(
                    axis=-1, create_scale=True, name="ffn_RMS_norm", eps=rms_norm_eps
                )
        else:
            raise ValueError(f"unrecognized norm_type : {norm_type}")

        # Get ffn activation function
        self._ffn_activation_fn = get_activation_fn(activation_name=ffn_activation_name)
        self._use_glu_in_fnn = use_glu_in_ffn

        # Define layers
        if use_glu_in_ffn:
            # user should multiply ffn_embed_dim by 2/3 when using GLU
            # to keep total number of parameters equal
            # see https://arxiv.org/pdf/2002.05202.pdf. for more details
            # we multiply by 2 here as the output will be split in 2 for GLU
            ffn_embed_dim = int(2 * ffn_embed_dim)

        self.fc1_linear = hk.Linear(
            output_size=ffn_embed_dim,
            with_bias=add_bias_ffn,
            name="fc1_linear_glu" if use_glu_in_ffn else "fc1_linear",
        )
        self.fc2_linear = hk.Linear(
            output_size=embed_dim, with_bias=add_bias_ffn, name="fc2_linear"
        )

    @hk.transparent
    def mlp(self, x: Embedding) -> Embedding:
        """
        Applies one linear layer, a ReLU activation, dropout, then a final linear layer.

        Args:
            x: Embeddings of shape (batch_size, seq_len, embed_dim).

        Returns:
            The transformed sequence embedding.
        """
        if self._use_glu_in_fnn:
            x1, x2 = jnp.split(self.fc1_linear(x), indices_or_sections=2, axis=-1)
            x = self._ffn_activation_fn(x1) * x2
        else:
            x = self._ffn_activation_fn(self.fc1_linear(x))

        x = self.fc2_linear(x)
        return x

    def __call__(
        self,
        embeddings: Embedding,
        attention_mask: AttentionMask,
    ) -> Embedding:
        """
        Computes the output embeddings of the encoder layer.
        if self.parallel_attention_ff, the model uses parallel MLP and attention

        Args:
            embeddings: Decoder layer input embeddings of shape
                (batch_size,seq_len,embed_dim).
            attention_mask: Mask to be applied in the attention layers.
                Triangular for autoregressive models.
                shape = (1, 1, seq_len, seq_len)

        Returns:
            The output embeddings that result from the application of the layer
        """
        if self.parallel_attention_ff:
            residuals = embeddings

            embeddings = self.attn_norm(embeddings)
            attn_outputs = self.sa_layer(
                query_inputs=embeddings,
                key_inputs=embeddings,
                value_inputs=embeddings,
                attention_mask=attention_mask,
            )
            mlp_ouputs = self.mlp(embeddings)
            return residuals + attn_outputs + mlp_ouputs
        else:
            normed_embeddings = self.attn_norm(embeddings)
            attn_outputs = embeddings + self.sa_layer(
                query_inputs=normed_embeddings,
                key_inputs=normed_embeddings,
                value_inputs=normed_embeddings,
                attention_mask=attention_mask,
            )
            mlp_ouputs = attn_outputs + self.mlp(self.ffn_norm(attn_outputs))

            return mlp_ouputs


class GptDecoder(hk.Module):
    """
    Creates the Gpt model ( decoder only ).
    """

    def __init__(
        self,
        config: GptConfig,
        name: Optional[str] = None,
    ):
        """
        Initializes the Decoder stack.

        Args:
            config: Configuration data class
            name: haiku module name
        """

        self._config = config

        super().__init__(name=name)

        self.token_embed = hk.Embed(
            vocab_size=config.vocab_size, embed_dim=config.embed_dim, name="token_embed"
        )
        self.lm_head = SimpleLMHead(
            embed_dim=config.embed_dim,
            alphabet_size=config.vocab_size,
            add_bias_lm_head=self._config.add_bias_lm_head,
        )
        if self._config.norm_type == "layer_norm":
            self.final_norm = hk.LayerNorm(
                axis=-1, create_scale=True, create_offset=True, name="final_layer_norm"
            )
        elif self._config.norm_type == "RMS_norm":
            self.final_norm = hk.RMSNorm(
                axis=-1, create_scale=True, name="final_RMS_norm"
            )
        else:
            raise ValueError(
                f"unrecognized norm_type in config {self._config.norm_type}"
            )

    @hk.experimental.name_like("__call__")
    def decoder_layer(self, layer_idx: int) -> GptDecoderLayer:
        """
        Returns the Gpt encoder layer.

        Args:
            layer_idx: the layer index

        Returns:
            The named GptDecoderLayer
        """
        return GptDecoderLayer(
            embed_dim=self._config.embed_dim,
            ffn_embed_dim=self._config.ffn_embed_dim,
            num_heads=self._config.num_heads,
            rope_config=self._config.rope_config,
            norm_type=self._config.norm_type,
            parallel_attention_ff=self._config.parallel_attention_ff,
            add_bias_ffn=self._config.add_bias_ffn,
            ffn_activation_name=self._config.ffn_activation_name,
            use_glu_in_ffn=self._config.use_glu_in_ffn,
            name=f"gpt_decoder_layer_{layer_idx}",
            num_kv_heads=self._config.num_kv_heads,
            add_bias_attn=self._config.add_bias_attn,
            rms_norm_eps=self._config.rms_norm_eps,
        )

    @hk.transparent
    def apply_transformer_layers(
        self,
        tokens_embeddings: Embedding,
        attention_mask: Optional[AttentionMask] = None,
    ) -> Embedding:
        """
        Takes as inputs the tokens embeddings and apply successively an attention mask
        and the transformer layers to obtain final embeddings ready to be decoded.

        Args:
            tokens_embeddings: Embeddings of the input token ids.
            attention_mask: Optional attention mask to allow for stacked examples or
                other non-standard behavior. Otherwise, will create default causal mask.

        Returns:
            Embeddings transformed through successive transformer layers.
        """

        # go through the transformer layer_stack
        layer_stack = [self.decoder_layer(i) for i in range(self._config.num_layers)]

        # use gradient checkpointing if required
        if self._config.use_gradient_checkpointing:
            # the remat-ed function cannot take control flow arguments
            layer_stack = [hk.remat(layer) for layer in layer_stack]

        if attention_mask is None:
            attention_mask = build_causal_attention_mask(1, tokens_embeddings.shape[1])

        # cast mask to current mixed precision policy
        attention_mask = attention_mask.astype(tokens_embeddings.dtype)

        embeddings = tokens_embeddings
        for i in range(len(layer_stack)):
            embeddings = layer_stack[i](
                embeddings=embeddings,
                attention_mask=attention_mask,
            )

        return embeddings

    def __call__(
        self, token_ids: Tokens, attention_mask: Optional[jnp.ndarray] = None
    ) -> TransformerOutput:
        """
        Compute the logits and embeddings from a sequence of tokens

        Args:
            token_ids: Sequence of token ids delivered to the decoder of shape
                (batch_size, seq_len)
            attention_mask: Optional attention mask to allow for stacked examples or
                other non-standard behavior. Otherwise, will create default causal mask.

        Returns:
             The logits over the token vocabulary at each time step.
        """

        # (batch_size,seq_len) -> (batch_size,seq_len,embed_dim)
        tokens_embeddings = self.token_embed(token_ids)

        # (batch_size,seq_len,embed_dim) -> (batch_size,seq_len,embed_dim)
        embeddings = self.apply_transformer_layers(tokens_embeddings, attention_mask)
        embeddings = self.final_norm(embeddings)

        # Get outputs
        outs = {}
        outs["embeddings"] = embeddings

        # Compute logits
        logits = self.lm_head(embeddings)["logits"]

        outs["logits"] = logits

        return outs  # type: ignore


SUPPORTED_FFN_ACTIVATIONS = ["gelu", "gelu-no-approx", "relu", "swish", "silu", "sin"]


def get_activation_fn(activation_name: str) -> Callable:
    """
    Return activation fn given its name.
    Args:
        activation_name: Activation name.

    Returns:
        activation function.
    """
    if activation_name not in SUPPORTED_FFN_ACTIVATIONS:
        raise NotImplementedError(
            f"Activation {activation_name} not supported yet. "
            f"Supported activations for feed forward "
            f"block are {SUPPORTED_FFN_ACTIVATIONS}"
        )
    if activation_name == "gelu-no-approx":
        activation_fn = lambda x: jax.nn.gelu(x, approximate=False)  # noqa: E731
    elif activation_name == "sin":
        activation_fn = lambda x: jnp.sin(x)  # noqa: E731
    else:
        activation_fn = getattr(jax.nn, activation_name)
    return activation_fn


def build_causal_attention_mask(batch_size: int, seq_len: int) -> AttentionMask:
    """
    Builds a batch of causal masks of shape (batch_size, 1, seq_len, seq_len) to feed
    to an attention layer.

    Args:
        batch_size: Batch size.
        seq_len: Length of the sequences.

    Returns:
        Batch of causal masks.
    """
    mask = jnp.ones((batch_size, 1, seq_len, seq_len))
    causal_mask = jnp.tril(mask)
    return causal_mask
