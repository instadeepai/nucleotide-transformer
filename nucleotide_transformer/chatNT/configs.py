from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass
class RotaryEmbeddingConfig:
    """
    Rotary Positional Embedding configuration
        max_seq_len: The number of positions to encode and cache.
        dim: Dimension of RoPE.
        theta: Rotation angle.
    """

    max_seq_len: int
    dim: int
    theta: float


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


@dataclass
class ESMTransformerConfig:
    """
    Parameters to initialize an ESM model. While the ESM architecture is an encoder-only
    model, different choices have been made for each version and this configuration aims
    to cover most of them.

    Args:
        alphabet_size: Token vocabulary.
        pad_token_id: ID of pad token.
        mask_token_id: ID of mask token.
        max_positions: Maximum sequence length.
        embed_scale: Correction ratio applied to the embeddings to make up for the
            norm difference between the input during training and inference.
        emb_layer_norm_before: Whether to use layer norm before the first attention
            layer.
        attention_heads: Number of attention heads.
        key_size: The dimension of the query, key, and values within each attention
            head, if not specified, it is set to attention_heads//embed_dim.
            It can be useful to set a custom key size if we want to impose the size of
            the query, key and value tensor ( for example, tensors shaped with
            power of 2 are more efficiently handled on TPUs ).
            Note: Parametrizing the model with a custom key size has been done in :
            Brown, Tom, et al. "Language models are few-shot learners."
            Advances in neural information processing systems 33 (2020): 1877-1901.
        embed_dim: Embedding dimension.
        ffn_embed_dim: Feed forward embedding dimension.
        num_layers: Number of attention blocks.
        positional_embedding: Type of positional embedding to use before the first
            attention layer. Options: "learned", "learned_standard" "sinusoidal" or
            None.
            NOTE: "learned" is the positional embeding of ESM, and "learned_standard" is
             a more standard one, used for example in DNAbert.
        lm_head: type of language model head. Options: "simple", "roberta" or None.
        add_bias_kv: Add bias in attention layer.
        add_bias_ffn: Add bias in feed forward network block.
        use_rotary_embedding: Whether to use rotary embeddings (for ESM2). Requires:
            positional_embeddings = None.
        rescaling_factor: Scaling factor to use for rotary embeddings.
        ffn_activation_name: Activation function to be used in FFN block. Supported
            names are "gelu", "relu", "swish".
        use_glu_in_ffn: Whether to use Gated Linear Unit (GLU) in Feed
            Forward Network (FFN) block. To do a swiGLU (gated-swish) put this arg
            to True and use swish as ffn_activation_name.
            Same principle for a gated-relu. To keep the same number of parameters in
            the FFN block, one should multiply by 2/3 the ffn_embed_dim when using GLU.
            See https://arxiv.org/pdf/2002.05202.pdf for more details.
        mask_before_attention: Use mask before attention layers (for EMS1b and ESM2).
        layer_norm_eps: the eps factor in the different layer norms of the model (refer
            to layer norm implementation)
        token_dropout: Token dropout.
        masking_ratio: Masking ratio (used if token dropout is enabled).
        masking_prob: Masking probability (used if token dropout is enabled).
        use_gradient_checkpointing: Whether to use gradient checkpointing (checkpoint
            gradients in the forward pass to reduce the computation in the backward).
    """

    alphabet_size: int
    pad_token_id: int
    mask_token_id: int

    max_positions: int = 1024
    embed_scale: float = 1.0

    # architecture
    emb_layer_norm_before: bool = False
    attention_heads: int = 20
    key_size: Optional[int] = None
    embed_dim: int = 1280
    ffn_embed_dim: int = 5120
    num_layers: int = 24
    positional_embedding: Optional[str] = "learned"
    lm_head: Optional[str] = "simple"
    add_bias_kv: bool = False
    add_bias_ffn: bool = True
    use_rotary_embedding: bool = False
    rescaling_factor: Optional[float] = None
    ffn_activation_name: str = "gelu-no-approx"
    use_glu_in_ffn: bool = False
    mask_before_attention: bool = False
    layer_norm_eps: float = 1e-5
    pre_layer_norm: bool = True
    bias_word_embedding: bool = False

    # dropout
    token_dropout: bool = False
    masking_ratio: float = 0.1
    masking_prob: float = 0.8

    # logging
    use_gradient_checkpointing: bool = False

    # return
    embeddings_layers_to_save: Tuple[int, ...] = ()
    attention_maps_to_save: List[Tuple[int, int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        """
        Checks that the given values are compatible.
        """

        if self.key_size is None:
            if not self.embed_dim % self.attention_heads == 0:
                raise ValueError(
                    f"When no key size is provided, the embedding dimension should be "
                    f"divisible by the number of heads, however provided embedding "
                    f"dimension is {self.embed_dim} and the number of heads is "
                    f"{self.attention_heads}."
                )
            self.key_size = self.embed_dim // self.attention_heads
        if self.positional_embedding is not None:
            if type(self.positional_embedding) != str:
                raise TypeError

            if self.positional_embedding not in [
                "learned",
                "sinusoidal",
                "learned_standard",
                "alibi_dnabert_2",
            ]:
                raise ValueError(
                    "The positional_embedding argument should either be None,"
                    "`learned`, `sinusoidal`, 'learned_standard' or 'alibi_dnabert_2'."
                )
        if self.lm_head is not None:
            if type(self.lm_head) != str:
                raise TypeError

            if self.lm_head not in ["simple", "roberta"]:
                raise ValueError(
                    "The lm_head argument should either be None,"
                    "`simple` or `roberta`."
                )

        if self.use_rotary_embedding and self.positional_embedding is not None:
            raise ValueError(
                "When using rotary embedding, positional_embedding must be set to none"
            )

        if self.add_bias_kv and self.use_rotary_embedding:
            raise ValueError(
                "Biases on key and values are not compatible with Rotary embeddings."
            )

        if self.positional_embedding == "alibi_dnabert_2":
            assert not self.add_bias_kv


@dataclass
class PerceiverResamplerConfig:
    """
    Parameters to initialize an PerceiverResampler model. Based on the ESM architecture.

    Args:
        emb_layer_norm_before: Whether to use layer norm before the first attention
            layer.
        attention_heads: Number of attention heads.
        key_size: The dimension of the query, key, and values within each attention
            head, if not specified, it is set to attention_heads//embed_dim.
            It can be useful to set a custom key size if we want to impose the size of
            the query, key and value tensor ( for example, tensors shaped with
            power of 2 are more efficiently handled on TPUs ).
            Note: Parametrizing the model with a custom key size has been done in :
            Brown, Tom, et al. "Language models are few-shot learners."
            Advances in neural information processing systems 33 (2020): 1877-1901.
        embed_dim: Embedding dimension.
        ffn_embed_dim: Feed forward embedding dimension.
        num_layers: Number of attention blocks.
        ffn_activation_name: Activation function to be used in FFN block. Supported
            names are "gelu", "relu", "swish".
        use_glu_in_ffn: Whether to use Gated Linear Unit (GLU) in Feed
            Forward Network (FFN) block. To do a swiGLU (gated-swish) put this arg
            to True and use swish as ffn_activation_name.
            Same principle for a gated-relu. To keep the same number of parameters in
            the FFN block, one should multiply by 2/3 the ffn_embed_dim when using GLU.
            See https://arxiv.org/pdf/2002.05202.pdf for more details.
        resampled_length: length of the resampled output of the module
        use_gradient_checkpointing: Whether to use gradient checkpointing (checkpoint
            gradients in the forward pass to reduce the computation in the backward).
    """

    # architecture
    emb_layer_norm_before: bool = False
    attention_heads: int = 20
    key_size: Optional[int] = None
    embed_dim: int = 1280
    ffn_embed_dim: int = 5120
    num_layers: int = 24
    add_bias_kv: bool = False
    add_bias_ffn: bool = True
    ffn_activation_name: str = "gelu-no-approx"
    use_glu_in_ffn: bool = False
    resampled_length: int = 64

    # performance
    use_gradient_checkpointing: bool = False

    def __post_init__(self) -> None:
        """
        Checks that the given values are compatible.
        """

        if self.key_size is None:
            if not self.embed_dim % self.attention_heads == 0:
                raise ValueError(
                    f"When no key size is provided, the embedding dimension should be "
                    f"divisible by the number of heads, however provided embedding "
                    f"dimension is {self.embed_dim} and the number of heads is "
                    f"{self.attention_heads}."
                )
            self.key_size = self.embed_dim // self.attention_heads
