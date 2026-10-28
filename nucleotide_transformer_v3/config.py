# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""NTv3 model configuration."""

import dataclasses
from typing import Literal

import numpy as np
from flax.typing import Dtype


@dataclasses.dataclass(kw_only=True)
class NTv3Config:
    """Configuration for the NTv3 model.

    U-Net encoder/decoder with a transformer tower at the bottleneck. Uses
    grouped-query attention (GQA) with per-head Q/K/V LayerNorms, genomic
    rotary position embeddings, GeGLU FFN, ``nnx.RMSNorm`` around each
    attention / MLP block, and a final RMSNorm after the trunk. Attention is
    computed via ``jax.nn.dot_product_attention``; the backend is selected by
    ``attention_implementation`` (both backends require K and V head dims to
    be equal).

    All parameters and activations default to float32. Per-submodule
    ``*_compute_dtype``/``*_param_dtype`` fields (mirroring
    ``LongRangeNTConfig``) can lower the *compute* dtype of the matmul/conv
    heavy paths — stem, down/up convolutions, transformer QKVO and transformer
    FFN — while ``*_param_dtype`` still defaults to float32, so stored weights
    are unaffected either way. Embeddings, all norms, the FiLM modulation
    projections and the LM head are deliberately not configurable here and
    always run in float32 (Adam moments and small-vocab MLM logits don't
    tolerate lower precision cleanly there). The one other exception is the cuDNN
    attention path (``attention_implementation="cudnn"``), which casts Q/K/V
    to bfloat16 for the kernel only (cuDNN has no fp32 FlashAttention).

    Training runs in float32; the compute dtypes are an *inference* knob. They
    are baked into each submodule at construction time, so switching an
    already-loaded model over means rebuilding it - see
    :func:`nucleotide_transformer_v3.dtypes.build_bf16_compute_model`.

    With the compute dtypes lowered, ``logits``, ``embedding`` and
    ``after_transformer_embedding`` still come back float32 (a norm or the LM
    head sits between them and the last low-precision matmul), but the
    encoder-side outputs ``before_transformer_embedding`` and
    ``embeddings_conv_{i}`` are the raw conv stream and come back in the
    convolution's compute dtype.

    Args:
        alphabet_size: number of possible tokens.
        pad_token_id: id of pad token.
        mask_token_id: id of mask token.
        token_embed_dim: token embedding dimension.
        conv_init_embed_dim: embedding dimension at the encoder start (after
            the stem). The encoder progressively widens up to ``embed_dim``
            at the trunk.
        embed_dim: embedding dimension in the Transformer torso (trunk).
        num_downsamples: number of downsampling steps in the encoder. Each
            step halves the sequence length. Total downsampling factor is
            ``2 ** (num_downsamples + 1)`` (includes an initial 2x pool).
            The decoder has ``num_downsamples + 1`` upsampling blocks.
        num_layers: number of Transformer layers.
        attention_heads: number of query attention heads.
        num_kv_heads: number of K/V heads (grouped-query attention).
            ``attention_heads`` must be divisible by ``num_kv_heads``. With
            ``num_kv_heads=2`` and ``attention_heads=8`` each K/V head is
            shared across a group of 4 Q heads.
        key_size: per-head dimension for Q/K projections. Must equal
            ``value_size`` because ``jax.nn.dot_product_attention`` requires
            K and V head dims to be equal.
        value_size: per-head dimension for V projections (must equal
            ``key_size``).
        attention_implementation: backend for
            ``jax.nn.dot_product_attention``. ``"xla"`` (default) is the
            portable fp32 path that materializes the S x S scores (O(S^2)
            memory). ``"cudnn"`` uses the cuDNN FlashAttention kernel
            (O(S) memory) and casts Q/K/V to bfloat16 for the kernel
            (cuDNN has no fp32 support); GPU-only and experimental.
        use_remat_in_transformer: whether to use gradient checkpointing in
            the transformer tower.
        use_remat_in_convs: whether to use gradient checkpointing in the
            conv encoder/decoder.
        conv_layers_to_save: 1-based indices of encoder blocks whose output
            should be returned under ``embeddings_conv_{i}`` keys.
        embeddings_layers_to_save: 1-based indices of transformer layers
            whose output should be returned under ``embeddings_{i}`` keys.
        deconv_layers_to_save: 1-based indices of decoder blocks whose
            output should be returned under ``embeddings_deconv_{i}`` keys.
        conditions_dims: dimensions of the per-example condition embeddings
            used for FiLM trunk conditioning (e.g. species). An empty tuple
            (the default) builds a fully unconditioned model whose parameter
            tree and forward pass are byte-identical to the base NTv3. When
            non-empty, each element is the embedding dim of one condition;
            FiLM (adaptive) norms in the transformer pre-norms (``norm_in`` /
            MLP ``norm1``) and in the decoder conv blocks consume them, while
            the stem, encoder, QKV norms, attention/MLP post-norms and the
            final trunk norm stay unconditioned.
        stem_compute_dtype: dtype for the stem (``TokenDnaEmbedder``)
            computation.
        stem_param_dtype: dtype for the stem parameters.
        down_convolution_compute_dtype: dtype for the encoder
            (``DownResBlock``) convolution computation.
        down_convolution_param_dtype: dtype for the encoder convolution
            parameters.
        up_convolution_compute_dtype: dtype for the decoder (``UpResBlock``)
            convolution computation.
        up_convolution_param_dtype: dtype for the decoder convolution
            parameters.
        transformer_qkvo_compute_dtype: dtype for the attention QKVO
            projections' computation.
        transformer_qkvo_param_dtype: dtype for the attention QKVO
            projections' parameters.
        transformer_ffn_compute_dtype: dtype for the FFN (``MLPBlock``)
            computation.
        transformer_ffn_param_dtype: dtype for the FFN parameters.
    """

    # Tokenizer info
    alphabet_size: int = 11
    pad_token_id: int = 1
    mask_token_id: int = 2
    token_embed_dim: int = 16

    # Architecture
    conv_init_embed_dim: int = 1024
    embed_dim: int = 2048
    num_downsamples: int = 6
    num_layers: int = 14
    attention_heads: int = 8
    num_kv_heads: int = 2
    key_size: int = 192
    value_size: int = 192
    attention_implementation: Literal["xla", "cudnn"] = "xla"

    # Optimization
    use_remat_in_transformer: bool = False
    use_remat_in_convs: bool = False

    # Intermediate outputs to save
    conv_layers_to_save: tuple[int, ...] = ()
    embeddings_layers_to_save: tuple[int, ...] = ()
    deconv_layers_to_save: tuple[int, ...] = ()

    conditions_dims: tuple[int, ...] = ()

    # dtypes (default float32 everywhere; see class docstring)
    stem_param_dtype: Dtype = "float32"
    stem_compute_dtype: Dtype = "float32"
    down_convolution_param_dtype: Dtype = "float32"
    down_convolution_compute_dtype: Dtype = "float32"
    up_convolution_param_dtype: Dtype = "float32"
    up_convolution_compute_dtype: Dtype = "float32"
    transformer_qkvo_param_dtype: Dtype = "float32"
    transformer_qkvo_compute_dtype: Dtype = "float32"
    transformer_ffn_param_dtype: Dtype = "float32"
    transformer_ffn_compute_dtype: Dtype = "float32"

    def __post_init__(self) -> None:
        if self.attention_heads % self.num_kv_heads != 0:
            raise ValueError(
                f"attention_heads ({self.attention_heads}) must be divisible "
                f"by num_kv_heads ({self.num_kv_heads})"
            )
        if self.key_size != self.value_size:
            raise ValueError(
                "jax.nn.dot_product_attention requires K and V head dims to "
                f"be equal; got key_size={self.key_size}, "
                f"value_size={self.value_size}."
            )
        if self.attention_implementation not in ("xla", "cudnn"):
            raise ValueError(
                "attention_implementation must be 'xla' or 'cudnn'; got "
                f"{self.attention_implementation!r}."
            )
        if self.conv_init_embed_dim > self.embed_dim:
            # The encoder progressively widens channels via np.linspace, and
            # DownResBlock's residual path only supports out_channels >=
            # in_channels (it zero-pads along the feature axis to match).
            raise ValueError(
                f"conv_init_embed_dim ({self.conv_init_embed_dim}) must be "
                f"<= embed_dim ({self.embed_dim}); the encoder widens "
                "monotonically and the residual cannot shrink channels."
            )

    @property
    def total_downsampling(self) -> int:
        """Total downsampling factor (initial pool + encoder blocks)."""
        result: int = 2 ** (self.num_downsamples + 1)
        return result

    @property
    def encoder_channels(self) -> list[int]:
        """Channel progression through the encoder.

        Returns ``num_downsamples + 1`` values linearly interpolating from
        ``conv_init_embed_dim`` to ``embed_dim``. The linspace is cast to
        ``int``, so intermediate stages are not constrained to powers of 2
        (default 1024 -> 2048 with 6 downsamples yields
        ``[1024, 1194, 1365, 1536, 1706, 1877, 2048]``).
        """
        channels: list[int] = (
            np.linspace(
                self.conv_init_embed_dim, self.embed_dim, self.num_downsamples + 1
            )
            .astype(int)
            .tolist()
        )
        return channels

    @property
    def decoder_x_channels(self) -> list[int]:
        """Input channel sizes for each decoder UpResBlock."""
        ch = self.encoder_channels
        n = self.num_downsamples
        # First block takes trunk (ch[-1]), rest take previous decoder output.
        return [ch[n]] + [ch[n]] + [ch[i] for i in range(n - 1, 0, -1)]

    @property
    def decoder_skip_channels(self) -> list[int]:
        """Skip-connection channel sizes for each decoder UpResBlock."""
        ch = self.encoder_channels
        n = self.num_downsamples
        return [ch[i] for i in range(n, -1, -1)]
