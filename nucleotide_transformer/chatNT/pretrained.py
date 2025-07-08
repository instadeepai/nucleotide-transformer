from transformers import AutoTokenizer

from nucleotide_transformer.chatNT.gpt_decoder import GptConfig, RotaryEmbeddingConfig
from nucleotide_transformer.chatNT.model import build_chat_nt_fn
from nucleotide_transformer.chatNT.multi_modal_perceiver_projection import (
    PerceiverResamplerConfig,
)
from nucleotide_transformer.chatNT.params import download_ckpt
from nucleotide_transformer.model import NucleotideTransformerConfig


def get_chatNT():
    gpt_config = GptConfig(
        vocab_size=32000,
        eos_token_id=2,
        embed_dim=4096,
        ffn_embed_dim=11008,
        num_heads=32,
        num_kv_heads=32,
        num_layers=32,
        rope_config=RotaryEmbeddingConfig(max_seq_len=2048, dim=128, theta=10000.0),
        add_bias_ffn=False,
        ffn_activation_name="silu",
        use_glu_in_ffn=True,
        add_bias_lm_head=False,
        norm_type="RMS_norm",
        rms_norm_eps=1e-06,
        parallel_attention_ff=False,
        use_gradient_checkpointing=False,
        add_bias_attn=False,
    )
    nt_config = NucleotideTransformerConfig(
        alphabet_size=4107,
        pad_token_id=1,
        mask_token_id=2,
        max_positions=2048,
        embed_scale=1.0,
        emb_layer_norm_before=False,
        attention_heads=16,
        key_size=64,
        embed_dim=1024,
        ffn_embed_dim=4096,
        num_layers=29,
        positional_embedding=None,
        lm_head="roberta",
        add_bias_kv=False,
        add_bias_ffn=False,
        use_rotary_embedding=True,
        rescaling_factor=None,
        ffn_activation_name="swish",
        use_glu_in_ffn=True,
        mask_before_attention=False,
        layer_norm_eps=1e-05,
        pre_layer_norm=True,
        bias_word_embedding=False,
        token_dropout=False,
        masking_ratio=0.0,
        masking_prob=0.0,
        use_gradient_checkpointing=False,
        embeddings_layers_to_save=(21,),
        attention_maps_to_save=[],
    )
    perceiver_resampler_config = PerceiverResamplerConfig(
        emb_layer_norm_before=False,
        attention_heads=32,
        key_size=128,
        embed_dim=4096,
        ffn_embed_dim=11008,
        num_layers=3,
        add_bias_kv=False,
        add_bias_ffn=True,
        ffn_activation_name="gelu-no-approx",
        use_glu_in_ffn=False,
        resampled_length=64,
        use_gradient_checkpointing=False,
    )

    bio_tokenizer = AutoTokenizer.from_pretrained(
        "InstaDeepAI/ChatNT",
        subfolder="bio_tokenizer",
    )
    english_tokenizer = AutoTokenizer.from_pretrained(
        "InstaDeepAI/ChatNT",
        subfolder="english_tokenizer",
    )
    seq_token_id = 32000

    forward_fn = build_chat_nt_fn(
        nt_config=nt_config,
        gpt_config=gpt_config,
        seq_token_id=seq_token_id,
        bio_pad_token_id=bio_tokenizer.pad_token_id,
        english_pad_token_id=english_tokenizer.pad_token_id,
        perceiver_resampler_config=perceiver_resampler_config,
        nt_name="dcnuc_v2_500M_multi_species",
        gpt_name="llama_decoder",
    )

    params = download_ckpt()
    params = {
        key.replace("bio_brain_decoder", "chat_nt_decoder").replace(
            "bio_brain_encoder", "chat_nt_encoder"
        ): value
        for key, value in params.items()
    }

    return forward_fn, params, english_tokenizer, bio_tokenizer
