import logging
import math
from typing import Any, Literal, Optional

import numpy as np
from pydantic import BaseModel, Field, root_validator


class MOJOConfig(BaseModel):
    name: Literal["MOJO"]
    alphabet_size: dict[str, int] = Field(default_factory=dict)
    token_embed_dim: int
    conv_init_embed_dim: int
    embed_dim: int
    num_downsamples: int
    filter_list: list = Field(default_factory=list)
    init_gene_embed_dim: int = 200
    project_gene_embedding: bool = False
    use_gene_embedding: bool = True

    sequence_length: int  # n_genes
    fixed_sequence_length: int | None = None

    use_remat_in_convs: bool = False
    use_remat_in_transformer: bool = False
    use_skip_connection: bool = True

    embeddings_layers_to_save: tuple[int, ...] = ()
    attention_layers_to_save: tuple[int, ...] = ()

    gene2vec_weights_path: str

    num_attention_heads: int = 16
    key_size: Optional[int] = None
    ffn_embed_dim: int = 512
    num_layers: int = 4
    layer_norm_eps: float = 1e-5
    stem_kernel_shape: int = 15
    num_hidden_layers_head: int = 0

    n_expressions_bins: dict[str, int]
    min_omic_value: dict[str, float]
    max_omic_value: dict[str, float]
    use_log_normalization: dict[str, bool]
    use_max_normalization: dict[str, bool]
    normalization_factor: dict[str, float]

    @root_validator
    @classmethod
    def validate_key_size(cls, values: dict[str, Any]) -> dict[str, Any]:
        """
        Checks that the given values are compatible.
        """
        key_size = values.get("key_size")
        if key_size is None:
            embed_dim = values["embed_dim"]
            num_attention_heads = values["num_attention_heads"]
            if not embed_dim % num_attention_heads == 0:
                raise ValueError(
                    f"When no key size is provided, the embedding dimension should be "
                    f"divisible by the number of heads, however provided embedding "
                    f"dimension is {embed_dim} and the number of heads is "
                    f"{num_attention_heads}."
                )
            values["key_size"] = embed_dim // num_attention_heads
        return values

    @root_validator
    @classmethod
    def create_filter_list(cls, values: dict[str, Any]) -> dict[str, Any]:
        """
        Checks that the given values are compatible.
        """
        num_downsamples: int = values["num_downsamples"]
        filter_list = (
            np.linspace(
                values.get("conv_init_embed_dim"),
                values.get("embed_dim"),
                num_downsamples + 1,
            )
            .astype(int)
            .tolist()
        )

        values["filter_list"] = filter_list
        return values

    @root_validator
    @classmethod
    def validate_gene_embedding(cls, values: dict[str, Any]) -> dict[str, Any]:
        """
        Checks that the given values are compatible.
        """
        use_gene_embedding = values.get("use_gene_embedding")
        if use_gene_embedding:
            init_gene_embed_dim = values["init_gene_embed_dim"]
            token_embed_dim = values["token_embed_dim"]
            if init_gene_embed_dim != token_embed_dim:
                project_gene_embedding = values["project_gene_embedding"]
                if not project_gene_embedding:
                    logging.warning(
                        f"Init gene embedding dimension ({init_gene_embed_dim})"
                        f"different than token embedding dimension ({token_embed_dim})."
                        f"Setting `project_gene_embedding` to True"
                    )
                    values["project_gene_embedding"] = True
        return values

    @root_validator
    @classmethod
    def compute_fixed_sequence_length(cls, values: dict[str, Any]) -> dict[str, Any]:
        num_downsamples: int = values["num_downsamples"]
        sequence_length: int = values["sequence_length"]
        downsample_factor = 2**num_downsamples
        fixed_sequence_length = (
            math.ceil(sequence_length / downsample_factor) * downsample_factor
        )
        values["fixed_sequence_length"] = fixed_sequence_length
        return values
