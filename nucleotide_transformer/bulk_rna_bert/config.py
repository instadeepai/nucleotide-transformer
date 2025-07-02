import logging
from typing import Any

from pydantic import root_validator
from pydantic.main import BaseModel


class BulkRNABertConfig(BaseModel):
    n_genes: int
    n_expressions_bins: int
    embed_dim: int
    init_gene_embed_dim: int
    project_gene_embedding: bool = False
    use_gene_embedding: bool = True

    # architecture
    num_attention_heads: int
    key_size: int | None = None
    ffn_embed_dim: int
    num_layers: int
    use_memory_efficient_attention: bool = False

    use_gradient_checkpointing: bool = False

    gene2vec_weights_path: str

    embeddings_layers_to_save: tuple[int, ...] = ()
    attention_layers_to_save: tuple[int, ...] = ()

    # RNASeq data processing
    rnaseq_tokenizer_bins: list[float] | None = None
    use_log_normalization: bool = True
    use_max_normalization: bool = True
    normalization_factor: float | None = None

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
    def validate_gene_embedding(cls, values: dict[str, Any]) -> dict[str, Any]:
        """
        Checks that the given values are compatible.
        """
        use_gene_embedding = values.get("use_gene_embedding")
        if use_gene_embedding:
            init_gene_embed_dim = values["init_gene_embed_dim"]
            embed_dim = values["embed_dim"]
            if init_gene_embed_dim != embed_dim:
                project_gene_embedding = values["project_gene_embedding"]
                if not project_gene_embedding:
                    logging.warning(
                        f"Init gene embedding dimension ({init_gene_embed_dim})"
                        f"different than embedding dimension ({embed_dim})."
                        f"Setting `project_gene_embedding` to True"
                    )
                    values["project_gene_embedding"] = True
        return values
