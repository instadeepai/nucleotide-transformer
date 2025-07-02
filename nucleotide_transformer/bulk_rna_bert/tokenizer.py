from typing import List

import numpy as np


class BinnedOmicTokenizer:
    """
    Tokenizer that bins gene expressions or methylation to convert them to tokens.
    """

    def __init__(
        self,
        n_expressions_bins: int,
        min_omic_value: float = 0.0,
        max_omic_value: float = 1.0,
        use_max_normalization: bool = True,
        normalization_factor: float = 1.0,
        prepend_cls_token: bool = False,
        fixed_sequence_length: int | None = None,
        unpadded_length: int | None = None,
    ):
        self._n_expressions_bins = n_expressions_bins
        self._use_max_normalization = use_max_normalization
        self._normalization_factor = normalization_factor
        self._prepend_cls_token = prepend_cls_token

        self._gene_expression_bins = np.linspace(
            min_omic_value, max_omic_value, self._n_expressions_bins
        )

        self._fixed_sequence_length = fixed_sequence_length
        self._unpadded_length = unpadded_length

        standard_tokens = list(map(str, range(len(self._gene_expression_bins))))
        self._pad_token = "<pad>"
        self._mask_token = "<mask>"
        self._class_token = "<cls>"
        self._unk_token = "<unk>"
        self._eos_token = "<eos>"
        self._bos_token = "<bos>"
        self._missing_modality_token = "<mmo>"

        if prepend_cls_token:
            special_tokens = [
                self._class_token,
                self._pad_token,
                self._mask_token,
            ]
        else:
            special_tokens = [
                self._pad_token,
                self._mask_token,
            ]

        self._all_tokens = standard_tokens + special_tokens
        self._standard_tokens = standard_tokens
        self._special_tokens = special_tokens

        self._tokens_to_ids = {tok: i for i, tok in enumerate(self._all_tokens)}
        self._ids_to_tokens = {i: tok for tok, i in self._tokens_to_ids.items()}

    @property
    def gene_expression_bins(self) -> np.ndarray:
        return self._gene_expression_bins

    @property
    def pad_token(self) -> str:
        return self._pad_token

    @property
    def mask_token(self) -> str:
        return self._mask_token

    @property
    def mask_token_id(self) -> int:
        """
        Property that returns id (int representation) of the mask token.

        Returns:
            Id (int representation) of the mask token.
        """
        return self.token_to_id(self.mask_token)

    @property
    def class_token(self) -> str:
        return self._class_token

    @property
    def unk_token(self) -> str:
        return self.unk_token

    @property
    def eos_token(self) -> str:
        return self._eos_token

    @property
    def bos_token(self) -> str:
        return self._bos_token

    @property
    def pad_id(self) -> int:
        return self.token_to_id(self.pad_token)

    @property
    def mask_id(self) -> int:
        return self.token_to_id(self.mask_token)

    @property
    def class_id(self) -> int:
        return self.token_to_id(self.class_token)

    @property
    def vocabulary(self) -> List[str]:
        return self._all_tokens

    @property
    def standard_tokens(self) -> List[str]:
        return self._standard_tokens

    @property
    def special_tokens(self) -> List[str]:
        return self._special_tokens

    def id_to_token(self, token_id: int) -> str:
        try:
            return self._ids_to_tokens[token_id]
        except KeyError:
            raise KeyError(f"Token id {token_id} not found in vocabulary")

    def token_to_id(self, token: str) -> int:
        try:
            return self._tokens_to_ids[token]
        except KeyError:
            raise KeyError(f"Token {token} not found in vocabulary")

    def tokenize(
        self,
        gene_expressions: np.ndarray | None,
        pad_to_fixed_length: bool = False,
        max_length: int | None = None,
    ) -> np.ndarray:
        """
        Tokenize a gene expression array and return an array of bin ids.

        Args:
            gene_expressions: Gene expressions sequence to be tokenized.
            pad_to_fixed_length: if True and fixed length is provided as attributed
            to the tokenizer, the sequence will be padded.
            max_length: allows to pass another max length than the one specified
            by self._fixed_sequence_length
        Returns:
            List of tokens ids.
        """
        if gene_expressions is None:
            assert self._unpadded_length is not None
            tokens_ids = np.array([self.mask_token_id] * self._unpadded_length)
        else:
            if self._use_max_normalization:
                gene_expressions /= self._normalization_factor
            tokens_ids = np.digitize(gene_expressions, self._gene_expression_bins)
            tokens_ids[gene_expressions == 0.0] = 0
        if self._prepend_cls_token:
            tokens_ids = np.concatenate([[self.class_id], tokens_ids])
        if pad_to_fixed_length:
            if self._fixed_sequence_length is not None:
                current_max_length = self._fixed_sequence_length
            else:
                assert max_length is not None
                current_max_length = max_length
            padded_tokens_ids = np.ones(
                current_max_length, dtype=tokens_ids.dtype
            ) * self.token_to_id(self._pad_token)
            padded_tokens_ids[: len(tokens_ids)] = tokens_ids
            return padded_tokens_ids
        return tokens_ids

    def batch_tokenize(
        self,
        gene_expressions: np.ndarray,
        pad_to_fixed_length: bool = False,
        max_length: int | None = None,
    ) -> np.ndarray:
        """
        Tokenizes a batch of gene expressions.

        Args:
            gene_expressions: gene expressions sequence to be tokenized.
            pad_to_fixed_length: if True and fixed length is provided as attributed
            to the tokenizer, the sequence will be padded.
            max_length: max length in the batch

        Returns:
            Tokenized gene expressions.
        """
        return np.vstack(
            [
                self.tokenize(
                    g, pad_to_fixed_length=pad_to_fixed_length, max_length=max_length
                )
                for g in gene_expressions
            ]
        )
