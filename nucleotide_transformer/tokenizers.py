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

from itertools import product
from typing import Dict, List, Optional, Tuple

import numpy as np
import regex as re

from nucleotide_transformer.constants import EXTRA_NUCLEOTIDES, NUCLEOTIDES


def _compute_k_mers(k: int) -> List[str]:
    """
    Generates all the different k-mers for nucleotides given a value of k.

    Args:
        k: The k parameter for k-mers.

    Returns:
        All the different k-mers.
    """
    return ["".join(elt) for elt in product(NUCLEOTIDES, repeat=k)]


def compute_tokens_to_ids_v2(k_mers: int) -> Tuple[Dict[str, int], List[str]]:
    """Compute the tokens to ids mapping that correspond to the tokenizer used to train
    the v2 models.

    Args:
        k_mers (int): k used for the kmers.

    Returns:
        Dict[str, int]: Corresponding tokens to ids mapping.
    """
    # Get tokenizer
    kmers_tokens = _compute_k_mers(k=k_mers)
    standard_tokens = kmers_tokens + NUCLEOTIDES + EXTRA_NUCLEOTIDES

    unk_token = "<unk>"
    pad_token = "<pad>"
    mask_token = "<mask>"
    class_token = "<cls>"
    eos_token = "<eos>"
    bos_token = "<bos>"

    special_tokens = [
        unk_token,
        pad_token,
        mask_token,
        class_token,
        eos_token,
        bos_token,
    ]
    all_tokens = special_tokens + standard_tokens
    tokens_to_ids = {tok: i for i, tok in enumerate(all_tokens)}

    return tokens_to_ids, standard_tokens


class StandardTokenizer:
    """
    Simple tokenizer that extracts pre-defined tokens from sequences using regex.
    """

    def __init__(
        self,
        standard_tokens: List[str],
        unk_token: str = "<unk>",
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        class_token: str = "<cls>",
        eos_token: str = "<eos>",
        bos_token: str = "<bos>",
        prepend_bos_token: bool = False,
        prepend_cls_token: bool = False,
        append_eos_token: bool = False,
        extra_special_tokens: Optional[List[str]] = None,
        tokens_to_ids: Optional[Dict[str, int]] = None,
    ):
        """
        Initializes a basic tokenizer instance.

        Args:
            standard_tokens: Standard tokens, where special tokens are omitted.
            unk_token: Unknown token.
            pad_token: Pad token.
            mask_token: Mask token.
            class_token: Class token.
            eos_token: End of speech tokens.
            bos_token: Beginning of sentence token.
            prepend_bos_token: Prepend beginning of sentence token.
            prepend_cls_token: Prepend class token.
            append_eos_token: Append end of speech token.
            extra_special_tokens: (Optional) Enable the user to define optionally
                additional special tokens. Since regex is used for tokenization, any
                special tokens that are also special tokens in regex must include
                a "\" escape seq. For instance "$" -> "\\$"
            tokens_to_ids: (Optional) Enable the user to optionally choose ids for
                the tokens. If you provide this argument the dictionary must include
                the following special tokens
                ["<unk>","<pad>","<mask>","<cls>","<eos>","<bos>"]
                or instantiation will fail. Additionally, if the ids in your dictionary
                do not start at 0 then an error will also be raised. If this argument is
                not specified, then ids are attributed automatically by the tokenizer
                during initialization.
        """

        # Define special tokens essential to masked language modelling
        special_tokens_1 = [unk_token, pad_token, mask_token, class_token]
        special_tokens_2 = [eos_token, bos_token]
        special_tokens = special_tokens_1 + special_tokens_2
        all_tokens = special_tokens_1 + standard_tokens + special_tokens_2

        if extra_special_tokens is not None:
            special_tokens.extend(extra_special_tokens)

        self._all_tokens = all_tokens
        self._standard_tokens = standard_tokens
        self._special_tokens = special_tokens

        self._unk_token = unk_token
        self._pad_token = pad_token
        self._mask_token = mask_token
        self._class_token = class_token
        self._eos_token = eos_token
        self._bos_token = bos_token
        self._prepend_bos_token = prepend_bos_token
        self._prepend_cls_token = prepend_cls_token
        self._append_eos_token = append_eos_token

        # Can only
        if self._prepend_bos_token and self._prepend_cls_token:
            raise ValueError(
                "Cannot prepend both BOS and CLS token, must choose only one"
            )

        # Matching between tokens and ids
        if tokens_to_ids is not None:
            if set(tokens_to_ids.keys()) != set(self._all_tokens):
                raise ValueError(
                    f"Specified matching between tokens and ids, "
                    f"but some tokens are missing or mismatch. "
                    f"Got specifications for tokens: {set(tokens_to_ids.keys())} "
                    f"and expected for {set(self._all_tokens)}"
                )
            sorted_tokens = np.sort(list(tokens_to_ids.values()))
            if np.any(sorted_tokens != np.arange(len(self._all_tokens))):
                raise ValueError(
                    f"Specified matching between tokens and ids, "
                    f"but some ids are missing or mismatch. "
                    f"Got specifications for ids: {sorted_tokens} "
                    f"and expected for {np.arange(len(self._all_tokens))}"
                )
            self._tokens_to_ids = tokens_to_ids
        else:
            self._tokens_to_ids = {tok: i for i, tok in enumerate(self._all_tokens)}

        self._ids_to_tokens = {i: tok for tok, i in self._tokens_to_ids.items()}
        self._compiled_regex = re.compile("|".join(self._all_tokens + [r"\S"]))  # noqa

    @property
    def vocabulary(self) -> List[str]:
        return self._all_tokens

    @property
    def standard_tokens(self) -> List[str]:
        return self._standard_tokens

    @property
    def vocabulary_size(self) -> int:
        """
        Property that returns the total number of tokens.

        Returns:
            Total number of tokens.
        """
        return len(self.vocabulary)

    @property
    def unk_token_id(self) -> int:
        """
        Property that returns id (int representation) of the unknown token.

        Returns:
            Id (int representation) of the unknown token.
        """
        return self.token_to_id(self.unk_token)

    @property
    def pad_token_id(self) -> int:
        """
        Property that returns id (int representation) of the pad token.

        Returns:
            Id (int representation) of the pad token.
        """
        return self.token_to_id(self.pad_token)

    @property
    def mask_token_id(self) -> int:
        """
        Property that returns id (int representation) of the mask token.

        Returns:
            Id (int representation) of the mask token.
        """
        return self.token_to_id(self.mask_token)

    @property
    def class_token_id(self) -> int:
        """
        Property that returns id (int representation) of the class token.

        Returns:
            Id (int representation) of the class token.
        """
        return self.token_to_id(self.class_token)

    @property
    def eos_token_id(self) -> int:
        """
        Property that returns id (int representation) of the eos token.

        Returns:
            Id (int representation) of the eos token.
        """
        return self.token_to_id(self.eos_token)

    @property
    def bos_token_id(self) -> int:
        """
        Property that returns id (int representation) of the bos token.

        Returns:
            Id (int representation) of the bos token.
        """
        return self.token_to_id(self.bos_token)

    @property
    def special_tokens(self) -> List[str]:
        return self._special_tokens

    @property
    def unk_token(self) -> str:
        return self._unk_token

    @property
    def pad_token(self) -> str:
        return self._pad_token

    @property
    def mask_token(self) -> str:
        return self._mask_token

    @property
    def class_token(self) -> str:
        return self._class_token

    @property
    def eos_token(self) -> str:
        return self._eos_token

    @property
    def bos_token(self) -> str:
        return self._bos_token

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

    def tokenize(self, sequence: str) -> Tuple[List[str], List[int]]:
        """
        Tokenizes a sequence and returns the list of tokens as well
        as the list of their IDs. Any character found in the sequence that does not
        correspond to any token in the vocabulary is replaced by the unk token.

        Args:
            sequence: Sequence to be tokenized.

        Returns:
            List of tokens.
            List of token ids.
        """
        tokens: List[str] = self._compiled_regex.findall(sequence)
        tokens = [
            tok if tok in self._tokens_to_ids.keys() else self._unk_token
            for tok in tokens
        ]
        if self._prepend_cls_token:
            tokens = [self._class_token] + tokens

        if self._prepend_bos_token:
            tokens = [self._bos_token] + tokens

        if self._append_eos_token:
            tokens.append(self._eos_token)

        tokens_ids = [self.token_to_id(tok) for tok in tokens]

        return tokens, tokens_ids

    def pad_tokens_batch(
        self, batch: List[Tuple[List[str], List[int]]]
    ) -> List[Tuple[List[str], List[int]]]:
        """
        Takes a batch of sequences tokens ids and returns a batch of padded sequences.

        Args:
            batch: List of tuples, each composed of a sequence's tokens and token ids.

        Returns:
            List of 2-elements tuple for each sequence in the input where the tuple is
            containing 1. the list of the str representations of the
            tokens for that sequence and 2. the list of the int representations of
            the tokens for that sequence. Pad Tokens are added so that each sequence
            of tokens in the batch has the same length (all sequences padded to the
            length of the longest sequence in the batch).
        """
        lengths = [len(t[0]) for t in batch]
        maximum_length = max(lengths)
        deltas = [maximum_length - length for length in lengths]
        padded_tokens = [
            t[0] + ([self.pad_token] * delta) for t, delta in zip(batch, deltas)
        ]
        padded_tokens_ids = [
            t[1] + ([self.pad_token_id] * delta) for t, delta in zip(batch, deltas)
        ]
        return [
            (toks, toks_ids) for toks, toks_ids in zip(padded_tokens, padded_tokens_ids)
        ]

    def batch_tokenize(self, sequences: List[str]) -> List[Tuple[List[str], List[int]]]:
        """
        Tokenizes a batch of sequences.
        Sequences are padded to the maximum length in the batch.

        Args:
            sequences: Batch of sequences to be tokenized.

        Returns:
            Batch of tokenized sequences as well as their token ids,
            where every sequence has been padded to the maximum length
            in the batch.
        """
        return self.pad_tokens_batch(  # type: ignore
            [self.tokenize(seq) for seq in sequences]
        )


class NucleotidesKmersTokenizer(StandardTokenizer):
    """
    This is a tokenizer specific for nucleotide sequences.
    It only considers sequence containing the tokens A, T, C, G and N.
    N is always considered as a special token and tokenized alone.
    """

    def __init__(
        self,
        k_mers: int,
        unk_token: str = "<unk>",
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        class_token: str = "<cls>",
        eos_token: str = "<eos>",
        bos_token: str = "<bos>",
        prepend_bos_token: bool = False,
        prepend_cls_token: bool = False,
        append_eos_token: bool = False,
        tokens_to_ids: Optional[Dict[str, int]] = None,
    ):
        """
        Instantiates a FixedSizeNucleotideKmersTokenizer.

        Args:
            k_mers: How many nucleotides to consider for generating vocabulary.
            unk_token: Unknown token.
            pad_token: Pad token.
            mask_token: Mask token.
            class_token: Class token.
            eos_token: End of speech tokens.
            bos_token: Beginning of sentence token.
            prepend_bos_token: Prepend beginning of sentence token.
            prepend_cls_token: Prepend class token.
            append_eos_token: Append end of speech token.
            tokens_to_ids: (Optional) Enable the user to optionally choose ids for
                the tokens. If you provide this argument the dictionary must include
                the following special tokens
                ["<unk>","<pad>","<mask>","<cls>","<eos>","<bos>"]
                or instantiation will fail. Additionally, if the ids in your dictionary
                do not start at 0 then an error will also be raised. If this argument is
                not specified, then ids are attributed automatically by the tokenizer
                during initialization.
        """
        kmers_tokens = _compute_k_mers(k_mers)
        standard_tokens = kmers_tokens + NUCLEOTIDES + EXTRA_NUCLEOTIDES

        StandardTokenizer.__init__(
            self,
            standard_tokens=standard_tokens,
            unk_token=unk_token,
            pad_token=pad_token,
            mask_token=mask_token,
            class_token=class_token,
            eos_token=eos_token,
            bos_token=bos_token,
            prepend_bos_token=prepend_bos_token,
            prepend_cls_token=prepend_cls_token,
            append_eos_token=append_eos_token,
            tokens_to_ids=tokens_to_ids,
        )

        self._k_mers = k_mers

    def tokenize(self, sequence: str) -> Tuple[List[str], List[int]]:
        """
        Tokenizes a sequence and returns the list of tokens as well
        as the list of their IDs. The tokenization algorithm first splits up the
        substrings of the input sequence in-between N characters.
        Then these substrings are split into pieces of length k, and if it
        is possible (edge cases) it adds up pieces of length 1.

        If a single character that does not correspond
        to any token is found, an error is raised.

        Args:
            sequence: Sequence to be tokenized.

        Returns:
            List of tokens.
            List of token ids.

        Example:
            Find below two tokenization examples when k_mers=5.

            ATCGAATGGCGATGCAC --> ATCGA ATGGC GATGC A C

            ATCGAATNGGCGATGCAC -> ATCGA A T N GGCGA TGCAC
        """
        splitted_seq = sequence.split("N")
        len_splitted = len(splitted_seq)
        tokens: List[str] = []

        for i, split in enumerate(splitted_seq):
            chunks = [
                split[i * self._k_mers : (i + 1) * self._k_mers]
                for i in range(len(split) // self._k_mers)
            ]
            if len(split) % self._k_mers != 0:
                chunks.append(split[(len(split) // self._k_mers) * self._k_mers :])

            for chunk in chunks:
                if len(chunk) == self._k_mers:
                    tokens.append(chunk)
                else:
                    for nucl in chunk:
                        tokens.append(nucl)
            if i < len_splitted - 1:
                tokens.append("N")

        if self._prepend_cls_token:
            tokens = [self._class_token] + tokens

        if self._prepend_bos_token:
            tokens = [self._bos_token] + tokens

        if self._append_eos_token:
            tokens.append(self._eos_token)

        tokens_ids = [self.token_to_id(tok) for tok in tokens]

        return tokens, tokens_ids


class FixedSizeNucleotidesKmersTokenizer(NucleotidesKmersTokenizer):
    """
    Simple tokenizer that naively extracts tokens. Used for amino-acids
    and nucleotides. This tokenizer also tokenizes batches to a
    fixed maximum length. If one of the sequences provided exceeds the maximum
    length, an exception is raised.
    """

    def __init__(
        self,
        k_mers: int,
        fixed_length: int,
        unk_token: str = "<unk>",
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        class_token: str = "<cls>",
        eos_token: str = "<eos>",
        bos_token: str = "<bos>",
        prepend_bos_token: bool = False,
        prepend_cls_token: bool = False,
        append_eos_token: bool = False,
        tokens_to_ids: Optional[Dict[str, int]] = None,
    ):
        """
        Instantiates a FixedSizeNucleotideKmersTokenizer.

        Args:
            k_mers: How many nucleotides to consider for generating vocabulary.
            unk_token: Unknown token.
            pad_token: Pad token.
            mask_token: Mask token.
            class_token: Class token.
            eos_token: End of speech tokens.
            bos_token: Beginning of sentence token.
            prepend_bos_token: Prepend beginning of sentence token.
            prepend_cls_token: Prepend class token.
            append_eos_token: Append end of speech token.
            fixed_length: Fixed length to pad all sequences in batches.
        """
        NucleotidesKmersTokenizer.__init__(
            self,
            unk_token=unk_token,
            pad_token=pad_token,
            mask_token=mask_token,
            class_token=class_token,
            eos_token=eos_token,
            bos_token=bos_token,
            prepend_bos_token=prepend_bos_token,
            prepend_cls_token=prepend_cls_token,
            append_eos_token=append_eos_token,
            k_mers=k_mers,
            tokens_to_ids=tokens_to_ids,
        )
        self._fixed_length = fixed_length

    @property
    def fixed_length(self) -> int:
        """
        Property that returns the pre-defined fixed sequence length.

        Returns:
            The pre-defined fixed sequence length.
        """
        return self._fixed_length

    def pad_tokens_batch(
        self, batch: List[Tuple[List[str], List[int]]]
    ) -> List[Tuple[List[str], List[int]]]:
        """
        Takes tokens and tokens ids of a batch of sequences, and returns a batch of
        padded sequences.

        Args:
            batch: List of tuples, each composed of a sequence's tokens and token ids.

        Returns:
            The padded list, where every sequence is padded to the fixed maximum length.
        """
        lengths = [len(t[0]) for t in batch]
        maximum_length = max(lengths)
        if maximum_length > self._fixed_length:
            raise ValueError(
                f"Found a sequence with length {maximum_length} that "
                f"exceeds the fixed length to tokenize ({self._fixed_length})."
            )
        deltas = [self._fixed_length - length for length in lengths]
        padded_tokens = [
            t[0] + ([self.pad_token] * delta) for t, delta in zip(batch, deltas)
        ]
        padded_tokens_ids = [
            t[1] + ([self.pad_token_id] * delta) for t, delta in zip(batch, deltas)
        ]
        return [
            (toks, toks_ids) for toks, toks_ids in zip(padded_tokens, padded_tokens_ids)
        ]
