from abc import ABC, abstractmethod
from itertools import product
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import regex as re

NUCLEOTIDES = ["A", "T", "C", "G"]
VALID_EXTRA_NUCLEOTIDES = ["N", "M", "Y", "B", "S", "W", "K", "H", "D", "V", "R"]
EXTRA_NUCLEOTIDES = ["N"]


def _compute_k_mers(k: int) -> List[str]:
    """
    Generates all the different k-mers for nucleotides given a value of k.

    Args:
        k: The k parameter for k-mers.

    Returns:
        All the different k-mers.
    """
    return ["".join(elt) for elt in product(NUCLEOTIDES, repeat=k)]


class BaseTokenizer(ABC):
    """
    This class represents a possible general abstraction for tokenizers in Trix.
    All the tokenizers used in examples will inherit from this class. However, users may
    themselves define other kinds of tokenizers as they see fit. In this abstraction, we
    distinguish standard tokens (e.g. words in english, or amino-acids in proteomics)
    from special tokens (e.g. pad token, mask token, etc ...). This abstraction
    introduces a set of special tokens that we found to be useful as defaults. However,
    classes that inherit from this abstraction can add other special tokens or not use
    any of the ones defined by default.
    """

    @property
    @abstractmethod
    def vocabulary(self) -> List[str]:
        """
        Property that returns the list of all tokens (in str representation)
        used by that tokenizer.

        Returns:
            The list of all tokens (in str representation) used by that tokenizer.
        """
        pass

    @property
    def vocabulary_size(self) -> int:
        """
        Property that returns the total number of tokens.

        Returns:
            Total number of tokens.
        """
        return len(self.vocabulary)

    @property
    @abstractmethod
    def standard_tokens(self) -> List[str]:
        """
        Property that returns the list of standards tokens (in str representation)
        used by that tokenizer.

        Returns:
            The list of standards tokens (in str representation) used by that tokenizer.
        """
        pass

    @property
    @abstractmethod
    def special_tokens(self) -> List[str]:
        """
        Property that returns the list of special tokens (in str representation)
        used by that tokenizer.

        Returns:
            The list of special tokens (in str representation) used by that tokenizer.
        """
        pass

    @property
    @abstractmethod
    def unk_token(self) -> str:
        """
        Property that returns the str representation of the unknown token.

        Returns:
            Str representation of the unknown token.
        """
        pass

    @property
    @abstractmethod
    def pad_token(self) -> str:
        """
        Property that returns the str representation of the pad token.

        Returns:
            Str representation of the pad token.
        """
        pass

    @property
    @abstractmethod
    def mask_token(self) -> str:
        """
        Property that returns the str representation of the mask token.

        Returns:
            Str representation of the mask token.
        """
        pass

    @property
    @abstractmethod
    def class_token(self) -> str:
        """
        Property that returns the str representation of the class token. Note that
        class and bos (begin of sequence) tokens can sometimes be confused. We introduce
        both  tokens for more granularity. For instance, one might append a class token
        at the beginning of a sequence that goes into an encoder and a bos token
        at the beginning of a sequence that goes into a decoder.

        Returns:
            Str representation of the class token.
        """
        pass

    @property
    @abstractmethod
    def eos_token(self) -> str:
        """
        Property that returns the str representation of the end of sequence token.

        Returns:
            Str representation of the end of sequence token.
        """
        pass

    @property
    @abstractmethod
    def bos_token(self) -> str:
        """
        Property that returns the str representation of the bos (beginning of sequence)
        token. Note that class and bos tokens can be sometimes be confused. We introduce
        both  tokens for more granularity. For instance, one might append a class token
        at the beginning of a sequence that goes into an encoder and a bos token
        at the beginning of a sequence that goes into a decoder.

        Returns:
            Str representation of the bos token.
        """
        pass

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

    @abstractmethod
    def id_to_token(self, token_id: int) -> str:
        """
        Retrieves the str representation of a token from its id (int representation).

        Args:
            token_id: Id of the token.

        Returns:
            The str representation of the token.
        """
        pass

    @abstractmethod
    def token_to_id(self, token: str) -> int:
        """
        Retrieves the id (int representation) of a token from its str representation.

        Args:
            token: The str representation of the token.

        Returns:
            The id of the token.
        """

    @abstractmethod
    def tokenize(self, sequence: str) -> Tuple[List[str], List[int]]:
        """
        Tokenizes a sequence and returns the list of tokens as well
        as the list of their IDs. Any character found in the sequence that does not
        correspond to any token in the vocabulary is replaced by the unk token.

        Args:
            sequence: The sequence to be tokenized.

        Returns:
            A 2-elements tuple containing 1. the list of the str representations of the
            tokens and 2. the list of the int representations of the tokens.
        """
        pass

    @abstractmethod
    def batch_tokenize(self, sequences: List[str]) -> List[Tuple[List[str], List[int]]]:
        """
        Tokenizes a batch of sequences.

        Args:
            sequences: List of sequences to be tokenized.

        Returns:
            List of 2-elements tuple for each sequence in the input where the tuple is
            containing 1. the list of the str representations of the
            tokens for that sequence and 2. the list of the int representations of
            the tokens for that sequence.
        """
        pass

    @abstractmethod
    def add_tokens(self, new_tokens: Union[List[str], str]) -> str:
        """
        Adds new tokens to existing vocabulary.

        Args:
            sequences: List of sequences to be tokenized.

        Returns:
            An string representing the new tokens added.
        """
        pass

    @abstractmethod
    def add_special_tokens(self, new_tokens: Union[List[str], str]) -> str:
        """
        Adds new special tokens to existing vocabulary.

        Args:
            sequences: List of sequences to be tokenized.

        Returns:
            An integer representing the number of new special tokens added.
        """
        pass

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


class StandardTokenizer(BaseTokenizer):
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
        special_tokens = [
            unk_token,
            pad_token,
            mask_token,
            class_token,
            eos_token,
            bos_token,
        ]

        if extra_special_tokens is not None:
            special_tokens.extend(extra_special_tokens)

        self._all_tokens = special_tokens + standard_tokens
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

    def np_tokenize(self, sequence: str) -> np.ndarray:
        """
        Version of the tokenize version that returns a numpy array of the tokens IDs.
        """
        tokens_ids = self.tokenize(sequence)[1]
        return np.array(tokens_ids, dtype=int)

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

    def batch_np_tokenize(self, sequences: List[str]) -> np.ndarray:
        """
        Version of the batch_tokenize method that returns a numpy array of the tokens
        IDs.
        """
        outs = self.batch_tokenize(sequences)
        tokens_ids = [np.array(out[1]) for out in outs]
        return np.stack(tokens_ids, axis=0, dtype=int)

    def add_tokens(self, new_tokens: Union[List[str], str]) -> str:
        """
        Adds new tokens to existing vocabulary.

        Args:
            sequences: List of sequences to be tokenized.

        Returns:
            An string representing the new special tokens added.
        """
        if not new_tokens:
            return "No new tokens to add"

        if not isinstance(new_tokens, list):
            new_tokens = [new_tokens]
        # remove any not unique tokens
        new_tokens = list(
            set(new_tokens) - set(self.vocabulary)
        )  # TODO : use self._all_tokens instead?

        if len(new_tokens) == 0:
            return "No new tokens to add"
        # update _tokens_to_ids, _ids_to_tokens, _compiled_regex,
        added_tokens_encoder = {
            token: self.vocabulary_size + i for i, token in enumerate(new_tokens)
        }
        added_tokens_decoder = {
            id_: token for token, id_ in added_tokens_encoder.items()
        }

        self._standard_tokens.extend(new_tokens)
        self._all_tokens.extend(new_tokens)
        self._tokens_to_ids.update(added_tokens_encoder)
        self._ids_to_tokens.update(added_tokens_decoder)
        self._compiled_regex = re.compile("|".join(self._all_tokens + [r"\S"]))

        return f'{len(new_tokens)} new tokens added : {", ".join(new_tokens)}'

    def add_special_tokens(self, new_tokens: Union[List[str], str]) -> str:
        """
        Adds new special tokens to existing vocabulary.

        Args:
            sequences: List of sequences to be tokenized.

        Returns:
            An string representing the new special tokens added.
        """
        if not new_tokens:
            return "No new tokens to add"

        if not isinstance(new_tokens, list):
            new_tokens = [new_tokens]
        # remove any not unique tokens
        new_tokens = list(
            set(new_tokens) - set(self.vocabulary)
        )  # TODO : use self._all_tokens instead?

        if len(new_tokens) == 0:
            return "No new tokens to add"
        # update _tokens_to_ids, _ids_to_tokens, _compiled_regex,
        added_tokens_encoder = {
            token: self.vocabulary_size + i for i, token in enumerate(new_tokens)
        }
        added_tokens_decoder = {
            id_: token for token, id_ in added_tokens_encoder.items()
        }

        self._special_tokens.extend(new_tokens)
        self._all_tokens.extend(new_tokens)
        self._tokens_to_ids.update(added_tokens_encoder)
        self._ids_to_tokens.update(added_tokens_decoder)
        self._compiled_regex = re.compile("|".join(self._all_tokens + [r"\S"]))

        return f'{len(new_tokens)} new tokens added : {", ".join(new_tokens)}'


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
