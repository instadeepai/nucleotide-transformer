import inspect
import re
from abc import ABC, abstractmethod
from typing import Any, Type

import numpy as np


class BaseTokenizer(ABC):
    """
    This class represents a possible general abstraction for tokenizers.
    All the tokenizers used in examples will inherit from this class. However, users may
    themselves define other kinds of tokenizers as they see fit. In this abstraction, we
    distinguish standard tokens (e.g. words in english, or amino-acids in proteomics)
    from special tokens (e.g. pad token, mask token, etc ...). This abstraction
    introduces a set of special tokens that we found to be useful as defaults. However,
    classes that inherit from this abstraction can add other special tokens or not use
    any of the ones defined by default.
    """

    # Registry of all tokenizer classes
    _tokenizer_registry: dict[str, Type["BaseTokenizer"]] = {}

    @property
    @abstractmethod
    def vocabulary(self) -> list[str]:
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
    def standard_tokens(self) -> list[str]:
        """
        Property that returns the list of standards tokens (in str representation)
        used by that tokenizer.

        Returns:
            The list of standards tokens (in str representation) used by that tokenizer.
        """
        pass

    @property
    @abstractmethod
    def special_tokens(self) -> list[str]:
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

    @classmethod
    def get_tokenizer_class(cls, name: str) -> Type["BaseTokenizer"]:
        """
        Get a tokenizer class by name from the registry.

        Args:
            name: The name of the tokenizer class to retrieve.

        Returns:
            The tokenizer class.

        Raises:
            ValueError: If the tokenizer class is not found in the registry.
        """
        tokenizer_class = cls._tokenizer_registry.get(name)
        if tokenizer_class is None:
            raise ValueError(
                f"Unknown tokenizer class: '{name}'. "
                f"Available tokenizers: {list(cls._tokenizer_registry.keys())}"
            )
        return tokenizer_class

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "BaseTokenizer":
        """Create a tokenizer instance from a configuration dictionary."""
        class_name = config.get("class_name")
        params = config.get("params")

        if not class_name or not params:
            raise ValueError(
                "Configuration must contain 'class_name' and 'params' keys."
            )

        # Get the tokenizer class from the registry
        tokenizer_class = cls.get_tokenizer_class(class_name)

        # Use dictionary unpacking to pass the parameters to the constructor
        return tokenizer_class(**params)

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
    def tokenize(self, sequence: str) -> tuple[list[str], list[int]]:
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
    def batch_tokenize(self, sequences: list[str]) -> list[tuple[list[str], list[int]]]:
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
    def add_tokens(self, new_tokens: list[str] | str) -> str:
        """
        Adds new tokens to existing vocabulary.

        Args:
            new_tokens: List of new tokens to add to vocabulary.

        Returns:
            An string representing the new tokens added.
        """
        pass

    @abstractmethod
    def add_special_tokens(self, new_tokens: list[str] | str) -> str:
        """
        Adds new special tokens to existing vocabulary.

        Args:
            new_tokens: List of new special tokens to add to vocabulary.
        Returns:
            An integer representing the number of new special tokens added.
        """
        pass

    def to_config(self) -> dict[str, Any]:
        """Convert a tokenizer instance to a configuration dictionary."""
        # Base parameters common to most tokenizers
        params = {
            "unk_token": self.unk_token,
            "pad_token": self.pad_token,
            "mask_token": self.mask_token,
            "class_token": self.class_token,
            "eos_token": self.eos_token,
            "bos_token": self.bos_token,
        }

        # Get constructor parameters
        sig = inspect.signature(self.__class__.__init__)
        for param_name in sig.parameters:
            if param_name != "self":
                # Try both public and private versions of the attribute
                value = None
                if hasattr(self, param_name):
                    value = getattr(self, param_name)
                elif hasattr(self, f"_{param_name}"):
                    value = getattr(self, f"_{param_name}")

                if value is not None:
                    params[param_name] = value

        return {
            "class_name": self.__class__.__name__,
            "params": params,
        }

    def pad_tokens_batch(
        self, batch: list[tuple[list[str], list[int]]]
    ) -> list[tuple[list[str], list[int]]]:
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

    def __init_subclass__(cls: Type["BaseTokenizer"], **kwargs: dict) -> None:
        """Automatically register tokenizer subclasses in the registry."""
        super().__init_subclass__(**kwargs)
        cls._tokenizer_registry[cls.__name__] = cls

class StandardTokenizer(BaseTokenizer):
    """
    Simple tokenizer that extracts pre-defined tokens from sequences using regex.
    """

    def __init__(
        self,
        standard_tokens: list[str],
        unk_token: str = "<unk>",
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        class_token: str = "<cls>",
        eos_token: str = "<eos>",
        bos_token: str = "<bos>",
        prepend_bos_token: bool = False,
        prepend_cls_token: bool = False,
        append_eos_token: bool = False,
        extra_special_tokens: list[str] | None = None,
        tokens_to_ids: dict[str, int] | None = None,
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
        self._extra_special_tokens = extra_special_tokens

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
        self._compiled_regex = re.compile(
            "|".join(rf"\b{re.escape(tok)}\b" for tok in self._all_tokens) + r"|\S"
        )

    @property
    def vocabulary(self) -> list[str]:
        return self._all_tokens

    @property
    def standard_tokens(self) -> list[str]:
        return self._standard_tokens

    @property
    def special_tokens(self) -> list[str]:
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

    def tokenize(self, sequence: str) -> tuple[list[str], list[int]]:
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
        tokens: list[str] = self._compiled_regex.findall(sequence)
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

    def batch_tokenize(self, sequences: list[str]) -> list[tuple[list[str], list[int]]]:
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

    def batch_np_tokenize(self, sequences: list[str]) -> np.ndarray:
        """
        Version of the batch_tokenize method that returns a numpy array of the tokens
        IDs.
        """
        outs = self.batch_tokenize(sequences)
        tokens_ids = [np.array(out[1]) for out in outs]
        return np.stack(tokens_ids, axis=0, dtype=int)

    def add_tokens(self, new_tokens: list[str] | str) -> str:
        """
        Adds new tokens to existing vocabulary.

        Args:
            new_tokens: List of sequences to be tokenized.

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
        self._compiled_regex = re.compile(
            "|".join(rf"\b{re.escape(tok)}\b" for tok in self._all_tokens) + r"|\S"
        )

        return f'{len(new_tokens)} new tokens added : {", ".join(new_tokens)}'

    def add_special_tokens(self, new_tokens: list[str] | str) -> str:
        """
        Adds new special tokens to existing vocabulary.

        Args:
            new_tokens: List of sequences to be tokenized.

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
        self._compiled_regex = re.compile(
            "|".join(rf"\b{re.escape(tok)}\b" for tok in self._all_tokens) + r"|\S"
        )

        return f'{len(new_tokens)} new tokens added : {", ".join(new_tokens)}'

class NucleotideTokenizer(StandardTokenizer):
    """
    Fast version of the 1 k-mer tokenizer that is used for the long_range_nt models.
    """

    def __init__(
        self,
        unk_token: str = "<unk>",
        pad_token: str = "<pad>",
        mask_token: str = "<mask>",
        class_token: str = "<cls>",
        eos_token: str = "<eos>",
        bos_token: str = "<bos>",
        prepend_bos_token: bool = False,
        prepend_cls_token: bool = False,
        append_eos_token: bool = False,
        extra_special_tokens: (
            list[str] | None
        ) = None,  # Additional special tokens for consistent serialization
        tokens_to_ids: dict[str, int] | None = None,
    ):
        StandardTokenizer.__init__(
            self,
            standard_tokens=["A", "T", "C", "G", "N"],
            unk_token=unk_token,
            pad_token=pad_token,
            mask_token=mask_token,
            class_token=class_token,
            eos_token=eos_token,
            bos_token=bos_token,
            prepend_bos_token=prepend_bos_token,
            prepend_cls_token=prepend_cls_token,
            append_eos_token=append_eos_token,
            extra_special_tokens=extra_special_tokens,
            tokens_to_ids=tokens_to_ids,
        )

    def np_tokenize(self, sequence: str, validate_input: bool = False) -> np.array:
        """
        Fast tokenization of DNA sequences. NOTE: It is intended to be used with
        sequences that are only made of standard tokens of length 1.

        Args:
            sequence: Biological sequence, should only include tokens that are within
                standard_tokens.
            validate_input: If True, it will check that the characters in the sequence
                are all within standard tokens. This adds an overhead and makes the
                tokenization slower so it defaults as False.
        """
        if validate_input:
            # check that all characters in the sequence are within standard tokens.
            if not set(sequence).issubset(self.standard_tokens):
                raise ValueError(
                    f"sequence must only include tokens among {self.standard_tokens}"
                )

        sequence = ", ".join(list(sequence))

        for token in self.standard_tokens:
            token_id = self._tokens_to_ids[token]
            sequence = sequence.replace(token, str(token_id))

        if self._prepend_cls_token:
            cls_tok_id = self._tokens_to_ids[self._class_token]
            sequence = f"{cls_tok_id}, " + sequence

        if self._prepend_bos_token:
            bos_tok_id = self._tokens_to_ids[self._bos_token]
            sequence = f"{bos_tok_id}, " + sequence

        if self._append_eos_token:
            eos_tok_id = self._tokens_to_ids[self._eos_token]
            sequence = sequence + f", {eos_tok_id}"

        tokens_ids = np.fromstring(sequence, dtype=int, sep=",")
        return tokens_ids

    def np_untokenize(self, tokens: np.ndarray) -> str:
        """
        Fast untokenization of a single DNA sequence.
        NOTE: It is approximately 5 times faster than the basic untokenization
        ([tokenizer.id_to_token(id) for id in tokens])

        Args:
            tokens (np.ndarray): Tokenized sequence with shape (num_tokens,)

        Returns:
            str: Decoded sequence

        Example:
            >>> tokens = np.array([6, 7, 8, 8, 7, 7, 1])
            >>> tokenizer = NucleotideTokenizer()
            >>> tokenizer.np_untokenize(tokens)
            'ATCCTT<pad>'
        """
        # Define ids_to_token in a np.ndarray way because it allows faster mapping of
        # ID to token for unokenizing
        ids_to_tokens = self._ids_to_tokens
        token_map = np.array(
            [ids_to_tokens[key] for key in sorted(ids_to_tokens.keys())]
        )

        # Efficiently map token IDs to their string equivalents
        string_tokens = token_map[tokens]

        return "".join(string_tokens)

    def batch_np_tokenize(self, sequences: list[str]) -> np.array:
        """
        Version of the batch_tokenize method that returns a numpy array of the tokens
        IDs.

        Args:
            sequences (list[str]): List of DNA sequences.

        Returns:
            np.array: Tokens with shape (num_sequences, max_tokenized_sequence_length)

        Example:
            >>> from nucleotide_transformer_v3.tokenizers import NucleotideTokenizer
            >>> tokenizer = NucleotideTokenizer()
            >>> sequences = ["ACGT", "AC"]
            >>> tokenizer.batch_np_tokenize(sequences)
            array([[6, 8, 9, 7],
                   [6, 8, 1, 1]], dtype=int32)

        Note for Developers:
            1) This overwrites the method `batch_np_tokenize` from StandardTokenizer
            because the current implementation is faster for a specific use case.
            The `for` loop on pre-defined array is faster than `np.stack` for big
            batches and long sequences

            2) This method does not use the method `pad_tokens_batch` because it expects
            list[tuple[list[str], list[int]]] as input. In the current method, since the
            method `np_tokenize` is used (instead of `tokenize`), we only have the
            tokens_ids (list[int]) and not the tokens themselves (list[str]), and thus
            cannot use `pad_tokens_batch` directly.
        """
        tokenized_sequences = [self.np_tokenize(sequence) for sequence in sequences]
        batch = np.full(
            shape=(len(sequences), max(len(tokens) for tokens in tokenized_sequences)),
            fill_value=self.pad_token_id,
            dtype=np.int32,
        )
        for i, tokenized_sequence in enumerate(tokenized_sequences):
            batch[i, : len(tokenized_sequence)] = tokenized_sequence

        return batch

    def batch_np_untokenize(self, tokens: np.ndarray) -> list[str]:
        """
        Fast untokenization of a batch of DNA sequences.

        Args:
            tokens (np.ndarray): Tokenized sequence with shape
                (num_sequences, num_tokens)

        Returns:
            list[str]: List of decoded sequences

        Example:
            >>> tokens = np.array(
                    [
                        [6, 6, 7, 8, 9, 7],
                        [8, 7, 9, 8, 7, 7]
                    ]
                )
            >>> tokenizer = NucleotideTokenizer()
            >>> tokenizer.batch_np_untokenize(tokens)
            ['AATCGT', 'CTGCTT']

        Note for Developers:
            This method could also be implemented by calling the method np_untokenize
            on each row (ex: [self.np_untokenize(seq) for seq in tokens]).
            However, the current implementation yields slightly faster execution,
            because it prevents from redefining the variable `token_map` at each call
            of the method np_untokenize
        """
        # Define ids_to_token in a np.ndarray way because it allows faster mapping of
        # ID to token for unokenizing
        ids_to_tokens = self._ids_to_tokens
        token_map = np.array(
            [ids_to_tokens[key] for key in sorted(ids_to_tokens.keys())]
        )

        # Efficiently map token IDs to their string equivalents
        string_tokens = token_map[tokens]

        return ["".join(row) for row in string_tokens]


def get_ntv3_tokenizer() -> NucleotideTokenizer:
    """
    Returns tokenizer for DNA sequences when using the NTv3 models.
    """
    tokenizer = NucleotideTokenizer(
        unk_token="<unk>",
        pad_token="<pad>",
        mask_token="<mask>",
        prepend_bos_token=False,
        prepend_cls_token=False,
        append_eos_token=False,
    )
    return tokenizer
