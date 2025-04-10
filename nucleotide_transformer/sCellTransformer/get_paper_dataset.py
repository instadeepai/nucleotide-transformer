import dataclasses
import gzip
import itertools
import json
import math
import random
from functools import partial
from typing import Any, Callable, Generator, Literal, TextIO

import numpy as np
import torch
from huggingface_hub import hf_hub_download
from torch.utils.data import IterableDataset

file_path = hf_hub_download(
    repo_id="InstaDeepAI/sCellTransformer", filename="data/protein_gene_map.json"
)
with open(file_path, "r") as f:
    PROTEIN_GENE_MAP = json.load(f)

TISSUES = [
    "blood",
    "brain",
    "heart",
    "intestine",
    "kidney",
    "lung",
    "pancreas",
    "others",
]


@dataclasses.dataclass
class DatasetConfig:
    # parameterizes the huggingface dataset to fetch CELLxGENES samples
    # Parameters:
    #     min_gene_expression: used to filter out genes with low expression.
    #     gene_expression_num_bins: number of bins to use for digitizing gene
    #         expression values.
    #     apply_log_trick: whether to apply the log trick to the gene
    #         expression values.
    min_gene_expression: float = 0.0
    gene_expression_num_bins: int = 51
    apply_log_trick: bool = False


class ScRNABagofCellsDataset(IterableDataset):
    """
    Dataset class for single cell data that returns samples
    containing a fixed number of cells.
    """

    def __init__(
        self,
        num_hosts: int,
        worker_id: int,
        epoch_num: int,
        num_cells_per_sample: int,
        num_files_per_chunk: int,
        fixed_sequence_length: int,
        expression_pad_token_id: int,
        id_pad_token_id: int,
        split: str = "train",
        preprocess_fn: Callable | None = None,
        resume_training: bool = False,
        seed: int = 0,
        add_all_metadata: bool = False,
        merge_fn: Callable | None = None,
    ):
        """
        Iterable dataset using pytorch utilities for fast loading
        of cell x gene data. Data currently exists on GCS.
        Args:
            num_hosts: Number of hosts used for distributed training
            worker_id: Current worker id
            epoch_num: Current epoch
            num_cells_per_sample: Number of cells to return per sample
            num_files_per_chunk: Number of files to download at a time.
                Helps handle the limited persistent memory issue in VMs.
            fixed_sequence_length: maximum sequence length.
            expression_pad_token_id: Padding token id for gene expressions
            id_pad_token_id: Padding token id for gene ids
            split: train or test split
            preprocess_fn: Optional preprocessing function.
            resume_training: Check if we are resuming training or not.
                If True, the file list is shuffled so that the training
                does not see the same order of samples.
            seed: Random seed used for shuffling.
            add_all_metadata: Include all metadata to records when
                de-serializing form JSON
            merge_fn: Callable to merge several sequences into
                a batch or single sequence. The default mere function (
                used for CellNT) merges multiple cells into a single
                sequence.
        """
        self.num_hosts = num_hosts
        self.worker_id = worker_id
        self.epoch_num = epoch_num
        self.split = split
        self.preprocess_fn = preprocess_fn
        assert (
            self.worker_id < self.num_hosts
        ), f"{self.worker_id} should be less than {self.num_hosts}"
        self.resume_training = resume_training
        self.seed = seed
        self.num_cells_per_sample = num_cells_per_sample
        self.num_files_per_chunk = num_files_per_chunk
        self.expression_pad_token_id = expression_pad_token_id
        self.id_pad_token_id = id_pad_token_id
        self.fixed_sequence_length = fixed_sequence_length
        self.add_all_metadata = add_all_metadata
        self.merge_fn = merge_fn

        # Set seed for shuffling
        random.seed(self.seed)
        dataset_file_path = hf_hub_download(
            repo_id="InstaDeepAI/sCellTransformer", filename="data/dataset.txt.gz"
        )
        self.file_paths = [dataset_file_path]

    def get_study_file_iterators(self, file_paths) -> tuple[list[TextIO], list[str]]:
        """
        Returns a list of file iterators for the study partitions and a list
        of failed downloads (here always empty).
        """

        res = []
        for file_path in file_paths:
            fp = gzip.open(file_path, "rt", encoding="utf-8")
            res.append(itertools.cycle(fp))

        return res, []

    def _pad_sequence(self, sample: dict[str, Any]) -> dict[str, Any]:
        """
        Pads gene expressions and ids to get a sequence length divisible
        by 2.
        """
        gene_ids = sample["gene_ids"]
        gene_expressions = sample["gene_expressions"]
        seq_len = gene_expressions.shape[0]
        padded_gene_expressions = (
            np.ones((self.fixed_sequence_length), dtype=gene_expressions.dtype)
            * self.expression_pad_token_id
        )
        padded_gene_ids = (
            np.ones((self.fixed_sequence_length), dtype=gene_ids.dtype)
            * self.id_pad_token_id
        )
        if "raw_gene_expressions" in sample:
            padded_raw_gene_expressions = (
                np.ones((self.fixed_sequence_length), dtype=np.float32)
                * self.expression_pad_token_id
            )
            padded_raw_gene_expressions[:seq_len] = sample["raw_gene_expressions"]
            sample["raw_gene_expressions"] = padded_raw_gene_expressions

        padded_gene_expressions[:seq_len] = gene_expressions
        padded_gene_ids[:seq_len] = gene_ids
        sample["gene_ids"] = padded_gene_ids
        sample["gene_expressions"] = padded_gene_expressions
        return sample

    def merge_records_into_sample(
        self, record_list: list[dict[str, Any]]
    ) -> dict[str, Any]:
        """
        Merges a list of cells to get a sequence of gene expressions
        and ids for all cells.
        Args:
            record_list: List of dictionaries
        Returns:
            A dictionary with all gene ids, gene expressions,
            and other metadata for a set of cells.
        """
        output_record = {}
        keys = list(record_list[0].keys())
        for key in keys:
            output_record[key] = [i[key] for i in record_list]
        for key in keys:
            if isinstance(output_record[key][0], np.ndarray):
                output_record[key] = np.concatenate(output_record[key])
            else:
                output_record[key] = list(itertools.chain(output_record[key]))
        return output_record

    def __iter__(self) -> Generator:
        file_paths_subset = self.file_paths
        random.shuffle(file_paths_subset)

        current_file_chunk = file_paths_subset[: self.num_files_per_chunk]
        file_paths_subset = file_paths_subset[self.num_files_per_chunk :]
        file_iters, failed_downloads = self.get_study_file_iterators(current_file_chunk)
        file_paths_subset = file_paths_subset + failed_downloads
        while True:
            # Stop if all files have been read
            if len(file_iters) == 0:
                break
            file_iter_to_use = random.sample(file_iters, 1)[0]
            sample_records = []
            idx = 0
            while idx < self.num_cells_per_sample:
                try:
                    record_text = next(file_iter_to_use)
                except StopIteration:
                    file_iter_to_use.close()
                    file_iters.remove(file_iter_to_use)  # type: ignore
                    if len(file_paths_subset) == 0:
                        break
                    new_file_path = file_paths_subset.pop()
                    new_file_ptrs, failed_file = self.get_study_file_iterators(
                        [new_file_path]
                    )
                    if len(new_file_ptrs) != 0:
                        new_file_ptr = new_file_ptrs[0]
                    else:
                        file_paths_subset = file_paths_subset + failed_file
                        break
                    file_iters.append(new_file_ptr)
                    # We drop the last batch from every study if it is not
                    # exactly equal to num cells per sequence.
                    break
                try:
                    record = json.loads(record_text)
                except json.JSONDecodeError:
                    continue
                dataset_id: str = record["dataset_id"]  # type: ignore
                tissue_type = record.get("tissue_general", "others")  # type: ignore

                record_dict = {
                    "gene_expressions": record[
                        "gene_expressions"  # type:ignore
                    ],
                    "gene_ids": record["gene_ids"],  # type: ignore
                    "tissue": [tissue_type],
                    "dataset_id": [dataset_id],
                    "source": ["scrna"],
                }
                # Keep other data if present, let preprocess_fn deal with it
                if self.add_all_metadata:
                    record_dict.update(
                        {
                            key: record[key]
                            for key in record.keys()
                            if key
                            not in [
                                "gene_ids",
                                "gene_expressions",
                                "tissue",
                                "dataset_id",
                                "source",
                            ]
                        }
                    )

                if self.preprocess_fn:
                    record_dict = self.preprocess_fn(record_dict)
                if record_dict is None:
                    continue
                idx += 1
                sample_records.append(record_dict)

            # Merge all cells into a single sequence
            if len(sample_records) == self.num_cells_per_sample:
                if not self.merge_fn:
                    sample = self.merge_records_into_sample(sample_records)
                else:
                    sample = self.merge_fn(sample_records)
                # sample = self._pad_sequence(sample)

            else:
                continue

            yield sample

    def __len__(self) -> int:
        # Calculated from the training dataset
        if self.split == "train":
            return 4500000 // self.num_cells_per_sample
        else:
            return 1000000 // self.num_cells_per_sample


def get_protein_gene_map() -> tuple[dict[int, int], np.ndarray, np.ndarray]:
    """
    All protein coding genes are used.
    Returns:
        gene_map: a map from the gene id in data to
            a position in the full sequence.
        gene_ids: a list of gene ids in the same order
        as the full sequence.
        gene_map_array: a numpy array of size 70000 with
            the gene ids in the same order as the full sequence.
    """
    # Fix type casting error
    gene_map = {int(k): PROTEIN_GENE_MAP[k] for k in PROTEIN_GENE_MAP}

    # Compute inverse gene_map and gene_ids
    inverse_gene_map = {v: k for k, v in gene_map.items()}
    gene_ids = np.array([inverse_gene_map[i] for i in range(len(gene_map.keys()))])

    # Convert gene_map to an indexable numpy array
    gene_map_array = np.zeros(70000, dtype=np.int32) - 1
    for k, v in gene_map.items():
        gene_map_array[k] = v

    return gene_map, gene_ids, gene_map_array


def binning_non_zero_gene_expressions(
    non_zero_gene_expressions: np.ndarray, binning_style: str, num_bins: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Combines bin edge calculation and applies digitizing by bins
    N.B convention is to apply cell by cell
    """
    if binning_style == "uniform":
        bins = np.linspace(
            min(non_zero_gene_expressions), max(non_zero_gene_expressions), num_bins
        )
        bins[-1] += 0.01
    elif binning_style == "quantile":
        bins = np.quantile(non_zero_gene_expressions, np.linspace(0, 1, num_bins))
        bins[-1] += 0.01
    else:
        raise ValueError(
            "Possible values for binning_style are 'uniform' or 'quantile'"
        )
    binned_gene_expressions = np.digitize(non_zero_gene_expressions, bins)
    assert binned_gene_expressions.min() >= 1
    assert binned_gene_expressions.max() <= num_bins - 1
    return binned_gene_expressions, bins


def construct_full_gene_sequence(
    example: dict[str, Any],
    gene_map_array: np.ndarray,
    gene_ids: np.ndarray,
    expression_pad_token_id: int,
    id_pad_token_id: int,
    num_bins: int = 51,
    apply_log_trick: bool = False,
    fixed_sequence_length: int | None = None,
    binning_style: Literal["quantile", "uniform"] = "uniform",
    include_raw_expressions: bool = False,
    include_bins: bool = False,
    normalize_by_sum: bool = False,
    use_spatial_dataset: bool = False,
) -> dict[str, Any] | None:
    """
    Process the current non-zero gene sequence
    to construct a full gene sequence with all
    possible genes from a gene_map.
    Args:
        gene_ids:
        example:
        gene_map_array:
        examples: dictionary containing non-zero
            gene_ids and gene expressions along
            with metadata.
        gene_map: a map from the gene id in data to
            a position in the full sequence.
        num_bins: number of bins to discretize the
            gene expressions.
        apply_log_trick: Optional. Applying logarithm on the raw expression of genes.
        expression_pad_token_id: Pad token used for expression sequences.
        id_pad_token_id: Pad token used for id sequences.
        fixed_sequence_length: Optional parameter to
            force the sequence to be a certain length.
        binning_style: Choice of binning algorithm.
            Currently supports quantile or uniform.
        include_raw_expressions: Whether to include
            unbinned expressions.
        include_bins: Whether to include bin edges.
        normalize_by_sum: Normalize non-zero gene expressions
            by total number of counts in a cell.
        use_spatial_dataset: Whether to include spatial info.
    Returns:
        processed sample with full gene id sequence
        and all gene expressions including 0.
    """
    other_keys = [
        key
        for key in example.keys()
        if key not in ["gene_ids", "gene_expressions", "spatial_x", "spatial_y"]
    ]

    expressed_gene_indexes = np.asarray(example["gene_ids"], dtype=np.int32)
    non_zero_gene_expressions = np.asarray(
        example["gene_expressions"], dtype=np.float32
    )

    if apply_log_trick:
        non_zero_gene_expressions = np.log1p(non_zero_gene_expressions)

    # Map the gene indexes to the full gene sequence
    # -1 indicates that the gene is not in the full sequence
    mapped_gene_indexes = gene_map_array[expressed_gene_indexes]
    new_position_indexes = mapped_gene_indexes[mapped_gene_indexes != -1]
    non_zero_gene_expressions = non_zero_gene_expressions[mapped_gene_indexes != -1]
    if normalize_by_sum:
        non_zero_gene_expressions /= non_zero_gene_expressions.sum()

    # Expand this to a full gene id and expression sequence now.
    gene_expressions = np.zeros(len(gene_ids), dtype=np.uint16)
    if len(non_zero_gene_expressions) == 0:
        return None
    binned_gene_expressions, bins = binning_non_zero_gene_expressions(
        non_zero_gene_expressions, binning_style, num_bins
    )
    gene_expressions[new_position_indexes] = binned_gene_expressions

    if include_raw_expressions:
        raw_gene_expressions = np.zeros(len(gene_ids))
        raw_gene_expressions[new_position_indexes] = non_zero_gene_expressions

    if fixed_sequence_length and fixed_sequence_length > len(gene_expressions):
        padded_gene_expressions = (
            np.ones(fixed_sequence_length) * expression_pad_token_id
        )
        padded_gene_ids = np.ones(fixed_sequence_length) * id_pad_token_id

        padded_gene_expressions[: len(gene_expressions)] = gene_expressions
        padded_gene_ids[: len(gene_ids)] = gene_ids
        if include_raw_expressions:
            padded_raw_gene_expressions = (
                np.ones(fixed_sequence_length) * expression_pad_token_id
            )
            padded_raw_gene_expressions[: len(raw_gene_expressions)] = (
                raw_gene_expressions
            )

    else:
        padded_gene_expressions = gene_expressions
        padded_gene_ids = gene_ids
        if include_raw_expressions:
            padded_raw_gene_expressions = raw_gene_expressions

    sample = {
        "gene_ids": np.asarray(padded_gene_ids, dtype=np.int32),
        "gene_expressions": np.asarray(padded_gene_expressions, dtype=np.int32),
    }
    if use_spatial_dataset:
        sample["spatial_x"] = np.asarray(example["spatial_x"], dtype=np.float32)
        sample["spatial_y"] = np.asarray(example["spatial_y"], dtype=np.float32)

    if include_bins:
        sample["bins"] = np.asarray(bins, dtype=np.float32)

    if include_raw_expressions:
        sample["raw_gene_expressions"] = np.asarray(
            padded_raw_gene_expressions, dtype=np.float32
        )

    for key in other_keys:
        sample[key] = example[key]

    return sample


def get_scrna_bag_of_cells_dataset(
    model_config: Any, dataset_config: DatasetConfig
) -> ScRNABagofCellsDataset:
    downsample_factor = 2**model_config.num_downsamples
    fixed_sequence_length = (
        math.ceil(len(PROTEIN_GENE_MAP) / downsample_factor) * downsample_factor
    )
    gene_map, gene_ids, gene_map_array = get_protein_gene_map()
    preprocess_fn = partial(
        construct_full_gene_sequence,
        num_bins=dataset_config.gene_expression_num_bins,
        expression_pad_token_id=model_config.pad_token_id,
        id_pad_token_id=len(gene_map),
        gene_map_array=gene_map_array,
        gene_ids=gene_ids,
        include_raw_expressions=True,
        include_bins=True,
        fixed_sequence_length=fixed_sequence_length,
        apply_log_trick=dataset_config.apply_log_trick,
    )

    dataset = ScRNABagofCellsDataset(
        epoch_num=0,
        worker_id=0,
        num_hosts=1,
        num_files_per_chunk=10,
        seed=0,
        preprocess_fn=preprocess_fn,
        resume_training=False,
        fixed_sequence_length=fixed_sequence_length,
        num_cells_per_sample=model_config.num_cells,
        expression_pad_token_id=model_config.pad_token_id,
        id_pad_token_id=len(gene_map),
        add_all_metadata=True,
        merge_fn=None,
    )
    return dataset


def get_dataset_dataloader(model_config, batch_size=1):
    dataset_config = DatasetConfig(
        min_gene_expression=0,
        gene_expression_num_bins=5,
        apply_log_trick=True,
    )
    dataset = get_scrna_bag_of_cells_dataset(model_config, dataset_config)
    dataloader = torch.utils.data.DataLoader(
        dataset,
        drop_last=True,
        batch_size=batch_size,
    )

    return dataloader
