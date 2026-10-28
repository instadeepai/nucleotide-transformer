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

"""In-memory helpers for the production genome-annotation entry point.

The research pipeline reads its sequence through a multiprocess grain dataloader
over the S3 bucket layout, and builds ground-truth targets alongside it. Production
has neither: the caller hands over FASTA files, so the whole contig is tokenized
once and windows are slices of it.

Window *placement* is not reimplemented here. ``select_windows`` from the annotate
pipeline is imported and used as-is, so production tiles a contig exactly as
``annotate/`` does, head pad and shift-back tail included. What this module adds is
the in-memory half of the research pipeline's ``annotate_one_chromosome``
bookkeeping — the N padding a boundary window needs, and the slice of each decoded
row that the window is responsible for — computed from the same geometry its
``_Tokenize`` derives it from.
"""

import gzip
import hashlib
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO
from urllib.parse import quote

import numpy as np

from nucleotide_transformer_v3.gene_annotation.annotate_utils import select_windows
from nucleotide_transformer_v3.gene_annotation.formats import (
    FORMAT_GFF3,
    FORMAT_GTF,
    gff3_lines_to_gtf,
)
from nucleotide_transformer_v3.logging import get_logger
from nucleotide_transformer_v3.tokenizers import NucleotideTokenizer

logger = get_logger()

# FASTA suffixes a directory input is globbed for, plus their gzipped forms.
# Matched case-insensitively: an assembly staged as ``chr2.FA`` is a FASTA, and
# skipping it would annotate fewer contigs than the caller staged without ever
# tripping the "directory holds no FASTA" guard.
FASTA_SUFFIXES: tuple[str, ...] = (".fa", ".fasta", ".fna")

# The character class GFF3 allows unescaped in column 1. The spec requires every
# other character to be percent-escaped, and a seqid may not be empty.
GFF3_SEQID_SAFE_CHARACTERS: frozenset[str] = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|"
)

# Stem of the name given to a record whose FASTA header carries no id at all.
UNNAMED_RECORD_STEM = "unnamed_record"

# Lowercase ACGTN — the RepeatMasker soft-mask convention.
SOFTMASKED_BASES = "acgtn"

# Bytes the nucleotide alphabet knows; anything else is folded to N and counted.
_KNOWN_BYTES = np.zeros(256, dtype=bool)
for _base in "ACGTN" + SOFTMASKED_BASES:
    _KNOWN_BYTES[ord(_base)] = True


# ---------------------------------------------------------------------------
# FASTA input resolution
# ---------------------------------------------------------------------------


def resolve_fasta_paths(fasta_paths: list[str]) -> list[Path]:
    """Expand the configured inputs into an ordered list of FASTA files.

    Each entry may be:

    - a local FASTA file;
    - a local directory — globbed for ``FASTA_SUFFIXES`` and their ``.gz`` forms
      in any case, and sorted by name, since our own bucket stores a genome as
      one FASTA per contig and a staged assembly should need no concatenation;
    - a remote FASTA (``s3://``, ``gs://``), downloaded to the local cache via
      ``download_and_cache_abstract_file`` — retried, size-verified, and durable
      for a ``datasets/`` key, so a second job reuses the bytes.

    A remote *directory* is rejected: listing a bucket prefix reliably is a
    different problem from reading a file, and naming the contigs explicitly is
    unambiguous.

    Args:
        fasta_paths: Files and/or directories, in the order they should be read.

    Returns:
        The resolved local files, directory entries sorted by name, duplicates
        removed while preserving first-seen order.

    Raises:
        ValueError: If no inputs are configured, an entry does not exist, a
            directory holds no FASTA, or a remote entry is not a FASTA file.
    """
    if not fasta_paths:
        raise ValueError(
            "fasta_paths is empty; pass at least one FASTA file or a directory "
            "of per-contig FASTAs."
        )
    resolved: list[Path] = []
    for entry in fasta_paths:
        if "://" in str(entry):
            if not is_fasta_name(str(entry).rsplit("/", 1)[-1]):
                raise ValueError(
                    f"Remote entry {entry} is not a FASTA file. A bucket prefix "
                    "cannot be listed here — name each contig's FASTA explicitly."
                )
            raise ValueError(
                f"Remote path {entry} is not supported. Download the FASTA "
                "first and pass the local path."
            )
        path = Path(entry)
        if not path.exists():
            raise ValueError(f"fasta_paths entry does not exist: {path}")
        if path.is_dir():
            found = sorted(
                child
                for child in path.iterdir()
                if is_fasta_name(child.name) and not child.is_dir()
            )
            if not found:
                raise ValueError(
                    f"Directory {path} holds no FASTA (looked for "
                    f"{', '.join(FASTA_SUFFIXES)} and their .gz forms, in any "
                    "case)."
                )
            resolved.extend(found)
        else:
            resolved.append(path)

    seen: set[str] = set()
    unique: list[Path] = []
    for path in resolved:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def is_fasta_name(name: str) -> bool:
    """Whether a file name looks like a FASTA, gzipped or not.

    Case-insensitive, because a file system is not: a directory holding
    ``chr1.fa`` beside ``chr2.FA`` would otherwise contribute one contig and
    drop the other silently, since the directory still holds *some* FASTA and
    the empty-directory guard never fires.

    Args:
        name: The file's base name.

    Returns:
        Whether the name ends in a FASTA suffix, optionally gzipped, in any
        case.
    """
    lowered = name.lower()
    stem = lowered[: -len(".gz")] if lowered.endswith(".gz") else lowered
    return stem.endswith(FASTA_SUFFIXES)


def escape_seqid(record_id: str) -> str:
    """Percent-escape a FASTA record id into a legal GFF3 ``seqid``.

    The GFF3 specification allows column 1 to hold any characters, but requires
    every character outside ``[a-zA-Z0-9.:^*$@!+_?-|]`` to be percent-escaped as
    its UTF-8 bytes. A raw FASTA header is under no such constraint, and three
    unescaped characters actively corrupt the file: a leading ``#`` turns every
    feature line for that contig into a comment, so the contig silently vanishes
    from any consumer, while ``;`` and ``=`` are the attribute-column
    separators.

    Escaping rather than rejecting is deliberate. Escaping is lossless and
    reversible (a consumer unescapes column 1 and gets the FASTA id back, which
    is what the spec's own escaping rule is for), whereas rejecting would refuse
    an assembly nobody can re-header without rewriting their FASTA, for ids the
    format is defined to carry. It also matches what :func:`part_path_for`
    already does with the same id, and keeps GFF3 the one place the answer is
    written down rather than adding a second, stricter notion of a legal contig
    name.

    Args:
        record_id: The record id as it appears in the FASTA header.

    Returns:
        The id with every character outside the GFF3 ``seqid`` class replaced by
        its percent-escaped UTF-8 bytes. An id that is already legal is returned
        unchanged, so the common case is a no-op.
    """
    return "".join(
        character
        if character in GFF3_SEQID_SAFE_CHARACTERS
        else "".join(f"%{byte:02X}" for byte in character.encode("utf-8"))
        for character in record_id
    )


def iter_fasta_records(paths: list[Path]) -> Iterator[tuple[str, str]]:
    """Yield (record id, sequence) for every record across every input file.

    One record is held in memory at a time. The id yielded here is the seqid
    every downstream output is keyed on (GFF3 column 1, the part file name, the
    FASTA descriptions), so it is normalised once, here, rather than at each
    consumer:

    - a header carrying characters GFF3 does not allow unescaped in column 1 is
      percent-escaped by :func:`escape_seqid`;
    - a header with no id at all (a bare ``>``) is named
      ``unnamed_record_<n>``, numbered across the whole input set, because an
      empty column 1 is not a GFF3 any parser accepts.

    Both are logged, since the emitted annotation then names a contig the input
    FASTA does not.

    Record ids are checked for collisions across the whole input set, after
    normalisation: annotating one contig twice would emit overlapping features
    under the same seqid, so it raises rather than silently corrupting the GFF3.
    Two headers that differ only in a character that escapes to the same text
    collide for exactly the same reason, and are caught by the same check.

    Args:
        paths: Resolved FASTA files, read in order.

    Yields:
        The normalised record id and its sequence, verbatim.

    Raises:
        ValueError: If a record id appears twice, or a file holds no record.
    """
    seen: dict[str, Path] = {}
    num_unnamed = 0
    for path in paths:
        num_records = 0
        for header_id, sequence in _iter_one_fasta(path):
            if header_id:
                record_id = escape_seqid(header_id)
                if record_id != header_id:
                    logger.warning(
                        f"{path}: record id {header_id!r} holds characters GFF3 "
                        f"does not allow unescaped in a seqid; annotating it as "
                        f"{record_id!r}."
                    )
            else:
                num_unnamed += 1
                record_id = f"{UNNAMED_RECORD_STEM}_{num_unnamed}"
                logger.warning(
                    f"{path}: a FASTA header carries no record id; annotating it "
                    f"as {record_id!r}, since an empty GFF3 seqid is unreadable."
                )
            if record_id in seen:
                raise ValueError(
                    f"Record {record_id!r} appears in both {seen[record_id]} and "
                    f"{path}; annotating one contig twice would emit duplicate "
                    "features under the same seqid."
                )
            seen[record_id] = path
            num_records += 1
            yield record_id, sequence
        if num_records == 0:
            raise ValueError(f"{path} holds no FASTA record.")


def _iter_one_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (record id, sequence) from one FASTA, gzipped or plain.

    Hand-rolled rather than ``Bio.SeqIO`` so the sequence is assembled with a
    single ``"".join`` per record and no ``Seq`` object is built. Our own bucket
    FASTAs are written unwrapped (one line per record), which this reads without
    copying every line twice.
    """
    header: str | None = None
    chunks: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                # A header may legally be bare ">" or ">" plus only whitespace,
                # so take the first field only when there is one rather than
                # indexing an empty split and killing the whole job. The empty
                # id is given a name by ``iter_fasta_records``, which is where
                # every id is made legal for GFF3 column 1.
                fields = line[1:].split()
                header = fields[0] if fields else ""
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)
        _reject_a_partial_read(path, handle)


def _file_size(path: Path) -> int:
    """Bytes the file holds."""
    return path.stat().st_size


def _reject_a_partial_read(path: Path, handle: IO[str]) -> None:
    """Refuse a FASTA whose final byte was never reached.

    A reader that stops early is indistinguishable, downstream, from an
    assembly that simply holds fewer contigs: the run annotates what it saw,
    writes a well-formed GFF3 and reports success. That is not hypothetical. A
    2.38 GB, 135-record assembly was once read as 72 records and 1.60 GB, and
    every line of the output was correct for the 68% of the genome it covered,
    so nothing anywhere in the pipeline objected. A user hits the same thing
    with an interrupted download, which leaves a plausible FASTA truncated at
    an arbitrary offset.

    So the bytes consumed are compared against the bytes the file holds rather
    than assumed equal. Reaching EOF is what the reader already does; the point
    here is that EOF was the file's *end* and not merely where the data stopped.

    Uncompressed input only: a text handle over gzip reports an uncompressed
    offset, which is not comparable to a compressed size, and a truncated gzip
    member already raises :class:`EOFError` of its own accord.

    Args:
        path: The file just read, used for its size and in the message.
        handle: The handle that read it, positioned wherever it stopped.

    Raises:
        ValueError: If the handle stopped before the end of the file.
    """
    if str(path).lower().endswith(".gz"):
        return
    consumed = handle.tell()
    total = _file_size(path)
    if consumed >= total:
        return
    raise ValueError(
        f"{path}: read {consumed:,} of {total:,} bytes before hitting end of "
        f"input, so {total - consumed:,} bytes were never annotated. The file "
        "is truncated, or it changed while it was being read. Annotating part "
        "of an assembly and reporting success is the one outcome this refuses."
    )


def _open_text(path: Path) -> IO[str]:
    """Open a FASTA for text reading, transparently decompressing ``.gz``.

    The suffix is matched case-insensitively, for the same reason
    :func:`is_fasta_name` does: ``chr4.fna.GZ`` is accepted as an input, so it has
    to be read as gzip rather than handed to the text reader as bytes.
    """
    if str(path).lower().endswith(".gz"):
        return gzip.open(path.open("rb"), mode="rt")
    return path.open("r")


# ---------------------------------------------------------------------------
# Tokenization
# ---------------------------------------------------------------------------


def build_token_lut(
    tokenizer: NucleotideTokenizer, use_softmasked: bool
) -> tuple[np.ndarray, int]:
    """Build the 256-row byte -> token-id table used to tokenize a whole contig.

    ``NucleotideTokenizer.np_tokenize`` joins the sequence into a comma-separated
    string and rewrites it once per vocabulary entry, which costs ~12 ms per 131 kb
    window; a table lookup over the raw bytes is ~15x faster and, applied to the
    contig once, makes per-window tokenization free (a window becomes a slice).

    The table also folds what ``sanitize_iupac_codes`` folds, with one deliberate
    difference. Without a soft-masked vocabulary, lowercase ``acgtn`` is mapped to
    the *uppercase* token rather than to ``N``: the regex class ``[^ACGTN]`` treats
    ``a`` as unknown, which is harmless in the research pipeline (it selects the
    uppercase FASTA variant for such a checkpoint) but would erase every
    repeat-masked base of a soft-masked assembly handed to production. Everything
    outside the alphabet still becomes ``N``, as in training.

    Args:
        tokenizer: The checkpoint's restored sequence tokenizer.
        use_softmasked: Whether the checkpoint's vocabulary carries the lowercase
            soft-mask tokens.

    Returns:
        The (256,) uint8 table indexed by byte value, and the ``N`` token id used
        for padding and for out-of-alphabet bytes.
    """
    n_token_id = tokenizer.token_to_id("N")
    lut = np.full(256, n_token_id, dtype=np.uint8)
    for base in "ACGTN":
        lut[ord(base)] = tokenizer.token_to_id(base)
    for base in SOFTMASKED_BASES:
        source = base if use_softmasked else base.upper()
        lut[ord(base)] = tokenizer.token_to_id(source)
    return lut, int(n_token_id)


def tokenize_contig(sequence: str, lut: np.ndarray) -> tuple[np.ndarray, int]:
    """Tokenize a whole contig through the byte table.

    Args:
        sequence: One FASTA record's sequence.
        lut: The table from :func:`build_token_lut`.

    Returns:
        The (len(sequence),) uint8 token ids — never shorter than the input, so
        genomic coordinates cannot shift — and the number of bytes that fell
        outside the alphabet and were folded to ``N``.
    """
    raw = np.frombuffer(sequence.encode("ascii", errors="replace"), dtype=np.uint8)
    tokens = lut[raw]
    num_folded = int((~_KNOWN_BYTES[raw]).sum())
    return tokens, num_folded


# ---------------------------------------------------------------------------
# Window planning
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ProductionWindow:
    """One window's model input bounds and the slice of its decode that is kept.

    The research pipeline's ``annotate_one_chromosome`` groups windows into passes
    because ``left_pad``
    applies to a whole grain loader; here padding is per-row array slicing, so
    every window carries its own and any set of them can share a batch.

    Attributes:
        left_pad: N token ids prepended to the row (non-zero only for a window
            reaching ahead of bp 0).
        element_start: First contig bp the row reads, i.e. the window start after
            the left pad is applied. Never negative.
        available: Real bp the row reads from ``element_start``.
        room: Row length after the left pad; ``left_pad + room`` is the window
            length, and ``room - available`` is the trailing N pad.
        keep_start: Genomic bp the retained slice begins at — the key this
            window's states are stored under.
        slice_start: First position of the decoded row that is retained; non-zero
            only for the shifted tail window, whose lead was already covered.
        slice_stop: One past the last retained position, cut at the contig end so
            padding never reaches the stitched state path.
    """

    left_pad: int
    element_start: int
    available: int
    room: int
    keep_start: int
    slice_start: int
    slice_stop: int

    @property
    def is_empty(self) -> bool:
        """Whether this window contributes no bp (fully covered by earlier ones)."""
        return self.slice_stop <= self.slice_start


def plan_contig_windows(
    *,
    sequence_length: int,
    keep_target_center_fraction: float,
    contig_length: int,
    chrom: str,
) -> tuple[list[ProductionWindow], int]:
    """Tile one contig and resolve each window's padding and retained slice.

    Placement comes from ``select_windows``, so the tiling — a stride grid anchored
    at bp 0, a first window reaching ``crop`` bp ahead of the contig, and a last
    window shifted back onto the final base — is the annotate pipeline's, not a
    second implementation of it. The per-window values mirror what the research
    pipeline's ``_Tokenize`` and ``annotate_one_chromosome`` compute for the same
    window: ``available`` and ``room`` reproduce its front-pad/cut/overhang-pad
    sequence, and ``slice_stop`` reproduces the ``valid_length`` its target mask
    encodes.

    Args:
        sequence_length: Window length fed to the model (bp).
        keep_target_center_fraction: Fraction of the window the head predicts.
        contig_length: True length of the contig (bp).
        chrom: Contig name, for the error messages ``select_windows`` raises.

    Returns:
        The planned windows in ascending genomic order (empty ones dropped) and
        the center offset, which is the whole per-side crop.
    """
    windows, center_offset = select_windows(
        sequence_length=sequence_length,
        keep_target_center_fraction=keep_target_center_fraction,
        chrom_lengths={chrom: contig_length},
        chrom=chrom,
    )
    if not windows:
        return [], center_offset

    crop = center_offset
    stride = sequence_length - 2 * crop
    *inner, last = windows

    # The shifted tail window starts inside its predecessor's kept region; those bp
    # belong to the earlier window, so the tail drops them from its own row.
    inner_kept_end = inner[-1].start + crop + stride if inner else 0
    overlaps = {id(w): 0 for w in inner}
    overlaps[id(last)] = max(0, inner_kept_end - (last.start + crop))

    planned: list[ProductionWindow] = []
    for window in windows:
        left_pad = max(0, -window.start)
        element_start = window.start + left_pad
        keep_start = window.start + crop
        keep_end = keep_start + stride
        room = sequence_length - left_pad
        existing = max(0, min(sequence_length, contig_length - element_start))
        available = min(room, existing)
        valid_length = max(0, min(keep_end, element_start + available) - keep_start)
        slice_start = overlaps[id(window)]
        candidate = ProductionWindow(
            left_pad=left_pad,
            element_start=element_start,
            available=available,
            room=room,
            keep_start=keep_start + slice_start,
            slice_start=slice_start,
            slice_stop=valid_length,
        )
        if not candidate.is_empty:
            planned.append(candidate)
    return planned, center_offset


def build_window_row(
    tokens: np.ndarray, window: ProductionWindow, n_token_id: int
) -> np.ndarray:
    """Assemble one window's token row from the tokenized contig.

    Reproduces the research pipeline's ``_Tokenize`` front-pad /
    cut-at-contig-end / pad-the-overhang sequence, in token space rather than
    string space.

    Args:
        tokens: The whole contig's token ids from :func:`tokenize_contig`.
        window: The planned window.
        n_token_id: Token id for ``N``, used for both pads.

    Returns:
        A (left_pad + room,) uint8 row, i.e. exactly one window length.
    """
    row = np.full(window.left_pad + window.room, n_token_id, dtype=tokens.dtype)
    stop = window.element_start + window.available
    row[window.left_pad : window.left_pad + window.available] = tokens[
        window.element_start : stop
    ]
    return row


# ---------------------------------------------------------------------------
# Cross-window stitching
# ---------------------------------------------------------------------------


def assemble_runs(
    states_plus: dict[int, np.ndarray],
    states_minus: dict[int, np.ndarray],
    base_codes_by_start: dict[int, np.ndarray] | None,
    expected_starts: list[int],
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray | None, int]]:
    """Group per-window decoded states into contiguous genomic runs.

    Adapted from ``AnnotationMetrics.assemble_runs`` without the metrics object.
    Row lengths are measured individually rather than assumed uniform: the shifted
    tail window contributes only its uncovered suffix, and a contig shorter than
    one kept region is cut at its end, so a single shared length would read the
    short final row as a gap and split the run.

    Args:
        states_plus: Plus-strand states keyed on genomic start.
        states_minus: Minus-strand states, same keys.
        base_codes_by_start: Optional per-base codes, same keys.
        expected_starts: Every start the plan produced.

    Returns:
        One (plus, minus, base codes or None, genomic start) tuple per contiguous
        run, ascending.

    Raises:
        RuntimeError: If any expected window is missing from the decode output.
    """
    if not expected_starts:
        return []
    missing = set(expected_starts) - set(states_plus)
    if missing:
        raise RuntimeError(
            f"{len(missing)} windows missing from decode output "
            f"(first 5: {sorted(missing)[:5]})"
        )
    sorted_starts = sorted(expected_starts)
    have_codes = bool(base_codes_by_start)

    runs: list[list[int]] = [[sorted_starts[0]]]
    for start in sorted_starts[1:]:
        previous = runs[-1][-1]
        previous_end = previous + states_plus[previous].shape[0]
        if start == previous_end:
            runs[-1].append(start)
        elif start < previous_end:
            # Two windows claiming the same bp: concatenating them would double
            # every position between them and shift the rest of the contig, so
            # one contig would emit two gene sets over the same coordinates.
            # AnnotationMetrics.assemble_runs in annotate/metrics.py raises here
            # and this did not, which left plan_contig_windows as the only thing
            # standing between an overlap and a silently doubled annotation. A
            # 129,776-geometry sweep says it is unreachable today; that makes
            # this cheap, not unnecessary.
            raise RuntimeError(
                f"window {start} overlaps the one ending at {previous_end} by "
                f"{previous_end - start} bp; each bp must be owned by exactly "
                "one window."
            )
        else:
            runs.append([start])

    return [
        (
            np.concatenate([states_plus[s] for s in run]),
            np.concatenate([states_minus[s] for s in run]),
            (
                np.concatenate([base_codes_by_start[s] for s in run])
                if have_codes and base_codes_by_start is not None
                else None
            ),
            run[0],
        )
        for run in runs
    ]


# ---------------------------------------------------------------------------
# Output: per-record part files, merged into one GFF3
# ---------------------------------------------------------------------------

# Suffix of a finished per-record part. Its presence is what a resumed run reads
# to decide a record is already done, so it is only ever created by an atomic
# rename from a ``.tmp`` sibling: a process killed mid-write leaves the staging
# file, which both the resume check and the merge ignore.
PART_SUFFIX = ".gff3part"

# First line of a part, carrying what the part was produced from so a resumed run
# can tell whether it is still valid. It is written as a GFF3 comment but never
# reaches the output: parts are read back through :func:`_split_part`, which
# takes it off again, so a merged file is exactly the features.
PART_HEADER_PREFIX = "#ntv3-part "

# Bounded sample a sequence fingerprint is taken over: this many blocks of this
# many characters, evenly spaced from the first character to the last. The cost
# is therefore the same for a 250 Mb contig as for a 250 kb one.
FINGERPRINT_NUM_BLOCKS = 16
FINGERPRINT_BLOCK_CHARS = 65_536


def part_path_for(parts_dir: Path, record_id: str) -> Path:
    """Path of one record's part file.

    The record id is percent-escaped again on top of the GFF3 escaping
    :func:`iter_fasta_records` applies, because a legal GFF3 seqid may still
    contain ``|``, ``*`` or ``?``, which are a glob or an unusable file name on
    some file systems. Escaping is injective, so distinct records keep distinct
    parts either way.

    Args:
        parts_dir: Directory the part files live in.
        record_id: The record's seqid.

    Returns:
        The part file's path, which exists only once the record is finished.
    """
    return parts_dir / f"{quote(record_id, safe='')}{PART_SUFFIX}"


def fingerprint_sequence(sequence: str) -> str:
    """Digest a bounded, evenly spaced sample of one record's sequence.

    Deliberately not a digest of the whole contig. A resumed run fingerprints
    every record it is handed, *including the ones it is about to skip*, which
    is the whole point of resuming, so the check has to cost the same on a
    250 Mb contig as on a 250 kb one. ``FINGERPRINT_NUM_BLOCKS`` blocks of
    ``FINGERPRINT_BLOCK_CHARS`` characters, spaced from the first character to
    the last, cap it at 1 MiB hashed per record; a contig that fits inside that
    budget is hashed whole.

    What it catches, together with the length recorded beside it in
    :class:`PartIdentity`:

    * every change in length, so every insertion, deletion, truncation or
      extension, and any different contig reused under the same id;
    * any length-preserving substitution landing in one of the sampled blocks,
      which for a contig up to 1 MiB is all of them.

    What it does **not** catch: a length-preserving edit that falls entirely
    between two sampled blocks of a contig larger than 1 MiB, for instance a
    single corrected base or a soft-mask change in the middle of a chromosome.
    Detecting those needs a full digest of the assembly, and the honest place
    for that is the caller's own checksum of the FASTA, not a per-resume scan.

    Args:
        sequence: One FASTA record's sequence, exactly as it will be annotated.

    Returns:
        The digest, as hex.
    """
    digest = hashlib.blake2b(digest_size=16)
    budget = FINGERPRINT_NUM_BLOCKS * FINGERPRINT_BLOCK_CHARS
    if len(sequence) <= budget:
        digest.update(sequence.encode("utf-8", errors="replace"))
        return digest.hexdigest()
    span = len(sequence) - FINGERPRINT_BLOCK_CHARS
    for index in range(FINGERPRINT_NUM_BLOCKS):
        start = round(index * span / (FINGERPRINT_NUM_BLOCKS - 1))
        block = sequence[start : start + FINGERPRINT_BLOCK_CHARS]
        digest.update(block.encode("utf-8", errors="replace"))
    return digest.hexdigest()


def build_settings_digest(
    provenance_lines: list[str],
    *,
    species: str,
    mask_species_condition: bool,
) -> str:
    """Digest everything about a run that decides what its features are.

    Built from the provenance block, which already records the checkpoint, the
    geometry, the decoder and the filters, plus the species conditioning, which
    changes the model's predictions without appearing in the output header.

    The output *format* is deliberately absent: parts always hold GFF3, so a
    GFF3 job and a GTF job over the same input produce the same parts and may
    resume from each other, which is the behaviour the writer documents.

    Args:
        provenance_lines: The run's ``##ntv3-*`` provenance lines.
        species: The species condition, empty when there is none.
        mask_species_condition: Whether the species conditioning was dropped.

    Returns:
        The digest, as hex.
    """
    payload = "\n".join(
        [
            *provenance_lines,
            f"species {species}",
            f"mask-species-condition {str(mask_species_condition).lower()}",
        ]
    )
    return hashlib.blake2b(payload.encode("utf-8"), digest_size=16).hexdigest()


@dataclass(frozen=True)
class PartIdentity:
    """What one part file was produced from.

    A part's file name records only the record id, which says nothing about the
    sequence that was under that id or the settings it was annotated with. This
    is what the part carries so a resumed run can tell that it is looking at
    finished work for the job in front of it rather than for a different one.

    Attributes:
        sequence_length: Length in bp of the record the part annotates. Kept
            beside the digest rather than folded into it so a mismatch can name
            both lengths.
        sequence_digest: :func:`fingerprint_sequence` of that record.
        settings_digest: :func:`build_settings_digest` of the run that wrote it.
    """

    sequence_length: int
    sequence_digest: str
    settings_digest: str

    def render(self) -> str:
        """The header line a part carries this identity in."""
        return (
            f"{PART_HEADER_PREFIX}length={self.sequence_length} "
            f"sequence={self.sequence_digest} settings={self.settings_digest}"
        )

    @classmethod
    def parse(cls, line: str) -> "PartIdentity | None":
        """Read an identity back from a header line.

        Args:
            line: The candidate header line.

        Returns:
            The identity, or None if the line is not a header or is malformed.
            Both mean the same thing to a resumed run: the part cannot vouch for
            itself, so it is not reused.
        """
        if not line.startswith(PART_HEADER_PREFIX):
            return None
        fields: dict[str, str] = {}
        for token in line[len(PART_HEADER_PREFIX) :].split():
            key, separator, value = token.partition("=")
            if separator:
                fields[key] = value
        if not {"length", "sequence", "settings"} <= fields.keys():
            return None
        try:
            length = int(fields["length"])
        except ValueError:
            return None
        return cls(
            sequence_length=length,
            sequence_digest=fields["sequence"],
            settings_digest=fields["settings"],
        )


def build_part_identity(sequence: str, settings_digest: str) -> PartIdentity:
    """The identity a part written for ``sequence`` by this run must carry.

    Args:
        sequence: The record's sequence, exactly as it will be annotated.
        settings_digest: This run's :func:`build_settings_digest`.

    Returns:
        The identity.
    """
    return PartIdentity(
        sequence_length=len(sequence),
        sequence_digest=fingerprint_sequence(sequence),
        settings_digest=settings_digest,
    )


def describe_part_mismatch(stored: PartIdentity | None, current: PartIdentity) -> str:
    """Say why a part cannot be reused for the record in front of it.

    Args:
        stored: The identity read out of the part, or None when it carries none
            (a part written before parts carried one, or a corrupted header).
        current: The identity this run would write for the same record.

    Returns:
        A sentence naming what differs, or an empty string when the part is
        valid for this run.
    """
    if stored is None:
        return "it carries no identity header, so what produced it is unknown"
    if stored.sequence_length != current.sequence_length:
        return (
            f"it was written for a {stored.sequence_length:,} bp sequence and the "
            f"record now under that id is {current.sequence_length:,} bp"
        )
    if stored.sequence_digest != current.sequence_digest:
        return (
            "the sequence under that record id has changed "
            f"({stored.sequence_digest} -> {current.sequence_digest})"
        )
    if stored.settings_digest != current.settings_digest:
        return (
            "it was produced under different settings "
            f"({stored.settings_digest} -> {current.settings_digest})"
        )
    return ""


def write_record_part(
    parts_dir: Path,
    record_id: str,
    lines: list[str],
    identity: PartIdentity | None = None,
) -> Path:
    """Write one record's GFF3 lines to its part file, sorted by start.

    Sorting per record and merging in record-id order reproduces ``write_gff3``'s
    global ``(seqid, start)`` order exactly, while holding only one record in
    memory. The sort is stable, so a gene row keeps preceding its own
    sub-features as ``find_genes`` emitted them even though they share a start.

    A record with no genes still gets an empty part, or a resumed run would
    re-annotate it every time.

    Args:
        parts_dir: Directory the part files live in.
        record_id: The record's seqid.
        lines: The record's GFF3 feature lines.
        identity: What the part was produced from, written as its header line so
            a resumed run can check it. None writes a part with no header, which
            a resumed run refuses to reuse.

    Returns:
        The path written.
    """
    parts_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(lines, key=lambda line: int(line.split("\t")[3]))
    if identity is not None:
        ordered = [identity.render(), *ordered]
    final = part_path_for(parts_dir, record_id)
    staging = final.parent / f"{final.name}.tmp"
    staging.write_text("".join(f"{line}\n" for line in ordered))
    staging.replace(final)
    return final


def _split_part(text: str) -> tuple[PartIdentity | None, list[str]]:
    """Separate a part's identity header from its feature lines."""
    lines = [line for line in text.splitlines() if line]
    if lines and lines[0].startswith(PART_HEADER_PREFIX):
        return PartIdentity.parse(lines[0]), lines[1:]
    return None, lines


def read_part(
    parts_dir: Path, record_id: str
) -> tuple[PartIdentity | None, list[str]]:
    """Read back one finished record's identity and GFF3 feature lines.

    Args:
        parts_dir: Directory the part files live in.
        record_id: The record's seqid.

    Returns:
        What the part says it was produced from (None when it says nothing) and
        its GFF3 feature lines, without trailing newlines.
    """
    return _split_part(part_path_for(parts_dir, record_id).read_text())


def max_gene_id_in_parts(parts_dir: Path) -> int:
    """Highest gene number any part already staged in ``parts_dir`` carries.

    Gene ids are ``gene_<n>_<strand>`` from a counter threaded across records, so
    a resumed run that restarted that counter at zero would hand ids it already
    used in the parts it skipped to the records it still has to do. GFF3 requires
    ``ID`` to be unique within a file, and our own Tiberius-style evaluator keys
    genes on it: on a collision it drops the second ``gene`` row and appends its
    CDS rows to the first, so two genes on two contigs are scored as one gene
    whose CDS spans both. Seeding the counter past everything on disk is what
    makes that impossible.

    Every part in the directory is read, not only the ones this job will reuse,
    because which records the job covers is not known until its FASTAs have been
    read, and a part it never touches costs it nothing but a skipped id. The
    price of being conservative is that gene numbering may have gaps after a
    resume; the price of not being is a spec violation.

    Args:
        parts_dir: Directory the part files live in; it need not exist.

    Returns:
        The highest gene number found, or 0 for an empty or absent directory.
    """
    if not parts_dir.exists():
        return 0
    highest = 0
    for path in parts_dir.iterdir():
        if not path.name.endswith(PART_SUFFIX):
            continue
        for line in path.read_text().splitlines():
            columns = line.split("\t")
            if len(columns) != 9 or columns[2] != "gene":
                continue
            fields = columns[8].removeprefix("ID=").split("_")
            if len(fields) == 3 and fields[0] == "gene" and fields[1].isdigit():
                highest = max(highest, int(fields[1]))
    return highest


def merge_parts(
    output_path: str | Path,
    parts_dir: Path,
    provenance: list[str],
    record_ids: list[str],
    *,
    output_format: str = FORMAT_GFF3,
    record_lengths: dict[str, int] | None = None,
) -> int:
    """Concatenate this job's parts into one annotation file, gzipping a ``.gz`` path.

    Parts are ordered by their record id, which is the key ``write_gff3`` sorts
    its ``seqid`` column on, so a merged GFF3 is byte-identical to what
    buffering the whole genome and calling ``write_gff3`` would have produced,
    comment lines aside.

    Only ``record_ids`` are merged, never whatever the directory happens to
    hold. Parts deliberately survive a killed job so a restart can skip finished
    records, which means a later job writing to the same output path meets the
    earlier one's parts; globbing the directory would silently emit records the
    current job never annotated.

    Parts always hold GFF3, whatever the requested output format: it is the
    canonical intermediate, so a job started as one format can be resumed as the
    other, and the GTF is by construction a re-expression of the same intervals.

    Args:
        output_path: Destination file; gzipped when it ends in ``.gz``.
        parts_dir: Directory holding the part files.
        provenance: ``##`` comment lines inserted after the version pragma
            (GFF3) or at the top of the file (GTF, which has no pragma).
        record_ids: Every record this job covered, annotated or resumed. Order
            is irrelevant; they are sorted here.
        output_format: ``"gff3"`` or ``"gtf"``, as
            ``formats.resolve_output_format`` reads it off ``output_path``.

    Returns:
        The number of feature lines written.

    Raises:
        RuntimeError: If a named record has no part, which means it was never
            finished and the output would be silently incomplete.
        ValueError: If ``output_format`` is not a known format.
    """
    if output_format not in {FORMAT_GFF3, FORMAT_GTF}:
        raise ValueError(
            f"output_format={output_format!r} is not one of "
            f"{FORMAT_GFF3!r}, {FORMAT_GTF!r}."
        )
    missing = [
        record_id
        for record_id in record_ids
        if not part_path_for(parts_dir, record_id).exists()
    ]
    if missing:
        raise RuntimeError(
            f"{len(missing)} record(s) have no part file in {parts_dir} and would "
            f"be missing from the output: {sorted(missing)[:5]}. A part is written "
            "only once a record is fully annotated, so this means the job did not "
            "finish them."
        )
    parts = sorted(
        (record_id, part_path_for(parts_dir, record_id)) for record_id in record_ids
    )
    # Read through _split_part so a part's identity header stays an internal
    # bookkeeping line and never reaches the merged annotation.
    body = [
        "".join(f"{line}\n" for line in _split_part(path.read_text())[1])
        for _, path in parts
    ]
    if output_format == FORMAT_GTF:
        body = [
            "".join(f"{line}\n" for line in gff3_lines_to_gtf(chunk.splitlines()))
            for chunk in body
        ]
    comments = "".join(f"{line}\n" for line in provenance)
    # ``##sequence-region`` per record, which is spec-conformant GFF3 and what a
    # downstream tool reads to know a contig's extent.
    regions = "".join(
        f"##sequence-region {record_id} 1 {record_lengths[record_id]}\n"
        for record_id in record_ids
        if record_lengths and record_id in record_lengths
    )
    # GTF has no version pragma and no sequence-region directive; GFF3 requires
    # ``##gff-version 3`` to be first.
    header = (
        comments
        if output_format == FORMAT_GTF
        else f"##gff-version 3\n{regions}{comments}"
    )
    text = header + "".join(body)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Case-insensitive to match resolve_output_format and the input side: a
    # path accepted as gzipped must also be WRITTEN gzipped, or the reader
    # gets a plain file named .GZ.
    if str(out).lower().endswith(".gz"):
        out.write_bytes(gzip.compress(text.encode("utf-8")))
    else:
        out.write_text(text)

    num_features = sum(1 for chunk in body for line in chunk.splitlines() if line)
    logger.info(
        f"Wrote {num_features} {output_format.upper()} features from "
        f"{len(parts)} record(s) to {out}"
    )
    return num_features


def build_provenance_lines(
    *,
    checkpoint_path: str,
    model_name: str,
    sequence_length: int,
    keep_target_center_fraction: float,
    decoder: str,
    apply_biological_filters: bool,
    min_cds_length: int,
    use_bf16_compute: bool,
    species: str | None,
    mask_species_condition: bool,
    species_conditioning_consumed: bool,
    backend: str,
) -> list[str]:
    """Build the ``##`` comments describing how a GFF3 was produced.

    GFF3 comment lines are inert for tabix and IGV, so this costs consumers
    nothing while making a shared annotation reproducible: without it, nobody
    receiving the file can tell which checkpoint, geometry or filters produced
    it, and the filter settings in particular change what the numbers mean.

    Args:
        checkpoint_path: The weights this run used.
        model_name: The model name saved in the checkpoint's training config.
        sequence_length: Input window in bp.
        keep_target_center_fraction: The decoded centre fraction.
        decoder: Decoder name.
        apply_biological_filters: Whether the in-frame-stop filter ran.
        min_cds_length: Minimum total CDS length kept.
        use_bf16_compute: Whether the trunk ran in bf16.
        species: The species the caller asked for, if any.
        mask_species_condition: Whether the caller asked to mask it.
        backend: The **resolved** platform the run executed on ("gpu"/"cpu"),
            not the requested ``backend`` string, which may be empty for
            "whatever JAX picks" and so cannot distinguish two hosts.
        species_conditioning_consumed: Whether the model actually reads the
            species condition. When it does not, the species lines are left out
            rather than recorded: writing ``##ntv3-species cow`` for a model
            that cannot see it invites a reader to attribute a difference to
            conditioning that conditioning did not cause. That mistake has
            already been made once from these very headers.

    Returns:
        One ``##ntv3-*`` line per recorded setting, in a stable order.
    """
    lines = [
        f"##ntv3-checkpoint {checkpoint_path}",
        f"##ntv3-model-name {model_name}",
        f"##ntv3-sequence-length {sequence_length}",
        f"##ntv3-keep-target-center-fraction {keep_target_center_fraction}",
        f"##ntv3-decoder {decoder}",
        f"##ntv3-apply-biological-filters {str(apply_biological_filters).lower()}",
        f"##ntv3-min-cds-length {min_cds_length}",
        f"##ntv3-use-bf16-compute {str(use_bf16_compute).lower()}",
        # The platform is an output-changing setting, so it belongs here and,
        # through here, in the settings digest. fp32 is reproducible across CPU
        # processes and NOT across GPU ones, and CPU bf16 is emulated rather
        # than native, so GPU-produced and CPU-produced parts are not
        # interchangeable. Without this a job started on GPU and resumed with
        # backend=cpu reuses the GPU parts beside fresh CPU ones and the merged
        # file records nothing that tells them apart. use_bf16_compute is
        # digested for precisely this reason, which made the omission an
        # inconsistency rather than a judgement call.
        f"##ntv3-backend {backend}",
    ]
    if species_conditioning_consumed:
        lines.append(f"##ntv3-species {species or ''}")
        lines.append(
            f"##ntv3-mask-species-condition {str(mask_species_condition).lower()}"
        )
    else:
        lines.append("##ntv3-species-conditioning inert")
    return lines
