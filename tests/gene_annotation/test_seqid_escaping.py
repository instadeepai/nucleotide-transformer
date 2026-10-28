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

"""Tests that FASTA ids are escaped into legal GFF3 seqids.

Spec: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

import string
from pathlib import Path
from urllib.parse import unquote

import pytest

from nucleotide_transformer_v3.gene_annotation.annotate_genome import run
from nucleotide_transformer_v3.gene_annotation.utils import (
    GFF3_SEQID_SAFE_CHARACTERS,
    escape_seqid,
    iter_fasta_records,
)

# The character class the GFF3 specification names for column 1, transcribed
# from the spec rather than from the implementation: "IDs may contain any
# characters, but must escape any characters not in the set
# [a-zA-Z0-9.:^*$@!+_?-|]".
SPEC_SEQID_SAFE_CHARACTERS: frozenset[str] = frozenset(
    string.ascii_letters + string.digits + ".:^*$@!+_?-|"
)

# Every printable ASCII character the spec reserves, plus tab.
RESERVED_ASCII: tuple[str, ...] = tuple(
    sorted(
        {chr(code) for code in range(0x20, 0x7F)} - SPEC_SEQID_SAFE_CHARACTERS | {"\t"}
    )
)

# The three whose damage is silent, called out so a regression names itself.
CORRUPTING: tuple[tuple[str, str], ...] = (("#", "%23"), (";", "%3B"), ("=", "%3D"))


def test_the_safe_class_is_the_one_the_specification_names() -> None:
    """The implementation's class is the spec's, character for character."""
    assert GFF3_SEQID_SAFE_CHARACTERS == SPEC_SEQID_SAFE_CHARACTERS


# ---------------------------------------------------------------------------
# Legal ids are left alone
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("character", sorted(SPEC_SEQID_SAFE_CHARACTERS))
def test_a_character_the_specification_allows_is_not_escaped(character: str) -> None:
    """Escaping a legal character would rename every contig in every assembly."""
    assert escape_seqid(character) == character


@pytest.mark.parametrize(
    "record_id",
    ["chr1", "NC_000001.11", "scaffold-12", "gi|1234|ref", "chr1:1-100", "Chr1.2"],
)
def test_an_already_legal_record_id_is_returned_unchanged(record_id: str) -> None:
    """The common case is a no-op, so ordinary assemblies keep their names."""
    assert escape_seqid(record_id) == record_id


# ---------------------------------------------------------------------------
# Everything the specification reserves is escaped
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("character", RESERVED_ASCII)
def test_a_character_the_specification_reserves_is_percent_escaped(
    character: str,
) -> None:
    """Swept over the whole reserved class, not the examples we thought of."""
    assert escape_seqid(character) == f"%{ord(character):02X}"


@pytest.mark.parametrize(("character", "escaped"), CORRUPTING)
def test_the_silently_corrupting_characters_are_escaped(
    character: str, escaped: str
) -> None:
    """``#``, ``;`` and ``=`` are the ones that damage the file, not the name."""
    assert escape_seqid(f"{character}chr1") == f"{escaped}chr1"


@pytest.mark.parametrize("whitespace", [" ", "\t"])
def test_whitespace_is_escaped(whitespace: str) -> None:
    """A GFF3 row is tab-separated and the spec forbids unescaped whitespace."""
    assert whitespace not in escape_seqid(f"chr{whitespace}1")


def test_a_leading_angle_bracket_is_escaped() -> None:
    """The spec singles it out: a seqid must not begin with an unescaped ``>``."""
    assert escape_seqid(">chr1") == "%3Echr1"


def test_the_percent_sign_is_itself_escaped() -> None:
    """Otherwise ``chr%41`` would unescape to ``chrA`` and rename the contig."""
    assert escape_seqid("chr%41") == "chr%2541"


@pytest.mark.parametrize(
    ("record_id", "escaped"),
    [
        ("chrü", "chr%C3%BC"),
        ("chrα", "chr%CE%B1"),
        ("染色体1", "%E6%9F%93%E8%89%B2%E4%BD%931"),
    ],
)
def test_a_non_ascii_record_id_is_escaped_as_its_utf8_bytes(
    record_id: str, escaped: str
) -> None:
    """The spec escapes bytes, so a multi-byte character becomes several."""
    assert escape_seqid(record_id) == escaped


def test_only_the_reserved_characters_are_touched_in_a_mixed_id() -> None:
    """A realistic id keeps its legal run and escapes the rest in place."""
    assert escape_seqid("#chr1;a=b 2") == "%23chr1%3Ba%3Db%202"


# ---------------------------------------------------------------------------
# Escaping is lossless
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "record_id",
    ["chr1", "#chr1", "chr 1", "a;b=c", "chr%41", "chrü", "gi|1234|ref", ">chr1"],
)
def test_escaping_round_trips_back_to_the_fasta_header(record_id: str) -> None:
    """A consumer unescaping column 1 gets the FASTA id back."""
    assert unquote(escape_seqid(record_id)) == record_id


def test_two_record_ids_that_differ_stay_different_once_escaped() -> None:
    """Injective, or two distinct contigs would collide as one seqid."""
    ids = ["#chr1", "chr1", "chr 1", "chr%201", "chr_1", "chrü", "chr%C3%BC"]
    assert len({escape_seqid(record_id) for record_id in ids}) == len(ids)


# ---------------------------------------------------------------------------
# The escaped id is what the pipeline carries
# ---------------------------------------------------------------------------


def test_a_reserved_header_reaches_the_pipeline_already_escaped(
    tmp_path: Path,
) -> None:
    """``iter_fasta_records`` escapes once, so no consumer sees the raw id."""
    fasta = tmp_path / "odd.fa"
    fasta.write_text(">#chr1\nACGT\n>chr;2\nACGT\n")

    records = list(iter_fasta_records([fasta]))

    assert [record_id for record_id, _ in records] == ["%23chr1", "chr%3B2"]


def _feature_rows(gff3: str) -> list[list[str]]:
    """The feature rows of a GFF3, split into columns."""
    return [line.split("\t") for line in gff3.splitlines() if not line.startswith("#")]


@pytest.mark.parametrize("record_id", ["#chrHead", "chr;Head", "a=b", "chrü", "chr%41"])
def test_a_reserved_contig_name_keeps_its_features_in_the_output(
    record_id: str,
    synthetic_records: dict[str, str],
    stub_model_load: None,
    production_config,
    tmp_path: Path,
) -> None:
    """``run()`` writes the same features under the escaped seqid as a plain one."""
    sequence = synthetic_records["chrHead"]
    rows = {}
    for name in ("chrHead", record_id):
        fasta = tmp_path / f"{len(rows)}.fa"
        fasta.write_text(f">{name}\n{sequence}\n")
        output = tmp_path / f"{len(rows)}.gff3"
        run(production_config(fasta_paths=[str(fasta)], output_path=str(output)))
        rows[name] = _feature_rows(output.read_text())

    plain, escaped = rows["chrHead"], rows[record_id]
    assert plain, "the synthetic contig produced no features"
    assert {row[0] for row in escaped} == {escape_seqid(record_id)}
    assert [row[1:] for row in escaped] == [row[1:] for row in plain]
