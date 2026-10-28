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

"""Tests that a resumed run only reuses parts written under the same settings."""

from pathlib import Path

import pytest

from nucleotide_transformer_v3.gene_annotation import annotate_genome
from nucleotide_transformer_v3.gene_annotation.utils import (
    PartIdentity,
    build_part_identity,
    build_provenance_lines,
    build_settings_digest,
    describe_part_mismatch,
    read_part,
    write_record_part,
)

# One complete, realistic run description. Every test below starts from this
# and changes exactly one thing, so a digest difference has one cause.
_BASE_RUN: dict[str, object] = {
    "checkpoint_path": "/weights/NTv3_annot_mammals",
    "model_name": "ntv3_annot_mammals",
    "sequence_length": 131_072,
    "keep_target_center_fraction": 0.75,
    "decoder": "hmm_constrained",
    "apply_biological_filters": True,
    "min_cds_length": 201,
    "use_bf16_compute": False,
    "species": "homo_sapiens",
    "mask_species_condition": False,
    "species_conditioning_consumed": True,
    "backend": "gpu",
}

_SEQUENCE = "ACGT" * 256


def _provenance(**overrides: object) -> list[str]:
    """The provenance block of a run described by ``_BASE_RUN`` plus overrides."""
    settings = {**_BASE_RUN, **overrides}
    return build_provenance_lines(**settings)  # type: ignore[arg-type]


def _digest(
    provenance_lines: list[str] | None = None,
    *,
    species: str = "homo_sapiens",
    mask_species_condition: bool = False,
) -> str:
    """The settings digest of a run, defaulting to the baseline everywhere."""
    lines = _provenance() if provenance_lines is None else provenance_lines
    return build_settings_digest(
        lines,
        species=species,
        mask_species_condition=mask_species_condition,
    )


# ---------------------------------------------------------------------------
# The digest covers the provenance block
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("index", range(len(_provenance())))
def test_every_provenance_line_participates_in_the_digest(index: int) -> None:
    """Changing any one provenance line changes the digest, including future lines."""
    lines = _provenance()
    mutated = list(lines)
    mutated[index] = f"{mutated[index]}-changed"
    assert _digest(mutated) != _digest(lines)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("checkpoint_path", "/weights/NTv3_annot_plants"),
        ("model_name", "ntv3_annot_plants"),
        ("sequence_length", 65_536),
        ("keep_target_center_fraction", 1.0),
        ("decoder", "hmm"),
        ("apply_biological_filters", False),
        ("min_cds_length", 300),
        ("use_bf16_compute", True),
        ("backend", "cpu"),
        ("species_conditioning_consumed", False),
    ],
)
def test_a_run_that_predicts_differently_gets_a_different_digest(
    field: str, value: object
) -> None:
    """Each output-changing setting reaches the digest through provenance."""
    assert _digest(_provenance(**{field: value})) != _digest()


@pytest.mark.parametrize("species", ["mus_musculus", ""])
def test_the_species_condition_reaches_the_digest(species: str) -> None:
    """The species condition is digested, as it is absent from the provenance."""
    assert _digest(species=species) != _digest()


def test_masking_the_species_condition_reaches_the_digest() -> None:
    """Dropping the condition is a different model input, so a different run."""
    assert _digest(mask_species_condition=True) != _digest()


# ---------------------------------------------------------------------------
# The digest is stable
# ---------------------------------------------------------------------------


def test_the_digest_is_stable_for_identical_settings() -> None:
    """Two separately built provenance blocks of one job give the same digest."""
    assert _digest(_provenance()) == _digest(_provenance())


def test_the_digest_is_hex_of_a_stable_width() -> None:
    """It is written into a part header, so it has to be one plain token."""
    digest = _digest()
    assert len(digest) == 32
    assert all(character in "0123456789abcdef" for character in digest)


# ---------------------------------------------------------------------------
# What a mismatch tells the user
# ---------------------------------------------------------------------------


def test_a_part_this_run_could_have_written_reports_no_mismatch() -> None:
    """The empty string is what lets a valid part be reused."""
    identity = build_part_identity(_SEQUENCE, _digest())
    assert describe_part_mismatch(identity, identity) == ""


def test_a_part_from_a_different_checkpoint_reports_the_settings_change() -> None:
    """The message names the difference and both digests."""
    stored_digest = _digest(_provenance(checkpoint_path="/weights/other"))
    current = build_part_identity(_SEQUENCE, _digest())
    stored = build_part_identity(_SEQUENCE, stored_digest)

    message = describe_part_mismatch(stored, current)

    assert "different settings" in message
    assert stored_digest in message
    assert current.settings_digest in message


def test_a_part_with_no_identity_header_is_reported_as_unknown() -> None:
    """A part written before parts carried an identity cannot vouch for itself."""
    current = build_part_identity(_SEQUENCE, _digest())
    assert "no identity header" in describe_part_mismatch(None, current)


def test_a_resized_record_reports_both_lengths() -> None:
    """Naming both is what tells the user which FASTA is the odd one out."""
    stored = build_part_identity(_SEQUENCE, _digest())
    current = build_part_identity(_SEQUENCE + "ACGT" * 64, _digest())

    message = describe_part_mismatch(stored, current)

    assert f"{stored.sequence_length:,}" in message
    assert f"{current.sequence_length:,}" in message


def test_an_edited_record_of_the_same_length_reports_the_sequence() -> None:
    """A substitution keeps the length, so the fingerprint is what catches it."""
    edited = "T" + _SEQUENCE[1:]
    stored = build_part_identity(_SEQUENCE, _digest())
    current = build_part_identity(edited, _digest())

    message = describe_part_mismatch(stored, current)

    assert "sequence under that record id has changed" in message
    assert stored.sequence_digest in message


def test_a_malformed_header_is_reported_rather_than_trusted() -> None:
    """An unparseable header is not an identity, so the part is not reused."""
    stored = PartIdentity.parse("#ntv3-part length=notanumber sequence=a settings=b")
    assert stored is None
    current = build_part_identity(_SEQUENCE, _digest())
    assert describe_part_mismatch(stored, current) != ""


# ---------------------------------------------------------------------------
# The staged part carries all of it
# ---------------------------------------------------------------------------


def test_a_staged_part_is_refused_when_the_checkpoint_changed(
    parts_dir: Path,
) -> None:
    """A part written under one checkpoint and read back is refused by another."""
    first_run = _digest(_provenance(checkpoint_path="/weights/first"))
    second_run = _digest(_provenance(checkpoint_path="/weights/second"))
    write_record_part(
        parts_dir,
        "chr1",
        ["chr1\tntv3\tgene\t1\t99\t.\t+\t.\tID=gene_1_p"],
        build_part_identity(_SEQUENCE, first_run),
    )

    stored, lines = read_part(parts_dir, "chr1")

    assert lines == ["chr1\tntv3\tgene\t1\t99\t.\t+\t.\tID=gene_1_p"]
    assert (
        describe_part_mismatch(stored, build_part_identity(_SEQUENCE, second_run))
        != ""
    )


def test_a_staged_part_is_reused_when_nothing_changed(parts_dir: Path) -> None:
    """A part read back under unchanged settings is accepted."""
    digest = _digest()
    write_record_part(
        parts_dir,
        "chr1",
        ["chr1\tntv3\tgene\t1\t99\t.\t+\t.\tID=gene_1_p"],
        build_part_identity(_SEQUENCE, digest),
    )

    stored, _ = read_part(parts_dir, "chr1")

    assert describe_part_mismatch(stored, build_part_identity(_SEQUENCE, digest)) == ""


# ---------------------------------------------------------------------------
# The resume decision in run()
# ---------------------------------------------------------------------------


def _count_annotated_records(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record the id of every record ``run()`` annotates rather than resumes."""
    annotated: list[str] = []
    real = annotate_genome.annotate_record

    def _spy(*args: object, **kwargs: object) -> object:
        annotated.append(str(kwargs["record_id"]))
        return real(*args, **kwargs)

    monkeypatch.setattr(annotate_genome, "annotate_record", _spy)
    return annotated


@pytest.mark.parametrize(
    ("resumed_overrides", "expect_reuse"),
    [
        ({}, True),
        ({"decoder": "hmm"}, False),
        ({"min_cds_length": 300}, False),
        # The stub reads no species, so the flag cannot change its output.
        ({"species": "mus_musculus"}, True),
    ],
)
def test_run_resumes_only_parts_written_under_the_same_settings(
    resumed_overrides: dict[str, object],
    expect_reuse: bool,
    synthetic_genome: list[str],
    stub_model_load: None,
    production_config,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A resumed ``run()`` reuses a part iff its settings digest still matches."""
    output = tmp_path / "out.gff3"
    first = production_config(fasta_paths=synthetic_genome, output_path=str(output))
    annotated = _count_annotated_records(monkeypatch)
    annotate_genome.run(first)
    first_output = output.read_text()
    records = list(annotated)
    assert records, "the first run annotated nothing"

    annotated.clear()
    annotate_genome.run(
        production_config(
            fasta_paths=synthetic_genome,
            output_path=str(output),
            resume=True,
            **resumed_overrides,
        )
    )

    assert annotated == ([] if expect_reuse else records)
    if not resumed_overrides:
        assert output.read_text() == first_output
