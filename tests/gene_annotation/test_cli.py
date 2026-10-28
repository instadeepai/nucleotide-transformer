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

"""Tests for the ``ntv3-annotate`` command line, from flags to the composed config."""

import dataclasses
from pathlib import Path
from typing import Any

import pytest

from nucleotide_transformer_v3.gene_annotation.cli import (
    _coerce,
    apply_overrides,
    build_parser,
    compose_config,
    main,
)
from nucleotide_transformer_v3.gene_annotation.config import ProductionAnnotationConfig
from nucleotide_transformer_v3.gene_annotation.model_registry import (
    FIXED_KEEP_TARGET_CENTER_FRACTION,
    list_model_config_names,
)

GOLDEN_DIR = Path(__file__).parent / "data"

_FIELDS: dict[str, Any] = {
    field.name: field for field in dataclasses.fields(ProductionAnnotationConfig)
}


def _field(name: str) -> Any:
    """The dataclass field ``_coerce`` types a value against."""
    return _FIELDS[name]


def _feature_lines(gff3: str) -> list[str]:
    """Non-comment, non-blank lines, i.e. everything provenance-independent."""
    return [line for line in gff3.splitlines() if line and not line.startswith("#")]


def _window_length(golden_geometry: tuple[int, int]) -> int:
    """The stub head's window length, from the golden geometry fixture."""
    crop, stride = golden_geometry
    return stride + 2 * crop


# ---------------------------------------------------------------------------
# build_parser
# ---------------------------------------------------------------------------


def test_the_parser_accepts_every_documented_flag() -> None:
    """One command line using all of them, checked field by field."""
    args = build_parser().parse_args(
        [
            "--genome",
            "a.fa",
            "contigs/",
            "--model_cfg",
            "mammals",
            "--out",
            "out.gff3",
            "--species",
            "mouse",
            "--batch_size",
            "4",
            "--codingseq",
            "cds.fa",
            "--protseq",
            "prot.fa",
            "--resume",
            "--seq_len",
            "1048576",
            "--override",
            "decoder=hmm",
            "-O",
            "min_cds_length=0",
        ]
    )

    assert args.genome == ["a.fa", "contigs/"]
    assert args.model_cfg == "mammals"
    assert args.out == "out.gff3"
    assert args.species == "mouse"
    assert args.batch_size == 4
    assert args.codingseq == "cds.fa"
    assert args.protseq == "prot.fa"
    assert args.resume is True
    assert args.seq_len == 1048576
    assert args.override == ["decoder=hmm", "min_cds_length=0"]


def test_an_unpassed_flag_parses_as_none_rather_than_a_default() -> None:
    """``None`` is how ``_apply_flags`` tells "unset" from "set to the default"."""
    args = build_parser().parse_args(["--genome", "a.fa"])

    assert args.resume is None
    assert args.batch_size is None
    assert args.seq_len is None
    assert args.model_cfg is None
    assert args.override == []


@pytest.mark.parametrize(
    "argv",
    [
        ["--genomes", "a.fa"],
        ["--model-cfg", "mammals"],
        ["--threads", "8"],
        ["--genome", "a.fa", "--out", "o.gff3", "--jit"],
    ],
)
def test_the_parser_rejects_an_unknown_flag(argv: list[str]) -> None:
    """A misspelled flag must not be swallowed as a positional or ignored."""
    with pytest.raises(SystemExit) as exit_info:
        build_parser().parse_args(argv)

    assert exit_info.value.code == 2


@pytest.mark.parametrize("flag", ["--batch_size", "--seq_len"])
def test_the_parser_rejects_a_non_integer_where_it_wants_an_integer(
    flag: str,
) -> None:
    """Typed at the parser, so it fails before anything is loaded."""
    with pytest.raises(SystemExit) as exit_info:
        build_parser().parse_args([flag, "lots"])

    assert exit_info.value.code == 2


# ---------------------------------------------------------------------------
# _coerce
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("field_name", "raw", "expected"),
    [
        ("batch_size", "4", 4),
        ("sequence_length", "1048576", 1048576),
        ("min_cds_length", "0", 0),
        ("keep_target_center_fraction", "0.75", 0.75),
        ("use_bf16_compute", "true", True),
        ("use_bf16_compute", "True", True),
        ("use_bf16_compute", "1", True),
        ("use_bf16_compute", "false", False),
        ("use_bf16_compute", "FALSE", False),
        ("use_bf16_compute", "0", False),
        ("decoder", "hmm", "hmm"),
        ("weights_dir", "s3://bucket/key", "s3://bucket/key"),
        ("species", "mouse", "mouse"),
    ],
)
def test_coerce_converts_a_value_to_the_fields_declared_type(
    field_name: str, raw: str, expected: object
) -> None:
    """Every type the config actually declares, read off the field itself."""
    coerced = _coerce(_field(field_name), raw)

    assert coerced == expected
    assert isinstance(coerced, type(expected))


@pytest.mark.parametrize("raw", ["null", "none", "None", "~"])
def test_coerce_spells_unset_the_way_hydra_did(raw: str) -> None:
    """Every Optional field defaults to None, so there has to be a way back."""
    assert _coerce(_field("species"), raw) is None
    assert _coerce(_field("keep_target_center_fraction"), raw) is None


def test_coerce_leaves_null_alone_on_a_field_that_cannot_be_none() -> None:
    """``decoder=null`` is a bad decoder name, not an unset decoder."""
    assert _coerce(_field("decoder"), "null") == "null"


def test_coerce_wraps_a_list_field_without_splitting_on_commas() -> None:
    """``fasta_paths=/my,dir/a.fa`` used to become two paths, neither real."""
    assert _coerce(_field("fasta_paths"), "/my,dir/a.fa") == ["/my,dir/a.fa"]


def test_coerce_types_the_items_of_a_list_field() -> None:
    """A non-``str`` list field gets typed items, not a scalar."""
    hypothetical = dataclasses.make_dataclass(
        "_HypotheticalConfig",
        [("counts", list[int], dataclasses.field(default_factory=list))],
    )

    assert _coerce(dataclasses.fields(hypothetical)[0], "3") == [3]


def test_coerce_resolves_an_annotation_left_as_a_string() -> None:
    """A string annotation such as ``"int"`` is resolved, not coerced to ``str``."""
    stringly = dataclasses.make_dataclass("_Stringly", [("batch_size", "int")])
    field = dataclasses.fields(stringly)[0]
    assert field.type == "int"

    assert _coerce(field, "4") == 4


@pytest.mark.parametrize(
    ("field_name", "raw", "message"),
    [
        ("batch_size", "many", "batch_size expects an int, got 'many'"),
        ("sequence_length", "1e6", "sequence_length expects an int"),
        (
            "keep_target_center_fraction",
            "most",
            "keep_target_center_fraction expects a float, got 'most'",
        ),
        ("use_bf16_compute", "yes", "use_bf16_compute expects true/false"),
        ("use_bf16_compute", "", "use_bf16_compute expects true/false"),
    ],
)
def test_coerce_says_what_it_expected_on_a_bad_value(
    field_name: str, raw: str, message: str
) -> None:
    """The message names the field and the type, not just "invalid"."""
    with pytest.raises(ValueError, match=message):
        _coerce(_field(field_name), raw)


# ---------------------------------------------------------------------------
# apply_overrides
# ---------------------------------------------------------------------------


def test_an_override_reaches_the_config() -> None:
    """The whole point of the flag, and the thing whose absence is silent."""
    config = apply_overrides(
        ProductionAnnotationConfig(),
        ["batch_size=8", "decoder=hmm", "use_bf16_compute=false"],
    )

    assert config.batch_size == 8
    assert config.decoder == "hmm"
    assert config.use_bf16_compute is False


def test_an_override_does_not_mutate_the_config_it_was_given() -> None:
    """``dataclasses.replace`` copies, and callers rely on that."""
    original = ProductionAnnotationConfig()

    updated = apply_overrides(original, ["batch_size=8"])

    assert original.batch_size == ProductionAnnotationConfig().batch_size
    assert updated.batch_size == 8


def test_overrides_left_unmentioned_keep_their_defaults() -> None:
    """An override sets one field, not the whole config."""
    config = apply_overrides(ProductionAnnotationConfig(), ["batch_size=8"])

    assert config.min_cds_length == ProductionAnnotationConfig().min_cds_length
    assert config.decoder == ProductionAnnotationConfig().decoder


def test_a_repeated_list_override_appends_rather_than_replaces() -> None:
    """Repeating the flag is the documented way to pass several paths."""
    config = apply_overrides(
        ProductionAnnotationConfig(),
        ["fasta_paths=a.fa", "fasta_paths=b.fa", "fasta_paths=c.fa"],
    )

    assert config.fasta_paths == ["a.fa", "b.fa", "c.fa"]


def test_a_repeated_scalar_override_takes_the_last_one() -> None:
    """A scalar cannot accumulate, so the later mention has to win."""
    config = apply_overrides(
        ProductionAnnotationConfig(), ["batch_size=2", "batch_size=16"]
    )

    assert config.batch_size == 16


def test_whitespace_around_a_key_is_not_a_different_key() -> None:
    """``--override ' batch_size=8'`` is a quoting accident, not a typo."""
    config = apply_overrides(ProductionAnnotationConfig(), ["  batch_size=8"])

    assert config.batch_size == 8


@pytest.mark.parametrize("key", ["batchsize", "batch-size", "decoder_type", ""])
def test_an_unknown_override_key_is_an_error_not_a_dropped_setting(
    key: str,
) -> None:
    """Accepting it would run the whole genome with the setting ignored."""
    with pytest.raises(ValueError, match="Unknown config field"):
        apply_overrides(ProductionAnnotationConfig(), [f"{key}=8"])


def test_the_unknown_key_error_lists_what_would_have_worked() -> None:
    """A user who misspells a field needs the spelling, not a refusal."""
    with pytest.raises(ValueError, match="batch_size") as error_info:
        apply_overrides(ProductionAnnotationConfig(), ["batchsize=8"])

    assert "'batchsize'" in str(error_info.value)


@pytest.mark.parametrize("override", ["batch_size", "", "decoder hmm"])
def test_an_override_without_an_equals_sign_is_an_error(override: str) -> None:
    """Otherwise it parses as a key with an empty value, or is ignored."""
    with pytest.raises(ValueError, match="is not of the form key=value"):
        apply_overrides(ProductionAnnotationConfig(), [override])


def test_a_bad_value_is_rejected_by_the_field_it_was_meant_for() -> None:
    """The typing in ``_coerce`` is reached through the public function."""
    with pytest.raises(ValueError, match="batch_size expects an int"):
        apply_overrides(ProductionAnnotationConfig(), ["batch_size=some"])


# ---------------------------------------------------------------------------
# compose_config
# ---------------------------------------------------------------------------


def _compose(argv: list[str]) -> ProductionAnnotationConfig:
    """Parse a command line and compose the config it describes."""
    return compose_config(build_parser().parse_args(argv))


def test_compose_config_maps_every_flag_onto_its_field() -> None:
    """The flag names are Tiberius', the field names are ours; this is the map."""
    config = _compose(
        [
            "--genome",
            "a.fa",
            "contigs/",
            "--model_cfg",
            "mammals",
            "--out",
            "out.gff3",
            "--species",
            "mouse",
            "--batch_size",
            "4",
            "--codingseq",
            "cds.fa",
            "--protseq",
            "prot.fa",
            "--resume",
            "--seq_len",
            "1048576",
        ]
    )

    assert config.fasta_paths == ["a.fa", "contigs/"]
    assert config.model_cfg == "mammals"
    assert config.output_path == "out.gff3"
    assert config.species == "mouse"
    assert config.batch_size == 4
    assert config.codingseq_path == "cds.fa"
    assert config.protseq_path == "prot.fa"
    assert config.resume is True
    assert config.sequence_length == 1048576


def test_a_flag_left_off_does_not_clobber_the_dataclass_default() -> None:
    """Every unset flag is None, and writing None would erase the default."""
    defaults = ProductionAnnotationConfig()

    config = _compose(["--genome", "a.fa", "--out", "out.gff3"])

    assert config.batch_size == defaults.batch_size
    assert config.sequence_length == defaults.sequence_length
    assert config.resume == defaults.resume
    assert config.species == defaults.species
    assert config.codingseq_path == defaults.codingseq_path
    assert config.use_bf16_compute == defaults.use_bf16_compute


def test_compose_config_carries_an_override_for_a_field_with_no_flag() -> None:
    """``--override`` is the only way to these, so it has to reach them."""
    config = _compose(
        [
            "--genome",
            "a.fa",
            "--out",
            "out.gff3",
            "--override",
            "decoder=hmm",
            "--override",
            "min_cds_length=0",
            "--override",
            "apply_biological_filters=false",
        ]
    )

    assert config.decoder == "hmm"
    assert config.min_cds_length == 0
    assert config.apply_biological_filters is False


def test_a_named_flag_wins_over_an_override_of_the_same_field() -> None:
    """The order the module docstring promises: overrides first, flags second."""
    config = _compose(
        [
            "--genome",
            "a.fa",
            "--out",
            "out.gff3",
            "--override",
            "batch_size=2",
            "--batch_size",
            "9",
        ]
    )

    assert config.batch_size == 9


def test_an_override_still_wins_where_no_flag_was_passed() -> None:
    """The precedence above must not make ``--override`` conditional on luck."""
    config = _compose(
        ["--genome", "a.fa", "--out", "out.gff3", "--override", "batch_size=2"]
    )

    assert config.batch_size == 2


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def test_help_exits_zero() -> None:
    """``--help`` is a success, not a usage error."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--help"])

    assert exit_info.value.code == 0


def test_list_cfg_exits_zero_and_names_the_registered_models(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The flag exists to answer "what may I pass to --model_cfg"."""
    assert main(["--list_cfg"]) == 0

    printed = capsys.readouterr().out
    for name in list_model_config_names():
        assert name in printed


def test_show_cfg_exits_zero_and_reports_the_crop_it_will_run_at(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The one question whose wrong answer is silent, answered without weights."""
    assert main(["--show_cfg", "mammals"]) == 0

    printed = capsys.readouterr().out
    assert "mammals" in printed
    assert str(FIXED_KEEP_TARGET_CENTER_FRACTION) in printed


def test_an_unknown_show_cfg_name_is_a_usage_error() -> None:
    """Reported as a usage error rather than an unhandled ValueError."""
    with pytest.raises(SystemExit) as exit_info:
        main(["--show_cfg", "tardigrades"])

    assert exit_info.value.code == 2


@pytest.mark.parametrize(
    ("argv", "expected_in_stderr"),
    [
        ([], "--genome is required"),
        (["--out", "out.gff3"], "--genome is required"),
        (["--genome", "a.fa"], "--out is required"),
        (["--genome", "a.fa", "--out", "out.gff3"], "--model_cfg is required"),
        (
            [
                "--genome",
                "a.fa",
                "--out",
                "out.gff3",
                "--model_cfg",
                "mammals",
                "--override",
                "batch_size=plenty",
            ],
            "batch_size expects an int",
        ),
        (
            [
                "--genome",
                "a.fa",
                "--out",
                "out.gff3",
                "--model_cfg",
                "mammals",
                "--override",
                "batchsize=2",
            ],
            "Unknown config field",
        ),
    ],
)
def test_main_refuses_an_incomplete_or_wrong_command_line(
    argv: list[str],
    expected_in_stderr: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Each refusal happens before any weight is fetched, and says why."""
    with pytest.raises(SystemExit) as exit_info:
        main(argv)

    assert exit_info.value.code == 2
    assert expected_in_stderr in capsys.readouterr().err


def test_main_annotates_the_synthetic_genome_end_to_end(
    synthetic_genome: list[str],
    stub_model_load: None,
    golden_geometry: tuple[int, int],
    tmp_path: Path,
) -> None:
    """The default CLI run reproduces the ``hmm_constrained_filtered`` golden.

    No ``--model_cfg``, so this also takes the ``weights_dir`` escape hatch.
    """
    output = tmp_path / "cli.gff3"

    status = main(
        [
            "--genome",
            *synthetic_genome,
            "--out",
            str(output),
            "--seq_len",
            str(_window_length(golden_geometry)),
            "--batch_size",
            "2",
            "--override",
            "weights_dir=/stub/checkpoint",
            "--override",
            "use_bf16_compute=false",
        ]
    )

    assert status == 0
    golden = GOLDEN_DIR / "golden_hmm_constrained_filtered.gff3"
    assert _feature_lines(output.read_text()) == _feature_lines(golden.read_text())
