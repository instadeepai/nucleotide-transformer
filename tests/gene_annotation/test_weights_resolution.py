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

"""Where a registered clade's weights come from, and what overrides that.

The shipped promise is that `pip install` plus one command annotates a genome,
so the entries must resolve to something loadable with nothing staged. None of
this downloads: it checks what the registry resolves to, not what the Hub
returns.
"""

import logging
from pathlib import Path

import pytest

from nucleotide_transformer_v3.gene_annotation.annotate_genome import (
    resolve_sequence_outputs,
)
from nucleotide_transformer_v3.gene_annotation.config import ProductionAnnotationConfig
from nucleotide_transformer_v3.gene_annotation.model_registry import (
    _WEIGHTS_DIR_ENV,
    apply_model_config,
    list_model_config_names,
    load_model_config,
    resolve_weights_uri,
)
from nucleotide_transformer_v3.pretrained import RELEASED_MODELS


@pytest.fixture(autouse=True)
def _no_override(monkeypatch: pytest.MonkeyPatch) -> None:
    """Start every test from the shipped default, whatever the shell has set."""
    monkeypatch.delenv(_WEIGHTS_DIR_ENV, raising=False)


@pytest.mark.parametrize("name", sorted(list_model_config_names()))
def test_every_entry_names_a_released_model(name: str) -> None:
    """A clade resolves to something the loader can fetch on its own.

    An entry pointing at a path that happens to exist on one machine would
    pass every other test here and fail for everyone else.
    """
    assert load_model_config(name).weights_uri in RELEASED_MODELS


def test_the_insect_entry_uses_the_published_spelling() -> None:
    """The clade is `insects`, like its repository `NTv3_annot_insects`.

    The exported weights keep their training-time name, which says `insecta`:
    the load-time identity check compares against it, so it must not follow.
    """
    config = load_model_config("insects")
    assert config.weights_uri == "NTv3_annot_insects"
    assert config.model_name == "new_ntv3_700M_annotate_insecta_1mb"


@pytest.mark.parametrize("name", sorted(list_model_config_names()))
def test_the_default_needs_nothing_staged(name: str) -> None:
    """With no override, the resolved location is the released name itself."""
    config = apply_model_config(
        ProductionAnnotationConfig(model_cfg=name, output_path="out.gff3")
    )
    assert config.weights_dir == load_model_config(name).weights_uri


@pytest.mark.parametrize("name", sorted(list_model_config_names()))
def test_the_environment_override_wins(
    name: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Someone who sets the variable means it, so it takes precedence."""
    monkeypatch.setenv(_WEIGHTS_DIR_ENV, str(tmp_path))
    config = apply_model_config(
        ProductionAnnotationConfig(model_cfg=name, output_path="out.gff3")
    )
    assert config.weights_dir == str(tmp_path / name)


def test_an_empty_override_is_not_an_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`VAR=""` once yielded `/mammals`, an absolute path at the root."""
    monkeypatch.setenv(_WEIGHTS_DIR_ENV, "")
    config = apply_model_config(
        ProductionAnnotationConfig(model_cfg="mammals", output_path="out.gff3")
    )
    assert config.weights_dir == "NTv3_annot_mammals"


def test_the_override_expands_a_home_relative_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`~` has to expand, or the run creates a literal `~` directory."""
    monkeypatch.setenv(_WEIGHTS_DIR_ENV, "~/staged")
    resolved = Path(
        apply_model_config(
            ProductionAnnotationConfig(model_cfg="plants", output_path="out.gff3")
        ).weights_dir
    )
    # `.absolute()` also yields an absolute path named "plants", with a literal
    # `~` directory inside it, so anchor on the real home instead.
    assert resolved == Path.home() / "staged" / "plants"
    assert "~" not in str(resolved)


def test_a_hand_written_entry_still_reports_an_unset_variable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The expansion path survives for entries that are not released models."""
    monkeypatch.delenv("NTV3_TEST_MISSING", raising=False)
    with pytest.raises(ValueError, match="NTV3_TEST_MISSING is not set"):
        resolve_weights_uri("mammals", "$NTV3_TEST_MISSING/x", Path("entry.yaml"))


def test_a_cpu_fallback_on_a_gpu_machine_is_reported(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The failure this catches is silence, so the warning is the behaviour.

    A plain ``pip install`` gives a CUDA-less jaxlib. On a GPU host JAX then
    falls back without erroring and a genome takes ~7.5 days instead of ~45
    minutes.
    """
    from nucleotide_transformer_v3.gene_annotation import annotate_genome

    class _CpuDevice:
        platform = "cpu"

    monkeypatch.setattr(annotate_genome.jax, "devices", lambda *_: [_CpuDevice()])
    monkeypatch.setattr(annotate_genome.Path, "exists", lambda _: True)

    with caplog.at_level(logging.WARNING):
        annotate_genome.resolve_devices("")

    assert "Annotating on the CPU" in caplog.text
    assert "nucleotide_transformer[gpu]" in caplog.text


def test_an_explicit_cpu_backend_is_not_second_guessed(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """``backend=cpu`` is how the tests run; warning there would be noise."""
    from nucleotide_transformer_v3.gene_annotation import annotate_genome

    class _CpuDevice:
        platform = "cpu"

    monkeypatch.setattr(annotate_genome.jax, "devices", lambda *_: [_CpuDevice()])
    monkeypatch.setattr(annotate_genome.Path, "exists", lambda _: True)

    with caplog.at_level(logging.WARNING):
        annotate_genome.resolve_devices("cpu")

    assert "Annotating on the CPU" not in caplog.text


@pytest.mark.parametrize("field", ["codingseq_path", "protseq_path"])
def test_a_sequence_output_may_not_share_the_annotation_path(field: str) -> None:
    """A FASTA output pointed at ``output_path`` would overwrite the annotation.

    The sequence merge runs after the annotation merge, so the collision is
    silent: exit code 0, both writes logged, and the GFF3 replaced by sequences.
    """
    config = ProductionAnnotationConfig(
        model_cfg="mammals", output_path="out.gff3", **{field: "out.gff3"}
    )
    with pytest.raises(ValueError, match="output_path"):
        resolve_sequence_outputs(config)


def _sequence_outputs(tmp_path: Path, fasta: str, **paths: str) -> list[object]:
    """Resolve the FASTA outputs of a run reading ``fasta`` and writing ``paths``."""
    config = ProductionAnnotationConfig(
        model_cfg="mammals",
        fasta_paths=[fasta],
        output_path=str(tmp_path / "out.gff3"),
        **paths,
    )
    return resolve_sequence_outputs(config)


@pytest.mark.parametrize(
    ("fasta", "paths", "message"),
    [
        ("genome.fa", {"protseq_path": "{tmp}/./out.gff3"}, "output_path"),
        (
            "genome.fa",
            {"codingseq_path": "{tmp}/s.fa", "protseq_path": "{tmp}//s.fa"},
            "codingseq",
        ),
        ("genome.fa", {"codingseq_path": "{tmp}/genome.fa"}, "input FASTA"),
        ("genome.fa", {"protseq_path": "{tmp}/link.fa"}, "input FASTA"),
        ("genome.fa", {"protseq_path": "{tmp}/hard.fa"}, "input FASTA"),
        ("contigs", {"codingseq_path": "{tmp}/contigs/cds.fa.gz"}, "input directory"),
    ],
)
def test_a_sequence_output_may_not_overwrite_an_output_or_an_input(
    tmp_path: Path, fasta: str, paths: dict[str, str], message: str
) -> None:
    """Collisions are found on resolved paths, inputs and input directories included."""
    (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
    (tmp_path / "link.fa").symlink_to(tmp_path / "genome.fa")
    (tmp_path / "hard.fa").hardlink_to(tmp_path / "genome.fa")
    (tmp_path / "contigs").mkdir()
    (tmp_path / "contigs" / "chr1.fa").write_text(">chr1\nACGT\n")
    resolved = {key: value.format(tmp=tmp_path) for key, value in paths.items()}

    with pytest.raises(ValueError, match=message):
        _sequence_outputs(tmp_path, str(tmp_path / fasta), **resolved)


def test_a_sequence_output_beside_the_inputs_is_allowed(tmp_path: Path) -> None:
    """Only a FASTA name inside an input directory is read back as genome."""
    (tmp_path / "contigs").mkdir()
    (tmp_path / "contigs" / "chr1.fa").write_text(">chr1\nACGT\n")

    outputs = _sequence_outputs(
        tmp_path,
        str(tmp_path / "contigs"),
        codingseq_path=str(tmp_path / "cds.fa"),
        protseq_path=str(tmp_path / "contigs" / "proteins.txt"),
    )

    assert len(outputs) == 2
