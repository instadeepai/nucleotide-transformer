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

"""The package's dependency surface is declared, not accidental.

A dependency surface is a decision a reviewer cannot enforce by
reading 7,000 lines, so it is enforced here: adding a dependency is allowed,
adding one silently is not.
"""

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest

import nucleotide_transformer_v3 as package

# Everything the package may import beyond the standard library.
ALLOWED_THIRD_PARTY: frozenset[str] = frozenset(
    {"jax", "jaxlib", "flax", "numpy", "yaml", "typing_extensions"}
)

# Imports that would undo the port. orbax/tensorstore/etils/s3fs come back with
# the Orbax checkpoint format; hydra/omegaconf with the research config layer;
# grain/mlflow/pandas/pyranges/pysam with the research data layer.
FORBIDDEN: frozenset[str] = frozenset(
    {
        "scripts",
        "orbax",
        "tensorstore",
        "etils",
        "s3fs",
        "gcsfs",
        "fsspec",
        "hydra",
        "omegaconf",
        "grain",
        "mlflow",
        "pandas",
        "pyranges",
        "pysam",
        "tensorflow",
        "ray",
    }
)

# The ported code, not the whole package. Since gene_annotation/ lives inside
# nucleotide_transformer_v3, scanning the package root would police modules
# this port does not own, and would fail the moment anyone adds a dependency
# to one of them. What the port owes is that ITS OWN sources stay within the
# declared surface; what it inherits is covered by the runtime check below.
_ROOT = Path(package.__file__).parent / "gene_annotation"


def _sources() -> list[Path]:
    return sorted(p for p in _ROOT.rglob("*.py") if "__pycache__" not in p.parts)


def _imported_roots(path: Path) -> set[str]:
    """Top-level module names one source file imports, absolute imports only."""
    roots: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
        if isinstance(node, ast.Import):
            roots |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            roots.add(node.module.split(".")[0])
    return roots


def test_sources_exist() -> None:
    """Guard the guard: an empty file list would make everything below vacuous."""
    assert len(_sources()) >= 10


def test_no_forbidden_imports() -> None:
    """Nothing reaches back into the research codebase or its heavy dependencies."""
    offenders: dict[str, set[str]] = {}
    for source in _sources():
        outside = _imported_roots(source) & FORBIDDEN
        if outside:
            offenders[str(source.relative_to(_ROOT))] = outside
    assert not offenders, (
        f"These modules import dependencies the port exists to shed: {offenders}. "
        "The package must install and run without the research codebase."
    )


def test_third_party_surface_is_declared() -> None:
    """No undeclared third-party dependency creeps in."""
    offenders: dict[str, set[str]] = {}
    for source in _sources():
        outside = {
            root
            for root in _imported_roots(source)
            if root not in sys.stdlib_module_names
            and root != "nucleotide_transformer_v3"
            and root not in ALLOWED_THIRD_PARTY
        }
        if outside:
            offenders[str(source.relative_to(_ROOT))] = outside
    assert not offenders, (
        f"Undeclared third-party imports: {offenders}. If intended, add the "
        "package to ALLOWED_THIRD_PARTY here and to setup.py's install_requires."
    )


@pytest.mark.parametrize("forbidden", sorted(FORBIDDEN - {"fsspec", "etils"}))
def test_forbidden_modules_are_not_imported_at_runtime(forbidden: str) -> None:
    """Importing the pipeline must not pull a forbidden module in transitively.

    Run in a fresh interpreter, because ``sys.modules`` in this one is ambient:
    anything else in the session can populate it. Under pytest's plugin
    autoload, hydra-core's plugin imports ``hydra`` and ``omegaconf`` before
    any test runs, and this failed while blaming the package. It passed only
    because CI installs nothing but pytest -- one extra plugin, or one
    co-resident test that imports pandas, and it went red for the wrong reason.

    ``etils`` and ``fsspec`` are exempt: flax imports them itself, so their
    presence says nothing about this package. The static checks above still
    forbid importing them here directly.
    """
    probe = (
        "import sys;"
        " import nucleotide_transformer_v3.gene_annotation.annotate_genome;"
        " import nucleotide_transformer_v3.gene_annotation.cli;"
        " import nucleotide_transformer_v3.gene_annotation.check_equivalence;"
        f" loaded = {{n.split('.')[0] for n in sys.modules}};"
        f" assert {forbidden!r} not in loaded, {forbidden!r}"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        env={**os.environ, "JAX_PLATFORMS": "cpu"},
    )
    assert result.returncode == 0, (
        f"importing the pipeline pulled in {forbidden!r}:\n{result.stderr[-2000:]}"
    )
