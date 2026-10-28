"""Generate deterministic local reproducibility artifacts for one NTv3 run."""

import hashlib
import json
import mimetypes
from pathlib import Path
from typing import Any

from ntv3_tools.errors import ArtifactError

RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
_GENERATED_PATHS = (
    "analysis.ipynb",
    "code/run_analysis.py",
    "browser.html",
    "ntv3_report.html",
    "report.json",
    "ro-crate-metadata.json",
)


def refuse_generated_overwrite(root: Path) -> None:
    """Refuse a report run when any generated output already exists."""
    existing = [relative for relative in _GENERATED_PATHS if (root / relative).exists()]
    if existing:
        raise ArtifactError(
            f"Generated report output already exists: {', '.join(existing)}. "
            "Use a new results/<run> directory."
        )


def write_replay_script(root: Path) -> Path:
    """Write the editable replay entry point for an approved plan."""
    code_dir = root / "code"
    code_dir.mkdir(exist_ok=True)
    path = code_dir / "run_analysis.py"
    source = '''"""Replay this approved NTv3 analysis into a new results directory."""

from pathlib import Path

from ntv3_tools.replay import main


if __name__ == "__main__":
    main(source_run=Path(__file__).resolve().parents[1])
'''
    path.write_text(source)
    return path


def write_analysis_notebook(root: Path, gene: str, cell: str) -> Path:
    """Write an editable artifact-analysis notebook without model execution."""
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                f"# {gene} in {cell}\n",
                "\n",
                (
                    "This notebook reads existing run artifacts with relative paths. "
                    "It never runs the model implicitly.\n"
                ),
                "\n",
                (
                    "**Caveats.** Predictions are not measurements. Attribution marks "
                    "model sensitivity and does not establish causality. Variant "
                    "effects need experimental validation. Missing optional artifacts "
                    "mean that workflow step was not run.\n"
                ),
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import json\n",
                "from pathlib import Path\n",
                "\n",
                "import yaml\n",
                "\n",
                'RUN = Path(".")\n',
                'plan = yaml.safe_load((RUN / "plan.yaml").read_text())\n',
                'report = json.loads((RUN / "report.json").read_text())\n',
                'manifest_path = RUN / "run_manifest.json"\n',
                "manifest = (\n",
                "    json.loads(manifest_path.read_text())\n",
                "    if manifest_path.exists()\n",
                "    else {}\n",
                ")\n",
                '{"plan": plan, "report_keys": sorted(report), "manifest": manifest}\n',
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## Load numerical and optional evidence artifacts\n",
                "\n",
                (
                    "Edit the paths or selections below to inspect a subset. Arrays "
                    "stay in their saved coordinate systems.\n"
                ),
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import json\n",
                "from pathlib import Path\n",
                "\n",
                "import numpy as np\n",
                "import pandas as pd\n",
                "\n",
                'RUN = Path(".")\n',
                'tracks = np.load(RUN / "tracks/tracks.npy")\n',
                "saliency = {\n",
                "    path.stem: np.load(path)\n",
                '    for path in sorted((RUN / "interpret").glob("*.npy"))\n',
                '    if (RUN / "interpret").exists()\n',
                "}\n",
                'disrupt_dir = RUN / "disrupt"\n',
                "disrupt = {\n",
                "    path.stem: np.load(path)\n",
                '    for path in sorted(disrupt_dir.glob("*.npy"))\n',
                "}\n",
                'disrupt_table_path = disrupt_dir / "disrupt_top.tsv"\n',
                "disrupt_table = (\n",
                '    pd.read_csv(disrupt_table_path, sep="\\t")\n',
                "    if disrupt_table_path.exists()\n",
                "    else None\n",
                ")\n",
                'critic_path = RUN / "critic/audit_log.json"\n',
                "critic = (\n",
                "    json.loads(critic_path.read_text())\n",
                "    if critic_path.exists()\n",
                "    else []\n",
                ")\n",
                'literature_path = RUN / "literature/search_manifest.json"\n',
                "literature = (\n",
                "    json.loads(literature_path.read_text())\n",
                "    if literature_path.exists()\n",
                "    else []\n",
                ")\n",
                "{\n",
                '    "tracks_shape": tracks.shape,\n',
                (
                    '    "saliency_shapes": {key: value.shape for key, value in '
                    "saliency.items()},\n"
                ),
                '    "disrupt_arrays": sorted(disrupt),\n',
                '    "critic_claims": len(critic),\n',
                '    "literature_records": len(literature),\n',
                "}\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## Editable analysis\n",
                "\n",
                (
                    "The following cell plots saved predictions only. Change the track "
                    "index, normalization, or genomic x-axis as needed.\n"
                ),
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import matplotlib.pyplot as plt\n",
                "import numpy as np\n",
                "\n",
                "track_index = 0\n",
                "values = np.asarray(tracks[track_index])\n",
                "figure, axis = plt.subplots(figsize=(10, 2.5))\n",
                'axis.plot(values, linewidth=0.8, color="#245b78")\n',
                (
                    'axis.set(title=f"Saved prediction, track {track_index}", '
                    'xlabel="Output bin", ylabel="Model output")\n'
                ),
                "figure.tight_layout()\n",
            ],
        },
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.11"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = root / "analysis.ipynb"
    path.write_text(json.dumps(notebook, indent=2, sort_keys=True) + "\n")
    return path


def _encoding_format(path: Path) -> str:
    """Return a stable encoding format for one artifact."""
    formats = {
        ".bib": "application/x-bibtex",
        ".ipynb": "application/x-ipynb+json",
        ".npy": "application/x-npy",
        ".tsv": "text/tab-separated-values",
        ".yaml": "application/yaml",
        ".yml": "application/yaml",
    }
    if path.suffix.lower() in formats:
        return formats[path.suffix.lower()]
    guessed, _ = mimetypes.guess_type(path.name)
    return guessed or "application/octet-stream"


def _artifact_paths(root: Path) -> list[Path]:
    """List safe regular artifact files in deterministic relative order."""
    files: list[Path] = []
    resolved_root = root.resolve()
    for path in root.rglob("*"):
        if not path.is_file() or path.name == "ro-crate-metadata.json":
            continue
        relative = path.relative_to(root)
        if ".env" in relative.parts or "secrets" in relative.parts:
            continue
        resolved = path.resolve()
        if resolved_root not in resolved.parents:
            continue
        files.append(path)
    return sorted(files, key=lambda path: path.relative_to(root).as_posix())


def _provenance_properties(root: Path) -> list[dict[str, Any]]:
    """Read selected checkpoint and repository provenance."""
    manifest_path = root / "run_manifest.json"
    if not manifest_path.is_file():
        return []
    payload = json.loads(manifest_path.read_text())
    if not isinstance(payload, dict):
        raise ArtifactError("run_manifest.json must contain a JSON object.")
    properties: list[dict[str, Any]] = []
    for source_name, output_name in (
        ("model", "checkpoint"),
        ("commit", "commit"),
        ("worktree_dirty", "worktree_dirty"),
    ):
        if source_name in payload:
            properties.append(
                {
                    "@type": "PropertyValue",
                    "name": output_name,
                    "value": payload[source_name],
                }
            )
    return properties


def write_ro_crate(root: Path, main_entity: str = "plan.yaml") -> Path:
    """Write deterministic RO-Crate JSON-LD for existing safe run artifacts."""
    file_entities: list[dict[str, Any]] = []
    for path in _artifact_paths(root=root):
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        file_entities.append(
            {
                "@id": relative,
                "@type": "File",
                "contentSize": len(content),
                "encodingFormat": _encoding_format(path=path),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
        )
    root_entity: dict[str, Any] = {
        "@id": "./",
        "@type": "Dataset",
        "hasPart": [{"@id": row["@id"]} for row in file_entities],
        "mainEntity": {"@id": main_entity},
        "name": "NTv3 locus-analysis run",
    }
    properties = _provenance_properties(root=root)
    if properties:
        root_entity["additionalProperty"] = properties
    descriptor = {
        "@id": "ro-crate-metadata.json",
        "@type": "CreativeWork",
        "about": {"@id": "./"},
        "conformsTo": {"@id": RO_CRATE_CONTEXT},
    }
    payload = {
        "@context": RO_CRATE_CONTEXT,
        "@graph": [root_entity, descriptor, *file_entities],
    }
    path = root / "ro-crate-metadata.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path
