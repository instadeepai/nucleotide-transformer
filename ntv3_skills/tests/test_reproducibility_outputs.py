"""Contracts for replay, notebook, and RO-Crate report outputs."""

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.plan import AnalysisPlan, write_plan
from ntv3_tools.provenance import build_manifest
from ntv3_tools.replay import (
    ReplayRequest,
    _reject_unsafe_input_path,
    replay_analysis,
)
from ntv3_tools.report import write_viewer_bundle


def make_plan(*, want_notebook: bool) -> AnalysisPlan:
    """Create a valid report plan."""
    return AnalysisPlan(
        species="human",
        assembly="hg38",
        gene="HBE1",
        chrom="chr11",
        start=5_260_000,
        end=5_276_384,
        tissue="K562",
        tracks=["K562 CAGE"],
        device="cpu",
        window_bp=16_384,
        want_report=True,
        want_notebook=want_notebook,
        want_critic=False,
        ucsc_tracks=[],
        analyses=["predict"],
    )


def prepare_run(root: Path, *, want_notebook: bool) -> bytes:
    """Write the minimum validated report inputs."""
    plan_path = write_plan(
        plan=make_plan(want_notebook=want_notebook),
        out_dir=str(root),
    )
    plan_bytes = plan_path.read_bytes()
    (root / "track_resolve.json").write_text(
        json.dumps(
            {
                "tracks": [
                    {
                        "query": "K562 CAGE",
                        "gt_accession": "ENCSR000CIO",
                        "model_head": "CNhs11250_P",
                    }
                ]
            }
        )
    )
    tracks = root / "tracks"
    tracks.mkdir()
    np.save(tracks / "tracks.npy", np.ones((1, 16)))
    (tracks / "accessions.json").write_text('["CNhs11250_P"]')
    (tracks / "metadata.json").write_text(
        '{"chrom":"chr11","start":5267168,"end":5269216,"bin_size":128}'
    )
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "model": "approved-checkpoint",
                "commit": "abc123",
                "worktree_dirty": True,
            }
        )
    )
    return plan_bytes


class ReproducibilityOutputTests(unittest.TestCase):
    def test_requested_outputs_are_complete_and_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plan_bytes = prepare_run(root=root, want_notebook=True)
            (root / ".env").write_text("API_TOKEN=do-not-package\n")
            (root / "secrets").mkdir()
            (root / "secrets" / "credential.txt").write_text("do-not-package\n")

            write_viewer_bundle(run_dir=str(root))

            replay = root / "code" / "run_analysis.py"
            notebook = root / "analysis.ipynb"
            crate_path = root / "ro-crate-metadata.json"
            self.assertTrue(replay.is_file())
            self.assertTrue(notebook.is_file())
            self.assertTrue(crate_path.is_file())
            self.assertEqual((root / "plan.yaml").read_bytes(), plan_bytes)

            notebook_payload = json.loads(notebook.read_text())
            notebook_text = json.dumps(notebook_payload)
            for cell in notebook_payload["cells"]:
                if cell["cell_type"] == "code":
                    compile("".join(cell["source"]), "analysis.ipynb", "exec")
            for artifact in (
                "plan.yaml",
                "report.json",
                "run_manifest.json",
                "tracks/tracks.npy",
                "interpret",
                "disrupt",
                "critic/audit_log.json",
                "literature/search_manifest.json",
            ):
                self.assertIn(artifact, notebook_text)
            self.assertIn("does not establish causality", notebook_text)

            report = json.loads((root / "report.json").read_text())
            statuses = {row["name"]: row["status"] for row in report["trace"]}
            self.assertEqual(statuses["replay"], "ok")
            self.assertEqual(statuses["notebook"], "ok")
            self.assertEqual(statuses["ro-crate"], "ok")

            crate = json.loads(crate_path.read_text())
            self.assertEqual(crate["@context"], "https://w3id.org/ro/crate/1.1/context")
            graph = crate["@graph"]
            self.assertEqual(graph[0]["@id"], "./")
            descriptor = next(
                row for row in graph if row["@id"] == "ro-crate-metadata.json"
            )
            self.assertEqual(descriptor["@type"], "CreativeWork")
            file_rows = [row for row in graph if row.get("@type") == "File"]
            paths = [row["@id"] for row in file_rows]
            self.assertEqual(paths, sorted(paths))
            self.assertNotIn("ro-crate-metadata.json", paths)
            self.assertIn("code/run_analysis.py", paths)
            self.assertIn("analysis.ipynb", paths)
            self.assertIn("report.json", paths)
            self.assertIn("plan.yaml", paths)
            for row in file_rows:
                artifact = root / row["@id"]
                self.assertFalse(Path(row["@id"]).is_absolute())
                self.assertEqual(row["contentSize"], artifact.stat().st_size)
                self.assertEqual(
                    row["sha256"],
                    hashlib.sha256(artifact.read_bytes()).hexdigest(),
                )
            provenance = {
                row["name"]: row["value"] for row in graph[0]["additionalProperty"]
            }
            self.assertEqual(provenance["checkpoint"], "approved-checkpoint")
            self.assertEqual(provenance["commit"], "abc123")
            self.assertIs(provenance["worktree_dirty"], True)
            serialized = crate_path.read_text()
            self.assertNotIn(str(root), serialized)
            self.assertNotIn("credential", serialized.lower())
            self.assertNotIn("do-not-package", serialized)
            for generated_path in (replay, notebook, root / "report.json", crate_path):
                self.assertNotIn(str(root), generated_path.read_text())

    def test_replay_copies_plan_exactly_and_skips_optional_workflows(self) -> None:
        worktree = Path(__file__).resolve().parents[2]
        results = worktree / "results"
        results.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(
            dir=results,
            prefix="replay-copy-source-",
        ) as source:
            source_path = Path(source)
            approved_plan = prepare_run(root=source_path, want_notebook=False)
            destination = results / f"{source_path.name}-destination"
            self.addCleanup(shutil.rmtree, destination, True)
            commands: list[list[str]] = []

            with patch(
                "ntv3_tools.replay._run_command",
                side_effect=lambda command, environment: commands.append(command),
            ):
                replay_analysis(
                    request=ReplayRequest(
                        source_run=source_path,
                        output_run=str(destination),
                        research_python=sys.executable,
                        claims=None,
                        literature_records=None,
                        variant=None,
                    )
                )

            self.assertEqual((destination / "plan.yaml").read_bytes(), approved_plan)
            command_names = [command[command.index("-m") + 2] for command in commands]
            self.assertEqual(command_names, ["lookup", "predict", "report"])
            self.assertFalse((destination / "claims.json").exists())
            self.assertFalse((destination / "literature").exists())
            self.assertFalse((destination / "disrupt").exists())

    def test_replay_keeps_a_venv_python_symlink(self) -> None:
        worktree = Path(__file__).resolve().parents[2]
        results = worktree / "results"
        results.mkdir(exist_ok=True)
        with (
            tempfile.TemporaryDirectory(dir=results, prefix="replay-venv-") as source,
            tempfile.TemporaryDirectory() as venv_root,
        ):
            source_path = Path(source)
            prepare_run(root=source_path, want_notebook=False)
            destination = results / f"{source_path.name}-destination"
            self.addCleanup(shutil.rmtree, destination, True)
            venv_python = Path(venv_root) / "venv" / "bin" / "python"
            venv_python.parent.mkdir(parents=True)
            venv_python.symlink_to(sys.executable)
            commands: list[list[str]] = []

            with patch(
                "ntv3_tools.replay._run_command",
                side_effect=lambda command, environment: commands.append(command),
            ):
                replay_analysis(
                    request=ReplayRequest(
                        source_run=source_path,
                        output_run=str(destination),
                        research_python=str(venv_python),
                        claims=None,
                        literature_records=None,
                        variant=None,
                    )
                )

            predict = next(command for command in commands if "predict" in command)
            self.assertEqual(predict[0], str(venv_python.absolute()))

    def test_replay_refuses_a_python_that_links_into_secrets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            secret = Path(tmp) / "secrets" / "python"
            secret.parent.mkdir()
            secret.write_text("")
            link = Path(tmp) / "python"
            link.symlink_to(secret)
            with self.assertRaisesRegex(RequestError, "secrets"):
                _reject_unsafe_input_path(path=str(link), label="--research-python")

    def test_manifest_outside_a_git_checkout_records_unknown(self) -> None:
        for error in (
            subprocess.CalledProcessError(128, "git"),
            FileNotFoundError("git"),
        ):
            with (
                self.subTest(error=type(error).__name__),
                patch("ntv3_tools.provenance.subprocess.run", side_effect=error),
            ):
                manifest = build_manifest(
                    model_name="NTv3_tracks", genome_build="hg38", objective="test"
                )
            self.assertEqual(manifest["commit"], "unknown")
            self.assertIsNone(manifest["worktree_dirty"])

    def test_notebook_is_absent_when_not_requested(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prepare_run(root=root, want_notebook=False)

            write_viewer_bundle(run_dir=str(root))

            self.assertFalse((root / "analysis.ipynb").exists())
            crate = json.loads((root / "ro-crate-metadata.json").read_text())
            self.assertNotIn(
                "analysis.ipynb",
                [row["@id"] for row in crate["@graph"] if row.get("@type") == "File"],
            )

    def test_report_refuses_partial_generated_output_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prepare_run(root=root, want_notebook=True)
            replay = root / "code" / "run_analysis.py"
            replay.parent.mkdir()
            replay.write_text("do not overwrite\n")

            with self.assertRaisesRegex(ArtifactError, "already exists"):
                write_viewer_bundle(run_dir=str(root))

            self.assertEqual(replay.read_text(), "do not overwrite\n")
            self.assertFalse((root / "report.json").exists())

    def test_generated_replay_rejects_existing_and_outside_destinations(self) -> None:
        worktree = Path(__file__).resolve().parents[2]
        results = worktree / "results"
        results.mkdir(exist_ok=True)
        with (
            tempfile.TemporaryDirectory(dir=results, prefix="replay-source-") as source,
            tempfile.TemporaryDirectory(
                dir=results, prefix="replay-existing-"
            ) as existing,
            tempfile.TemporaryDirectory(prefix="replay-outside-") as outside_root,
        ):
            source_path = Path(source)
            approved_plan = prepare_run(root=source_path, want_notebook=False)
            write_viewer_bundle(run_dir=str(source_path))
            replay = source_path / "code" / "run_analysis.py"
            destinations = (Path(existing), Path(outside_root) / "new-run")
            for output in destinations:
                completed = subprocess.run(
                    [
                        sys.executable,
                        str(replay),
                        "--research-python",
                        "/usr/bin/false",
                        "--output-run",
                        str(output),
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn("--output-run", completed.stderr)
            self.assertEqual((source_path / "plan.yaml").read_bytes(), approved_plan)


if __name__ == "__main__":
    unittest.main()
