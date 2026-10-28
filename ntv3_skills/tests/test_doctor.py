"""Find the JAX Python in this repo before asking for one. No model, no checkpoint."""

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from ntv3_tools.doctor import (
    GPU_INSTALL_COMMAND,
    INSTALL_COMMAND,
    DoctorReport,
    device_for,
    doctor_payload,
    gpu_line,
    jax_python_candidates,
    jax_python_line,
    run_doctor,
)

_IMPORTS_BOTH = "3.11.13\n0.4.35\nnucleotide_transformer_v3\n"


def _python(root: Path, venv: str) -> Path:
    """Create an empty stand-in interpreter at <root>/<venv>/bin/python."""
    path = root / venv / "bin" / "python"
    path.parent.mkdir(parents=True)
    path.write_text("")
    return path


def _ok(command: list[str], **_: object) -> subprocess.CompletedProcess:
    """Fake probe: the interpreter imports jax and the package."""
    return subprocess.CompletedProcess(command, 0, _IMPORTS_BOTH, "")


class CandidateTests(unittest.TestCase):
    def test_env_var_then_this_repo_venv_and_nothing_else(self) -> None:
        candidates = jax_python_candidates(
            environ={"NTV3_JAX_PYTHON": "/opt/jax/bin/python"},
            repo_root=Path("/work/nucleotide-transformer-private"),
        )

        self.assertEqual(
            candidates,
            [
                ("/opt/jax/bin/python", "$NTV3_JAX_PYTHON"),
                (
                    "/work/nucleotide-transformer-private/.venv/bin/python",
                    ".venv in this repo (pip install -e .)",
                ),
            ],
        )

    def test_sibling_checkouts_are_never_candidates_even_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            repo = root / "nucleotide-transformer-private"
            repo.mkdir()
            for sibling in (
                "other-checkout",
            ):
                _python(root, f"{sibling}/.venv")
            candidates = jax_python_candidates(environ={}, repo_root=repo)

        self.assertEqual(
            candidates,
            [
                (
                    str(repo / ".venv" / "bin" / "python"),
                    ".venv in this repo (pip install -e .)",
                )
            ],
        )

    def test_planner_venv_is_not_a_jax_candidate(self) -> None:
        candidates = jax_python_candidates(environ={}, repo_root=Path("/work/repo"))

        self.assertNotIn("ntv3_skills", json.dumps(candidates))

    def test_env_var_naming_a_venv_uses_its_python(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            candidates = jax_python_candidates(
                environ={"NTV3_JAX_PYTHON": tmp}, repo_root=Path("/work/repo")
            )

        self.assertEqual(candidates[0][0], str(Path(tmp, "bin", "python")))


class RunDoctorTests(unittest.TestCase):
    def test_a_python_must_import_jax_and_the_package(self) -> None:
        def jax_only(command: list[str], **_: object) -> subprocess.CompletedProcess:
            if Path(command[0]).parents[1].name == "env":
                return subprocess.CompletedProcess(
                    command,
                    1,
                    "3.11.13\n0.4.35\n",
                    "ModuleNotFoundError: No module named 'nucleotide_transformer_v3'\n",
                )
            return _ok(command)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            configured = _python(root, "env")
            repo_venv = _python(root, "repo/.venv")
            report = run_doctor(
                candidates=[
                    (str(configured), "$NTV3_JAX_PYTHON"),
                    (str(repo_venv), ".venv in this repo (pip install -e .)"),
                ],
                run=jax_only,
            )

        self.assertEqual(report.jax_python, str(repo_venv))
        failed, passed = report.candidates
        self.assertFalse(failed.passed)
        self.assertEqual(failed.jax_version, "0.4.35")
        self.assertFalse(failed.imports_package)
        self.assertIn("nucleotide_transformer_v3", failed.detail)
        self.assertTrue(passed.passed)
        self.assertTrue(passed.imports_package)
        self.assertEqual(
            (passed.python_version, passed.jax_version), ("3.11.13", "0.4.35")
        )

    def test_missing_repo_venv_names_the_install_and_borrows_nothing(self) -> None:
        def refuse_run(command: list[str], **_: object) -> subprocess.CompletedProcess:
            raise AssertionError(f"Ran {command}")

        report = run_doctor(
            candidates=[("/nowhere/.venv/bin/python", ".venv in this repo")],
            run=refuse_run,
        )
        line = jax_python_line(report=report)

        self.assertIsNone(report.jax_python)
        self.assertEqual(report.candidates[0].detail, "not found")
        self.assertIn(INSTALL_COMMAND, line)
        self.assertIn("Predict cannot start", line)
        self.assertIn(".venv is missing", line)
        self.assertNotIn("PYTHONPATH", line)

    def test_install_command_is_cpu_only(self) -> None:
        self.assertEqual(
            INSTALL_COMMAND,
            "python3.11 -m venv .venv && .venv/bin/pip install -e .",
        )

    def test_selected_line_names_both_imports(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            python = _python(Path(tmp), "repo/.venv")
            report = run_doctor(
                candidates=[(str(python), ".venv in this repo (pip install -e .)")],
                run=_ok,
            )
        line = jax_python_line(report=report)

        self.assertIn(str(python), line)
        self.assertIn("jax 0.4.35", line)
        self.assertIn("nucleotide_transformer_v3", line)
        self.assertNotIn("cannot start", line)
        self.assertNotIn("PYTHONPATH", line)

    def test_probe_imports_only_jax_and_this_package(self) -> None:
        probes: list[str] = []

        def record(command: list[str], **_: object) -> subprocess.CompletedProcess:
            probes.append(command[-1])
            return _ok(command)

        with tempfile.TemporaryDirectory() as tmp:
            python = _python(Path(tmp), "repo/.venv")
            run_doctor(candidates=[(str(python), "repo")], run=record)

        (probe,) = probes
        self.assertIn("import jax", probe)
        self.assertIn("import nucleotide_transformer_v3", probe)

    def test_payload_holds_paths_and_versions_only(self) -> None:
        def fake_run(
            command: list[str], **kwargs: object
        ) -> subprocess.CompletedProcess:
            self.assertNotIn("SECRET_TOKEN", json.dumps(kwargs.get("env") or {}))
            return _ok(command)

        with tempfile.TemporaryDirectory() as tmp:
            python = _python(Path(tmp), "jax-env")
            candidates = jax_python_candidates(
                environ={"NTV3_JAX_PYTHON": str(python), "SECRET_TOKEN": "hf_s3cr3t"},
                repo_root=Path(tmp, "repo"),
            )
            payload = doctor_payload(report=run_doctor(candidates, run=fake_run))

        self.assertNotIn("hf_s3cr3t", json.dumps(payload))
        self.assertEqual(payload["jax_python"], str(python))
        self.assertNotIn("research_root", payload)
        self.assertEqual(
            set(payload["candidates"][0]),
            {
                "python",
                "rule",
                "passed",
                "python_version",
                "jax_version",
                "imports_package",
                "detail",
                "jax_backend",
            },
        )

    def _report_on(self, backend: str) -> DoctorReport:
        def probe(command: list[str], **_: object) -> subprocess.CompletedProcess:
            return subprocess.CompletedProcess(
                command, 0, _IMPORTS_BOTH + backend + "\n", ""
            )

        with tempfile.TemporaryDirectory() as tmp:
            python = _python(Path(tmp), "repo/.venv")
            return run_doctor(candidates=[(str(python), "repo")], run=probe)

    def test_cpu_backend_names_what_needs_a_gpu(self) -> None:
        report = self._report_on("cpu")
        line = gpu_line(report=report)

        self.assertIsNotNone(line)
        self.assertIn("runs on cpu", line)
        self.assertIn("Annotate", line)
        self.assertIn("16,384 bp", line)
        self.assertIn(GPU_INSTALL_COMMAND, line)
        self.assertEqual(doctor_payload(report=report)["gpu"], line)
        self.assertEqual(device_for(report=report), "cpu")

    def test_gpu_backend_runs_everything_on_cuda(self) -> None:
        report = self._report_on("gpu")

        self.assertIsNone(gpu_line(report=report))
        self.assertNotIn("gpu", doctor_payload(report=report))
        self.assertEqual(device_for(report=report), "cuda")

    def test_no_python_has_no_gpu_line_and_plans_cpu(self) -> None:
        report = run_doctor(candidates=[("/nowhere/.venv/bin/python", "repo")])

        self.assertIsNone(gpu_line(report=report))
        self.assertEqual(device_for(report=report), "cpu")

    def test_probe_sees_the_gpu_visibility_variables_and_no_secret(self) -> None:
        envs: list[dict[str, str]] = []

        def record(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
            envs.append(dict(kwargs["env"]))
            return _ok(command)

        visibility = {
            "CUDA_VISIBLE_DEVICES": "",
            "JAX_PLATFORMS": "cpu",
            "LD_LIBRARY_PATH": "/opt/cuda/lib64",
        }
        with (
            tempfile.TemporaryDirectory() as tmp,
            mock.patch.dict("os.environ", {**visibility, "HF_TOKEN": "hf_s3cr3t"}),
        ):
            python = _python(Path(tmp), "repo/.venv")
            run_doctor(candidates=[(str(python), "repo")], run=record)

        (env,) = envs
        self.assertEqual({key: env[key] for key in visibility}, visibility)
        self.assertNotIn("HF_TOKEN", env)

    def test_payload_without_a_python_carries_the_install_command(self) -> None:
        payload = doctor_payload(
            report=run_doctor(candidates=[("/nowhere/.venv/bin/python", "repo")])
        )

        self.assertIsNone(payload["jax_python"])
        self.assertEqual(payload["install"], INSTALL_COMMAND)


if __name__ == "__main__":
    unittest.main()
