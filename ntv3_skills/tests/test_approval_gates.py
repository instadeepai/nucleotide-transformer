"""Gates before any model load: the run name, the questions, approval, critic.

No model, no network, and no reused example locus: values are placeholders.
"""

import argparse
import contextlib
import dataclasses
import io
import json
import re
import tempfile
import unittest
from pathlib import Path
from typing import ClassVar
from unittest import mock

import yaml
from ntv3_tools.approval import (
    RESOLVED_BLOCK,
    approve_run,
    require_approval,
)
from ntv3_tools.cli import (
    ASSEMBLY_QUESTION,
    GPU_WINDOW_QUESTION,
    LANES_QUESTION,
    SPECIES_QUESTION,
    TISSUE_QUESTION,
    TRACKS_QUESTION,
    WINDOW_QUESTION,
    _cmd_predict,
)
from ntv3_tools.cli import main as cli_main
from ntv3_tools.doctor import DoctorReport
from ntv3_tools.errors import RequestError
from ntv3_tools.gt import search_encode_experiment
from ntv3_tools.plan import read_plan, run_dir_name
from ntv3_tools.report import _run_claims, _run_critic, lane_label

from tests.test_skill_resolution import _JAX


def _on(backend: str) -> DoctorReport:
    """The doctor fixture with its selected JAX Python on ``backend``."""
    chosen = dataclasses.replace(_JAX.candidates[-1], jax_backend=backend)
    return dataclasses.replace(_JAX, candidates=[*_JAX.candidates[:-1], chosen])


_PLAN_ARGS = [
    "plan",
    "--species", "human",
    "--assembly", "hg38",
    "--gene", "GENE7",
    "--chrom", "chr4",
    "--start", "1000",
    "--end", "5096",
    "--device", "cpu",
    "--window-bp", "4096",
    "--ucsc-tracks", "none",
    "--analyses", "predict", "saliency",
    "--want-report",
]  # fmt: skip
# Everything a gene-lookup plan needs except the window and the browser lanes.
_NO_WINDOW_NO_LANES = [
    "--species", "human", "--assembly", "hg38", "--gene", "GENE7",
    "--tissue", "tissue-x", "--tracks", "assay-y", "--device", "cpu",
    "--analyses", "predict", "--want-report",
]  # fmt: skip


class GateTests(unittest.TestCase):
    def setUp(self) -> None:
        root = tempfile.TemporaryDirectory()
        self.addCleanup(root.cleanup)
        self.results = Path(root.name).resolve()
        for target, value in (
            ("ntv3_tools.cli._RESULTS_ROOT", self.results),
            ("ntv3_tools.cli.run_doctor", mock.Mock(return_value=_JAX)),
            (
                "ntv3_tools.cli.resolve_gene_locus",
                mock.Mock(side_effect=AssertionError),
            ),
            (
                "ntv3_tools.cli.search_encode_experiment",
                mock.Mock(side_effect=AssertionError("ENCODE searched")),
            ),
        ):
            patcher = mock.patch(target, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def _plan(self, *extra: str) -> str:
        with (
            mock.patch("sys.argv", ["ntv3", *_PLAN_ARGS, *extra]),
            contextlib.redirect_stdout(io.StringIO()) as out,
        ):
            cli_main()
        return out.getvalue()

    def _refused(self, *extra: str) -> str:
        with (
            mock.patch("sys.argv", ["ntv3", *_PLAN_ARGS, *extra]),
            self.assertLogs("ntv3_tools", level="ERROR") as logs,
            self.assertRaises(SystemExit) as exit_,
        ):
            cli_main()
        self.assertEqual(exit_.exception.code, 2)
        return "\n".join(logs.output)

    def _bare(self, *argv: str) -> tuple[int, str]:
        """Run `ntv3 plan` with exactly these flags; return exit code and all output."""
        out, err = io.StringIO(), io.StringIO()
        with (
            mock.patch("sys.argv", ["ntv3", "plan", *argv]),
            contextlib.redirect_stdout(out),
            contextlib.redirect_stderr(err),
            self.assertLogs("ntv3_tools", level="ERROR") as logs,
            self.assertRaises(SystemExit) as exit_,
        ):
            cli_main()
        return int(exit_.exception.code), "\n".join(
            [out.getvalue(), err.getvalue(), *logs.output]
        )

    def test_a_gene_alone_asks_for_the_species_not_for_usage(self) -> None:
        code, output = self._bare(
            "--gene", "BRCA1", "--device", "cpu", "--window-bp", "4096",
            "--analyses", "predict", "--resolve-only",
        )  # fmt: skip

        self.assertEqual(code, 2)
        self.assertIn("Which species?", output)
        self.assertNotIn("the following arguments are required", output)
        self.assertNotIn("usage:", output)
        self.assertNotIn("Which tissue", output)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_a_species_alone_asks_for_the_assembly(self) -> None:
        code, output = self._bare(
            "--species", "human", "--gene", "BRCA1", "--device", "cpu",
            "--window-bp", "4096", "--analyses", "predict", "--resolve-only",
        )  # fmt: skip

        self.assertEqual(code, 2)
        self.assertIn(
            "Which assembly? Pass --assembly. No assembly is assumed.", output
        )
        self.assertNotIn("usage:", output)
        self.assertNotIn("Which tissue", output)
        self.assertEqual(list(self.results.iterdir()), [])

    def _asked(self, argv: list[str]) -> tuple[str, str]:
        """Run ``ntv3 <argv>``; return the logged question and stderr."""
        stderr = io.StringIO()
        with (
            mock.patch("sys.argv", ["ntv3", *argv]),
            contextlib.redirect_stderr(stderr),
            self.assertLogs("ntv3_tools", level="ERROR") as logs,
            self.assertRaises(SystemExit) as exit_,
        ):
            cli_main()
        self.assertEqual(exit_.exception.code, 2)
        return "\n".join(logs.output), stderr.getvalue()

    def test_a_gene_with_no_species_asks_for_the_species_first(self) -> None:
        message, stderr = self._asked(
            [
                "plan",
                "--gene",
                "BRCA1",
                "--device",
                "cpu",
                "--window-bp",
                "4096",
                "--analyses",
                "predict",
                "--want-report",
            ]
        )

        self.assertIn(SPECIES_QUESTION, message)
        self.assertNotIn(TISSUE_QUESTION, message)
        self.assertNotIn("usage:", stderr)
        self.assertIsNone(
            re.search(r"human|hg38|mouse|mm10", message.split(":", 2)[-1])
        )
        self.assertEqual(list(self.results.iterdir()), [])

    def test_a_species_with_no_assembly_asks_for_the_assembly(self) -> None:
        message, stderr = self._asked(
            [
                "plan",
                "--species",
                "mouse",
                "--gene",
                "BRCA1",
                "--device",
                "cpu",
                "--window-bp",
                "4096",
                "--analyses",
                "predict",
                "--want-report",
            ]
        )

        self.assertIn(ASSEMBLY_QUESTION, message)
        self.assertNotIn(TISSUE_QUESTION, message)
        self.assertNotIn("usage:", stderr)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_no_tissue_exits_2_with_a_question_that_names_no_track(self) -> None:
        message = self._refused("--tracks", "ENCSR100AAA")

        self.assertIn(TISSUE_QUESTION, message)
        self.assertIsNone(re.search(r"ENCSR|CNhs|K562|keratinocyte|retina", message))
        self.assertEqual(list(self.results.iterdir()), [])

    def test_an_assay_encode_lacks_for_the_tissue_is_one_question(self) -> None:
        def no_experiments(tissue: str, assay: str, species: str) -> str:
            return search_encode_experiment(
                tissue=tissue,
                assay=assay,
                get_json=lambda url, params=None: {"@graph": [], "total": 0},
                species=species,
            )

        with mock.patch("ntv3_tools.cli.search_encode_experiment", no_experiments):
            message = self._refused(
                "--tissue", "tissue-x", "--tracks", "assay-y", "--resolve-only"
            )  # fmt: skip

        self.assertIn(
            "ENCODE has no released assay-y experiment for tissue-x. Which ENCSR id "
            "should be predicted for assay-y in tissue-x? Pass it in --tracks. No id "
            "is assumed.",
            message,
        )
        for agent_talk in ("Do not invent", "Ask, or pass"):
            self.assertNotIn(agent_talk, message)
        self.assertIsNone(re.search(r"ENCSR\d{3}[A-Z]{3}", message))
        self.assertEqual(list(self.results.iterdir()), [])

    def test_several_encode_experiments_are_one_question_that_picks_none(
        self,
    ) -> None:
        hits = {
            "@graph": [{"accession": "ENCSR900AAA"}, {"accession": "ENCSR900BBB"}],
            "total": 2,
        }

        def two_experiments(tissue: str, assay: str, species: str) -> str:
            return search_encode_experiment(
                tissue=tissue,
                assay=assay,
                get_json=lambda url, params=None: hits,
                species=species,
            )

        with mock.patch("ntv3_tools.cli.search_encode_experiment", two_experiments):
            message = self._refused(
                "--tissue", "tissue-x", "--tracks", "assay-y", "--resolve-only"
            )  # fmt: skip

        self.assertIn(
            "ENCODE has several released assay-y experiments for tissue-x: "
            "ENCSR900AAA, ENCSR900BBB. Which one should be used? Pass its ENCSR id "
            "in --tracks. None is picked for you.",
            message,
        )
        self.assertNotIn("Ask.", message)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_no_tracks_exits_2_with_a_question(self) -> None:
        message = self._refused("--tissue", "tissue-x")

        self.assertIn(TRACKS_QUESTION, message)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_no_window_is_one_question_before_lookup_encode_or_any_write(
        self,
    ) -> None:
        # Gene lookup and ENCODE are patched to fail the test if reached.
        code, output = self._bare(*_NO_WINDOW_NO_LANES, "--ucsc-tracks", "none")

        self.assertEqual(code, 2)
        self.assertIn(WINDOW_QUESTION, output)
        self.assertIn("16384 or more runs with a warning", WINDOW_QUESTION)
        self.assertNotIn(LANES_QUESTION, output)
        self.assertEqual(output.count("?"), 1)
        self.assertNotIn("usage:", output)
        self.assertEqual(list(self.results.iterdir()), [])

    def _cpu_block(self, window_bp: int, *analyses: str, device: str = "cpu") -> str:
        """Print the Resolved block of a plan over chr4 from 1000; write nothing."""
        argv = [
            "ntv3", "plan", "--resolve-only", "--species", "human",
            "--assembly", "hg38", "--gene", "GENE7", "--chrom", "chr4",
            "--start", "1000", "--end", str(1000 + window_bp),
            "--window-bp", str(window_bp), "--device", device,
            "--ucsc-tracks", "none", "--tissue", "tissue-x",
            "--tracks", "ENCSR100AAA", "--analyses", *analyses, "--want-report",
        ]  # fmt: skip
        with (
            mock.patch("sys.argv", argv),
            contextlib.redirect_stdout(io.StringIO()) as out,
        ):
            cli_main()
        return out.getvalue()

    def test_a_short_cpu_window_has_no_warning(self) -> None:
        self.assertNotIn("CPU warning", self._cpu_block(4096, "predict", "saliency"))

    def test_a_cpu_window_from_16384_plans_with_a_warning(self) -> None:
        predict = self._cpu_block(16384, "predict")
        short_saliency = self._cpu_block(131072, "predict", "saliency")
        saliency = self._cpu_block(1048576, "predict", "saliency")

        self.assertIn("CPU warning: a 16,384 bp window on CPU is slow", predict)
        self.assertIn("needs a lot of memory", predict)
        self.assertNotIn("Saliency needs", predict)
        self.assertIn("CPU warning", short_saliency)
        self.assertNotIn("Saliency needs", short_saliency)
        self.assertIn("CPU warning: a 1,048,576 bp window on CPU is slow", saliency)
        self.assertIn("Saliency needs", saliency)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_gpu_saliency_above_131072_plans_with_a_warning(self) -> None:
        long_saliency = self._cpu_block(262144, "predict", "saliency", device="cuda")
        tested = self._cpu_block(131072, "predict", "saliency", device="cuda")
        long_predict = self._cpu_block(1048576, "predict", device="cuda")

        self.assertIn(
            "GPU warning: Saliency over a 262,144 bp window may run out of GPU memory.",
            long_saliency,
        )
        self.assertNotIn("80 GB", long_saliency)
        self.assertNotIn("CPU warning", long_saliency)
        self.assertNotIn("GPU warning", tested)
        self.assertNotIn("GPU warning", long_predict)

    def test_a_bad_window_is_refused_before_the_lanes_question_or_any_lookup(
        self,
    ) -> None:
        # Gene lookup and ENCODE are patched to fail the test if reached.
        for window, message in (
            ("5000", "must be divisible by 128"),
            ("2097152", "exceeds the released 1,048,576 bp"),
        ):
            with self.subTest(window=window):
                code, output = self._bare(*_NO_WINDOW_NO_LANES, "--window-bp", window)

                self.assertEqual(code, 2)
                self.assertIn(message, output)
                self.assertNotIn(LANES_QUESTION, output)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_an_unknown_device_is_refused(self) -> None:
        typo = ["cdua" if arg == "cpu" else arg for arg in _NO_WINDOW_NO_LANES]
        argv = ["ntv3", "plan", *typo, "--window-bp", "4096", "--ucsc-tracks", "none"]
        with (
            mock.patch("sys.argv", argv),
            contextlib.redirect_stderr(io.StringIO()) as err,
            self.assertRaises(SystemExit) as exit_,
        ):
            cli_main()

        self.assertEqual(exit_.exception.code, 2)
        self.assertIn("invalid choice: 'cdua'", err.getvalue())
        self.assertEqual(list(self.results.iterdir()), [])

    def test_without_a_device_the_window_question_follows_the_jax_backend(
        self,
    ) -> None:
        no_device = [
            arg for arg in _NO_WINDOW_NO_LANES if arg not in ("--device", "cpu")
        ]
        for backend, question in (
            ("gpu", GPU_WINDOW_QUESTION),
            ("cpu", WINDOW_QUESTION),
        ):
            with (
                self.subTest(backend=backend),
                mock.patch("ntv3_tools.cli.run_doctor", return_value=_on(backend)),
            ):
                code, output = self._bare(*no_device, "--ucsc-tracks", "none")

                self.assertEqual(code, 2)
                self.assertIn(question, output)
                self.assertEqual(output.count("?"), 1)
        self.assertIn("up to 1048576", GPU_WINDOW_QUESTION)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_without_a_device_a_gpu_plans_cuda_and_a_long_window(self) -> None:
        argv = [
            "ntv3", "plan", "--species", "human", "--assembly", "hg38",
            "--gene", "GENE7", "--chrom", "chr4", "--start", "1000",
            "--end", "132072", "--window-bp", "131072", "--ucsc-tracks", "none",
            "--tissue", "tissue-x", "--tracks", "ENCSR100AAA",
            "--analyses", "predict", "--want-report",
        ]  # fmt: skip
        with (
            mock.patch("ntv3_tools.cli.run_doctor", return_value=_on("gpu")),
            mock.patch("sys.argv", argv),
            contextlib.redirect_stdout(io.StringIO()) as out,
        ):
            cli_main()
        (run_dir,) = self.results.iterdir()

        self.assertEqual(read_plan(path=run_dir / "plan.yaml").device, "cuda")
        self.assertIn(
            "chr4:1,000-132,072 (131,072 bp) on cuda — your span equals the cuda "
            "131,072 bp window; cuda because this JAX python sees a GPU",
            out.getvalue(),
        )

    def test_a_window_above_the_released_length_is_one_question(self) -> None:
        code, output = self._bare(
            "--species", "human", "--assembly", "hg38", "--gene", "GENE7",
            "--chrom", "chr4", "--start", "0", "--end", "2097152",
            "--window-bp", "2097152", "--device", "cuda", "--ucsc-tracks", "none",
            "--tissue", "tissue-x", "--tracks", "ENCSR100AAA",
            "--analyses", "predict", "--want-report",
        )  # fmt: skip

        self.assertEqual(code, 2)
        self.assertIn("exceeds the released 1,048,576 bp", output)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_no_browser_lane_choice_is_one_question_and_writes_nothing(
        self,
    ) -> None:
        code, output = self._bare(*_NO_WINDOW_NO_LANES, "--window-bp", "4096")

        self.assertEqual(code, 2)
        self.assertIn(
            "Include gene and ClinVar browser lanes? Pass --ucsc-tracks with "
            "names, or pass none. No lane is added for you.",
            output,
        )
        self.assertNotIn(WINDOW_QUESTION, output)
        self.assertEqual(output.count("?"), 1)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_none_adds_no_lane_and_named_lanes_are_kept_as_named(self) -> None:
        block = self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")
        (run_dir,) = self.results.iterdir()

        self.assertEqual(read_plan(path=run_dir / "plan.yaml").ucsc_tracks, [])
        self.assertIn("Browser lanes: none — you passed none", block)
        for absent in ("knownGene", "clinvarMain", "ncbiRefSeq"):
            self.assertNotIn(absent, block)

        block = self._plan(
            "--tissue", "tissue-x", "--tracks", "ENCSR100AAA", "--start", "2000",
            "--end", "6096", "--ucsc-tracks", "knownGene", "clinvarMain",
        )  # fmt: skip
        plan = read_plan(
            path=self.results / "hg38_chr4_2000_6096_predict+saliency" / "plan.yaml"
        )
        self.assertEqual(plan.ucsc_tracks, ["knownGene", "clinvarMain"])
        self.assertIn(
            "Browser lanes: knownGene (GENCODE); clinvarMain (ClinVar) — named in "
            "the request",
            block,
        )

    def test_a_run_directory_made_by_a_second_writer_is_refused_unwritten(
        self,
    ) -> None:
        run_dir = self.results / "hg38_chr4_1000_5096_predict+saliency"

        def other_plan_wins(**_: object) -> object:
            # Another `ntv3 plan` creates the directory after our exists check.
            run_dir.mkdir()
            return _JAX

        with mock.patch("ntv3_tools.cli.run_doctor", side_effect=other_plan_wins):
            message = self._refused("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")

        self.assertIn(
            "Which run directory? results/hg38_chr4_1000_5096_predict+saliency "
            "already exists and is never overwritten.",
            message,
        )
        self.assertEqual(list(run_dir.iterdir()), [])

    def test_the_resolved_block_never_says_the_model_loaded(self) -> None:
        block = self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")

        self.assertNotIn("loaded", block)
        self.assertIn("loads only after approve", block)
        self.assertTrue(block.rstrip().endswith("nothing loads before it."))

    def test_the_run_directory_is_the_formula_and_is_never_reused(self) -> None:
        block = self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")

        (run_dir,) = self.results.iterdir()
        plan = read_plan(path=run_dir / "plan.yaml")
        self.assertEqual(run_dir.name, "hg38_chr4_1000_5096_predict+saliency")
        self.assertEqual(run_dir.name, run_dir_name(plan=plan))
        self.assertFalse(plan.approved)
        self.assertTrue(plan.want_critic)
        self.assertEqual((run_dir / RESOLVED_BLOCK).read_text().strip(), block.strip())
        for field in (
            "Assembly: hg38",
            "Locus: chr4 start 1000 end 5096",
            "Gene: none looked up",
            "Tissue: tissue-x",
            "Tracks: ENCSR100AAA",
            "Analyses: predict, saliency",
            "Critic: on",
            f"Run directory: results/{run_dir.name}",
            f"ntv3 approve --run-dir results/{run_dir.name}",
        ):
            self.assertIn(field, block)

        message = self._refused("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")
        self.assertIn(
            f"Which run directory? results/{run_dir.name} already exists", message
        )

    def test_predict_without_approval_never_starts_the_jax_python(self) -> None:
        self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")
        (run_dir,) = self.results.iterdir()

        with (
            mock.patch(
                "ntv3_tools.cli.importlib.util.find_spec", return_value=object()
            ),
            mock.patch(
                "ntv3_tools.cli.subprocess.run", side_effect=AssertionError("started")
            ),
            self.assertRaisesRegex(RequestError, "is not approved"),
        ):
            _cmd_predict(argparse.Namespace(plan=str(run_dir / "plan.yaml")))

    def test_only_approve_approves_and_an_edit_after_it_is_refused(self) -> None:
        self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")
        (run_dir,) = self.results.iterdir()
        plan_path = run_dir / "plan.yaml"

        # A hand-set flag is not an approval.
        payload = yaml.safe_load(plan_path.read_text())
        plan_path.write_text(yaml.safe_dump({**payload, "approved": True}))
        with self.assertRaisesRegex(RequestError, "only `ntv3 approve`"):
            require_approval(run_dir=run_dir, plan=read_plan(path=plan_path))
        plan_path.write_text(yaml.safe_dump(payload))

        plan = approve_run(run_dir=run_dir, critic=False)
        self.assertTrue(plan.approved)
        self.assertFalse(plan.want_critic)
        require_approval(run_dir=run_dir, plan=read_plan(path=plan_path))

        edited = yaml.safe_load(plan_path.read_text())
        plan_path.write_text(yaml.safe_dump({**edited, "tracks": ["ENCSR200BBB"]}))
        with self.assertRaisesRegex(RequestError, "changed after approval"):
            require_approval(run_dir=run_dir, plan=read_plan(path=plan_path))

    def test_approve_needs_the_resolved_block(self) -> None:
        self._plan("--tissue", "tissue-x", "--tracks", "ENCSR100AAA")
        (run_dir,) = self.results.iterdir()
        (run_dir / RESOLVED_BLOCK).rename(run_dir / "moved.txt")

        with self.assertRaisesRegex(RequestError, "No Resolved block"):
            approve_run(run_dir=run_dir)


class LaneLabelTests(unittest.TestCase):
    def test_a_typed_assay_keeps_its_words_and_gains_the_resolved_id(self) -> None:
        self.assertEqual(
            lane_label(query="assay-y", model_head="ENCSR100AAA"),
            "assay-y ENCSR100AAA",
        )

    def test_an_id_passed_alone_is_its_own_label(self) -> None:
        self.assertEqual(
            lane_label(query="ENCSR100AAA", model_head="ENCSR100AAA"),
            "ENCSR100AAA",
        )


def _critic(root: Path, plan: dict[str, object]) -> dict[str, object]:
    return _run_critic(
        root=root,
        enabled=bool(plan["want_critic"]),
        claims=_run_claims(root=root, plan=plan),
        log_path=root / "critic" / "audit_log.json",
    )


class CriticTests(unittest.TestCase):
    _PLAN: ClassVar[dict[str, object]] = {
        "model": "NTv3_tracks",
        "tracks": ["ENCSR100AAA"],
        "analyses": ["predict"],
        "want_critic": True,
    }

    def test_critic_off_is_reported_not_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status = _critic(Path(tmp), {**self._PLAN, "want_critic": False})
            self.assertFalse(Path(tmp, "critic").exists())

        self.assertEqual(status["enabled"], False)
        self.assertIn("critic was off", status["note"])

    def test_critic_on_runs_and_records_its_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for rel in (
                "tracks/tracks.npy",
                "tracks/metadata.json",
                "run_manifest.json",
            ):
                (root / rel).parent.mkdir(parents=True, exist_ok=True)
                (root / rel).write_text("x")
            (root / "tracks" / "accessions.json").write_text('["ENCSR100AAA"]')
            log = root / "critic" / "audit_log.json"

            status = _critic(root, self._PLAN)
            decisions = json.loads(log.read_text())

        self.assertEqual((status["kept"], status["suppressed"]), (1, 0))
        self.assertEqual(decisions[0]["claim_id"], "predict")

    def test_a_missing_artifact_is_suppressed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status = _critic(Path(tmp), self._PLAN)

        self.assertEqual((status["kept"], status["suppressed"]), (0, 1))


if __name__ == "__main__":
    unittest.main()
