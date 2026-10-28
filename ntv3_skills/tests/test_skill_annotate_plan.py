"""Isolated tests for the Annotate ask: intake, plan, locus, approval."""

import argparse
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from ntv3_tools.annotate import (
    ANNOTATION_INPUT,
    ANNOTATION_WINDOW_BP,
    AnnotateRequest,
    LocusWindow,
    annotation_resolved_block,
    clade_config,
    count_tiles,
    locus_window,
    plan_annotation,
    predict_run_window,
    read_fasta,
    shift_gff3,
)
from ntv3_tools.annotation_report import read_gff3
from ntv3_tools.approval import (
    ANNOTATION_APPROVAL,
    ANNOTATION_RESOLVED_BLOCK,
    approve_annotation,
    pending_annotation_plan,
    require_annotation_approval,
)
from ntv3_tools.cli import _cmd_annotate
from ntv3_tools.doctor import DoctorReport, JaxCandidate
from ntv3_tools.errors import ArtifactError, RequestError

_HALF = ANNOTATION_WINDOW_BP // 2


def _doctor(backend: str) -> DoctorReport:
    candidate = JaxCandidate(
        "/x/.venv/bin/python", "rule", True, "3.11.0", "0.5.3", True, "ok", backend
    )
    return DoctorReport(jax_python=candidate.python, candidates=[candidate])


def _fasta(tmp: str, text: str = ">chr1\n" + "ACGT" * 256 + "\n") -> Path:
    path = Path(tmp, "a.fa")
    path.write_text(text)
    return path


def _request(
    fasta_path: Path | None, out_dir: Path, **changes: object
) -> AnnotateRequest:
    values: dict[str, object] = {
        "clade": "insects",
        "decoder": "hmm_constrained",
        "out_dir": str(out_dir),
        "fasta": None if fasta_path is None else str(fasta_path),
    }
    values.update(changes)
    return AnnotateRequest(**values)  # type: ignore[arg-type]


class IntakeTests(unittest.TestCase):
    def test_contigs_are_read_in_windows_of_the_central_75_percent(self) -> None:
        stride = 786_432
        self.assertEqual(
            [count_tiles(n) for n in (16, stride, stride + 1, ANNOTATION_WINDOW_BP)],
            [1, 1, 2, 2],
        )

    def test_unknown_clade_is_refused(self) -> None:
        with self.assertRaisesRegex(
            RequestError, "Which clade\\? There is no NTv3 annotation model"
        ):
            clade_config("fungi")

    def test_clade_records_the_model_config_name(self) -> None:
        self.assertEqual(
            [clade_config(name) for name in ("mammals", "plants", "insects")],
            ["mammals", "plants", "insects"],
        )
        self.assertEqual(clade_config("insecta"), "insects")
        self.assertEqual(clade_config(" Insecta "), "insects")
        config_dir = Path(__file__).resolve().parents[2] / (
            "nucleotide_transformer_v3/gene_annotation/model_cfg"
        )
        self.assertEqual(
            sorted(path.stem for path in config_dir.glob("*.yaml")),
            ["insects", "mammals", "plants"],
        )

    def test_reads_soft_masked_and_gzipped_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            plain = Path(tmp, "a.fa")
            plain.write_text(">chr1 desc\nACGTacgt\nNN\n>chr2\nRYk\n")
            packed = Path(tmp, "a.fa.gz")
            packed.write_bytes(gzip.compress(plain.read_bytes()))

            for path in (plain, packed):
                records = read_fasta(path)
                self.assertEqual([r.seqid for r in records], ["chr1", "chr2"])
                self.assertEqual(records[0].sequence, "ACGTacgtNN")

    def test_rejects_malformed_fasta(self) -> None:
        cases = {
            "no_header.fa": ("ACGT\n", "'>' header"),
            "repeat.fa": (">a\nAC\n>a\nGT\n", "empty or repeated"),
            "bad_base.fa": (">a\nAC-GT\n", "non-IUPAC"),
            "empty_record.fa": (">a\n>b\nAC\n", "no bases"),
        }
        with tempfile.TemporaryDirectory() as tmp:
            for name, (text, message) in cases.items():
                path = Path(tmp, name)
                path.write_text(text)
                with (
                    self.subTest(name=name),
                    self.assertRaisesRegex(RequestError, message),
                ):
                    read_fasta(path)
            with self.assertRaisesRegex(RequestError, "does not exist"):
                read_fasta(Path(tmp, "missing.fa"))


class AnnotationResolvedBlockTests(unittest.TestCase):
    def test_block_states_model_geometry_filters_and_jax(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp, "results", "fly")
            lines = annotation_resolved_block(
                request=_request(_fasta(tmp), out_dir),
                run_dir="results/fly",
                doctor=_doctor("gpu"),
                weights_cached=True,
            )
            self.assertFalse(out_dir.exists())
        text = "\n".join(lines)

        self.assertEqual(lines[0], "Resolved")
        for expected in (
            "InstaDeepAI/NTv3_annot_insects is in the Hugging Face cache",
            "hmm_constrained",
            "1,048,576 bp, central 75% kept",
            "CDS shorter than 201 bp dropped",
            "Forward: on a GPU, after approve",
            "25 GB of GPU memory",
            "annotation.gff3 in results/fly",
            "JAX python: /x/.venv/bin/python",
        ):
            self.assertIn(expected, text)
        self.assertNotIn("will refuse", text)
        for line in lines[1:-1]:
            self.assertIn(" — ", line)
        self.assertEqual(lines[-1], "Approve to run the forward, or change a line.")

    def test_block_names_the_download_and_a_cpu_only_jax(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lines = annotation_resolved_block(
                request=_request(_fasta(tmp), Path(tmp, "results", "fly")),
                run_dir="results/fly",
                doctor=_doctor("cpu"),
                weights_cached=False,
            )
        text = "\n".join(lines)

        self.assertIn("the forward downloads 2.8 GB", text)
        self.assertIn("Hugging Face login with approved access", text)
        self.assertIn("this JAX python runs on cpu, so the forward will refuse", text)
        self.assertTrue(lines[-1].startswith("Approve to download"))

    def test_block_refuses_bad_clade_decoder_or_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = _fasta(tmp)
            for changes, message in (
                ({"clade": "fungi"}, "Which clade\\? There is no NTv3 annotation"),
                ({"fasta": str(Path(tmp, "no.fa"))}, "does not exist"),
                ({"decoder": "viterbi"}, "Unknown decoder"),
            ):
                with (
                    self.subTest(**{k: str(v) for k, v in changes.items()}),
                    self.assertRaisesRegex(RequestError, message),
                ):
                    annotation_resolved_block(
                        request=_request(fasta, Path(tmp), **changes),
                        run_dir="results/x",
                        doctor=_doctor("gpu"),
                        weights_cached=True,
                    )


class PlanAnnotationTests(unittest.TestCase):
    def test_fasta_plan_is_written_unapproved_with_its_block(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = _fasta(tmp, ">chr1\n" + "ACGT" * 256 + "\n>chr2\n" + "acgt" * 4)
            out_dir = Path(tmp, "results", "annotate-run")
            plan_path = plan_annotation(
                request=_request(fasta, out_dir), block=["Resolved", "Approve"]
            )
            plan = json.loads(plan_path.read_text())
            names = sorted(path.name for path in out_dir.iterdir())
            block = Path(out_dir, ANNOTATION_RESOLVED_BLOCK).read_text()

        self.assertEqual(names, ["annotation_plan.json", ANNOTATION_RESOLVED_BLOCK])
        self.assertEqual(block, "Resolved\nApprove\n")
        self.assertFalse(plan["approved"])
        self.assertIsNone(plan["locus"])
        self.assertEqual(plan["model"], "InstaDeepAI/NTv3_annot_insects")
        self.assertEqual(plan["output"], "annotation.gff3")
        self.assertEqual(plan["tiles"], 2)
        self.assertEqual(
            plan["records"],
            [
                {"seqid": "chr1", "length": 1024},
                {"seqid": "chr2", "length": 16},
            ],
        )

    def test_invalid_ask_or_nonempty_run_writes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp, "results", "bad")
            with self.assertRaisesRegex(RequestError, "must start with a '>'"):
                plan_annotation(
                    request=_request(_fasta(tmp, "ACGT\n"), out_dir),
                    block=["Resolved"],
                )
            self.assertFalse(out_dir.exists())

            used = Path(tmp, "results", "used")
            used.mkdir(parents=True)
            Path(used, "x.txt").write_text("x")
            with self.assertRaisesRegex(ArtifactError, "nonempty run directory"):
                plan_annotation(request=_request(_fasta(tmp), used), block=["Resolved"])

    def test_locus_plan_fetches_the_window_into_the_predict_run(self) -> None:
        window = LocusWindow("hg38", "chr4", 100, 356)
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp, "results", "hg38_chr4_predict")
            run.mkdir(parents=True)
            with mock.patch(
                "ntv3_tools.annotate.fetch_ucsc_sequence", return_value="acgT" * 64
            ) as fetch:
                plan_path = plan_annotation(
                    request=_request(None, run),
                    block=["Resolved"],
                    window=window,
                )
            plan = json.loads(plan_path.read_text())
            fasta = Path(run, ANNOTATION_INPUT).read_text()

        fetch.assert_called_once_with(
            assembly="hg38", chrom="chr4", start=100, end=356, keep_soft_mask=True
        )
        self.assertTrue(fasta.startswith(">chr4\nacgT"))
        self.assertEqual(
            plan["locus"],
            {"assembly": "hg38", "chrom": "chr4", "start": 100, "end": 356},
        )
        self.assertEqual(plan["records"], [{"seqid": "chr4", "length": 256}])
        self.assertEqual(plan["fasta"], str(run.resolve() / ANNOTATION_INPUT))


class LocusWindowTests(unittest.TestCase):
    def test_window_is_centred_then_shifted_to_fit_the_chromosome(self) -> None:
        size = 10_000_000
        self.assertEqual(
            locus_window(center=5_000_000, chrom_size=size),
            (5_000_000 - _HALF, 5_000_000 + _HALF),
        )
        self.assertEqual(
            locus_window(center=10, chrom_size=size), (0, ANNOTATION_WINDOW_BP)
        )
        self.assertEqual(
            locus_window(center=9_999_000, chrom_size=size),
            (size - ANNOTATION_WINDOW_BP, size),
        )
        self.assertEqual(locus_window(center=100, chrom_size=500), (0, 500))

    def test_predict_run_window_reads_the_plan_and_refuses_a_second_one(
        self,
    ) -> None:
        plan = SimpleNamespace(
            assembly="hg38", chrom="chr4", start=6_000_000, end=6_004_096
        )
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            with self.assertRaisesRegex(RequestError, "has no plan.yaml"):
                predict_run_window(run_dir=run, chrom_size=lambda *_: 10**8)
            Path(run, "plan.yaml").write_text("{}\n")
            with mock.patch("ntv3_tools.annotate.read_plan", return_value=plan):
                window = predict_run_window(run_dir=run, chrom_size=lambda *_: 10**8)
            Path(run, "annotation_plan.json").write_text("{}\n")
            with self.assertRaisesRegex(RequestError, "already has annotation_plan"):
                predict_run_window(run_dir=run, chrom_size=lambda *_: 10**8)

        self.assertEqual(
            window,
            LocusWindow("hg38", "chr4", 6_002_048 - _HALF, 6_002_048 + _HALF),
        )


_WINDOW_GFF3 = (
    "##gff-version 3\n"
    "##sequence-region chr4 1 256\n"
    "##ntv3-decoder hmm_constrained\n"
    "chr4\tNTv3\tgene\t11\t90\t.\t+\t.\tID=g1\n"
    "chr4\tNTv3\tCDS\t11\t40\t.\t+\t0\tParent=g1\n"
    "chr4\tNTv3\tCDS\t61\t90\t.\t+\t0\tParent=g1\n"
)


class ShiftGff3Tests(unittest.TestCase):
    def test_shift_writes_chromosome_coordinates_the_report_reads(self) -> None:
        locus = LocusWindow("hg38", "chr4", 1000, 1256)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp, "window.gff3")
            source.write_text(_WINDOW_GFF3)
            out = Path(tmp, "annotation.gff3")
            shift_gff3(source=source, destination=out, window=locus)
            text = out.read_text()
            gff3 = read_gff3(out)

        self.assertIn("##sequence-region chr4 1001 1256", text)
        self.assertEqual(gff3.provenance["decoder"], "hmm_constrained")
        self.assertEqual(
            [(g["chrom"], g["start"], g["end"]) for g in gff3.genes],
            [("chr4", 1010, 1090)],
        )
        self.assertEqual(
            gff3.genes[0]["exons"],
            [{"start": 1010, "end": 1040}, {"start": 1060, "end": 1090}],
        )

    def test_shift_refuses_another_seqid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp, "window.gff3")
            source.write_text(_WINDOW_GFF3.replace("chr4\tNTv3", "chr9\tNTv3"))
            with self.assertRaisesRegex(ArtifactError, "names seqid chr9, not chr4"):
                shift_gff3(
                    source=source,
                    destination=Path(tmp, "out.gff3"),
                    window=LocusWindow("hg38", "chr4", 0, 9),
                )


class AnnotationApprovalTests(unittest.TestCase):
    def _planned(self, tmp: str) -> Path:
        out_dir = Path(tmp, "results", "run")
        plan_annotation(
            request=_request(_fasta(tmp), out_dir), block=["Resolved", "Approve"]
        )
        return out_dir

    def test_approve_then_require_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = self._planned(tmp)
            self.assertTrue(pending_annotation_plan(run_dir=run))
            plan = approve_annotation(run_dir=run)
            self.assertFalse(pending_annotation_plan(run_dir=run))
            self.assertTrue(Path(run, ANNOTATION_APPROVAL).is_file())
            require_annotation_approval(run_dir=run, plan=plan)
            with self.assertRaisesRegex(RequestError, "already approved"):
                approve_annotation(run_dir=run)

    def test_no_critic_at_approval_turns_the_critic_off(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = self._planned(tmp)
            self.assertTrue(
                json.loads(Path(run, "annotation_plan.json").read_text())["want_critic"]
            )
            plan = approve_annotation(run_dir=run, critic=False)
            require_annotation_approval(run_dir=run, plan=plan)

        self.assertFalse(plan["want_critic"])

    def test_unapproved_or_edited_plan_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = self._planned(tmp)
            plan_path = Path(run, "annotation_plan.json")
            with self.assertRaisesRegex(RequestError, "is not approved"):
                require_annotation_approval(
                    run_dir=run, plan=json.loads(plan_path.read_text())
                )
            approve_annotation(run_dir=run)
            plan = json.loads(plan_path.read_text())
            plan["decoder"] = "argmax"
            with self.assertRaisesRegex(RequestError, "changed after approval"):
                require_annotation_approval(run_dir=run, plan=plan)

            plan = json.loads(plan_path.read_text())
            Path(run, ANNOTATION_RESOLVED_BLOCK).write_text("Resolved\nedited\n")
            with self.assertRaisesRegex(RequestError, "changed after approval"):
                require_annotation_approval(run_dir=run, plan=plan)


class AnnotateCliTests(unittest.TestCase):
    def _args(self, **changes: object) -> argparse.Namespace:
        values: dict[str, object] = {
            "clade": "mammals",
            "fasta": None,
            "locus_run": None,
            "decoder": "hmm_constrained",
            "out": None,
            "resolve_only": True,
        }
        values.update(changes)
        return argparse.Namespace(**values)

    def test_cli_needs_exactly_one_input_and_a_run_for_a_fasta(self) -> None:
        for changes, message in (
            ({}, "What should be annotated\\?"),
            (
                {"fasta": "a.fa", "locus_run": "results/x"},
                "What should be annotated\\?",
            ),
            ({"fasta": "a.fa"}, "Pass --out results/<run>"),
        ):
            with (
                self.subTest(**{k: str(v) for k, v in changes.items()}),
                self.assertRaisesRegex(RequestError, message),
            ):
                _cmd_annotate(self._args(**changes))

    def test_cli_confines_output_to_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = self._args(fasta=str(Path(tmp, "a.fa")), out=tmp)
            with self.assertRaisesRegex(RequestError, "named directory under"):
                _cmd_annotate(args)


if __name__ == "__main__":
    unittest.main()
