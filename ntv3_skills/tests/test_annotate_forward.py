"""The annotation forward in the JAX Python, with ``ntv3-annotate`` stubbed.

Run in the repo ``.venv``: it imports JAX and the gene_annotation package. No
weights are downloaded and no model is loaded.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from ntv3_tools.annotate import AnnotateRequest, LocusWindow, plan_annotation
from ntv3_tools.approval import approve_annotation
from ntv3_tools.errors import ArtifactError, RequestError

_HAS_JAX_PACKAGE = all(
    importlib.util.find_spec(name) is not None
    for name in ("jax", "nucleotide_transformer_v3")
)
if _HAS_JAX_PACKAGE:
    import jax
    from ntv3_tools.annotate_science import (
        ANNOTATION_MANIFEST,
        run_annotation_forward,
    )

    from nucleotide_transformer_v3.gene_annotation.cli import (
        build_parser,
        compose_config,
    )
    from nucleotide_transformer_v3.gene_annotation.model_registry import (
        apply_model_config,
    )

_GPU = SimpleNamespace(device_kind="Fake GPU")


class _FakeAnnotate:
    """Records the argv and writes a one-gene GFF3 where ``--out`` says."""

    def __init__(self, seqid: str) -> None:
        self.seqid = seqid
        self.argv: list[str] | None = None

    def __call__(self, argv: list[str]) -> int:
        self.argv = argv
        Path(argv[argv.index("--out") + 1]).write_text(
            "##gff-version 3\n"
            f"{self.seqid}\tNTv3\tgene\t11\t90\t.\t+\t.\tID=g1\n"
            f"{self.seqid}\tNTv3\tCDS\t11\t90\t.\t+\t0\tParent=g1\n"
        )
        return 0


def _plan(tmp: str, window: LocusWindow | None = None, approve: bool = True) -> Path:
    run = Path(tmp, "results", "run")
    if window is None:
        fasta = Path(tmp, "a.fa")
        fasta.write_text(">chr1\n" + "ACGT" * 64 + "\n")
        request = AnnotateRequest("mammals", "hmm", str(run), str(fasta))
    else:
        run.mkdir(parents=True)
        request = AnnotateRequest("mammals", "hmm", str(run))
    with mock.patch(
        "ntv3_tools.annotate.fetch_ucsc_sequence", return_value="acgt" * 64
    ):
        plan = plan_annotation(request=request, block=["Resolved"], window=window)
    if approve:
        approve_annotation(run_dir=run)
    return plan


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
@mock.patch("ntv3_tools.annotate_science.require_gpu", return_value=_GPU)
class AnnotationForwardTests(unittest.TestCase):
    def test_unapproved_plan_loads_nothing(self, _gpu: mock.Mock) -> None:
        fake = _FakeAnnotate("chr1")
        with tempfile.TemporaryDirectory() as tmp:
            plan = _plan(tmp, approve=False)
            with self.assertRaisesRegex(RequestError, "is not approved"):
                run_annotation_forward(plan_path=plan, annotate_main=fake)
        self.assertIsNone(fake.argv)

    def test_fasta_run_writes_the_gff3_and_manifest(self, _gpu: mock.Mock) -> None:
        fake = _FakeAnnotate("chr1")
        with tempfile.TemporaryDirectory() as tmp:
            plan = _plan(tmp)
            out = run_annotation_forward(plan_path=plan, annotate_main=fake)
            manifest = json.loads(Path(plan.parent, ANNOTATION_MANIFEST).read_text())
            self.assertEqual(out, plan.parent / "annotation.gff3")
            with self.assertRaisesRegex(ArtifactError, "holds one annotation"):
                run_annotation_forward(plan_path=plan, annotate_main=fake)

        self.assertEqual((manifest["genes"], manifest["device"]), (1, "Fake GPU"))
        assert fake.argv is not None
        self.assertEqual(fake.argv[fake.argv.index("--model_cfg") + 1], "mammals")
        self.assertIn("decoder=hmm", fake.argv)

    def test_locus_run_writes_chromosome_coordinates(self, _gpu: mock.Mock) -> None:
        fake = _FakeAnnotate("chr4")
        with tempfile.TemporaryDirectory() as tmp:
            plan = _plan(tmp, window=LocusWindow("hg38", "chr4", 5000, 5256))
            out = run_annotation_forward(plan_path=plan, annotate_main=fake)
            rows = [
                line.split("\t")
                for line in out.read_text().splitlines()
                if not line.startswith("#")
            ]
            files = sorted(path.name for path in plan.parent.iterdir())

        self.assertEqual((rows[0][0], rows[0][3], rows[0][4]), ("chr4", "5011", "5090"))
        self.assertEqual(
            files,
            [
                "annotation.gff3",
                "annotation_approval.json",
                "annotation_input.fa",
                "annotation_manifest.json",
                "annotation_plan.json",
                "annotation_resolved.txt",
            ],
        )


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class GpuTests(unittest.TestCase):
    def test_no_gpu_is_refused_before_any_load(self) -> None:
        if any(device.platform == "gpu" for device in jax.devices()):
            self.skipTest("a GPU is visible")
        fake = _FakeAnnotate("chr1")
        with tempfile.TemporaryDirectory() as tmp:
            plan = _plan(tmp)
            with self.assertRaisesRegex(RequestError, "needs a GPU and JAX sees none"):
                run_annotation_forward(plan_path=plan, annotate_main=fake)
        self.assertIsNone(fake.argv)


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class PackageContractTests(unittest.TestCase):
    def test_the_forward_argv_composes_in_the_real_package(self) -> None:
        argv = ["--genome", "a.fa", "--model_cfg", "mammals", "--out", "a.gff3"]
        args = build_parser().parse_args([*argv, "--override", "decoder=argmax"])
        config = apply_model_config(compose_config(args))

        self.assertEqual(config.decoder, "argmax")
        self.assertEqual(config.sequence_length, 1_048_576)
        self.assertEqual(config.keep_target_center_fraction, 0.75)
        self.assertTrue(config.apply_biological_filters)
        self.assertEqual(config.min_cds_length, 201)


if __name__ == "__main__":
    unittest.main()
