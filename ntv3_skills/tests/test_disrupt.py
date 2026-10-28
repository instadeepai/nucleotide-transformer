"""Isolated tests for single-SNV disruption policy and artifacts."""

import argparse
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from ntv3_tools.cli import _cmd_disrupt
from ntv3_tools.disrupt import (
    DisruptArtifacts,
    NullSamplingRequest,
    Snv,
    apply_snv,
    empirical_absolute_percentile,
    generate_null_snvs,
    signed_absmax,
    validate_snv,
    write_disrupt_artifacts,
)
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.report import _trace_from_tree


class SnvValidationTests(unittest.TestCase):
    def test_validates_one_based_snv_and_reference(self) -> None:
        variant = Snv(chrom="chr1", position=102, ref="C", alt="T")

        offset = validate_snv(
            variant=variant,
            window_chrom="chr1",
            window_start=100,
            sequence="ACGT",
        )

        self.assertEqual(offset, 1)

    def test_rejects_invalid_or_mismatched_snv(self) -> None:
        cases = [
            (Snv("chr2", 102, "C", "T"), "chromosome"),
            (Snv("chr1", 100, "C", "T"), "outside"),
            (Snv("chr1", 102, "CC", "T"), "single"),
            (Snv("chr1", 102, "C", "C"), "different"),
            (Snv("chr1", 102, "N", "T"), "A/C/G/T"),
            (Snv("chr1", 102, "A", "T"), "reference"),
        ]
        for variant, message in cases:
            with (
                self.subTest(variant=variant),
                self.assertRaisesRegex(RequestError, message),
            ):
                validate_snv(
                    variant=variant,
                    window_chrom="chr1",
                    window_start=100,
                    sequence="ACGT",
                )

    def test_apply_snv_does_not_mutate_reference_sequence(self) -> None:
        sequence = "ACGT"

        alternate = apply_snv(sequence=sequence, offset=1, alt="T")

        self.assertEqual(sequence, "ACGT")
        self.assertEqual(alternate, "ATGT")

    def test_apply_snv_keeps_a_soft_masked_base_lowercase(self) -> None:
        self.assertEqual(apply_snv(sequence="ACgt", offset=2, alt="t"), "ACtt")
        self.assertEqual(apply_snv(sequence="ACgt", offset=1, alt="t"), "ATgt")


class NullSamplingTests(unittest.TestCase):
    def test_null_is_deterministic_and_excludes_query_and_non_acgt(self) -> None:
        query = Snv(chrom="chr1", position=101, ref="A", alt="C")

        request = NullSamplingRequest(
            sequence="ANCT",
            chrom="chr1",
            window_start=100,
            query=query,
            null_n=20,
            seed=7,
        )
        first = generate_null_snvs(request=request)
        second = generate_null_snvs(request=request)

        self.assertEqual(first, second)
        self.assertNotIn(query, first)
        self.assertTrue(all(snv.position != 102 for snv in first))
        self.assertTrue(all(snv.ref in "ACGT" and snv.alt in "ACGT" for snv in first))


class EffectPolicyTests(unittest.TestCase):
    def test_signed_absmax_preserves_sign_per_sample_and_track(self) -> None:
        delta = np.array(
            [
                [[1.0, -2.0], [-4.0, 1.0], [3.0, 0.5]],
                [[-2.0, 1.0], [1.0, 5.0], [0.0, -3.0]],
            ]
        )

        reduced = signed_absmax(delta=delta, axis=1)

        np.testing.assert_array_equal(
            reduced,
            np.array([[-4.0, -2.0], [-2.0, 5.0]]),
        )

    def test_percentile_uses_finite_sample_pseudocount(self) -> None:
        percentile = empirical_absolute_percentile(
            query_effect=2.0,
            null_effects=np.array([1.0, -2.0, 3.0]),
        )

        self.assertEqual(percentile, 75.0)


class ArtifactSafetyTests(unittest.TestCase):
    def test_existing_artifact_refuses_all_writes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            disrupt_dir = Path(tmp, "disrupt")
            disrupt_dir.mkdir()
            existing = disrupt_dir / "variant.json"
            existing.write_text('{"keep": true}\n')
            artifacts = DisruptArtifacts(
                variant={"chrom": "chr1"},
                ref_tracks=np.zeros((1, 2)),
                alt_tracks=np.ones((1, 2)),
                delta_tracks=np.ones((1, 2)),
                top_rows=[{"track": "track", "effect": 1.0}],
                effect_percentile={"percentile": 50.0},
                log2fc_png=b"png",
                log2fc_pdf=b"pdf",
            )

            with self.assertRaises(ArtifactError):
                write_disrupt_artifacts(run_dir=tmp, artifacts=artifacts)

            self.assertEqual(json.loads(existing.read_text()), {"keep": True})
            self.assertEqual(
                {path.name for path in disrupt_dir.iterdir()}, {"variant.json"}
            )

    def test_writes_log2fc_figure_with_disruption_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_disrupt_artifacts(
                run_dir=tmp,
                artifacts=DisruptArtifacts(
                    variant={"chrom": "chr1"},
                    ref_tracks=np.ones((1, 2)),
                    alt_tracks=np.ones((1, 2)),
                    delta_tracks=np.zeros((1, 2)),
                    top_rows=[{"track": "track", "effect": 0.0}],
                    effect_percentile={"percentile": 50.0},
                    log2fc_png=b"png",
                    log2fc_pdf=b"pdf",
                ),
            )

            self.assertEqual(
                Path(tmp, "figures", "disrupt_log2fc.png").read_bytes(),
                b"png",
            )
            self.assertEqual(
                Path(tmp, "figures", "disrupt_log2fc.pdf").read_bytes(),
                b"pdf",
            )

    def test_report_trace_includes_disrupt_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            disrupt_dir = Path(tmp, "disrupt")
            disrupt_dir.mkdir()
            (disrupt_dir / "effect_percentile.json").write_text("{}\n")

            trace = _trace_from_tree(root=Path(tmp))

            disrupt = next(step for step in trace if step["name"] == "disrupt")
            self.assertEqual(disrupt["status"], "ok")
            self.assertEqual(disrupt["artifact"], "disrupt/effect_percentile.json")


class DisruptDispatchTests(unittest.TestCase):
    def test_dispatch_refuses_an_unapproved_plan_and_leaves_it_alone(self) -> None:
        results_root = Path(__file__).resolve().parents[2] / "results"
        results_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as tmp:
            plan_path = Path(tmp, "plan.yaml")
            approved_plan = "model: approved-checkpoint\ntracks:\n- approved-head\n"
            plan_path.write_text(approved_plan)
            args = argparse.Namespace(
                plan=str(plan_path),
                chrom="chr1",
                position=101,
                ref="A",
                alt="C",
                null_n=10,
                seed=7,
            )

            with (
                patch("ntv3_tools.cli.importlib.util.find_spec", return_value=object()),
                patch(
                    "ntv3_tools.cli.subprocess.run",
                    side_effect=AssertionError("dispatched"),
                ),
                self.assertRaises(RequestError),
            ):
                _cmd_disrupt(args=args)

            self.assertEqual(plan_path.read_text(), approved_plan)


if __name__ == "__main__":
    unittest.main()
