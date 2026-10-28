"""Tests for first-class optional analysis report payloads."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml
from ntv3_tools.errors import ArtifactError
from ntv3_tools.report import (
    REPORT_STAGE_LABELS,
    _analysis_stages,
    _explain,
    _load_disruption,
    _load_literature,
)


class AnalysisReportLoaderTests(unittest.TestCase):
    def test_marks_unselected_downstream_stage_as_not_requested(self) -> None:
        stages = _analysis_stages(
            requested_stages={"literature", "claims"},
            completed={"disruption": False, "literature": False, "claims": False},
        )

        by_id = {stage["id"]: stage for stage in stages}
        self.assertEqual(by_id["literature"]["status"], "pending")
        self.assertEqual(by_id["claims"]["status"], "pending")
        self.assertEqual(by_id["disruption"]["status"], "not_requested")

    def test_only_the_skills_analyses_are_stages(self) -> None:
        # Attention, interval shuffle, design, and motifs left the skill.
        self.assertEqual(
            set(REPORT_STAGE_LABELS), {"disruption", "literature", "claims"}
        )

    def test_explain_traces_compare_to_its_checkpoint_window_and_crop(self) -> None:
        plan = {
            "model": "NTv3_tracks",
            "species": "mouse",
            "assembly": "mm10",
            "device": "cpu",
            "chrom": "chr7",
            "start": 1_000,
            "end": 5_096,
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "disrupt").mkdir()
            (root / "plan.yaml").write_text(yaml.safe_dump(plan))
            (root / "track_resolve.json").write_text(
                json.dumps(
                    {"tracks": [{"query": "liver DNase", "gt_accession": "ENCSR1"}]}
                )
            )
            variant = {
                "chrom": "chr7",
                "position": 3_001,
                "ref": "C",
                "alt": "T",
                "checkpoint": "NTv3_tracks",
                "window": {"chrom": "chr7", "start": 1_000, "end": 5_096},
                "prediction": {
                    "chrom": "chr7",
                    "start": 2_280,
                    "end": 3_816,
                    "bin_size": 128,
                },
            }
            (root / "disrupt" / "variant.json").write_text(json.dumps(variant))

            explain = _explain(
                root=root,
                plan=plan,
                run_manifest={},
                prediction=("chr7", 2_280, 3_816, 128),
                track_rows=[
                    {
                        "query": "liver DNase",
                        "label": "liver DNase ENCSR1",
                        "modelHead": "ENCSR1",
                    }
                ],
                saliency_tracks=[],
            )

            self.assertEqual(explain["window"]["bp"], 4_096)
            self.assertEqual(explain["crop"]["bp"], 1_536)
            self.assertEqual(explain["nBins"], 12)
            self.assertEqual(explain["tracks"][0]["gtAccession"], "ENCSR1")
            self.assertIsNone(explain["manifestModel"])
            self.assertEqual(explain["letters"], [])
            self.assertEqual(
                explain["compare"],
                {
                    "label": "chr7:3001 C>T",
                    "checkpoint": "NTv3_tracks",
                    "window": {
                        "chrom": "chr7",
                        "start": 1_000,
                        "end": 5_096,
                        "bp": 4_096,
                    },
                    "crop": {
                        "chrom": "chr7",
                        "start": 2_280,
                        "end": 3_816,
                        "bp": 1_536,
                    },
                    "binSize": 128,
                },
            )
            self.assertEqual(
                explain["sources"],
                ["plan.yaml", "track_resolve.json", "disrupt/variant.json"],
            )

            del variant["prediction"]
            (root / "disrupt" / "variant.json").write_text(json.dumps(variant))
            with self.assertRaisesRegex(ArtifactError, "variant.json provenance"):
                _explain(
                    root=root,
                    plan=plan,
                    run_manifest={},
                    prediction=("chr7", 2_280, 3_816, 128),
                    track_rows=[],
                    saliency_tracks=[],
                )

    def test_loads_disruption_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "disrupt").mkdir()
            (root / "figures").mkdir()
            (root / "disrupt" / "variant.json").write_text(
                json.dumps(
                    {
                        "chrom": "chr11",
                        "position": 101,
                        "ref": "A",
                        "alt": "G",
                        "null_n": 100,
                        "seed": 7,
                        "prediction": {
                            "chrom": "chr11",
                            "start": 64,
                            "end": 128,
                            "bin_size": 1,
                        },
                    }
                )
            )
            np.save(root / "disrupt" / "delta_tracks.npy", np.full((1, 64), -0.5))
            (root / "disrupt" / "effect_percentile.json").write_text(
                json.dumps(
                    {
                        "signed_effect": -0.5,
                        "absolute_effect": 0.5,
                        "percentile": 91.0,
                        "null_n": 100,
                        "seed": 7,
                        "per_track": [
                            {
                                "track": "K562 CAGE",
                                "model_head": "head",
                                "signed_effect": -0.5,
                                "percentile": 91.0,
                            }
                        ],
                        "method": {"track_reduction": "signed abs-max"},
                    }
                )
            )
            (root / "figures" / "disrupt_log2fc.png").write_bytes(b"png")
            (root / "figures" / "disrupt_log2fc.pdf").write_bytes(b"pdf")

            payload = _load_disruption(root=root)

            assert payload is not None
            self.assertEqual(payload["variant"]["label"], "chr11:101 A>G")
            self.assertEqual(payload["percentile"], 91.0)
            self.assertEqual(payload["perTrack"][0]["modelHead"], "head")
            self.assertEqual(len(payload["perTrack"][0]["delta"]), 64)
            self.assertEqual(payload["prediction"]["binSize"], 1)
            self.assertEqual(payload["figure"]["png"], "figures/disrupt_log2fc.png")

            effect_path = root / "disrupt" / "effect_percentile.json"
            effect = json.loads(effect_path.read_text())
            effect.update(null_n=0, percentile=None)
            effect["per_track"][0]["percentile"] = None
            effect_path.write_text(json.dumps(effect))
            payload = _load_disruption(root=root)

            assert payload is not None
            self.assertIsNone(payload["percentile"])
            self.assertIsNone(payload["perTrack"][0]["percentile"])
            self.assertIn("no null was run", payload["caveat"])

            effect.update(percentile=50.0)
            effect_path.write_text(json.dumps(effect))
            with self.assertRaisesRegex(ArtifactError, "invalid values"):
                _load_disruption(root=root)

    def test_loads_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "literature").mkdir()
            (root / "literature" / "citations.bib").write_text("@article{x}\n")
            (root / "literature" / "search_manifest.json").write_text(
                json.dumps(
                    {
                        "records": [
                            {
                                "citekey": "Avsec2021Enformer",
                                "title": "Effective gene expression prediction",
                                "authors": ["Žiga Avsec"],
                                "year": 2021,
                                "source_url": "https://nature.example/paper",
                                "retrieval_date": "2026-09-08",
                                "search_query": "primary paper",
                                "relevance": "Sequence prediction context.",
                                "doi": "10.1038/example",
                            }
                        ]
                    }
                )
            )

            literature = _load_literature(root=root)

            assert literature is not None
            self.assertEqual(literature["records"][0]["year"], 2021)
            self.assertEqual(literature["bib"], "literature/citations.bib")


if __name__ == "__main__":
    unittest.main()
