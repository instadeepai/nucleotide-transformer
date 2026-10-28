"""JAX-environment tests for single-SNV Compare on token forwards."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import numpy as np

_JAX_AVAILABLE = importlib.util.find_spec("jax") is not None
if _JAX_AVAILABLE:
    from ntv3_tools.disrupt import Snv
    from ntv3_tools.disrupt_science import CompareInputs, run_compare
    from ntv3_tools.model import LoadedModel
    from ntv3_tools.plan import AnalysisPlan

_HEADS = ["ENCSR000AAA", "CNhs11250_P", "ENCSR000DDD"]


class _SequenceModel:
    """Fake loaded model: each track is the count of G/g tokens, per bin."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, list[int]]] = []
        self.config = SimpleNamespace(
            bigwigs_per_file_assembly={"hg38": _HEADS},
            keep_target_center_fraction=0.375,
        )


def _loaded(model: _SequenceModel) -> "LoadedModel":
    def forward(loaded: object, sequence: str, species: str, track_indices: list[int]):
        model.calls.append((sequence, list(track_indices)))
        n_out = int(0.375 * len(sequence))
        signal = float(sequence.count("G") + 2 * sequence.count("g"))
        return {
            "bigwig_tracks_logits": np.full(
                (1, n_out, len(track_indices)), signal, dtype=np.float32
            )
        }

    loaded = SimpleNamespace(model=model, config=model.config, checkpoint="NTv3_tracks")
    loaded.forward = forward
    return cast("LoadedModel", loaded)


def _plan(variant: "Snv") -> "AnalysisPlan":
    return AnalysisPlan(
        species="human",
        assembly="hg38",
        gene="HBE1",
        chrom="chr11",
        start=1_000,
        end=1_256,
        tissue="K562",
        tracks=["CAGE", "DNase"],
        device="cpu",
        window_bp=256,
        want_report=True,
        want_notebook=False,
        want_critic=False,
        ucsc_tracks=[],
        analyses=["predict", "compare"],
        variant=variant,
    )


@unittest.skipUnless(_JAX_AVAILABLE, "needs this repo's .venv (jax)")
class CompareTests(unittest.TestCase):
    def _run(self, null_n: int, reference: np.ndarray | None) -> tuple[Path, list]:
        from unittest.mock import patch

        sequence = "ACGTacgt" * 32
        # Offset 132 is "a" (lowercase): the alternate must stay soft-masked.
        variant = Snv(chrom="chr11", position=1_000 + 132 + 1, ref="A", alt="G")
        model = _SequenceModel()
        loaded = _loaded(model)
        tmp = tempfile.mkdtemp()
        Path(tmp, "track_resolve.json").write_text(
            json.dumps(
                {
                    "tracks": [
                        {"query": "CAGE", "model_head": "CNhs11250_P"},
                        {
                            "query": "DNase",
                            "gt_accession": "ENCSR000DDD",
                            "model_head": "ENCSR000DDD",
                        },
                    ]
                }
            )
        )
        with patch("ntv3_tools.disrupt_science.run_forward", loaded.forward):
            out = run_compare(
                inputs=CompareInputs(
                    plan=_plan(variant),
                    run_dir=Path(tmp),
                    loaded=loaded,
                    sequence=sequence,
                    reference=reference,
                ),
                variant=variant,
                null_n=null_n,
                seed=3,
            )
        return out, model.calls

    def test_reuses_the_reference_and_runs_one_alt_forward_on_approved_heads(
        self,
    ) -> None:
        reference = np.full((2, 96), 96.0, dtype=np.float32)

        out, calls = self._run(null_n=0, reference=reference)

        ((alt_sequence, indices),) = calls
        self.assertEqual(indices, [1, 2])
        self.assertEqual(alt_sequence[132], "g")
        delta = np.load(out / "delta_tracks.npy")
        self.assertEqual(delta.shape, (2, 96))
        self.assertTrue(np.allclose(delta, 2.0))
        effect = json.loads((out / "effect_percentile.json").read_text())
        self.assertIsNone(effect["percentile"])
        self.assertEqual(effect["per_track"][0]["model_head"], "CNhs11250_P")
        variant = json.loads((out / "variant.json").read_text())
        self.assertEqual(variant["prediction"]["bin_size"], 1)
        self.assertTrue(variant["in_prediction"])
        self.assertTrue((out.parent / "figures" / "disrupt_log2fc.png").is_file())

    def test_a_seeded_null_ranks_the_effect(self) -> None:
        out, calls = self._run(null_n=4, reference=None)

        self.assertEqual(len(calls), 1 + 1 + 4)
        effect = json.loads((out / "effect_percentile.json").read_text())
        self.assertEqual(effect["null_n"], 4)
        self.assertTrue(0.0 < effect["percentile"] <= 100.0)


if __name__ == "__main__":
    unittest.main()
