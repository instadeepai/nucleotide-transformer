"""Predict loads NTv3_tracks through this repo and passes track_indices.

Runs in this repo's .venv (jax and nucleotide_transformer_v3). A fake model, or
a tiny randomly initialised model of the released class, stands in for the
weights, so nothing is downloaded and no forward of the real model runs.
"""

import hashlib
import importlib.util
import inspect
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

_SKILLS_ROOT = Path(__file__).resolve().parents[1]
if str(_SKILLS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SKILLS_ROOT))
_HAS_JAX_PACKAGE = all(
    importlib.util.find_spec(name) is not None
    for name in ("jax", "nucleotide_transformer_v3")
)
_FORBIDDEN = (
    "papers",
    "scripts",
)
_HG38_HEADS = ["ENCSR000AAA", "ENCSR000BBB", "CNhs11250_P", "ENCSR000DDD"]

if _HAS_JAX_PACKAGE:
    import jax
    import jax.numpy as jnp
    from flax import nnx
    from ntv3_tools.approval import approve_run, write_resolved_block
    from ntv3_tools.disrupt import Snv
    from ntv3_tools.errors import ArtifactError, RequestError
    from ntv3_tools.locus import LocusRequest
    from ntv3_tools.model import (
        LoadedModel,
        load_model,
        resolve_device,
        token_gradients,
    )
    from ntv3_tools.plan import EXCLUDED_ANALYSES, AnalysisPlan, write_plan
    from ntv3_tools.predict import TrackPredictionRequest, predict_tracks
    from ntv3_tools.report import write_viewer_bundle
    from ntv3_tools.saliency import _track_saliency
    from ntv3_tools.science_cli import _compare_ground_truth, run_prediction

    from nucleotide_transformer_v3.tokenizers import get_ntv3_tokenizer
    from nucleotide_transformer_v3.with_tracks_heads import (
        DiscreteConditionedNTv3WithTracksHeads,
        DiscreteConditionedNTv3WithTracksHeadsConfig,
    )


class _FakeTokenizer:
    """Records what it tokenizes; one id per base, soft mask untouched."""

    def __init__(self) -> None:
        self.seen: list[str] = []

    def batch_np_tokenize(self, sequences: list[str]) -> np.ndarray:
        self.seen.extend(sequences)
        return np.zeros((len(sequences), len(sequences[0])), dtype=np.int32)


class _FakeSpeciesTokenizer:
    def token_to_id(self, token: str) -> int:
        return {"human": 3, "mouse": 4}[token]


class _FakeTracksModel:
    """Stands in for DiscreteConditionedNTv3WithTracksHeads: records each call."""

    def __init__(self) -> None:
        self.calls: list[tuple[tuple, dict]] = []
        self.config = SimpleNamespace(
            bigwigs_per_file_assembly={"hg38": _HG38_HEADS, "mm10": ["ENCSR000MMM"]},
            keep_target_center_fraction=0.375,
        )

    def __call__(self, *args: object, **kwargs: object) -> dict[str, np.ndarray]:
        self.calls.append((args, kwargs))
        tokens = np.asarray(args[0])
        track_indices = kwargs.get("track_indices")
        if track_indices is None:
            raise AssertionError("Forward without track_indices builds every head.")
        n_out = int(0.375 * tokens.shape[1])
        columns = np.asarray(track_indices, dtype=np.float32)
        return {
            "bigwig_tracks_logits": np.broadcast_to(
                columns, (tokens.shape[0], n_out, columns.shape[0])
            ).copy()
        }


def _loaded() -> "LoadedModel":
    model = _FakeTracksModel()
    return LoadedModel(
        model=model,
        tokenizer=_FakeTokenizer(),
        condition_tokenizers=[_FakeSpeciesTokenizer()],
        config=model.config,
        device=jax.devices("cpu")[0],
        checkpoint="NTv3_tracks",
    )


def _tiny_loaded() -> "LoadedModel":
    """The released tracks class, shrunk and randomly initialised, in eval mode."""
    config = DiscreteConditionedNTv3WithTracksHeadsConfig(
        alphabet_size=16,
        token_embed_dim=8,
        conv_init_embed_dim=16,
        embed_dim=32,
        num_downsamples=2,
        num_layers=1,
        attention_heads=2,
        num_kv_heads=1,
        key_size=8,
        value_size=8,
        conditions_vocab_size=[8],
        conditions_names=["species"],
        species_embed_dim=4,
        enable_trunk_conditioning=False,
        bigwigs_per_file_assembly={"hg38": _HG38_HEADS[:3]},
        number_bigwig_tracks=3,
        bed_elements_names=[],
        bed_num_classes=3,
        keep_target_center_fraction=0.375,
        condition_bigwig_head=True,
    )
    model = DiscreteConditionedNTv3WithTracksHeads(config, rngs=nnx.Rngs(0))
    leaves, treedef = jax.tree_util.tree_flatten(nnx.state(model, nnx.Param))
    keys = jax.random.split(jax.random.key(1), len(leaves))
    noisy = [
        0.2 * jax.random.normal(key, leaf.shape, leaf.dtype)
        for key, leaf in zip(keys, leaves, strict=True)
    ]
    nnx.update(model, jax.tree_util.tree_unflatten(treedef, noisy))
    model.eval()
    return LoadedModel(
        model=model,
        tokenizer=get_ntv3_tokenizer(include_softmasked_tokens=True),
        condition_tokenizers=[_FakeSpeciesTokenizer()],
        config=config,
        device=jax.devices("cpu")[0],
        checkpoint="NTv3_tracks",
    )


def _hbe1_plan(**overrides: object) -> "AnalysisPlan":
    values: dict[str, object] = {
        "species": "human",
        "assembly": "hg38",
        "gene": "HBE1",
        "chrom": "chr11",
        "start": 5_266_144,
        "end": 5_270_240,
        "tissue": "K562",
        "tracks": ["CAGE"],
        "device": "cpu",
        "window_bp": 4_096,
        "want_report": True,
        "want_notebook": False,
        "want_critic": False,
        "ucsc_tracks": [],
        "analyses": ["predict"],
    }
    values.update(overrides)
    return AnalysisPlan(**values)


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class ImportPathTests(unittest.TestCase):
    def test_predict_entry_imports_only_this_repo(self) -> None:
        script = f"""
import importlib.abc, sys
class Refuse(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in {_FORBIDDEN!r}:
            raise ImportError("refused " + name)
sys.meta_path.insert(0, Refuse())
sys.path.insert(0, {str(_SKILLS_ROOT)!r})
import ntv3_tools.model, ntv3_tools.predict, ntv3_tools.science_cli
loaded = sorted(m for m in sys.modules if m.split(".")[0] in {_FORBIDDEN!r})
print("forbidden:", loaded)
print("loader:", ntv3_tools.model.load_tracks_model.__module__)
"""
        env = {"PATH": "/usr/bin:/bin"}
        completed = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            env=env,
            check=False,
            timeout=300,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("forbidden: []", completed.stdout)
        self.assertIn("loader: nucleotide_transformer_v3.pretrained", completed.stdout)


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class LoadTests(unittest.TestCase):
    def test_load_model_calls_load_tracks_model_with_the_released_name(self) -> None:
        model = _FakeTracksModel()
        with (
            patch(
                "ntv3_tools.model.load_tracks_model",
                return_value=(model, _FakeTokenizer(), [_FakeSpeciesTokenizer()]),
            ) as loader,
            patch(
                "ntv3_tools.model.nnx",
                SimpleNamespace(state=lambda module: {}, update=lambda *args: None),
            ),
        ):
            loaded = load_model(model_name="InstaDeepAI/NTv3_tracks", device="cpu")

        loader.assert_called_once_with("NTv3_tracks")
        self.assertIs(loaded.model, model)
        self.assertEqual(loaded.checkpoint, "NTv3_tracks")
        self.assertEqual(loaded.device.platform, "cpu")

    def test_old_internal_key_is_refused_before_any_load(self) -> None:
        with (
            patch("ntv3_tools.model.load_tracks_model") as loader,
            self.assertRaisesRegex(KeyError, "NTv3_tracks"),
        ):
            load_model(model_name="ntv3_106M_7downsample_post_trained_131kb")
        loader.assert_not_called()


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class PredictForwardTests(unittest.TestCase):
    def test_predict_passes_only_the_requested_head_columns(self) -> None:
        loaded = _loaded()
        sequence = ("ACGTacgt" * 512)[:4096]
        prediction = predict_tracks(
            request=TrackPredictionRequest(
                model_name="NTv3_tracks",
                locus=LocusRequest(
                    chrom="chr11", start=5_266_144, end=5_270_240, species="human"
                ),
                track_ids=["CNhs11250_P", "ENCSR000AAA"],
                device="cpu",
                loaded=loaded,
                sequence=sequence,
            )
        )

        ((args, kwargs),) = loaded.model.calls
        inspect.signature(DiscreteConditionedNTv3WithTracksHeads.__call__).bind(
            None, *args, **kwargs
        )
        self.assertEqual(np.asarray(kwargs["track_indices"]).tolist(), [2, 0])
        self.assertTrue(kwargs["skip_lm_head"])
        tokens, conditions = args
        self.assertEqual(tokens.shape, (1, 4096))
        self.assertEqual([np.asarray(ids).tolist() for ids in conditions], [[3]])
        self.assertEqual(loaded.tokenizer.seen, [sequence])
        self.assertEqual(prediction.tracks.shape, (2, 1536))
        self.assertEqual(prediction.tracks[:, 0].tolist(), [2.0, 0.0])
        self.assertEqual(prediction.accessions, ["CNhs11250_P", "ENCSR000AAA"])
        self.assertEqual(
            (prediction.prediction_start, prediction.prediction_end),
            (5_267_424, 5_268_960),
        )

    def test_a_head_the_model_lacks_fails_before_the_forward(self) -> None:
        loaded = _loaded()
        with self.assertRaisesRegex(ValueError, "ENCSR999ZZZ"):
            predict_tracks(
                request=TrackPredictionRequest(
                    model_name="NTv3_tracks",
                    locus=LocusRequest(
                        chrom="chr11", start=0, end=4096, species="human"
                    ),
                    track_ids=["ENCSR999ZZZ"],
                    loaded=loaded,
                    sequence="A" * 4096,
                )
            )
        self.assertEqual(loaded.model.calls, [])


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class AttributionRefusalTests(unittest.TestCase):
    def test_a_plan_edited_to_name_an_excluded_analysis_fetches_and_loads_nothing(
        self,
    ) -> None:
        for method in EXCLUDED_ANALYSES:
            with self.subTest(method=method), tempfile.TemporaryDirectory() as tmp:
                plan_path = write_plan(plan=_hbe1_plan(), out_dir=tmp)
                plan_path.write_text(
                    plan_path.read_text().replace(
                        "analyses:\n- predict\n", f"analyses:\n- predict\n- {method}\n"
                    )
                )
                self.assertIn(f"- {method}\n", plan_path.read_text())
                Path(tmp, "track_resolve.json").write_text(
                    json.dumps(
                        {"tracks": [{"query": "CAGE", "model_head": "CNhs11250_P"}]}
                    )
                )
                with (
                    patch(
                        "ntv3_tools.science_cli.fetch_ucsc_sequence",
                        side_effect=AssertionError("fetched"),
                    ),
                    patch(
                        "ntv3_tools.science_cli.load_model",
                        side_effect=AssertionError("loaded"),
                    ),
                    self.assertRaisesRegex(
                        RequestError, f"{method} is outside this skill"
                    ),
                ):
                    run_prediction(plan_path=str(plan_path))


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class ApprovalGateTests(unittest.TestCase):
    def test_an_unapproved_plan_fetches_and_loads_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            plan_path = write_plan(plan=_hbe1_plan(), out_dir=tmp)
            Path(tmp, "track_resolve.json").write_text(
                json.dumps({"tracks": [{"query": "CAGE", "model_head": "CNhs11250_P"}]})
            )
            with (
                patch(
                    "ntv3_tools.science_cli.fetch_ucsc_sequence",
                    side_effect=AssertionError("fetched"),
                ),
                patch(
                    "ntv3_tools.science_cli.load_model",
                    side_effect=AssertionError("loaded"),
                ),
                self.assertRaisesRegex(RequestError, "is not approved"),
            ):
                run_prediction(plan_path=str(plan_path))


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class SaliencyTests(unittest.TestCase):
    def test_token_gradients_match_a_one_hot_input_and_leave_the_model_alone(
        self,
    ) -> None:
        loaded = _tiny_loaded()
        sequence = "ACGTacgtNnTTGGCCAA" * 7 + "AC"
        before = jax.tree_util.tree_map(np.asarray, nnx.state(loaded.model))

        gradients = token_gradients(
            loaded=loaded, sequence=sequence, species="human", track_index=1
        )

        tokens = jnp.asarray(loaded.tokenizer.batch_np_tokenize([sequence]))
        graphdef, state = nnx.split(loaded.model)

        def from_one_hot(one_hot: jax.Array) -> jax.Array:
            model = nnx.merge(graphdef, state)
            table = model.embed_layer.embedding.value
            model.embed_layer = lambda ids: one_hot @ table
            out = model(
                tokens,
                [jnp.asarray([3])],
                track_indices=jnp.asarray([1]),
                skip_bed_head=True,
                skip_lm_head=True,
            )
            return out["bigwig_tracks_logits"][0, :, 0].mean()

        expected = jax.grad(from_one_hot)(jax.nn.one_hot(tokens, 16))[0]
        after = jax.tree_util.tree_map(np.asarray, nnx.state(loaded.model))
        self.assertEqual(gradients.shape, (128, 16))
        np.testing.assert_allclose(gradients, expected, rtol=1e-5, atol=1e-7)
        self.assertGreater(np.abs(gradients).max(), 0.0)
        self.assertTrue(
            jax.tree_util.tree_all(
                jax.tree_util.tree_map(np.array_equal, before, after)
            )
        )

    def test_saliency_keeps_the_reference_base_in_its_own_case(self) -> None:
        loaded = _tiny_loaded()
        sequence = "ACGTacgtNn" * 12 + "ACGTACGT"

        scores = _track_saliency(
            loaded=loaded, sequence=sequence, species="human", track_index=0
        )

        gradients = token_gradients(
            loaded=loaded, sequence=sequence, species="human", track_index=0
        )
        ids = loaded.tokenizer.token_to_id
        self.assertEqual(scores.shape, (128, 4))
        for position, column in ((0, 0), (5, 1)):
            base = sequence[position]
            row = gradients[
                position, [ids(b if base.isupper() else b.lower()) for b in "ACGT"]
            ]
            self.assertAlmostEqual(
                float(scores[position, column]),
                float(row[column] - row.mean()),
                places=6,
            )
            self.assertEqual(np.count_nonzero(scores[position]), 1)
        self.assertFalse(scores[8].any() or scores[9].any())


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class DeviceTests(unittest.TestCase):
    def _no_cuda(self, kind: str | None = None) -> list:
        if kind == "cuda":
            raise RuntimeError("no cuda backend")
        return [self.cpu]

    def setUp(self) -> None:
        self.cpu = jax.devices("cpu")[0]

    def test_cuda_with_no_cuda_device_refuses_instead_of_using_cpu(self) -> None:
        with patch("ntv3_tools.model.jax.devices", side_effect=self._no_cuda):
            for device in ("cuda", "gpu"):
                with (
                    self.subTest(device=device),
                    self.assertRaisesRegex(RequestError, "sees no CUDA device"),
                ):
                    resolve_device(device=device)

    def test_auto_still_runs_on_cpu(self) -> None:
        with patch("ntv3_tools.model.jax.devices", side_effect=self._no_cuda):
            self.assertEqual(resolve_device(device="auto"), self.cpu)


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class DefaultPathTests(unittest.TestCase):
    def test_predict_saliency_and_compare_share_one_load_and_reach_the_browser(
        self,
    ) -> None:
        sequence = ("ACGTacgt" * 512)[:4096]
        # Offset 2000 is inside the 1,536 bp crop and holds "A".
        variant = Snv(chrom="chr11", position=5_266_144 + 2_000 + 1, ref="A", alt="G")
        loaded = _tiny_loaded()
        with tempfile.TemporaryDirectory() as tmp:
            plan_path = write_plan(
                plan=_hbe1_plan(
                    analyses=["predict", "saliency", "compare"], variant=variant
                ),
                out_dir=tmp,
            )
            Path(tmp, "track_resolve.json").write_text(
                json.dumps({"tracks": [{"query": "CAGE", "model_head": "CNhs11250_P"}]})
            )
            # The user approves the Resolved block before anything loads.
            write_resolved_block(run_dir=tmp, lines=["Resolved", "test block"])
            approve_run(run_dir=tmp)
            with (
                patch(
                    "ntv3_tools.science_cli.fetch_ucsc_sequence", return_value=sequence
                ),
                patch("ntv3_tools.science_cli.load_model", return_value=loaded) as load,
                patch(
                    "ntv3_tools.science_cli.bigwig_reader_problem",
                    return_value="test",
                ),
            ):
                run_prediction(plan_path=str(plan_path))
            write_viewer_bundle(run_dir=tmp)

            root = Path(tmp)
            load.assert_called_once()
            manifest = json.loads((root / "run_manifest.json").read_text())
            self.assertFalse(manifest["ground_truth_compared"])
            self.assertEqual(
                manifest["ground_truth_note"], "no ENCODE experiment for CAGE"
            )
            self.assertEqual(
                (manifest["device"], manifest["device_requested"]), ("cpu", "cpu")
            )
            self.assertEqual(
                manifest["input_sha256"],
                hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
            )
            self.assertEqual(
                np.load(root / "interpret" / "CNhs11250_P.npy").shape, (4096, 4)
            )
            self.assertEqual(
                np.load(root / "disrupt" / "delta_tracks.npy").shape, (1, 1536)
            )
            effect = json.loads(
                (root / "disrupt" / "effect_percentile.json").read_text()
            )
            self.assertEqual(effect["null_n"], 0)
            self.assertIsNone(effect["percentile"])
            self.assertTrue(
                json.loads((root / "disrupt" / "variant.json").read_text())[
                    "in_prediction"
                ]
            )
            html = (root / "browser.html").read_text()
            tag = '<script type="application/json" id="ntv3-run">'
            data = json.loads(html.split(tag, 1)[1].split("</script>", 1)[0])
            labels = [lane["label"] for lane in data["lanes"]]
            # A typed assay is labelled with the head it resolved to.
            self.assertEqual(
                labels,
                [
                    "Predicted CAGE CNhs11250_P",
                    "Saliency CAGE CNhs11250_P",
                    "Alt \u2212 ref CAGE CNhs11250_P",
                ],
            )
            self.assertEqual(
                data["marks"], [{"pos": 5_268_144, "label": "chr11:5268145 A>G"}]
            )
            self.assertEqual(
                data["compare"]["rows"][0][:2], ["CAGE CNhs11250_P", "CNhs11250_P"]
            )
            self.assertEqual(data["compare"]["rows"][0][3], "—")

    def test_a_ground_truth_failure_keeps_the_forward(self) -> None:
        sequence = ("ACGTacgt" * 512)[:4096]
        loaded = _tiny_loaded()
        with tempfile.TemporaryDirectory() as tmp:
            plan_path = write_plan(plan=_hbe1_plan(), out_dir=tmp)
            root = Path(tmp)
            (root / "track_resolve.json").write_text(
                json.dumps(
                    {
                        "tracks": [
                            {
                                "query": "CAGE",
                                "model_head": "CNhs11250_P",
                                "gt_accession": "ENCSR000CIO",
                            }
                        ]
                    }
                )
            )
            (root / "gt").mkdir()
            (root / "gt" / "manifest.json").write_text(
                json.dumps(
                    [
                        {
                            "accession": "ENCSR000CIO",
                            "url": "https://example.org/ENCFF000AAA.bigWig",
                        }
                    ]
                )
            )
            write_resolved_block(run_dir=tmp, lines=["Resolved", "test block"])
            approve_run(run_dir=tmp)
            with (
                patch(
                    "ntv3_tools.science_cli.fetch_ucsc_sequence", return_value=sequence
                ),
                patch("ntv3_tools.science_cli.load_model", return_value=loaded),
                patch(
                    "ntv3_tools.science_cli.bigwig_reader_problem", return_value=None
                ),
                patch(
                    "ntv3_tools.science_cli.write_ground_truth_comparison",
                    side_effect=RequestError("boom."),
                ),
            ):
                run_prediction(plan_path=str(plan_path))

            self.assertEqual(np.load(root / "tracks" / "tracks.npy").shape, (1, 1536))
            manifest = json.loads((root / "run_manifest.json").read_text())
            self.assertFalse(manifest["ground_truth_compared"])
            self.assertEqual(manifest["ground_truth_note"], "boom")
            self.assertFalse((root / "gt" / "observed.npy").exists())

    def test_an_existing_observed_array_is_a_note_not_an_exit(self) -> None:
        sequence = ("ACGTacgt" * 512)[:4096]
        loaded = _tiny_loaded()
        with tempfile.TemporaryDirectory() as tmp:
            plan_path = write_plan(plan=_hbe1_plan(), out_dir=tmp)
            root = Path(tmp)
            (root / "track_resolve.json").write_text(
                json.dumps(
                    {
                        "tracks": [
                            {
                                "query": "CAGE",
                                "model_head": "CNhs11250_P",
                                "gt_accession": "ENCSR000CIO",
                            }
                        ]
                    }
                )
            )
            (root / "gt").mkdir()
            (root / "gt" / "manifest.json").write_text(
                json.dumps(
                    [
                        {
                            "accession": "ENCSR000CIO",
                            "url": "https://example.org/ENCFF000AAA.bigWig",
                        }
                    ]
                )
            )
            np.save(root / "gt" / "observed.npy", np.zeros((1, 1536)))
            write_resolved_block(run_dir=tmp, lines=["Resolved", "test block"])
            approve_run(run_dir=tmp)
            with (
                patch(
                    "ntv3_tools.science_cli.fetch_ucsc_sequence", return_value=sequence
                ),
                patch("ntv3_tools.science_cli.load_model", return_value=loaded),
                patch(
                    "ntv3_tools.science_cli.bigwig_reader_problem", return_value=None
                ),
            ):
                run_prediction(plan_path=str(plan_path))

            self.assertEqual(np.load(root / "tracks" / "tracks.npy").shape, (1, 1536))
            manifest = json.loads((root / "run_manifest.json").read_text())
            self.assertFalse(manifest["ground_truth_compared"])
            self.assertIn("Refusing to overwrite", manifest["ground_truth_note"])


@unittest.skipUnless(_HAS_JAX_PACKAGE, "needs this repo's .venv (jax)")
class CompareGroundTruthTests(unittest.TestCase):
    def _compare(self, rows: list[dict], error: Exception | None = None) -> tuple:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "gt").mkdir()
            (root / "gt" / "manifest.json").write_text("[]")
            with (
                patch(
                    "ntv3_tools.science_cli.bigwig_reader_problem", return_value=None
                ),
                patch(
                    "ntv3_tools.science_cli.write_ground_truth_comparison",
                    side_effect=error,
                ) as write,
            ):
                result = _compare_ground_truth(
                    run_dir=root,
                    track_rows=rows,
                    comparison=SimpleNamespace(),
                    predicted=np.zeros((len(rows), 2)),
                )
            return result, write.call_count

    def test_one_track_without_an_accession_does_not_skip_the_others(self) -> None:
        rows = [
            {"query": "CAGE", "gt_accession": "ENCSR000CIO"},
            {"query": "DNase"},
        ]
        (compared, note), calls = self._compare(rows=rows)
        self.assertTrue(compared)
        self.assertEqual(note, "no ENCODE experiment for DNase")
        self.assertEqual(calls, 1)

    def test_no_track_with_an_accession_compares_nothing(self) -> None:
        (compared, note), calls = self._compare(rows=[{"query": "CAGE"}])
        self.assertFalse(compared)
        self.assertEqual(note, "no ENCODE experiment for CAGE")
        self.assertEqual(calls, 0)

    def test_every_ground_truth_error_becomes_a_note(self) -> None:
        rows = [{"query": "CAGE", "gt_accession": "ENCSR000CIO"}]
        for error in (
            ArtifactError("Refusing to overwrite gt/observed.npy."),
            ValueError("Predictions must be finite."),
            json.JSONDecodeError("Expecting value", "", 0),
        ):
            with self.subTest(error=type(error).__name__):
                (compared, note), _ = self._compare(rows=rows, error=error)
                self.assertFalse(compared)
                self.assertTrue(note)


if __name__ == "__main__":
    unittest.main()
