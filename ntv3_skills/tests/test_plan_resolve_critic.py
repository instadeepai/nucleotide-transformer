"""Cheap unit tests for plan / resolve / critic / gt. No GPU."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ntv3_tools.critic import Claim, audit_claims
from ntv3_tools.errors import ArtifactError, ExternalServiceError, RequestError
from ntv3_tools.gt import pick_encode_bigwig, resolve_gt, search_encode_files
from ntv3_tools.lookups import LookupRequest, run_lookup
from ntv3_tools.plan import AnalysisPlan, validate_plan, write_plan
from ntv3_tools.registry import resolve_track_ids
from ntv3_tools.report import write_viewer_bundle
from ntv3_tools.resolve import resolve_tracks
from ntv3_tools.ucsc import fetch_track


class PlanTests(unittest.TestCase):
    def test_missing_species_raises(self) -> None:
        plan = AnalysisPlan(
            species="",
            assembly="hg38",
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            tissue="K562",
            tracks=["K562 CAGE"],
            device="cpu",
            window_bp=16384,
            want_report=True,
            want_notebook=True,
            want_critic=True,
            ucsc_tracks=[],
            analyses=[],
        )
        with self.assertRaises(RequestError):
            validate_plan(plan=plan)

    def test_write_plan(self) -> None:
        plan = AnalysisPlan(
            species="human",
            assembly="hg38",
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            tissue="K562",
            tracks=["K562 CAGE"],
            device="cpu",
            window_bp=16384,
            want_report=True,
            want_notebook=False,
            want_critic=True,
            ucsc_tracks=["clinvarMain"],
            analyses=["predict"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = write_plan(plan=plan, out_dir=tmp)
            self.assertTrue(path.exists())

    def test_write_plan_refuses_nonempty_run_directory(self) -> None:
        plan = AnalysisPlan(
            species="human",
            assembly="hg38",
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            tissue="K562",
            tracks=["K562 DNase", "K562 CAGE", "K562 GATA1"],
            device="cpu",
            window_bp=16384,
            want_report=True,
            want_notebook=False,
            want_critic=True,
            ucsc_tracks=["knownGene"],
            analyses=["predict", "saliency"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "existing.txt").write_text("preserve me\n")

            with self.assertRaisesRegex(
                ArtifactError, "Refusing to use nonempty run directory"
            ):
                write_plan(plan=plan, out_dir=tmp)

            self.assertEqual(Path(tmp, "existing.txt").read_text(), "preserve me\n")
            self.assertFalse(Path(tmp, "plan.yaml").exists())

    def test_attribution_track_must_be_selected_for_prediction(self) -> None:
        plan = AnalysisPlan(
            species="human",
            assembly="hg38",
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            tissue="K562",
            tracks=["K562 DNase", "K562 CAGE"],
            device="cpu",
            window_bp=16384,
            want_report=True,
            want_notebook=False,
            want_critic=True,
            ucsc_tracks=["knownGene"],
            analyses=["predict", "saliency"],
            attribution_tracks=["K562 GATA1"],
        )

        with self.assertRaisesRegex(
            RequestError, "Attribution tracks must be selected prediction tracks"
        ):
            validate_plan(plan=plan)

    def test_writes_explicit_attribution_track(self) -> None:
        plan = AnalysisPlan(
            species="human",
            assembly="hg38",
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            tissue="K562",
            tracks=["K562 DNase", "K562 CAGE", "K562 GATA1"],
            device="cpu",
            window_bp=16384,
            want_report=True,
            want_notebook=False,
            want_critic=True,
            ucsc_tracks=["knownGene"],
            analyses=["predict", "saliency"],
            attribution_tracks=["K562 GATA1"],
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = write_plan(plan=plan, out_dir=tmp)

            payload = path.read_text()

        self.assertIn("attribution_tracks:\n- K562 GATA1\n", payload)


class ResolveTests(unittest.TestCase):
    def test_k562_cage_has_dual_ids(self) -> None:
        resolved = resolve_tracks(
            species="human",
            tissue="K562",
            names=["K562 CAGE", "K562 DNase"],
        )
        self.assertEqual(resolved.assembly, "hg38")
        cage = resolved.tracks[0]
        self.assertEqual(cage.gt_accession, "ENCSR000CIO")
        self.assertEqual(cage.model_head, "CNhs11250_P")
        self.assertEqual(resolved.tracks[1].gt_accession, "ENCSR921NMD")

    def test_cnhs_is_model_head_not_gt(self) -> None:
        ids = resolve_track_ids(name_or_accession="CNhs11250_P")
        self.assertEqual(ids.model_head, "CNhs11250_P")
        self.assertIsNone(ids.gt_accession)

    def test_unknown_name_uses_encode_search(self) -> None:
        resolved = resolve_tracks(
            species="human",
            tissue="HepG2",
            names=["DNase"],
            encode_search=lambda tissue, assay: "ENCSR999FAKE",
        )
        self.assertEqual(resolved.tracks[0].gt_accession, "ENCSR999FAKE")

    def test_unset_species_raises(self) -> None:
        with self.assertRaises(RequestError):
            resolve_tracks(species="", tissue="K562", names=["K562 CAGE"])

    def test_hepg2_vs_k562_track_raises(self) -> None:
        with self.assertRaises(RequestError):
            resolve_tracks(
                species="human",
                tissue="HepG2",
                names=["K562 CAGE"],
            )

    def test_unknown_track_raises(self) -> None:
        with self.assertRaises(KeyError):
            resolve_tracks(
                species="human",
                tissue="HepG2",
                names=["POLR2A"],
            )

    def test_mouse_is_mm10(self) -> None:
        resolved = resolve_tracks(
            species="mouse",
            tissue="mESC",
            names=["ENCSR000CIO"],
        )
        self.assertEqual(resolved.assembly, "mm10")


class CriticTests(unittest.TestCase):
    def test_missing_evidence_is_suppress(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            decisions = audit_claims(
                claims=[
                    Claim(
                        claim_id="loop",
                        text="This is a proven enhancer-promoter loop.",
                        evidence_paths=["missing.tsv"],
                        citekey=None,
                        kind="other",
                    )
                ],
                run_dir=tmp,
            )
        self.assertEqual(decisions[0].verdict, "SUPPRESS")

    def test_overclaim_is_suppress_even_with_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "ok.txt").write_text("x")
            decisions = audit_claims(
                claims=[
                    Claim(
                        claim_id="loop",
                        text="This is a proven enhancer-promoter loop.",
                        evidence_paths=["ok.txt"],
                        citekey=None,
                        kind="other",
                    )
                ],
                run_dir=tmp,
            )
        self.assertEqual(decisions[0].verdict, "SUPPRESS")

    def test_explicit_causal_kind_is_suppressed_with_present_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "motif_hits.tsv").write_text("matrix_id\n")
            decisions = audit_claims(
                claims=[
                    Claim(
                        claim_id="causal-activation",
                        text=(
                            "GATA1 binding at the detected motifs causally activates "
                            "HBE1 expression in K562 cells."
                        ),
                        evidence_paths=["motif_hits.tsv"],
                        citekey=None,
                        kind="causal_biological",
                    )
                ],
                run_dir=tmp,
            )
        self.assertEqual(decisions[0].verdict, "SUPPRESS")
        self.assertIn("causal", decisions[0].reason.lower())

    def _audit_one(self, run_dir: str, paths: list[str]) -> str:
        return audit_claims(
            claims=[
                Claim(
                    claim_id="c",
                    text="NTv3 predicted the track.",
                    evidence_paths=paths,
                    citekey=None,
                    kind="model",
                )
            ],
            run_dir=run_dir,
        )[0].verdict

    def test_claim_with_no_evidence_path_is_suppressed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self._audit_one(tmp, []), "SUPPRESS")

    def test_evidence_outside_the_run_is_suppressed(self) -> None:
        with (
            tempfile.TemporaryDirectory() as tmp,
            tempfile.TemporaryDirectory() as other,
        ):
            run = Path(tmp, "run")
            run.mkdir()
            outside = Path(other, "outside.txt")
            outside.write_text("x")
            Path(tmp, "sibling.txt").write_text("x")
            (run / "link.txt").symlink_to(outside)
            (run / "ok.txt").write_text("x")
            for paths in ([str(outside)], ["../sibling.txt"], ["link.txt"]):
                with self.subTest(paths=paths):
                    self.assertEqual(self._audit_one(str(run), paths), "SUPPRESS")
            self.assertEqual(self._audit_one(str(run), ["ok.txt"]), "KEEP")

    def test_the_run_directory_or_a_subdirectory_is_not_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "tracks").mkdir()
            for paths in ([""], ["."], ["tracks"], ["tracks/"]):
                with self.subTest(paths=paths):
                    self.assertEqual(self._audit_one(tmp, paths), "SUPPRESS")

    def test_lit_without_citekey_is_suppress(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "ok.txt").write_text("x")
            decisions = audit_claims(
                claims=[
                    Claim(
                        claim_id="lit",
                        text="A paper supports this.",
                        evidence_paths=["ok.txt"],
                        citekey=None,
                        kind="literature",
                    )
                ],
                run_dir=tmp,
                bib_path=str(Path(tmp, "empty.bib")),
            )
        self.assertEqual(decisions[0].verdict, "SUPPRESS")


class GtTests(unittest.TestCase):
    def test_pick_prefers_fold_change_on_matching_assembly(self) -> None:
        files = [
            {
                "accession": "ENCFFAAA",
                "assembly": "GRCh38",
                "output_type": "signal",
                "href": "/aaa.bigWig",
            },
            {
                "accession": "ENCFFBBB",
                "assembly": "GRCh38",
                "output_type": "fold change over control",
                "href": "/bbb.bigWig",
            },
        ]
        chosen = pick_encode_bigwig(files=files, assembly="hg38")
        self.assertEqual(chosen["accession"], "ENCFFBBB")

    def test_cnhs_is_not_encode_gt(self) -> None:
        with self.assertRaises(KeyError):
            resolve_gt(assembly="hg38", accession="CNhs11250_P")

    def test_search_uses_experiments_path(self) -> None:
        seen: list[tuple[str, dict[str, str] | None]] = []

        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            seen.append((url, params))
            return {"@graph": [], "files": []}

        search_encode_files(accession="ENCSR000CIO", get_json=fake_get)
        self.assertEqual(
            seen[0][1]["dataset"],
            "/experiments/ENCSR000CIO/",
        )

    def test_search_falls_back_when_encode_search_endpoint_is_unavailable(
        self,
    ) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            if url.endswith("/search/"):
                raise ExternalServiceError("ENCODE search returned 404")
            return {
                "files": [
                    {
                        "accession": "ENCFFFALLBACK",
                        "assembly": "GRCh38",
                        "file_format": "bigWig",
                        "status": "released",
                    }
                ]
            }

        files = search_encode_files(accession="ENCSR000CIO", get_json=fake_get)

        self.assertEqual(files[0]["accession"], "ENCFFFALLBACK")

    def test_search_uses_archived_bigwig_when_experiment_has_no_released_bigwig(
        self,
    ) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            if url.endswith("/search/"):
                raise ExternalServiceError("ENCODE search returned 404")
            return {
                "files": [
                    {
                        "accession": "ENCFFARCHIVED",
                        "assembly": "GRCh38",
                        "file_format": "bigWig",
                        "status": "archived",
                    }
                ]
            }

        files = search_encode_files(accession="ENCSR921NMD", get_json=fake_get)

        self.assertEqual(files[0]["accession"], "ENCFFARCHIVED")

    def test_resolve_gt_uses_injected_search(self) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            return {
                "@graph": [
                    {
                        "accession": "ENCFF079PRE",
                        "assembly": "GRCh38",
                        "output_type": "signal",
                        "cloud_metadata": {
                            "url": "https://example.com/ENCFF079PRE.bigWig"
                        },
                    }
                ]
            }

        source = resolve_gt(assembly="hg38", accession="ENCSR000CIO", get_json=fake_get)
        self.assertEqual(source.encff, "ENCFF079PRE")
        self.assertTrue(source.url.endswith("ENCFF079PRE.bigWig"))


class UcscTests(unittest.TestCase):
    def test_encsr_is_rejected(self) -> None:
        with self.assertRaises(RequestError):
            fetch_track(
                genome="hg38",
                track="ENCSR000CIO",
                chrom="chr11",
                start=1,
                end=10,
            )


class CheckpointGuardTests(unittest.TestCase):
    def test_missing_model_defaults_to_ntv3_tracks(self) -> None:
        from ntv3_tools.cli import _require_checkpoint

        self.assertEqual(_require_checkpoint(model=None), "NTv3_tracks")

    def test_hub_id_of_the_tracks_model_resolves_to_its_released_name(self) -> None:
        from ntv3_tools.registry import resolve_model

        self.assertEqual(resolve_model(name="InstaDeepAI/NTv3_tracks"), "NTv3_tracks")

    def test_other_public_models_are_refused(self) -> None:
        from ntv3_tools.cli import _require_checkpoint

        for name in ("NTv3_100M_post", "InstaDeepAI/NTv3_650M_post"):
            with self.subTest(name=name), self.assertRaises(RequestError):
                _require_checkpoint(model=name)

    def test_old_internal_keys_are_not_a_load_path(self) -> None:
        from ntv3_tools.registry import resolve_model

        for name in (
            "ntv3_106M_7downsample_post_trained_131kb",
            "ntv3_650M_7downsample_post_trained_1mb",
            "106M-1mb",
        ):
            with (
                self.subTest(name=name),
                self.assertRaisesRegex(KeyError, "NTv3_tracks"),
            ):
                resolve_model(name=name)


class LookupTests(unittest.TestCase):
    def test_lookup_writes_resolve_and_gt(self) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            return {
                "@graph": [
                    {
                        "accession": "ENCFF079PRE",
                        "assembly": "GRCh38",
                        "output_type": "signal",
                        "cloud_metadata": {
                            "url": "https://example.com/ENCFF079PRE.bigWig"
                        },
                    }
                ]
            }

        with tempfile.TemporaryDirectory() as tmp:
            resolved = run_lookup(
                request=LookupRequest(
                    species="human",
                    tissue="K562",
                    names=["K562 CAGE"],
                    out_dir=tmp,
                ),
                get_json=fake_get,
            )
            self.assertEqual(resolved.tracks[0].gt_accession, "ENCSR000CIO")
            self.assertEqual(resolved.tracks[0].model_head, "CNhs11250_P")
            self.assertTrue((Path(tmp) / "track_resolve.json").exists())
            self.assertTrue((Path(tmp) / "gt" / "manifest.json").exists())
            manifest = json.loads((Path(tmp) / "gt" / "manifest.json").read_text())
            self.assertEqual(manifest[0]["accession"], "ENCSR000CIO")
            self.assertTrue(manifest[0]["url"])


class ViewerBundleTests(unittest.TestCase):
    def test_report_json_writes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "plan.yaml").write_text(
                "species: human\nassembly: hg38\ngene: HBE1\n"
                "chrom: chr11\nstart: 5260000\nend: 5276384\n"
                "tissue: K562\ntracks: [K562 CAGE]\n"
                "device: cpu\nwindow_bp: 16384\n"
                "want_report: true\nwant_critic: false\n"
                "analyses: [predict, saliency]\n"
                "ucsc_tracks: [clinvarMain]\n"
                "model: NTv3_tracks\n"
                "want_notebook: false\n"
            )
            Path(tmp, "track_resolve.json").write_text(
                json.dumps(
                    {
                        "species": "human",
                        "assembly": "hg38",
                        "tissue": "K562",
                        "tracks": [
                            {
                                "query": "K562 CAGE",
                                "gt_accession": "ENCSR000CIO",
                                "model_head": "CNhs11250_P",
                            }
                        ],
                    }
                )
            )
            tracks_dir = Path(tmp, "tracks")
            tracks_dir.mkdir()
            np.save(tracks_dir / "tracks.npy", np.ones((1, 16)))
            (tracks_dir / "accessions.json").write_text('["CNhs11250_P"]')
            (tracks_dir / "metadata.json").write_text(
                '{"chrom":"chr11","start":5260000,"end":5276384,"bin_size":1024}'
            )
            path = write_viewer_bundle(run_dir=tmp)
            payload = json.loads(path.read_text())
            self.assertEqual(payload["gene"], "HBE1")
            self.assertEqual(payload["locus"]["chrom"], "chr11")
            self.assertEqual(payload["ucscTracks"], ["clinvarMain"])
            self.assertIn("gtSources", payload)
            self.assertIn("saliencyTracks", payload)
            self.assertFalse((Path(tmp) / "analysis.ipynb").exists())


if __name__ == "__main__":
    unittest.main()
