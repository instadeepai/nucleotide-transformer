"""Behavioral contracts for an NTv3 skill run."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from ntv3_tools.cli import _require_run_directory
from ntv3_tools.cli import main as cli_main
from ntv3_tools.disrupt import Snv
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.gt import pick_encode_bigwig, resolve_gt
from ntv3_tools.locus import resolve_gene_locus
from ntv3_tools.lookups import LookupRequest, run_lookup
from ntv3_tools.plan import (
    ALLOWED_ANALYSES,
    EXCLUDED_ANALYSES,
    AnalysisPlan,
    AttributionConfig,
    IntegratedGradientsConfig,
    IsmConfig,
    SaliencyLogoConfig,
    SmoothGradConfig,
    read_plan,
    validate_plan,
    write_plan,
)
from ntv3_tools.report import write_viewer_bundle


def make_plan(**overrides: object) -> AnalysisPlan:
    """Create a valid plan with explicit test overrides."""
    values: dict[str, object] = {
        "species": "human",
        "assembly": "hg38",
        "gene": "HBE1",
        "chrom": "chr11",
        "start": 5_260_000,
        "end": 5_276_384,
        "tissue": "K562",
        "tracks": ["K562 CAGE"],
        "device": "cpu",
        "window_bp": 16_384,
        "want_report": True,
        "want_notebook": False,
        "want_critic": True,
        "ucsc_tracks": ["clinvarMain"],
        "analyses": ["predict", "saliency"],
    }
    values.update(overrides)
    return AnalysisPlan(**values)


def read_plan_without_variant() -> AnalysisPlan:
    """Read a plan written before Compare existed: no variant key at all."""
    with tempfile.TemporaryDirectory() as tmp:
        path = write_plan(plan=make_plan(), out_dir=tmp)
        path.write_text(path.read_text().replace("variant: null\n", ""))
        return read_plan(path=path)


class PlanContractTests(unittest.TestCase):
    def test_species_and_assembly_must_agree(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(plan=make_plan(species="mouse", assembly="hg38"))

    def test_default_cpu_4096_profile_is_valid_for_both_species(self) -> None:
        for species, assembly in (("human", "hg38"), ("mouse", "mm10")):
            for analyses in (["predict"], ["predict", "saliency"]):
                with self.subTest(species=species, analyses=analyses):
                    validate_plan(
                        plan=make_plan(
                            species=species,
                            assembly=assembly,
                            start=5_266_144,
                            end=5_270_240,
                            window_bp=4_096,
                            analyses=analyses,
                            ucsc_tracks=[],
                        )
                    )

    def test_coordinates_must_be_valid(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(plan=make_plan(start=100, end=100))

    def test_window_length_must_match_coordinates(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(plan=make_plan(window_bp=1_024))

    def test_at_least_one_output_is_required(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(plan=make_plan(want_report=False, want_notebook=False))

    def test_logo_defaults_are_loaded_when_block_is_absent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_plan(plan=make_plan(), out_dir=tmp)
            payload = path.read_text().replace(
                "saliency_logo:\n  window_bp: 96\n  center_bp: null\n", ""
            )
            path.write_text(payload)
            loaded = read_plan(path=path)

        self.assertEqual(loaded.saliency_logo, SaliencyLogoConfig())

    def test_logo_window_must_fit_the_locus(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(
                plan=make_plan(saliency_logo=SaliencyLogoConfig(window_bp=20_000))
            )

    def test_logo_center_must_be_inside_the_locus(self) -> None:
        with self.assertRaises(RequestError):
            validate_plan(
                plan=make_plan(saliency_logo=SaliencyLogoConfig(center_bp=5_300_000))
            )

    def test_saliency_requires_predict(self) -> None:
        with self.assertRaisesRegex(RequestError, "Saliency requires predict"):
            validate_plan(plan=make_plan(analyses=["saliency"]))

    def test_the_allowed_set_is_predict_saliency_compare(self) -> None:
        self.assertEqual(ALLOWED_ANALYSES, ("predict", "saliency", "compare"))
        validate_plan(plan=make_plan(analyses=["predict", "saliency"]))
        # Annotate is its own command, so a plan never runs it as a forward.
        with self.assertRaisesRegex(RequestError, "ntv3 annotate"):
            validate_plan(plan=make_plan(analyses=["predict", "annotate"]))

    def test_excluded_analyses_never_validate(self) -> None:
        self.assertEqual(
            EXCLUDED_ANALYSES,
            (
                "integrated_gradients",
                "smoothgrad",
                "ism",
                "attention",
                "knockout",
                "design",
                "ledidi",
                "motif",
            ),
        )
        # An approved ISM interval does not make ism valid either.
        attribution = AttributionConfig(ism=IsmConfig(start=5_260_000, end=5_260_128))
        for method in EXCLUDED_ANALYSES:
            for analyses in (
                [method],
                ["predict", method],
                ["predict", "saliency", method],
            ):
                with (
                    self.subTest(analyses=analyses),
                    self.assertRaisesRegex(
                        RequestError, f"{method} is outside this skill"
                    ),
                ):
                    validate_plan(
                        plan=make_plan(analyses=analyses, attribution=attribution)
                    )
            with (
                tempfile.TemporaryDirectory() as tmp,
                self.assertRaisesRegex(RequestError, "outside this skill"),
            ):
                write_plan(plan=make_plan(analyses=["predict", method]), out_dir=tmp)

    def test_excluded_commands_are_refused_before_any_dispatch(self) -> None:
        commands = (
            ["knockout", "--plan", "p.yaml", "--chrom", "chr11", "--seed", "7"],
            ["design", "--plan", "p.yaml", "--track", "CAGE", "--target-delta", "1"],
            ["attention", "--plan", "p.yaml"],
            ["motif", "--plan", "p.yaml", "--method", "saliency", "--jaspar", "MA1"],
        )
        for argv in commands:
            with (
                self.subTest(command=argv[0]),
                mock.patch("sys.argv", ["ntv3", *argv]),
                mock.patch(
                    "ntv3_tools.cli.subprocess.run", side_effect=AssertionError("ran")
                ),
                self.assertLogs("ntv3_tools", level="ERROR") as logs,
                self.assertRaises(SystemExit) as exit_,
            ):
                cli_main()
            self.assertEqual(exit_.exception.code, 2)
            self.assertIn(f"{argv[0]} is outside this skill", logs.output[0])

    def test_help_offers_no_refused_command_or_setting(self) -> None:
        for argv in (["--help"], ["plan", "--help"]):
            stdout = io.StringIO()
            with (
                self.subTest(argv=argv),
                mock.patch("sys.argv", ["ntv3", *argv]),
                contextlib.redirect_stdout(stdout),
                self.assertRaises(SystemExit),
            ):
                cli_main()
            text = stdout.getvalue().lower()
            self.assertIn("usage: ntv3", text)
            for name in (*EXCLUDED_ANALYSES, "ig-n-steps", "smoothgrad", "ism-"):
                self.assertNotIn(name.replace("_", "-"), text.replace("_", "-"))
        top = self._help(["--help"])
        self.assertIn("registered only at 1,048,576 bp", top)
        self.assertIn(
            "{plan,approve,doctor,lookup,predict,compare,disrupt,annotate,"
            "annotate-forward,critic,literature,report}",
            top,
        )
        self.assertNotIn("motif", self._help(["--help"]))
        for alias in ("compare", "disrupt"):
            self.assertIn("is not a second Compare", self._help([alias, "--help"]))

    def _help(self, argv: list[str]) -> str:
        stdout = io.StringIO()
        with (
            mock.patch("sys.argv", ["ntv3", *argv]),
            contextlib.redirect_stdout(stdout),
            self.assertRaises(SystemExit),
        ):
            cli_main()
        return " ".join(stdout.getvalue().split())

    def test_plan_command_refuses_before_the_gene_lookup(self) -> None:
        with (
            mock.patch(
                "sys.argv",
                [
                    "ntv3",
                    "plan",
                    "--species",
                    "human",
                    "--assembly",
                    "hg38",
                    "--gene",
                    "HBE1",
                    "--tissue",
                    "K562",
                    "--tracks",
                    "CAGE",
                    "--device",
                    "cpu",
                    "--window-bp",
                    "4096",
                    "--want-report",
                    "--analyses",
                    "predict",
                    "ism",
                ],
            ),
            mock.patch(
                "ntv3_tools.cli.resolve_gene_locus", side_effect=AssertionError("ucsc")
            ),
            self.assertLogs("ntv3_tools", level="ERROR") as logs,
            self.assertRaises(SystemExit),
        ):
            cli_main()
        self.assertIn("ism is outside this skill", logs.output[0])

    def test_an_unknown_analysis_is_refused(self) -> None:
        with self.assertRaisesRegex(RequestError, "Unsupported analyses"):
            validate_plan(plan=make_plan(analyses=["predict", "deeplift"]))

    def test_compare_takes_one_in_window_snv_and_survives_the_round_trip(
        self,
    ) -> None:
        variant = Snv(chrom="chr11", position=5_268_001, ref="C", alt="T")
        plan = make_plan(analyses=["predict", "compare"], variant=variant)
        validate_plan(plan=plan)
        with tempfile.TemporaryDirectory() as tmp:
            loaded = read_plan(path=write_plan(plan=plan, out_dir=tmp))

        self.assertEqual(loaded.variant, variant)
        self.assertIsNone(read_plan_without_variant().variant)

    def test_compare_refuses_a_missing_outside_or_unpaired_variant(self) -> None:
        variant = Snv(chrom="chr11", position=5_268_001, ref="C", alt="T")
        cases = [
            ({"analyses": ["predict", "compare"]}, "exactly one variant"),
            ({"variant": variant}, "exactly one variant"),
            ({"analyses": ["compare"], "variant": variant}, "requires predict"),
            (
                {
                    "analyses": ["predict", "compare"],
                    "variant": Snv("chr11", 5_276_385, "C", "T"),
                },
                "outside the window",
            ),
            (
                {
                    "analyses": ["predict", "compare"],
                    "variant": Snv("chr11", 5_268_001, "C", "C"),
                },
                "two different",
            ),
        ]
        for overrides, message in cases:
            with (
                self.subTest(overrides=overrides),
                self.assertRaisesRegex(RequestError, message),
            ):
                validate_plan(plan=make_plan(**overrides))

    def test_attribution_runtime_parameters_are_bounded(self) -> None:
        invalid_configs = (
            AttributionConfig(
                integrated_gradients=IntegratedGradientsConfig(n_steps=1)
            ),
            AttributionConfig(smoothgrad=SmoothGradConfig(n_samples=0)),
            AttributionConfig(smoothgrad=SmoothGradConfig(noise_sigma=1.1)),
            AttributionConfig(ism=IsmConfig(chunk_size=0)),
        )
        for attribution in invalid_configs:
            with self.assertRaises(RequestError):
                validate_plan(plan=make_plan(attribution=attribution))

    def test_attribution_block_round_trips(self) -> None:
        attribution = AttributionConfig(
            integrated_gradients=IntegratedGradientsConfig(n_steps=32),
            smoothgrad=SmoothGradConfig(n_samples=12, noise_sigma=0.2, seed=7),
            ism=IsmConfig(chunk_size=64),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = write_plan(
                plan=make_plan(attribution=attribution),
                out_dir=tmp,
            )
            loaded = read_plan(path=path)
        self.assertEqual(loaded.attribution, attribution)


class SafetyContractTests(unittest.TestCase):
    def test_run_directory_must_be_under_worktree_results(self) -> None:
        with self.assertRaises(RequestError):
            _require_run_directory(path="/tmp/ntv3-outside-results")

    def test_worktree_results_directory_is_accepted(self) -> None:
        resolved = _require_run_directory(path="results/contract-test")
        self.assertTrue(resolved.endswith("/results/contract-test"))


def _search(symbol: str, *positions: str) -> dict:
    """A UCSC /search payload with one hgnc match per position."""
    return {
        "positionMatches": [
            {
                "name": "hgnc",
                "matches": [
                    {"position": position, "posName": symbol} for position in positions
                ],
            }
        ]
    }


def _gene_rows(*rows: tuple[str, str, int, int, str]) -> mock.Mock:
    """A UCSC knownGene fetch returning (name, geneName, start, end, strand) rows."""
    return mock.Mock(
        return_value={
            "knownGene": [
                {
                    "name": name,
                    "geneName": gene,
                    "chromStart": start,
                    "chromEnd": end,
                    "strand": strand,
                }
                for name, gene, start, end, strand in rows
            ]
        }
    )


class LocusContractTests(unittest.TestCase):
    def _resolve(self, gene: str, search: dict, rows: mock.Mock, window_bp: int):
        return resolve_gene_locus(
            gene=gene,
            assembly="hg38",
            window_bp=window_bp,
            get_json=lambda _url, _params=None: search,
            fetch_gene_rows=rows,
        )

    def test_plus_strand_gene_is_centred_on_its_tss(self) -> None:
        rows = _gene_rows(("ENST1", "SOX2", 181_711_924, 181_714_436, "+"))
        locus = self._resolve(
            "SOX2", _search("SOX2", "chr3:181711925-181714436"), rows, 16_384
        )

        self.assertEqual(
            (locus.chrom, locus.strand, locus.tss), ("chr3", "+", 181_711_924)
        )
        self.assertEqual(locus.end - locus.start, 16_384)
        self.assertEqual((locus.start + locus.end) // 2, 181_711_924)
        rows.assert_called_once_with(
            genome="hg38",
            track="knownGene",
            chrom="chr3",
            start=181_711_924,
            end=181_714_436,
        )

    def test_minus_strand_gene_is_centred_on_its_last_base(self) -> None:
        rows = _gene_rows(("ENST00000357654.9", "BRCA1", 43_044_294, 43_125_364, "-"))
        locus = self._resolve(
            "BRCA1", _search("BRCA1", "chr17:43044295-43125364"), rows, 4_096
        )

        self.assertEqual(locus.tss, 43_125_363)
        self.assertEqual((locus.start, locus.end), (43_123_315, 43_127_411))
        self.assertIn("ENST00000357654.9", locus.source)

    def test_a_long_gene_window_holds_the_promoter(self) -> None:
        rows = _gene_rows(("ENST00000275493.7", "EGFR", 55_019_016, 55_211_628, "+"))
        locus = self._resolve(
            "EGFR", _search("EGFR", "chr7:55019017-55211628"), rows, 4_096
        )

        self.assertEqual((locus.start, locus.end), (55_016_968, 55_021_064))

    def test_gene_search_refuses_multiple_distinct_exact_loci(self) -> None:
        for positions in (
            ("chr1:100-200", "chr2:100-200"),
            ("chr1:100-200", "chr1:5000-6000"),
        ):
            rows = _gene_rows()
            with (
                self.subTest(positions=positions),
                self.assertRaisesRegex(RequestError, "several exact loci"),
            ):
                self._resolve("GENE", _search("GENE", *positions), rows, 4_096)
            rows.assert_not_called()

    def test_overlapping_isoforms_are_one_locus(self) -> None:
        rows = _gene_rows(("T1", "GENE", 149, 300, "+"))
        locus = self._resolve(
            "GENE", _search("GENE", "chr1:100-200", "chr1:150-300"), rows, 4_096
        )

        self.assertEqual(locus.tss, 149)
        self.assertEqual(rows.call_args.kwargs["start"], 99)
        self.assertEqual(rows.call_args.kwargs["end"], 300)

    def test_no_transcript_of_that_name_is_a_question(self) -> None:
        rows = _gene_rows(("T1", "OTHER", 99, 200, "+"))
        with self.assertRaisesRegex(RequestError, "no transcript named 'GENE'"):
            self._resolve("GENE", _search("GENE", "chr1:100-200"), rows, 4_096)

    def test_mm10_prefers_a_curated_nm_transcript_over_a_shorter_xm_model(
        self,
    ) -> None:
        # ncbiRefSeq rows have no canonical tag; the shorter XM_ starts elsewhere.
        rows = mock.Mock(
            return_value={
                "ncbiRefSeq": [
                    {
                        "name": name,
                        "name2": "Sox2",
                        "txStart": start,
                        "txEnd": end,
                        "strand": "+",
                    }
                    for name, start, end in (
                        ("XM_006535149.1", 34_705_000, 34_706_000),
                        ("NR_000001.1", 34_704_800, 34_706_000),
                        ("NM_011443.4", 34_704_553, 34_706_610),
                    )
                ]
            }
        )
        locus = resolve_gene_locus(
            gene="Sox2",
            assembly="mm10",
            window_bp=4_096,
            get_json=lambda _url, _params=None: _search(
                "Sox2", "chr3:34704554-34706610"
            ),
            fetch_gene_rows=rows,
        )

        self.assertEqual((locus.transcript, locus.tss), ("NM_011443.4", 34_704_553))


class GroundTruthContractTests(unittest.TestCase):
    def test_resolve_hydrates_stable_cloud_url(self) -> None:
        def fake_get(
            url: str,
            _params: dict[str, str] | None = None,
        ) -> dict[str, object]:
            if "/files/ENCFFPORTAL/" in url:
                return {
                    "accession": "ENCFFPORTAL",
                    "cloud_metadata": {
                        "url": "https://encode-public.s3.amazonaws.com/object.bigWig"
                    },
                }
            return {
                "@graph": [
                    {
                        "accession": "ENCFFPORTAL",
                        "assembly": "GRCh38",
                        "output_type": "signal",
                        "href": (
                            "https://www.encodeproject.org/files/ENCFFPORTAL/"
                            "@@download/ENCFFPORTAL.bigWig"
                        ),
                    }
                ]
            }

        source = resolve_gt(
            assembly="hg38",
            accession="ENCSRTEST",
            get_json=fake_get,
        )

        self.assertEqual(
            source.url,
            "https://encode-public.s3.amazonaws.com/object.bigWig",
        )

    def test_wrong_assembly_is_not_silently_accepted(self) -> None:
        with self.assertRaises(KeyError):
            pick_encode_bigwig(
                files=[
                    {
                        "accession": "ENCFFMOUSE",
                        "assembly": "mm10",
                        "output_type": "signal",
                        "href": "/mouse.bigWig",
                    }
                ],
                assembly="hg38",
            )

    def test_experiment_object_is_used_when_file_search_is_empty(self) -> None:
        calls: list[str] = []

        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            calls.append(url)
            if url.endswith("/search/"):
                return {"@graph": []}
            if "/experiments/" in url:
                return {"files": ["/files/ENCFFFALLBACK/"]}
            if "/files/" in url:
                return {
                    "accession": "ENCFFFALLBACK",
                    "assembly": "GRCh38",
                    "file_format": "bigWig",
                    "status": "released",
                    "output_type": "read-depth normalized signal",
                    "cloud_metadata": {
                        "url": "https://example.org/ENCFFFALLBACK.bigWig"
                    },
                }
            raise AssertionError(url)

        source = resolve_gt(
            assembly="hg38",
            accession="ENCSRTEST",
            get_json=fake_get,
        )

        self.assertEqual(source.encff, "ENCFFFALLBACK")
        self.assertEqual(source.origin, "encode-public")
        self.assertTrue(any("/experiments/" in url for url in calls))


class LookupContractTests(unittest.TestCase):
    def test_failed_ground_truth_lookup_writes_no_partial_manifest(self) -> None:
        def empty_get(url: str, params: dict[str, str] | None = None) -> dict:
            if "/experiments/" in url:
                return {"files": []}
            return {"@graph": []}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(KeyError):
                run_lookup(
                    request=LookupRequest(
                        species="human",
                        tissue="K562",
                        names=["K562 CAGE"],
                        out_dir=tmp,
                    ),
                    get_json=empty_get,
                )
            self.assertFalse(Path(tmp, "track_resolve.json").exists())
            self.assertFalse(Path(tmp, "gt", "manifest.json").exists())


class ReportContractTests(unittest.TestCase):
    def test_report_requires_an_approved_plan(self) -> None:
        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaises(ArtifactError),
        ):
            write_viewer_bundle(run_dir=tmp)

    def test_report_requires_prediction_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_plan(plan=make_plan(), out_dir=tmp)
            with self.assertRaises(ArtifactError):
                write_viewer_bundle(run_dir=tmp)

    def test_report_rejects_empty_completed_critic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_plan(plan=make_plan(), out_dir=tmp)
            Path(tmp, "track_resolve.json").write_text(
                '{"species":"human","assembly":"hg38","tissue":"K562",'
                '"tracks":[{"query":"K562 CAGE","gt_accession":"ENCSR000CIO",'
                '"model_head":"CNhs11250_P"}]}'
            )
            tracks = Path(tmp, "tracks")
            tracks.mkdir()
            np.save(tracks / "tracks.npy", np.ones((1, 16)))
            (tracks / "accessions.json").write_text('["CNhs11250_P"]')
            (tracks / "metadata.json").write_text(
                '{"chrom":"chr11","start":5267168,"end":5269216,"bin_size":128}'
            )
            critic = Path(tmp, "critic")
            critic.mkdir()
            (critic / "audit_log.json").write_text("[]")

            with self.assertRaisesRegex(ArtifactError, "critic"):
                write_viewer_bundle(run_dir=tmp)

    def test_notebook_is_only_written_when_requested(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_plan(plan=make_plan(want_notebook=False), out_dir=tmp)
            Path(tmp, "track_resolve.json").write_text(
                '{"species":"human","assembly":"hg38","tissue":"K562",'
                '"tracks":[{"query":"K562 CAGE","gt_accession":"ENCSR000CIO",'
                '"model_head":"CNhs11250_P"}]}'
            )
            tracks = Path(tmp, "tracks")
            tracks.mkdir()
            np.save(tracks / "tracks.npy", np.ones((1, 16)))
            (tracks / "accessions.json").write_text('["CNhs11250_P"]')
            (tracks / "metadata.json").write_text(
                '{"chrom":"chr11","start":5267168,"end":5269216,"bin_size":128}'
            )
            ground_truth = Path(tmp, "gt")
            ground_truth.mkdir()
            np.save(
                ground_truth / "observed.npy",
                np.linspace(0.0, 1.0, 16, dtype=float)[None, :],
            )
            (ground_truth / "manifest.json").write_text(
                '[{"accession":"ENCSR000CIO","encff":"ENCFF079PRE",'
                '"url":"https://www.encodeproject.org/files/ENCFF079PRE/@@download/'
                'ENCFF079PRE.bigWig","origin":"encode-public",'
                '"label":"signal of unique reads"}]'
            )
            tables = Path(tmp, "tables")
            tables.mkdir()
            (tables / "track_agreement.csv").write_text(
                "track_query,model_head,gt_accession,file_accession,n_valid,"
                "pearson,spearman,caveat\n"
                "K562 CAGE,CNhs11250_P,ENCSR000CIO,ENCFF079PRE,16,0.5,0.4,"
                "Model logits and experimental signal differ; correlations compare "
                "shape.\n"
            )
            interpret = Path(tmp, "interpret")
            interpret.mkdir()
            figures = Path(tmp, "figures")
            figures.mkdir()
            (figures / "saliency_CNhs11250_P.png").touch()
            (figures / "saliency_CNhs11250_P.pdf").touch()
            np.save(interpret / "CNhs11250_P.npy", np.ones((16, 4)))
            (interpret / "CNhs11250_P.json").write_text(
                '{"accession":"CNhs11250_P","start":5260000,"end":5260016,'
                '"seq":"ACGTACGTACGTACGT","logo":{'
                '"selection":"strongest_peak","center_bp":5260008,'
                '"start":5260000,"end":5260016,'
                '"window_bp":16,"png":"figures/saliency_CNhs11250_P.png",'
                '"pdf":"figures/saliency_CNhs11250_P.pdf"}}'
            )
            report_path = write_viewer_bundle(run_dir=tmp)
            report = json.loads(report_path.read_text())
            self.assertFalse(Path(tmp, "analysis.ipynb").exists())
            self.assertEqual(report["tracks"]["start"], 5_267_168)
            self.assertEqual(report["tracks"]["end"], 5_269_216)
            self.assertEqual(report["locus"]["focusStart"], 5_267_168)
            self.assertEqual(report["locus"]["focusEnd"], 5_269_216)
            self.assertEqual(len(report["tracks"]["assays"][0]["observed"]), 16)
            self.assertEqual(report["trackAgreement"][0]["nValid"], 16)
            self.assertEqual(report["trackAgreement"][0]["pearson"], 0.5)
            self.assertEqual(
                report["saliencyTracks"][0]["name"], "K562 CAGE CNhs11250_P"
            )
            self.assertEqual(report["saliencyTracks"][0]["modelHead"], "CNhs11250_P")
            self.assertNotIn("weights", report["saliencyTracks"][0])
            self.assertEqual(report["saliencyTracks"][0]["values"], [4.0] * 16)
            self.assertEqual(
                report["saliencyTracks"][0]["logo"]["png"],
                "figures/saliency_CNhs11250_P.png",
            )
            stages = {stage["id"]: stage for stage in report["analysisStages"]}
            # The plan has the critic on, so the report ran it.
            self.assertEqual(stages.pop("claims")["status"], "complete")
            self.assertEqual(
                {stage["status"] for stage in stages.values()}, {"not_requested"}
            )
            explain = report["explain"]
            self.assertEqual(explain["model"], make_plan().model)
            self.assertEqual(explain["species"], "human")
            self.assertEqual(
                explain["window"],
                {"chrom": "chr11", "start": 5_260_000, "end": 5_276_384, "bp": 16_384},
            )
            self.assertEqual(
                explain["crop"],
                {"chrom": "chr11", "start": 5_267_168, "end": 5_269_216, "bp": 2_048},
            )
            self.assertEqual((explain["binSize"], explain["nBins"]), (128, 16))
            self.assertEqual(
                explain["tracks"],
                [
                    {
                        "label": "K562 CAGE CNhs11250_P",
                        "modelHead": "CNhs11250_P",
                        "gtAccession": "ENCSR000CIO",
                        "gtUnit": "signal of unique reads",
                    }
                ],
            )
            self.assertEqual(explain["letters"][0]["modelHead"], "CNhs11250_P")
            self.assertIsNone(explain["compare"])
            self.assertIn("tracks/metadata.json", explain["sources"])
            report_trace = next(
                step for step in report["trace"] if step["name"] == "report"
            )
            self.assertEqual(report_trace["status"], "ok")
            self.assertEqual(report_trace["artifact"], "report.json")
            trace_names = {step["name"] for step in report["trace"]}
            self.assertTrue(
                trace_names.isdisjoint(
                    {"knockout", "design", "motif", "attention", "literature"}
                )
            )
            self.assertIn("disrupt", trace_names)

            # Artifacts of the removed analyses are not loaded or traced.
            attention = Path(tmp, "attention")
            attention.mkdir()
            np.save(attention / "last_layer_mean.npy", np.ones((2, 2)))
            literature = Path(tmp, "literature")
            literature.mkdir()
            (literature / "citations.bib").write_text("@article{x}\n")
            (literature / "search_manifest.json").write_text('{"records": []}')

            refreshed_path = write_viewer_bundle(run_dir=tmp, refresh=True)
            refreshed = json.loads(refreshed_path.read_text())
            refreshed_stages = {
                stage["id"]: stage for stage in refreshed["analysisStages"]
            }
            self.assertEqual(refreshed_stages["literature"]["status"], "complete")
            self.assertIn(
                {
                    "name": "literature",
                    "artifact": "literature/search_manifest.json",
                    "status": "ok",
                },
                refreshed["trace"],
            )
            self.assertNotIn("attention", refreshed_stages)
            self.assertNotIn("attention", refreshed)
            self.assertNotIn("attention", {step["name"] for step in refreshed["trace"]})
            self.assertTrue(Path(tmp, "code", "run_analysis.py").is_file())


if __name__ == "__main__":
    unittest.main()
