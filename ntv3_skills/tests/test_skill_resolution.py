"""Resolve a sentence that already names region, tissue, and assay. No network."""

import argparse
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from ntv3_tools.cli import _cmd_plan
from ntv3_tools.doctor import INSTALL_COMMAND, DoctorReport, JaxCandidate
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.gt import search_encode_experiment
from ntv3_tools.lookups import LookupRequest, run_lookup
from ntv3_tools.plan import AnalysisPlan
from ntv3_tools.resolve import resolve_tracks
from ntv3_tools.resolved import _catalog_choice, gene_and_clinvar_lanes, resolved_block

_NINE_K562_CAGE_EXPERIMENTS = {
    "@graph": [{"accession": f"ENCSR{index:03d}CGE"} for index in range(9)]
}
_CAGE_BIGWIG = {
    "@graph": [
        {
            "accession": "ENCFF079PRE",
            "assembly": "GRCh38",
            "output_type": "signal",
            "cloud_metadata": {"url": "https://example.com/ENCFF079PRE.bigWig"},
        }
    ]
}


def _refuse_search(tissue: str, assay: str) -> str:
    """Fail the test if a hinted pair reaches the ENCODE experiment search."""
    raise AssertionError(f"Searched ENCODE for {tissue} {assay}.")


class HintResolutionTests(unittest.TestCase):
    def test_bare_assay_with_tissue_hits_the_documented_pair(self) -> None:
        resolved = resolve_tracks(
            species="human",
            tissue="K562",
            names=["CAGE"],
            encode_search=_refuse_search,
        )

        row = resolved.tracks[0]
        self.assertEqual(row.query, "CAGE")
        self.assertEqual(row.gt_accession, "ENCSR000CIO")
        self.assertEqual(row.model_head, "CNhs11250_P")
        self.assertIsNone(row.encode_hits)
        self.assertIn("K562 CAGE", row.source)

    def test_assay_and_tissue_match_without_case(self) -> None:
        resolved = resolve_tracks(
            species="human",
            tissue="k562",
            names=["cage", "dnase"],
            encode_search=_refuse_search,
        )

        self.assertEqual(
            [(row.gt_accession, row.model_head) for row in resolved.tracks],
            [("ENCSR000CIO", "CNhs11250_P"), ("ENCSR921NMD", "ENCSR921NMD")],
        )

    def test_assay_without_a_hint_records_the_live_search(self) -> None:
        resolved = resolve_tracks(
            species="human",
            tissue="HepG2",
            names=["CAGE"],
            encode_search=lambda tissue, assay: "ENCSR999FAKE",
        )

        row = resolved.tracks[0]
        self.assertEqual(row.gt_accession, "ENCSR999FAKE")
        self.assertEqual(row.model_head, "")
        self.assertEqual(row.encode_hits, ("ENCSR999FAKE",))


class HintLookupTests(unittest.TestCase):
    def test_k562_cage_lookup_writes_both_ids_without_asking(self) -> None:
        searched_experiments: list[dict[str, str]] = []

        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            if params and params.get("type") == "Experiment":
                searched_experiments.append(params)
                return _NINE_K562_CAGE_EXPERIMENTS
            return _CAGE_BIGWIG

        with tempfile.TemporaryDirectory() as tmp:
            run_lookup(
                request=LookupRequest(
                    species="human",
                    tissue="K562",
                    names=["CAGE"],
                    out_dir=tmp,
                ),
                get_json=fake_get,
            )
            payload = json.loads(Path(tmp, "track_resolve.json").read_text())

        self.assertEqual(searched_experiments, [])
        row = payload["tracks"][0]
        self.assertEqual(row["query"], "CAGE")
        self.assertEqual(row["gt_accession"], "ENCSR000CIO")
        self.assertEqual(row["model_head"], "CNhs11250_P")
        self.assertNotIn("encode_hits", row)

    def test_live_search_hits_are_written_when_the_search_ran(self) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            if params and params.get("type") == "Experiment":
                return {"@graph": [{"accession": "ENCSR999FAKE"}]}
            return _CAGE_BIGWIG

        with tempfile.TemporaryDirectory() as tmp:
            run_lookup(
                request=LookupRequest(
                    species="human",
                    tissue="HepG2",
                    names=["CAGE"],
                    out_dir=tmp,
                ),
                get_json=fake_get,
            )
            payload = json.loads(Path(tmp, "track_resolve.json").read_text())

        self.assertEqual(payload["tracks"][0]["encode_hits"], ["ENCSR999FAKE"])


class EncodeExperimentSearchTests(unittest.TestCase):
    def test_mouse_liver_rna_seq_lists_the_search_term_hits(self) -> None:
        queries: list[dict[str, str]] = []

        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            queries.append(dict(params or {}))
            if params and params.get("assay_title") == "RNA-seq":
                return {"@graph": [], "total": 0}
            return {
                "@graph": [
                    {"accession": "ENCSR000AAA"},
                    {"accession": "ENCSR000BBB"},
                ],
                "total": 22,
            }

        result = None
        with self.assertRaises(RequestError) as caught:
            result = search_encode_experiment(
                tissue="liver",
                assay="RNA-seq",
                get_json=fake_get,
                species="mouse",
            )

        self.assertIsNone(result)
        message = str(caught.exception)
        self.assertIn("ENCSR000AAA", message)
        self.assertIn("ENCSR000BBB", message)
        self.assertIn("22", message)
        self.assertEqual(len(queries), 2)
        self.assertEqual(queries[1]["searchTerm"], "RNA-seq")
        self.assertEqual(
            queries[1]["replicates.library.biosample.donor.organism.scientific_name"],
            "Mus musculus",
        )

    def test_both_queries_empty_asks_without_an_accession(self) -> None:
        def fake_get(url: str, params: dict[str, str] | None = None) -> dict:
            return {"@graph": [], "total": 0}

        with self.assertRaises(RequestError) as caught:
            search_encode_experiment(
                tissue="liver",
                assay="RNA-seq",
                get_json=fake_get,
                species="mouse",
            )

        message = str(caught.exception)
        self.assertEqual(
            message,
            "ENCODE has no released RNA-seq experiment for liver. Which ENCSR id "
            "should be predicted for RNA-seq in liver? Pass it in --tracks. No id "
            "is assumed.",
        )
        self.assertNotRegex(message, r"ENCSR\d{3}[A-Z]{3}")


_JAX = DoctorReport(
    jax_python="/work/repo/.venv/bin/python",
    candidates=[
        JaxCandidate(
            "/opt/missing/bin/python",
            "$NTV3_JAX_PYTHON",
            False,
            None,
            None,
            False,
            "not found",
        ),
        JaxCandidate(
            "/work/repo/.venv/bin/python",
            ".venv in this repo (pip install -e .)",
            True,
            "3.11.13",
            "0.4.35",
            True,
            "imports jax and nucleotide_transformer_v3",
        ),
    ],
)


def _hbe1_plan(**overrides: object) -> AnalysisPlan:
    """Return the plan for "HBE1 at chr11 5266144 5270240 on hg38 in K562"."""
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
        "want_critic": True,
        "ucsc_tracks": ["knownGene", "clinvarMain"],
        "analyses": ["predict"],
    }
    values.update(overrides)
    return AnalysisPlan(**values)


class BrowserCatalogTests(unittest.TestCase):
    def test_the_human_lane_question_names_gencode_and_clinvar(self) -> None:
        catalogs = gene_and_clinvar_lanes(assembly="hg38")

        self.assertEqual(
            [(choice.track, choice.label) for choice in catalogs],
            [("knownGene", "GENCODE"), ("clinvarMain", "ClinVar")],
        )
        self.assertEqual(catalogs[0].rule, "named in the request")

    def test_mouse_gene_lane_is_ncbi_refseq(self) -> None:
        catalogs = gene_and_clinvar_lanes(assembly="mm10")

        self.assertEqual(
            [(choice.track, choice.label) for choice in catalogs],
            [("ncbiRefSeq", "NCBI RefSeq")],
        )

    def test_no_regulatory_or_histone_lane_is_added_unless_the_plan_names_it(
        self,
    ) -> None:
        tracks = [choice.track for choice in gene_and_clinvar_lanes(assembly="hg38")]

        for absent in (
            "cCREregistry",
            "wgEncodeRegDnaseClustered",
            "phastCons100way",
            "wgEncodeRegMarkH3k27acK562",
            "wgEncodeRegMarkH3k4me3K562",
        ):
            self.assertNotIn(absent, tracks)

    def test_a_track_without_a_browser_lane_is_refused(self) -> None:
        with self.assertRaisesRegex(ArtifactError, "noSuchTrack"):
            _catalog_choice(track="noSuchTrack", assembly="hg38", rule="test")


def _block(doctor: DoctorReport = _JAX, weights_cached: bool = True) -> list[str]:
    """Return the HBE1 K562 CAGE Resolved block."""
    return resolved_block(
        plan=_hbe1_plan(),
        tracks=resolve_tracks(
            species="human",
            tissue="K562",
            names=["CAGE"],
            encode_search=_refuse_search,
        ),
        doctor=doctor,
        window_rule="your span equals the CPU 4,096 bp window",
        weights_cached=weights_cached,
    )


class ResolvedBlockTests(unittest.TestCase):
    def test_hbe1_k562_cage_is_one_block_that_ends_in_approval(self) -> None:
        lines = _block()
        text = "\n".join(lines)

        self.assertEqual(lines[0], "Resolved")
        for expected in (
            "CNhs11250_P",
            "ENCSR000CIO",
            "chr11:5,266,144-5,270,240 (4,096 bp)",
            "chr11:5,267,424-5,268,960 (1,536 bp)",
            "Model: NTv3_tracks (InstaDeepAI/NTv3_tracks)",
            "load_tracks_model",
            "in the Hugging Face cache — no download",
            "knownGene (GENCODE)",
            "clinvarMain (ClinVar)",
            "/work/repo/.venv/bin/python",
            "ntv3_report.html, published as the NTv3 report artifact",
            "Critic: on — default",
        ):
            self.assertIn(expected, text)
        # The published page is ntv3_report.html; no other report page is named.
        self.assertNotIn(" report.html", text.replace("ntv3_report.html", ""))
        for absent in (
            "artifact.html",
            "PYTHONPATH",
            "Research root",
            "ntv3_106M",
            "internal JAX",
        ):
            self.assertNotIn(absent, text)
        self.assertEqual(text.count("Resolved"), 1)
        for line in lines[1:-1]:
            self.assertIn(" — ", line)
        self.assertTrue(lines[-1].startswith("Approve to run lookup"))

    def test_uncached_weights_make_the_download_part_of_the_approval(self) -> None:
        lines = _block(weights_cached=False)

        weights = next(line for line in lines if line.startswith("Weights:"))
        self.assertIn("not in the Hugging Face cache", weights)
        self.assertTrue(
            lines[-1].startswith("Approve to download InstaDeepAI/NTv3_tracks")
        )

    def test_no_jax_python_names_the_candidates_that_failed(self) -> None:
        lines = _block(
            doctor=DoctorReport(jax_python=None, candidates=_JAX.candidates[:1])
        )

        jax_line = next(line for line in lines if line.startswith("JAX python"))
        self.assertIn("/opt/missing/bin/python (not found)", jax_line)
        self.assertIn(INSTALL_COMMAND, jax_line)
        self.assertIn("Predict cannot start", jax_line)
        self.assertEqual(
            lines[-1],
            "Predict cannot start from this install. Do not approve a forward.",
        )


class PlanResolveOnlyTests(unittest.TestCase):
    def _args(self, run: str, **overrides: object) -> argparse.Namespace:
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
            "ucsc_tracks": ["knownGene", "clinvarMain"],
            "analyses": ["predict"],
            "model": None,
            "attribution_tracks": [],
            "saliency_logo_window_bp": 96,
            "saliency_logo_center_bp": None,
            "variant": None,
            "resolve_only": True,
        }
        values.update(overrides)
        return argparse.Namespace(**values)

    def setUp(self) -> None:
        # The plan names its own directory; keep it out of the real results/.
        root = tempfile.TemporaryDirectory()
        self.addCleanup(root.cleanup)
        self.results = Path(root.name)
        patcher = mock.patch("ntv3_tools.cli._RESULTS_ROOT", self.results)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_resolve_only_prints_the_block_offline_and_writes_nothing(self) -> None:
        run = "test-resolve-only-hbe1"
        stdout = io.StringIO()
        with (
            mock.patch("ntv3_tools.cli.run_doctor", return_value=_JAX),
            mock.patch(
                "ntv3_tools.cli.search_encode_experiment", side_effect=AssertionError
            ),
            contextlib.redirect_stdout(stdout),
        ):
            _cmd_plan(self._args(run))

        self.assertEqual(list(self.results.iterdir()), [])
        block = stdout.getvalue()
        self.assertIn("CNhs11250_P", block)
        self.assertIn("clinvarMain", block)
        self.assertNotIn("cCREregistry", block)
        self.assertNotIn("TotalCounts", block)
        self.assertNotIn("wgEncodeRegMarkH3k27acK562", block)
        self.assertTrue(block.rstrip().splitlines()[-1].startswith("Approve"))

    def test_compare_and_saliency_are_one_block_centered_on_the_variant(self) -> None:
        run = "test-resolve-only-compare"
        stdout = io.StringIO()
        with (
            mock.patch("ntv3_tools.cli.run_doctor", return_value=_JAX),
            mock.patch(
                "ntv3_tools.cli.search_encode_experiment", side_effect=AssertionError
            ),
            mock.patch("ntv3_tools.cli.resolve_gene_locus", side_effect=AssertionError),
            contextlib.redirect_stdout(stdout),
        ):
            _cmd_plan(
                self._args(
                    run,
                    chrom=None,
                    start=None,
                    end=None,
                    analyses=["predict", "saliency", "compare"],
                    variant=["chr11", "5268001", "c", "t"],
                )
            )

        block = stdout.getvalue()
        self.assertEqual(list(self.results.iterdir()), [])
        self.assertIn("chr11:5,265,952-5,270,048 (4,096 bp)", block)
        self.assertIn("centered on the variant chr11:5,268,001", block)
        self.assertIn("Saliency: CAGE", block)
        self.assertIn("Compare: chr11:5,268,001 C>T (1-based) — inside the crop", block)


if __name__ == "__main__":
    unittest.main()
