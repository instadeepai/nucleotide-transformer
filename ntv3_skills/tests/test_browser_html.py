"""The browser page: a small static template that Python fills with one run."""

import argparse
import json
import re
import subprocess
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from ntv3_tools.browser import (
    BROWSER_HTML,
    RUN_DATA_TAG,
    TEMPLATE_PATH,
    _gene_item,
    render_browser_html,
)
from ntv3_tools.cli import _cmd_report
from ntv3_tools.report import write_viewer_bundle
from ntv3_tools.report_app import REPORT_APP_HTML

from tests.run_fixtures import (
    add_saliency_and_compare,
    predicted_run,
    write_annotation,
    write_gene_locus,
)

_CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
_DEFAULT_LANES = ("GENCODE", "ClinVar", "ENCODE cCREs", "DNase clusters", "phastCons")


def _run_data(html: str) -> dict | None:
    """Return the run data embedded in a page."""
    match = re.search(
        r'<script type="application/json" id="ntv3-run">(.*?)</script>', html, re.DOTALL
    )
    assert match is not None
    return json.loads(match.group(1))


class GeneRowTests(unittest.TestCase):
    def test_mouse_genepred_rows_use_tx_coordinates(self) -> None:
        # mm10 knownGene comes back as genePred: txStart/txEnd, no chromStart.
        row = {
            "name": "ENSMUST00000000001.1",
            "chrom": "chr2",
            "strand": "-",
            "txStart": 1_000,
            "txEnd": 5_000,
            "exonStarts": "1000,4000,",
            "exonEnds": "1200,5000,",
        }

        gene = _gene_item(row)

        self.assertEqual((gene["s"], gene["e"]), (1_000, 5_000))
        self.assertEqual(gene["b"], [[1_000, 1_200], [4_000, 5_000]])
        self.assertEqual(gene["d"], "-")
        # No symbol column in the row: the transcript id stays; nothing is guessed.
        self.assertEqual(gene["n"], "ENSMUST00000000001.1")

    def test_mouse_genepred_rows_label_the_gene_symbol_when_present(self) -> None:
        # mm10 genePred tracks such as ncbiRefSeq carry the symbol in name2.
        row = {
            "name": "NM_000000001.1",
            "name2": "Gene1",
            "chrom": "chr2",
            "strand": "+",
            "txStart": 2_000,
            "txEnd": 3_000,
            "exonStarts": "2000,",
            "exonEnds": "3000,",
        }

        self.assertEqual(_gene_item(row)["n"], "Gene1")
        self.assertEqual(
            _gene_item({**row, "name2": "", "geneSymbol": "Gene2"})["n"], "Gene2"
        )

    def test_human_biggenepred_rows_keep_gene_name(self) -> None:
        row = {
            "name": "ENST00000000001.1",
            "geneName": "GENE3",
            "chrom": "chr3",
            "strand": "+",
            "chromStart": 10,
            "chromEnd": 90,
            "blockSizes": "20,30,",
            "chromStarts": "0,50,",
        }

        gene = _gene_item(row)

        self.assertEqual(gene["n"], "GENE3")
        self.assertEqual(gene["b"], [[10, 30], [60, 90]])


class TemplateTests(unittest.TestCase):
    def test_template_is_small_self_contained_and_has_one_data_tag(self) -> None:
        html = TEMPLATE_PATH.read_text()

        self.assertEqual(html.count(RUN_DATA_TAG), 1)
        self.assertLess(len(html.encode()), 250_000)
        self.assertIsNone(re.search(r"<script[^>]*\ssrc=", html))
        self.assertIsNone(re.search(r"<link[^>]*\shref=", html))
        self.assertNotIn("fetch(", html)
        self.assertNotIn("XMLHttpRequest", html)
        self.assertNotIn("http://", html.replace("http://www.w3.org/2000/svg", ""))
        self.assertNotIn("https://", html)
        self.assertIn('name="viewport"', html)
        self.assertIn("<title>NTv3 browser</title>", html)

    def test_empty_state_names_the_default_lanes_and_no_locus(self) -> None:
        html = render_browser_html(run_data=None)

        self.assertIsNone(_run_data(html))
        for lane in _DEFAULT_LANES:
            self.assertIn(f'"{lane}"', html)
        self.assertIn("Predicted tracks appear when a run has them.", html)
        self.assertIn("An NTv3 annotation lane appears when the run has a GFF3.", html)
        for absent in ("chr11", "HBE1", "K562", "Agarwal", "agarwal"):
            self.assertNotIn(absent, html)

    def test_run_data_cannot_close_the_script(self) -> None:
        html = render_browser_html(run_data={"title": "HBE1</script><b>"})

        self.assertNotIn(RUN_DATA_TAG, html)
        self.assertNotIn("HBE1</script><b>", html)
        self.assertEqual(_run_data(html), {"title": "HBE1</script><b>"})


class ReportWritesBrowserTests(unittest.TestCase):
    def test_page_holds_the_ucsc_lanes_and_the_predicted_track(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            data = _run_data(Path(tmp, BROWSER_HTML).read_text())
            crate = Path(tmp, "ro-crate-metadata.json").read_text()
            self.assertFalse(Path(tmp, "report.html").exists())

        assert data is not None
        labels = [lane["label"] for lane in data["lanes"]]
        self.assertEqual(labels, ["GENCODE", "ClinVar", "Predicted CAGE CNhs11250_P"])
        gene = data["lanes"][0]["items"][0]
        self.assertEqual(gene["n"], "HBE1")
        self.assertEqual(len(gene["b"]), 3)
        predicted = data["lanes"][2]
        self.assertEqual(predicted["kind"], "signal")
        self.assertEqual(predicted["items"][0], [5267424, 5267936, 0.2])
        self.assertEqual(data["focus"], [5267424, 5268960])
        self.assertIn(
            "An NTv3 annotation lane appears when the run has a GFF3.", data["notes"]
        )
        self.assertIn("NTv3_tracks", dict(data["explain"])["Model"])
        self.assertIn('"browser.html"', crate)

    def test_saliency_and_compare_add_lanes_a_variant_mark_and_a_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            add_saliency_and_compare(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            data = _run_data(Path(tmp, BROWSER_HTML).read_text())

        assert data is not None
        lanes = {lane["label"]: lane for lane in data["lanes"]}
        self.assertEqual(
            list(lanes),
            [
                "GENCODE",
                "ClinVar",
                "Predicted CAGE CNhs11250_P",
                "Saliency CAGE CNhs11250_P",
                "Alt − ref CAGE CNhs11250_P",
            ],
        )
        self.assertEqual(
            lanes["Saliency CAGE CNhs11250_P"]["items"][3], [5268403, 5268404, 0.5]
        )
        self.assertEqual(lanes["Saliency CAGE CNhs11250_P"]["color"], "saliency")
        self.assertEqual(
            lanes["Alt − ref CAGE CNhs11250_P"]["items"][1], [5267936, 5268448, -0.3]
        )
        self.assertEqual(
            data["marks"], [{"pos": 5268500, "label": "chr11:5268501 A>G"}]
        )
        self.assertEqual(
            data["compare"]["rows"],
            [["CAGE CNhs11250_P", "CNhs11250_P", "-0.3", "—"]],
        )
        self.assertIn("no null was run", data["compare"]["summary"])
        self.assertIn("chr11:5268501 A>G", dict(data["explain"])["Compare"])

    def test_compare_names_the_saved_clinvar_record_only_on_an_exact_match(
        self,
    ) -> None:
        for hgvs, expected in (("g.5268501A>G", True), ("g.5268501A>T", False)):
            with self.subTest(hgvs=hgvs), tempfile.TemporaryDirectory() as tmp:
                predicted_run(Path(tmp))
                add_saliency_and_compare(Path(tmp))
                clinvar = Path(tmp, "ucsc", "clinvarMain.json")
                payload = json.loads(clinvar.read_text())
                payload["clinvarMain"][0].update(
                    name="VCV000000001", _jsonHgvsTable=json.dumps([[hgvs, "-"]])
                )
                clinvar.write_text(json.dumps(payload))
                write_viewer_bundle(run_dir=tmp)
                data = _run_data(Path(tmp, BROWSER_HTML).read_text())

            assert data is not None
            compare = dict(data["explain"])["Compare"]
            mark = data["marks"][0]["label"]
            if expected:
                self.assertEqual(
                    compare, "chr11:5268501 A>G (ClinVar VCV000000001) with NTv3_tracks"
                )
                self.assertEqual(mark, "chr11:5268501 A>G (ClinVar VCV000000001)")
            else:
                self.assertNotIn("ClinVar", compare)
                self.assertEqual(mark, "chr11:5268501 A>G")

    def test_a_gene_run_without_a_knowngene_lane_gives_the_page_a_null_source(
        self,
    ) -> None:
        # Lookup writes locus_metadata.json only with a knownGene lane, so the
        # report reads locus.json, whose source is text, not {track, url}.
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_gene_locus(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            locus = json.loads(Path(tmp, "report.json").read_text())["locus"]

        self.assertIsNone(locus["source"])
        self.assertEqual((locus["tss"], locus["strand"]), (5269944, "-"))

    def test_refresh_adds_the_annotation_lane_in_the_same_view(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            write_annotation(Path(tmp))
            write_viewer_bundle(run_dir=tmp, refresh=True)
            data = _run_data(Path(tmp, BROWSER_HTML).read_text())

        assert data is not None
        annotation = data["lanes"][-1]
        self.assertEqual(annotation["label"], "NTv3 annotation")
        self.assertEqual(
            [gene["n"] for gene in annotation["items"]], ["gene_1_p", "gene_2_m"]
        )
        self.assertEqual(
            annotation["items"][1]["b"], [[5268000, 5268300], [5269200, 5269500]]
        )
        self.assertEqual(data["view"], {"start": 5266144, "end": 5270240})
        self.assertIn("hmm_constrained", dict(data["explain"])["Decoder"])

    def test_annotate_run_frames_the_predicted_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            data = _run_data(Path(tmp, BROWSER_HTML).read_text())

        assert data is not None
        self.assertEqual(data["view"], {"start": 5265500, "end": 5270000})
        self.assertIn("Predicted tracks appear when a run has them.", data["notes"])

    def test_report_command_says_publish_and_names_no_path_to_open(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            args = argparse.Namespace(run_dir=tmp, refresh=False, stages=None, contig=None)
            with (
                mock.patch("ntv3_tools.cli._require_run_directory", return_value=tmp),
                self.assertLogs("ntv3_tools", level="INFO") as logs,
            ):
                _cmd_report(args)
            self.assertTrue(Path(tmp, BROWSER_HTML).is_file())

        message = "\n".join(logs.output)
        self.assertIn("Artifact tool", message)
        self.assertIn("do not show its path", message)
        self.assertIn("Publish ntv3_report.html", message)
        self.assertIsNone(re.search(r"(?<!ntv3_)report\.html", message))
        self.assertNotIn("browser.html", message)
        self.assertNotIn("pnpm", message)
        self.assertNotIn("artifact limit", message)

    def test_report_command_warns_when_the_page_is_over_the_artifact_limit(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            args = argparse.Namespace(
                run_dir=tmp, refresh=False, stages=None, contig=None
            )
            with (
                mock.patch("ntv3_tools.cli._require_run_directory", return_value=tmp),
                mock.patch("ntv3_tools.cli.ARTIFACT_PAGE_LIMIT_BYTES", 1_000),
                self.assertLogs("ntv3_tools", level="WARNING") as logs,
            ):
                _cmd_report(args)

        message = "\n".join(logs.output)
        self.assertIn("over the 16 MB artifact limit", message)
        self.assertIn("Do not publish it", message)


def _chrome_dom(html_path: Path, width: int) -> str:
    """Render a page in headless Chrome with every host unresolvable.

    Returns the rendered document without its scripts, so the embedded run data
    cannot satisfy an assertion meant for the page.
    """
    with tempfile.TemporaryDirectory() as tmp:
        dom_path = Path(tmp, "dom.html")
        with dom_path.open("w") as dom:
            chrome = subprocess.Popen(
                [
                    str(_CHROME),
                    "--headless=new",
                    "--disable-gpu",
                    "--no-first-run",
                    "--no-default-browser-check",
                    f"--user-data-dir={Path(tmp, 'profile')}",
                    "--host-resolver-rules=MAP * ~NOTFOUND",
                    f"--window-size={width},900",
                    "--virtual-time-budget=5000",
                    "--dump-dom",
                    html_path.as_uri(),
                ],
                stdout=dom,
                stderr=subprocess.DEVNULL,
            )
            deadline = time.monotonic() + 60
            try:
                while "</html>" not in dom_path.read_text()[-100:]:
                    if chrome.poll() is not None or time.monotonic() > deadline:
                        break
                    time.sleep(0.5)
            finally:
                chrome.kill()
                chrome.wait()
        text = dom_path.read_text()
    return re.sub(r"<script\b[^>]*>.*?</script>", "", text, flags=re.DOTALL)


@unittest.skipUnless(_CHROME.is_file(), "Google Chrome is not installed")
class OfflineRenderTests(unittest.TestCase):
    def test_empty_state_renders_the_default_lanes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, BROWSER_HTML)
            path.write_text(render_browser_html(run_data=None))
            dom = _chrome_dom(path, width=400)

        self.assertIn("No run loaded", dom)
        for lane in _DEFAULT_LANES:
            self.assertIn(f"<b>{lane}</b>", dom)

    def test_run_renders_tracks_gencode_and_the_annotation_lane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_annotation(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, BROWSER_HTML), width=400)

        self.assertNotIn("No run loaded", dom)
        for expected in (
            "<b>GENCODE</b>",
            "<b>Predicted CAGE CNhs11250_P</b>",
            "<b>NTv3 annotation</b>",
            "gene_1_p",
            "HBE1",
            "InstaDeepAI/NTv3_annot_mammals",
            "<svg",
        ):
            self.assertIn(expected, dom)

    def test_compare_run_renders_color_pickers_the_mark_and_the_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            add_saliency_and_compare(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, BROWSER_HTML), width=900)

        for expected in (
            "<b>Saliency CAGE CNhs11250_P</b>",
            "<b>Alt − ref CAGE CNhs11250_P</b>",
            'aria-label="Color for Predicted CAGE CNhs11250_P"',
            "<title>chr11:5268501 A&gt;G</title>",
            "<th>Alt − ref (largest)</th>",
            "<td>CNhs11250_P</td>",
        ):
            self.assertIn(expected, dom)
        self.assertEqual(dom.count('type="color"'), 5)

    def test_report_page_renders_a_gene_run_without_a_knowngene_lane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_gene_locus(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, REPORT_APP_HTML), width=900)

        self.assertNotIn("No run loaded", dom)
        self.assertIn("HBE1", dom)

    def test_report_page_renders_a_compare_run_without_a_null(self) -> None:
        # The default Compare runs no null SNVs: percentile null, nullN 0.
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            add_saliency_and_compare(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, REPORT_APP_HTML), width=900)

        self.assertNotIn("No run loaded", dom)
        self.assertIn("No null SNVs were run, so the effect is not ranked.", dom)

    def test_report_page_lanes_name_their_source_and_skip_an_empty_rail(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, REPORT_APP_HTML), width=1200)
            write_annotation(Path(tmp))
            write_viewer_bundle(run_dir=tmp, refresh=True)
            annotated = _chrome_dom(Path(tmp, REPORT_APP_HTML), width=1200)

        self.assertIn("NTv3 · CAGE CNhs11250_P", dom)
        # 1.5 is the fixture's largest prediction; the old scale read
        # "0–1.50 model logit". The unit is in the Explain table instead.
        self.assertIn("0–1.5<", dom)
        self.assertNotIn("0–1.5 model logit", dom)
        self.assertNotIn("NTv3 annotation", dom)
        self.assertIn("NTv3 annotation", annotated)

    def test_report_page_explain_table_gives_each_lane_its_unit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            add_saliency_and_compare(Path(tmp))
            write_viewer_bundle(run_dir=tmp)
            dom = _chrome_dom(Path(tmp, REPORT_APP_HTML), width=1200)

        self.assertIn("<th>Unit</th>", dom)
        for lane, unit in (
            ("NTv3 · CAGE CNhs11250_P", "model logit"),
            ("Saliency · CAGE CNhs11250_P", "contribution"),
        ):
            self.assertRegex(dom, f"<td>{lane}</td><td>.*?</td><td>{unit}</td>")


if __name__ == "__main__":
    unittest.main()
