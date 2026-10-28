"""The report app page: the Vite UI in one file, with one run embedded.

Uses a neutral placeholder report, not a reused example locus. No network.
"""

import json
import re
import tempfile
import unittest
from pathlib import Path

from ntv3_tools.report_app import (
    DATA_TAG,
    REPORT_APP_HTML,
    SHELL_PATH,
    render_report_app,
    write_report_app,
)

_REPORT = {
    "gene": "GENE5",
    "cell": "tissue-x",
    "locus": {"chrom": "chr9", "start": 1_000, "end": 5_096},
    "ucscData": {"knownGene": {"chrom": "chr9", "start": 1_000, "end": 5_096}},
    "note": "</script><b>",
}


def _embedded(html: str) -> dict:
    match = re.search(
        r'<script type="application/json" id="ntv3-report">(.*?)</script>',
        html,
        re.DOTALL,
    )
    assert match is not None
    return json.loads(match.group(1))


class ShellTests(unittest.TestCase):
    def test_shell_is_one_self_contained_file_with_one_empty_slot(self) -> None:
        shell = SHELL_PATH.read_text()

        self.assertEqual(shell.count(DATA_TAG), 1)
        self.assertEqual(shell.count('<script type="module">'), 1)
        self.assertIsNone(re.search(r"<script[^>]*\ssrc=", shell))
        self.assertIsNone(re.search(r"<link[^>]*\shref=", shell))
        self.assertNotIn("HBE1", shell)


class RenderTests(unittest.TestCase):
    def test_the_run_is_embedded_and_cannot_close_its_tag(self) -> None:
        html = render_report_app(report=_REPORT)

        self.assertNotIn(DATA_TAG, html)
        self.assertNotIn("</script><b>", html)
        self.assertEqual(_embedded(html), _REPORT)

    def test_write_report_app_reads_report_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "report.json").write_text(json.dumps(_REPORT))

            path = write_report_app(run_dir=tmp)

            self.assertEqual(path.name, REPORT_APP_HTML)
            self.assertEqual(_embedded(path.read_text())["gene"], "GENE5")


if __name__ == "__main__":
    unittest.main()
