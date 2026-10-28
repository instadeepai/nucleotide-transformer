"""The report embeds only what the page draws, so a long run still publishes."""

import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ntv3_tools.errors import ArtifactError
from ntv3_tools.report import UCSC_ROW_FIELDS, slim_ucsc, write_viewer_bundle

from tests.run_fixtures import add_saliency_and_compare, predicted_run

_UCSC_TS = (
    Path(__file__).resolve().parents[1] / "report" / "src" / "browser" / "ucsc.ts"
)


def _four_digits(value: float) -> bool:
    return value == float(f"{value:.4g}")


class UcscRowTests(unittest.TestCase):
    def test_rows_keep_the_drawn_fields_and_the_envelope(self) -> None:
        payload = {
            "chrom": "chr11",
            "start": 10,
            "end": 20,
            "_request_url": "https://api.genome.ucsc.edu/getData/track",
            "clinvarMain": [
                {
                    "chrom": "chr11",
                    "chromStart": 12,
                    "chromEnd": 13,
                    "name": "NM_005330.5(HBE1):c.1A>G",
                    "_mouseOver": "x" * 500,
                    "_jsonHgvsTable": "[]",
                }
            ],
        }

        slim = slim_ucsc(payload)

        self.assertEqual(
            {key: value for key, value in slim.items() if key != "clinvarMain"},
            {key: value for key, value in payload.items() if key != "clinvarMain"},
        )
        self.assertEqual(
            slim["clinvarMain"],
            [
                {
                    "chrom": "chr11",
                    "chromStart": 12,
                    "chromEnd": 13,
                    "name": "NM_005330.5(HBE1):c.1A>G",
                }
            ],
        )

    def test_every_row_field_the_app_reads_is_kept(self) -> None:
        read = set(re.findall(r"\brow\.(\w+)", _UCSC_TS.read_text()))

        self.assertTrue(read)
        self.assertEqual(sorted(read - UCSC_ROW_FIELDS), [])


class EmbeddedArrayTests(unittest.TestCase):
    def test_report_rounds_arrays_drops_weights_and_slims_ucsc(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            predicted_run(run)
            add_saliency_and_compare(run)
            np.save(run / "tracks" / "tracks.npy", np.array([[0.123456789, 1.5, 0.4]]))
            clinvar = json.loads((run / "ucsc" / "clinvarMain.json").read_text())
            clinvar["clinvarMain"][0]["_mouseOver"] = "long tooltip"
            (run / "ucsc" / "clinvarMain.json").write_text(json.dumps(clinvar))

            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())
            on_disk = json.loads((run / "ucsc" / "clinvarMain.json").read_text())

        saliency = report["saliencyTracks"][0]
        assay = report["tracks"]["assays"][0]
        delta = report["disruption"]["perTrack"][0]["delta"]
        self.assertNotIn("weights", saliency)
        self.assertEqual(saliency["values"][3], 0.5)
        self.assertEqual(saliency["values"][9], -0.25)
        self.assertIn(0.1235, assay["predicted"])
        for values in (assay["predicted"], saliency["values"], delta):
            self.assertTrue(all(_four_digits(v) for v in values if v is not None))
        row = report["ucscData"]["clinvarMain"]["clinvarMain"][0]
        self.assertNotIn("_mouseOver", row)
        self.assertEqual(row["name"], "NM_005330.5(HBE1):c.1A>G")
        # The run directory keeps the full rows; only the embedded copy is slim.
        self.assertIn("_mouseOver", on_disk["clinvarMain"][0])

    def test_non_finite_saliency_or_compare_arrays_are_refused(self) -> None:
        for name in ("interpret/CNhs11250_P.npy", "disrupt/delta_tracks.npy"):
            with self.subTest(array=name), tempfile.TemporaryDirectory() as tmp:
                run = Path(tmp)
                predicted_run(run)
                add_saliency_and_compare(run)
                values = np.load(run / name).astype(float)
                values.flat[1] = np.nan
                np.save(run / name, values)

                with self.assertRaisesRegex(ArtifactError, "finite"):
                    write_viewer_bundle(run_dir=tmp)


if __name__ == "__main__":
    unittest.main()
