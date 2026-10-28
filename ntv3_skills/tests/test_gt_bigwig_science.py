"""Science-environment GT persistence tests with a mocked BigWig boundary."""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ntv3_tools.gt_bigwig import bigwig_reader_problem, can_read_bigwig

_HAS_PYBIGWIG = importlib.util.find_spec("pyBigWig") is not None
if _HAS_PYBIGWIG:
    from ntv3_tools.gt_bigwig import (
        GroundTruthComparison,
        write_ground_truth_comparison,
    )


class FakeBigWig:
    """Minimal remote BigWig fake."""

    def __init__(self) -> None:
        self.closed = False

    def values(
        self,
        _chrom: str,
        start: int,
        end: int,
        numpy: bool,
    ) -> np.ndarray:
        """Return deterministic base values over the requested span."""
        if not numpy:
            raise AssertionError("NumPy mode is required.")
        return np.arange(end - start, dtype=float)

    def close(self) -> None:
        """Record boundary cleanup."""
        self.closed = True


class BigWigReaderTests(unittest.TestCase):
    def test_missing_module_is_named(self) -> None:
        with patch("ntv3_tools.gt_bigwig.pyBigWig", None):
            self.assertIn("not installed", bigwig_reader_problem())
            self.assertFalse(can_read_bigwig())

    def test_wheel_without_libcurl_cannot_read(self) -> None:
        build = SimpleNamespace(remote=0, numpy=1)
        with patch("ntv3_tools.gt_bigwig.pyBigWig", build):
            self.assertIn("libcurl", bigwig_reader_problem())
            self.assertFalse(can_read_bigwig())

    def test_build_without_numpy_cannot_read(self) -> None:
        build = SimpleNamespace(remote=1, numpy=0)
        with patch("ntv3_tools.gt_bigwig.pyBigWig", build):
            self.assertIn("NumPy", bigwig_reader_problem())
            self.assertFalse(can_read_bigwig())

    def test_remote_numpy_build_can_read(self) -> None:
        build = SimpleNamespace(remote=1, numpy=1)
        with patch("ntv3_tools.gt_bigwig.pyBigWig", build):
            self.assertIsNone(bigwig_reader_problem())
            self.assertTrue(can_read_bigwig())


@unittest.skipUnless(_HAS_PYBIGWIG, "pyBigWig is science-environment only")
class GroundTruthPersistenceTests(unittest.TestCase):
    def test_writes_observed_array_and_sourced_agreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "gt").mkdir()
            (root / "track_resolve.json").write_text(
                json.dumps(
                    {
                        "tracks": [
                            {
                                "query": "K562 CAGE",
                                "model_head": "CNhs11250_P",
                                "gt_accession": "ENCSR000CIO",
                            }
                        ]
                    }
                )
            )
            (root / "gt" / "manifest.json").write_text(
                json.dumps(
                    [
                        {
                            "accession": "ENCSR000CIO",
                            "encff": "ENCFF079PRE",
                            "url": (
                                "https://www.encodeproject.org/files/ENCFF079PRE/"
                                "@@download/ENCFF079PRE.bigWig"
                            ),
                        }
                    ]
                )
            )
            reader = FakeBigWig()
            with patch("ntv3_tools.gt_bigwig.pyBigWig.open", return_value=reader):
                observed_path, agreement_path = write_ground_truth_comparison(
                    comparison=GroundTruthComparison(
                        run_dir=root,
                        chrom="chr11",
                        start=100,
                        end=108,
                        bin_size=4,
                    ),
                    predicted=np.array([[1.5, 5.5]]),
                )

            np.testing.assert_allclose(
                np.load(observed_path),
                np.array([[1.5, 5.5]]),
            )
            with agreement_path.open(newline="") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["gt_accession"], "ENCSR000CIO")
            self.assertEqual(row["file_accession"], "ENCFF079PRE")
            self.assertEqual(row["n_valid"], "2")
            self.assertAlmostEqual(float(row["pearson"]), 1.0)
            self.assertTrue(reader.closed)

    def test_a_track_without_an_accession_gets_a_nan_row_and_no_agreement(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "gt").mkdir()
            (root / "track_resolve.json").write_text(
                json.dumps(
                    {
                        "tracks": [
                            {"query": "CAGE", "gt_accession": "ENCSR000CIO"},
                            {"query": "DNase", "gt_accession": None},
                        ]
                    }
                )
            )
            (root / "gt" / "manifest.json").write_text(
                json.dumps(
                    [
                        {
                            "accession": "ENCSR000CIO",
                            "url": "https://example.org/ENCFF079PRE.bigWig",
                        }
                    ]
                )
            )
            with patch("ntv3_tools.gt_bigwig.pyBigWig.open", return_value=FakeBigWig()):
                observed_path, agreement_path = write_ground_truth_comparison(
                    comparison=GroundTruthComparison(
                        run_dir=root, chrom="chr11", start=100, end=108, bin_size=4
                    ),
                    predicted=np.array([[1.5, 5.5], [0.0, 1.0]]),
                )

            observed = np.load(observed_path)
            np.testing.assert_allclose(observed[0], [1.5, 5.5])
            self.assertTrue(np.isnan(observed[1]).all())
            with agreement_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["track_query"] for row in rows], ["CAGE"])


if __name__ == "__main__":
    unittest.main()
