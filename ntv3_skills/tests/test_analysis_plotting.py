"""Tests for report-ready disruption figures."""

import unittest

import numpy as np
from ntv3_tools.disrupt_plotting import DisruptionPlotData, plot_log2fc_tracks


class DisruptionPlotTests(unittest.TestCase):
    def test_writes_paper_style_log2fc_png_and_pdf(self) -> None:
        artifacts = plot_log2fc_tracks(
            data=DisruptionPlotData(
                reference=np.array([[1.0, 2.0, 1.0]], dtype=np.float32),
                alternate=np.array([[1.0, 4.0, 0.5]], dtype=np.float32),
                labels=["K562 CAGE"],
                variant_bin=1,
                bin_size_bp=1,
                window_bp=None,
            )
        )

        self.assertTrue(artifacts.png.startswith(b"\x89PNG"))
        self.assertTrue(artifacts.pdf.startswith(b"%PDF"))

    def test_rejects_negative_tracks_for_log_ratio(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            plot_log2fc_tracks(
                data=DisruptionPlotData(
                    reference=np.array([[-0.1, 1.0]], dtype=np.float32),
                    alternate=np.ones((1, 2), dtype=np.float32),
                    labels=["track"],
                    variant_bin=1,
                    bin_size_bp=1,
                    window_bp=None,
                )
            )


if __name__ == "__main__":
    unittest.main()
