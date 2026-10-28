"""Pure tests for saliency sequence-logo window selection."""

import unittest

import numpy as np
from ntv3_tools.plotting import select_saliency_window


class SaliencyWindowTests(unittest.TestCase):
    def test_default_uses_l1_strongest_peak_without_sign_cancellation(self) -> None:
        scores = np.zeros((128, 4), dtype=np.float32)
        scores[20] = [10.0, -10.0, 0.0, 0.0]
        scores[80] = [2.0, 0.0, 0.0, 0.0]

        window = select_saliency_window(
            scores=scores,
            genomic_start=1_000,
            window_bp=32,
        )

        self.assertEqual(window.selection, "strongest_peak")
        self.assertEqual(window.center_bp, 1_020)
        self.assertEqual((window.left, window.right), (4, 36))

    def test_explicit_genomic_center_overrides_peak(self) -> None:
        scores = np.zeros((128, 4), dtype=np.float32)
        scores[20, 0] = 10.0

        window = select_saliency_window(
            scores=scores,
            genomic_start=1_000,
            window_bp=32,
            center_bp=1_080,
        )

        self.assertEqual(window.selection, "explicit_center")
        self.assertEqual((window.start, window.end), (1_064, 1_096))

    def test_start_boundary_shifts_to_keep_requested_width(self) -> None:
        scores = np.zeros((128, 4), dtype=np.float32)
        scores[0, 0] = 10.0

        window = select_saliency_window(
            scores=scores,
            genomic_start=1_000,
            window_bp=32,
        )

        self.assertEqual((window.left, window.right), (0, 32))

    def test_end_boundary_shifts_to_keep_requested_width(self) -> None:
        scores = np.zeros((128, 4), dtype=np.float32)
        scores[-1, 0] = 10.0

        window = select_saliency_window(
            scores=scores,
            genomic_start=1_000,
            window_bp=32,
        )

        self.assertEqual((window.left, window.right), (96, 128))

    def test_short_sequence_uses_all_bases(self) -> None:
        scores = np.zeros((12, 4), dtype=np.float32)

        window = select_saliency_window(
            scores=scores,
            genomic_start=1_000,
            window_bp=96,
        )

        self.assertEqual((window.left, window.right), (0, 12))

    def test_center_outside_scores_raises(self) -> None:
        scores = np.zeros((12, 4), dtype=np.float32)

        with self.assertRaises(ValueError):
            select_saliency_window(
                scores=scores,
                genomic_start=1_000,
                window_bp=8,
                center_bp=999,
            )


if __name__ == "__main__":
    unittest.main()
