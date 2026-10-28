"""Ground-truth alignment and agreement tests without network access."""

import unittest

import numpy as np
from ntv3_tools.gt_compare import (
    aggregate_base_values,
    finite_pair_correlations,
    read_binned_bigwig,
)


class FakeBigWig:
    """Minimal pyBigWig-compatible boundary fake."""

    def __init__(self, values: np.ndarray) -> None:
        self.values_array = values
        self.calls: list[tuple[str, int, int, bool]] = []

    def values(
        self,
        chrom: str,
        start: int,
        end: int,
        numpy: bool,
    ) -> np.ndarray:
        """Return configured values and record the exact interval."""
        self.calls.append((chrom, start, end, numpy))
        return self.values_array


class GroundTruthAlignmentTests(unittest.TestCase):
    def test_aggregate_preserves_unavailable_bins_as_nan(self) -> None:
        values = np.array([1.0, 3.0, np.nan, np.nan, 2.0, 4.0])

        binned = aggregate_base_values(values=values, n_bins=3)

        np.testing.assert_allclose(binned[[0, 2]], np.array([2.0, 3.0]))
        self.assertTrue(np.isnan(binned[1]))

    def test_bigwig_boundary_reads_exact_prediction_span(self) -> None:
        bigwig = FakeBigWig(np.arange(8, dtype=float))

        observed = read_binned_bigwig(
            bigwig=bigwig,
            chrom="chr11",
            start=5_267_168,
            end=5_267_176,
            n_bins=2,
        )

        self.assertEqual(
            bigwig.calls,
            [("chr11", 5_267_168, 5_267_176, True)],
        )
        np.testing.assert_allclose(observed, np.array([1.5, 5.5]))


class CorrelationTests(unittest.TestCase):
    def test_correlations_use_only_finite_pairs(self) -> None:
        predicted = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
        observed = np.array([2.0, 4.0, 8.0, np.nan, 10.0])

        agreement = finite_pair_correlations(
            predicted=predicted,
            observed=observed,
        )

        self.assertEqual(agreement.n_valid, 3)
        self.assertAlmostEqual(agreement.pearson, 1.0)
        self.assertAlmostEqual(agreement.spearman, 1.0)
        self.assertEqual(agreement.caveat, "")

    def test_constant_signal_has_explicit_undefined_reason(self) -> None:
        agreement = finite_pair_correlations(
            predicted=np.ones(4),
            observed=np.arange(4, dtype=float),
        )

        self.assertTrue(np.isnan(agreement.pearson))
        self.assertTrue(np.isnan(agreement.spearman))
        self.assertIn("constant predicted", agreement.caveat)


if __name__ == "__main__":
    unittest.main()
