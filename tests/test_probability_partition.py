import unittest
from unittest.mock import MagicMock, patch
import polars as pl
import numpy as np
from pathlib import Path
from RIsearch_pipeline.services.probability import ProbabilityService, RT


class TestProbabilityPartition(unittest.TestCase):
    def setUp(self):
        self.service = ProbabilityService(None)

    def test_z_score_off_targets_only(self):
        """Test P(OT) calculation with only off-targets (Z = sum of off-targets)."""
        # Data: 2 off-targets
        # 1. dG = -10.0, Expr = 100
        # 2. dG = -5.0, Expr = 100
        df = pl.DataFrame(
            {
                "energy": [-10.0, -5.0],
                "exp_value": [100.0, 100.0],
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [120, 220],
            }
        )

        # Expected Weights
        # W1 = 100 * exp(-(-10)/RT) = 100 * exp(10/RT)
        # W2 = 100 * exp(-(-5)/RT) = 100 * exp(5/RT)
        w1 = 100.0 * np.exp(10.0 / RT)
        w2 = 100.0 * np.exp(5.0 / RT)
        z = w1 + w2

        expected_p1 = w1 / z
        expected_p2 = w2 / z

        res = self.service.calculate_probabilities(df)

        p1 = res["P_off_target"][0]
        p2 = res["P_off_target"][1]

        # Assert sum is 1.0 (since Z covers all sites in this closed universe)
        self.assertAlmostEqual(p1 + p2, 1.0, places=4)
        self.assertAlmostEqual(p1, expected_p1, places=4)
        self.assertAlmostEqual(p2, expected_p2, places=4)

    @patch(
        "RIsearch_pipeline.services.probability.ProbabilityService._run_risearch_binary"
    )
    def test_z_score_with_on_target(self, mock_risearch):
        """Test P(OT) including On-Target in Z (Open Universe)."""
        # Mock RIsearch to return -30.0 kcal/mol for on-target (very stable)
        mock_risearch.return_value = -30.0

        # One Off-Target: dG = -10.0, Expr = 100
        df = pl.DataFrame(
            {
                "energy": [-10.0],
                "exp_value": [100.0],
                "chrom": ["chr1"],
                "start": [100],
                "end": [120],
            }
        )

        on_target_path = Path("on.fa")
        query_path = Path("siRNA.fa")
        on_target_expr = 1000.0

        # Expected Weights
        # W_off = 100 * exp(10/RT)
        # W_on = 1000 * exp(30/RT)
        w_off = 100.0 * np.exp(10.0 / RT)
        w_on = 1000.0 * np.exp(30.0 / RT)
        z = w_off + w_on

        expected_p_off = w_off / z

        # Calculate
        res = self.service.calculate_probabilities(
            df,
            on_target_path=on_target_path,
            query_path=query_path,
            on_target_expression=on_target_expr,
        )

        p_off = res["P_off_target"][0]

        # W_on should be huge compared to W_off because exp(30/RT) >>> exp(10/RT)
        # So P_off should be very small.
        self.assertLess(p_off, 0.001)
        self.assertAlmostEqual(p_off, expected_p_off, places=10)

        # Verify Mock was called
        mock_risearch.assert_called_once_with(query_path, on_target_path)


if __name__ == "__main__":
    unittest.main()
