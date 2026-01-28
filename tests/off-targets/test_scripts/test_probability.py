import unittest
from unittest.mock import MagicMock
import polars as pl
import numpy as np
from math import log
from RIsearch_pipeline.services.probability import ProbabilityService, RT


class TestProbabilityService(unittest.TestCase):
    def test_calculate_simple(self):
        """Test calculation without accessibility."""
        service = ProbabilityService(None)

        df = pl.DataFrame({"energy": [-10.0, -5.0, 0.0]})

        res, _ = service.calculate_probabilities(df)

        self.assertIn("dG_total", res.columns)
        self.assertIn("P_off_target", res.columns)

        pass

    def test_formula_correction(self):
        """Verify the formula correction."""
        service = ProbabilityService(None)
        df = pl.DataFrame({"energy": [-10.0]})

        res, _ = service.calculate_probabilities(df)
        p = res["P_off_target"][0]

        self.assertGreater(p, 0.9)

    def test_with_accessibility(self):
        """Test with accessibility service mocked."""
        acc_service = MagicMock()
        acc_service.query.return_value = np.array([0.5, 0.5], dtype=np.float32)

        service = ProbabilityService(acc_service)

        df = pl.DataFrame(
            {
                "energy": [-10.0],
                "chrom": ["chr1"],
                "start": [10],
                "end": [12],
                "strand": ["+"],
            }
        )

        res, _ = service.calculate_probabilities(df)

        self.assertIn("opening_energy", res.columns)
        dG_open = res["opening_energy"][0]

        self.assertAlmostEqual(dG_open, 0.5, places=2)


if __name__ == "__main__":
    unittest.main()
