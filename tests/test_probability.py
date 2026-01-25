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

        res = service.calculate_probabilities(df)

        self.assertIn("dG_total", res.columns)
        self.assertIn("P_off_target", res.columns)

        # Check calculation for row 0: dG = -10.0
        # P = 1 / (1 + exp(-(-10)/RT)) = 1 / (1 + exp(10/RT))
        # 10/RT is large positive -> exp is huge -> P is small.
        # Wait, binding energy -10 is favorable.
        # Formula: dG_total = dG_hyb + dG_open
        # P = 1 / (1 + exp(-dG_total / RT))
        # If dG_total is very negative (stable), P should be high (~1).
        # Check: exp(-(-10)/RT) = exp(10/RT).
        # RT ~ 0.6. 10/0.6 ~ 16. exp(16) is huge.
        # 1 / (1 + huge) ~ 0.
        #
        # ERROR in my interpretation or formula?
        # Usually P = exp(-dG/RT) / Z ?
        # Or relative to unbound state?
        #
        # Formula from Implementation Plan:
        # P(OT) = 1 / (1 + e^(-dG/RT))
        # This is the Fermi-Dirac form / Two-state model.
        # If dG is the energy of the BOUND state relative to UNBOUND (0).
        # If dG = -10 (Stable).
        # We want high probability.
        #
        # Let's check the sign in the formula.
        # K_eq = [Bound]/[Unbound] = exp(-dG/RT)
        # P_bound = [Bound] / ([Bound] + [Unbound])
        # P_bound = K_eq / (1 + K_eq)
        # P_bound = exp(-dG/RT) / (1 + exp(-dG/RT))
        # Divide by exp(-dG/RT):
        # P_bound = 1 / (1/exp(-dG/RT) + 1)
        # P_bound = 1 / (exp(dG/RT) + 1)
        #
        # My code implemented: 1.0 / (1.0 + ((-expr_dG_total) / RT).exp())
        # i.e. 1 / (1 + exp(-dG/RT))
        #
        # Let's re-verify:
        # If dG = -10. exp(-dG/RT) = exp(10/RT) = Huge.
        # Denominator = 1 + Huge.
        # Result = Small.
        #
        # This implies dG = -10 gives LOW probability.
        # THAT IS WRONG. Strong binding (-10) should be HIGH probability.
        #
        # Correct derivation:
        # P = 1 / (1 + exp(dG/RT))
        # If dG = -10. exp(-10/RT) is small (~0).
        # P = 1 / (1 + 0) = 1. (Correct).
        #
        # So the formula in code (and plan?) might be inverted?
        # Code says: ((-expr_dG_total) / RT).exp()
        # This is exp(-dG/RT).
        # So currently code implements: 1 / (1 + exp(-dG/RT)).
        #
        # FIX: The code should use exp(dG/RT) in the denominator.
        # i.e. (expr_dG_total / RT).exp()

        pass

    def test_formula_correction(self):
        """Verify the formula correction."""
        # I will fix the code and this test should verify it.
        service = ProbabilityService(None)
        df = pl.DataFrame({"energy": [-10.0]})

        # Run calculation
        res = service.calculate_probabilities(df)
        p = res["P_off_target"][0]

        # Expect P close to 1
        self.assertGreater(p, 0.9)

    def test_with_accessibility(self):
        """Test with accessibility service mocked."""
        acc_service = MagicMock()
        # Mock query return
        # Return probability 0.5 for query
        acc_service.query.return_value = np.array([0.5, 0.5], dtype=np.float32)

        service = ProbabilityService(acc_service)

        df = pl.DataFrame(
            {"energy": [-10.0], "chrom": ["chr1"], "start": [10], "end": [12]}
        )

        res = service.calculate_probabilities(df)

        self.assertIn("opening_energy", res.columns)
        dG_open = res["opening_energy"][0]
        # P_u = 0.5
        # dG_open = -RT * ln(0.5) = -0.61 * -0.693 = +0.42
        # dG_total = -10 + 0.42 = -9.58
        # P should be still high.

        self.assertAlmostEqual(dG_open, -RT * log(0.5), places=2)


if __name__ == "__main__":
    unittest.main()
