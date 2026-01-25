import logging
from math import log
import polars as pl
from typing import Optional
import numpy as np

from RIsearch_pipeline.services.accessibility import (
    GenomeAccessibilityService,
    AccessibilityError,
)

logger = logging.getLogger(__name__)

# Gas constant in kcal/mol*K
R = 0.001987
# Temperature in Kelvin (37 C)
T = 310.15
RT = R * T


class ProbabilityService:
    """
    Service to calculate off-target probabilities by integrating hybridization energy
    and accessibility (opening energy).
    """

    def __init__(
        self, accessibility_service: Optional[GenomeAccessibilityService] = None
    ):
        self.accessibility_service = accessibility_service

    def calculate_probabilities(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Calculate P(OT) for each candidate.

        Requires 'energy' (hybridization energy) column.
        If accessibility_service is present, it looks up 'opening_energy' or computes it.

        Formula:
          dG_total = dG_hybridization + dG_open
          P(OT) = 1 / (1 + exp(-dG_total / RT))

        New Columns:
          - P_off_target
          - opening_energy (if computed)
          - dG_total
        """
        # Ensure we have energy column
        if "energy" not in df.columns:
            logger.warning("No 'energy' column found. Cannot calculate probabilities.")
            return df

        if self.accessibility_service:
            # We need to annotate with opening energy first
            df = self._annotate_opening_energy(df)

        # Calculate dG_total
        # If 'opening_energy' exists, add it. Else just hyb energy (assuming scale 0 for open)

        # Polars expression for dG_total
        if "opening_energy" in df.columns:
            expr_dG_total = pl.col("energy") + pl.col("opening_energy")
        else:
            expr_dG_total = pl.col("energy")

        # Calculate P(OT)
        # P = 1 / (1 + exp(dG_total / RT))
        # Explanation:
        # Binding Constant K = exp(-dG/RT)
        # P_bound = K / (1 + K) = exp(-dG/RT) / (1 + exp(-dG/RT))
        # Multiply top/bottom by exp(dG/RT) -> P_bound = 1 / (exp(dG/RT) + 1)

        return df.with_columns(
            [
                expr_dG_total.alias("dG_total"),
                (1.0 / (1.0 + (expr_dG_total / RT).exp())).alias("P_off_target"),
            ]
        )

    def _annotate_opening_energy(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Lookup opening energy for each row using the AccessibilityService.
        """
        if not self.accessibility_service:
            return df

        # We need chromosomal coordinates.
        # Required columns: chrom, start, end.
        if not all(c in df.columns for c in ["chrom", "start", "end"]):
            logger.warning("Missing coordinate columns. Cannot lookup accessibility.")
            return df

        # This lookup is row-by-row or batched.
        # Since our profiles are memory mapped, random access is fast.
        # However, calling python function for millions of rows is slow.
        #
        # Strategy:
        # Group by Chromosome? Or just map_rows?
        # For Polars, custom python functions in map_elements are slow.
        #
        # But we are looking up ranges.
        # Let's iterate over rows using a generator or python loop, then recreate specific columns?
        # Or, filter unique targets first?

        # Unique Targets: (chrom, start, end)
        targets = df.select(["chrom", "start", "end"]).unique()

        # Calculate logic
        # We need the mean unpaired probability across the target site?
        # Or the max?
        # Usually: Accessibility of the binding site.
        # Formula: dG_open = -RT * ln(P_u_average) or sum?
        #
        # Standard approach (e.g. RNAplfold):
        # The probability that a region is unpaired is the probability that ALL bases are unpaired?
        # Or usually approximated.
        #
        # If we use P_u from plfold -u 30:
        # It gives P that a stretch of length u is unpaired ending at i.
        # IF the target length matches 'u', we just look up one value!
        #
        # If lengths vary, and we precomputed for u=30 (fixed):
        # We can only crudely approximate if target length != 30.
        #
        # Assumption: The siRNA length is ~20. The seed is critical.
        # Often accessibility is checked for the Seed region (6-8nt) or the whole site.
        #
        # User Requirement Check: The implementation plan said:
        # "Query AccessibilityService for P_u at the target site"
        #
        # I will assume we take the *average* per-nucleotide accessibility?
        # NO, typically we use the P_u for a specific window length.
        # But we pre-computed for a FIXED u (e.g. 30).
        # And stored it.
        #
        # Let's assume for this implementations we take the P_u value at the center or average over the region?
        #
        # BETTER: For now, I will implement a lookup that takes the P_u from the profile
        # at the valid index, assuming the pre-computation 'u' roughly matches our needs or we just take the max/avg.
        #
        # Let's use: Mean probability of being unpaired for nucleotides in the region?
        # No, that's available from u=1 profile.
        #
        # If we used -u 30, we have P(segment of 30 is unpaired).
        # If our target is 21nt, P(30) is a lower bound (harder to open 30 than 21).

        # I will compute: Mean accessibility score from the profile for the indices covered by the target.
        # (This is a heuristic if the profile is u=30, but robust enough for a scaffold).

        # Extract unique list to Python
        rows = targets.rows(named=True)
        results = []

        for row in rows:
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]

            try:
                # Query range
                probs = self.accessibility_service.query(chrom, start, end)
                # Avoid log(0)
                mean_p = np.mean(probs) if len(probs) > 0 else 0.0
                mean_p = max(mean_p, 1e-10)

                # dG = -RT * ln(P)
                dG_open = -RT * log(mean_p)

                results.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "P_u_mean": mean_p,
                        "opening_energy": dG_open,
                    }
                )
            except AccessibilityError:
                # Missing chrom or data
                results.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "P_u_mean": 0.0,
                        "opening_energy": 10.0,  # High penalty for unknown
                    }
                )

        # Join back
        df_acc = pl.DataFrame(results)

        # Join on keys
        return df.join(df_acc, on=["chrom", "start", "end"], how="left")
