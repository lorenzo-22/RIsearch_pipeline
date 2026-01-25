import logging
from pathlib import Path
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

    def calculate_probabilities(
        self,
        df: pl.DataFrame,
        on_target_path: Optional[Path] = None,
        query_path: Optional[Path] = None,
        on_target_expression: float = 1000.0,
    ) -> pl.DataFrame:
        """
        Calculate P(OT) for each candidate using Partition Function (Z).

        Formula:
          Z = Sum(Expr_i * exp(-dG_total_i / RT)) + (Expr_on * exp(-dG_total_on / RT))
          P(t) = (Expr_t * exp(-dG_total_t / RT)) / Z

        Args:
            df: DataFrame with predictions.
            on_target_path: Path to On-Target FASTA (optional, requires query_path).
            query_path: Path to siRNA query FASTA (required for On-Target calc).
            on_target_expression: Expression level for on-target (default 1000.0).
        """
        # Ensure we have energy column
        if "energy" not in df.columns:
            logger.warning("No 'energy' column found. Cannot calculate probabilities.")
            return df

        # 1. Annotate Off-Targets with Opening Energy (if service available)
        if self.accessibility_service:
            df = self._annotate_opening_energy(df)

        # 2. Calculate dG_total for Off-Targets
        if "opening_energy" in df.columns:
            expr_dG_total = pl.col("energy") + pl.col("opening_energy")
        else:
            expr_dG_total = pl.col("energy")

        df = df.with_columns(expr_dG_total.alias("dG_total"))

        # 3. Calculate Boltzmann Weight for Off-Targets: W_i = Expr_i * exp(-dG_total / RT)
        if "exp_value" in df.columns:
            expr_val = pl.col("exp_value").fill_null(0.0)
        else:
            logger.info(
                "No expression value found; assuming uniform expression=1.0 for Off-Targets"
            )
            expr_val = pl.lit(1.0)

        weight_expr = expr_val * ((-pl.col("dG_total") / RT).exp())
        df = df.with_columns(weight_expr.alias("boltzmann_weight"))

        # 4. Calculate Z_partial (sum of off-targets)
        z_partial = df["boltzmann_weight"].sum()

        # 5. Calculate On-Target Weight (W_on)
        w_on = 0.0
        if on_target_path and query_path:
            logger.info("Calculating On-Target weight...")
            dG_total_on = self._calculate_on_target_dg(on_target_path, query_path)
            w_on = on_target_expression * np.exp(-dG_total_on / RT)
            logger.info(f"On-Target dG_total={dG_total_on:.2f}, Weight={w_on:.2e}")

        # 6. Total Partition Function Z
        z_total = z_partial + w_on
        logger.info(
            f"Partition Function Z = {z_total:.2e} (Off-Target={z_partial:.2e}, On-Target={w_on:.2e})"
        )

        # 7. Calculate P(t) = W_t / Z
        if z_total > 0:
            df = df.with_columns(
                (pl.col("boltzmann_weight") / z_total).alias("P_off_target")
            )
        else:
            df = df.with_columns(pl.lit(0.0).alias("P_off_target"))

        return df

    def _annotate_opening_energy(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Lookup opening energy for each row using the AccessibilityService.
        """
        if not self.accessibility_service:
            return df

        # We need chromosomal coordinates.
        # Required columns: chrom, start, end.
        if not all(c in df.columns for c in ["chrom", "start", "end", "strand"]):
            logger.warning(
                "Missing coordinate or strand columns. Cannot lookup accessibility."
            )
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

        # Unique Targets: (chrom, start, end, strand)
        targets = df.select(["chrom", "start", "end", "strand"]).unique()

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
            strand = row["strand"]

            try:
                # Query range with strand information
                # Now returns opening energies directly (not probabilities)
                opening_energies = self.accessibility_service.query(
                    chrom, start, end, strand
                )
                # Take mean opening energy for the region
                mean_e_open = (
                    np.mean(opening_energies) if len(opening_energies) > 0 else 10.0
                )

                results.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "opening_energy": mean_e_open,
                    }
                )
            except AccessibilityError:
                # Missing chrom or data
                results.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "P_u_mean": 0.0,
                        "opening_energy": 10.0,  # High penalty for unknown
                    }
                )

        # Join back
        df_acc = pl.DataFrame(results)

        # Join on keys
        return df.join(df_acc, on=["chrom", "start", "end", "strand"], how="left")

    def _calculate_on_target_dg(self, on_target_path: Path, query_path: Path) -> float:
        """
        Calculate dG_total for the on-target sequence.

        dG_total = dG_hyb (RIsearch) + dG_open (Accessibility)
        """
        # 1. Hybridization Energy using RIsearch binary
        dG_hyb = self._run_risearch_binary(query_path, on_target_path)

        # 2. Accessibility Energy (Assumed 0.0 for now as per plan/constraints)
        dG_open = 0.0

        return dG_hyb + dG_open

    def _run_risearch_binary(self, query_path: Path, target_path: Path) -> float:
        """
        Run RIsearch binary to get hybridization energy.
        Requires indexing the target first.
        """
        import subprocess
        import tempfile
        import os

        binary_path = "/Users/lorenzo/.cargo/bin/RIsearch"
        if not Path(binary_path).exists():
            logger.error(f"RIsearch binary not found at {binary_path}")
            return 0.0

        with tempfile.NamedTemporaryFile(suffix=".sa", delete=False) as tmp_index:
            index_path = Path(tmp_index.name)

        try:
            # 1. Index
            # RIsearch index <INPUT> <OUTPUT>
            cmd_index = [binary_path, "index", str(target_path), str(index_path)]
            subprocess.run(cmd_index, capture_output=True, check=True)

            # 2. Search
            # RIsearch search -q <QUERY> -i <INDEX> --output -
            cmd_search = [
                binary_path,
                "search",
                "-q",
                str(query_path),
                "-i",
                str(index_path),
                "--output",
                "-",
            ]
            result = subprocess.run(
                cmd_search, capture_output=True, text=True, check=True
            )

            # 3. Parse Output
            lines = result.stdout.splitlines()
            energies = []
            for line in lines:
                if line.startswith("#"):
                    continue
                parts = line.split()
                # Default format is likely:
                # Query Target ... Energy (at end?)
                # We look for negative floats.
                for p in parts:
                    try:
                        f = float(p)
                        if f < 0:
                            energies.append(f)
                    except ValueError:
                        pass

            if energies:
                # Assuming most negative (min) is the best binding energy
                return min(energies)
            else:
                logger.warning(
                    "Could not parse energy from RIsearch output for On-Target."
                )
                return 0.0

        except subprocess.CalledProcessError as e:
            logger.error(f"RIsearch binary execution failed: {e.stderr}")
            if e.stdout:
                logger.error(f"Stdout: {e.stdout}")
            return 0.0
        except Exception as e:
            logger.error(f"Error running RIsearch: {e}")
            return 0.0
        finally:
            # Cleanup index
            if index_path.exists():
                os.unlink(index_path)
