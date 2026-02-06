from loguru import logger
from pathlib import Path
import polars as pl
from typing import Optional, Tuple, Dict
import numpy as np

from RIsearch_pipeline.services.accessibility import (
    GenomeAccessibilityService,
    AccessibilityError,
)

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
        self,
        accessibility_service: Optional[GenomeAccessibilityService] = None,
        use_rnaplfold_cli: bool = False,
    ):
        self.accessibility_service = accessibility_service
        self.use_rnaplfold_cli = use_rnaplfold_cli

    def calculate_probabilities(
        self,
        df: pl.DataFrame,
        on_target_path: Optional[Path] = None,
        query_path: Optional[Path] = None,
        on_target_expression: float = 1000.0,
        on_target_accessibility_path: Optional[Path] = None,
        on_target_risearch_path: Optional[Path] = None,
    ) -> Tuple[pl.DataFrame, Dict]:
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
        # 5. Calculate On-Target Weight (W_on)
        w_on = 0.0
        dG_total_on = 0.0
        dG_open_on = 0.0
        dG_hyb_on = 0.0

        if on_target_path and query_path:
            logger.info("Calculating On-Target weight...")
            dG_hyb_on, dG_open_on = self._calculate_on_target_components(
                on_target_path,
                query_path,
                on_target_accessibility_path,
                on_target_risearch_path,
            )
            dG_total_on = dG_hyb_on + dG_open_on
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

        # Append On-Target candidate to DataFrame (to match legacy output behavior)
        if on_target_path and query_path:
            # Re-run to get metadata (or fetch from parsed file again - slightly inefficient but consistent)
            if on_target_risearch_path:
                _hyb, _start, _end, _strand = self._parse_risearch_file(
                    on_target_risearch_path
                )
            else:
                _hyb, _start, _end, _strand = self._run_risearch_binary(
                    query_path, on_target_path
                )

            # P_on = W_on / Z
            p_on_val = w_on / z_total if z_total > 0 else 0.0

            # Construct row
            row_dict = {
                "chrom": "onTarget",
                "start": _start,
                "end": _end,
                "strand": _strand,
                "transcript_id": "onTarget",
                "gene_id": "onTarget",
                "energy": _hyb,
                "opening_energy": dG_open_on,
                "dG_total": dG_total_on,
                "exp_value": on_target_expression,
                "P_off_target": p_on_val,
            }

            ont_df = pl.DataFrame([row_dict])

            # Align schema with main df
            for col in df.columns:
                if col not in ont_df.columns:
                    # Determine type from df
                    dtype = df.schema[col]
                    ont_df = ont_df.with_columns(pl.lit(None, dtype=dtype).alias(col))
                else:
                    # Cast existing column to match (e.g. Int64 -> Int32)
                    target_dtype = df.schema[col]
                    current_dtype = ont_df.schema[col]
                    if target_dtype != current_dtype:
                        ont_df = ont_df.with_columns(pl.col(col).cast(target_dtype))

            ont_df = ont_df.select(df.columns)
            df = pl.concat([df, ont_df], how="vertical")

        # Metadata collection
        metadata = {
            "z_total": float(z_total),
            "z_off_target": float(z_partial),
            "w_on_target": float(w_on),
            "dG_on_target": float(dG_total_on),
            "has_on_target": bool(on_target_path and query_path),
        }

        # Add P_on to metadata if applicable
        if on_target_path and query_path:
            metadata["p_on_target"] = float(w_on / z_total) if z_total > 0 else 0.0

        return df, metadata

    def calculate_probabilities_per_sirna(
        self,
        df: pl.DataFrame,
        alpha_gamma_pairs: Optional[list[tuple[float, float]]] = None,
        theta_values: Optional[list[float]] = None,
        on_target_map: Optional[dict[str, str]] = None,
        on_target_expression: float = 1000.0,
    ) -> Tuple[pl.DataFrame, Dict]:
        """
        Calculate P(OT) per siRNA using per-siRNA partition functions.

        When processing multiple siRNAs, each siRNA should have its own partition
        function Z. This method uses native Polars group_by operations for
        efficient parallel computation.

        Args:
            df: DataFrame with predictions.
            alpha_gamma_pairs: List of (alpha, gamma) tuples for clamping.
            theta_values: List of theta scaling values.
            on_target_map: Mapping of siRNA_id -> on_target_transcript_id.
            on_target_expression: Expression level for on-targets (default 1000.0).

        Formula per siRNA:
          Z_s = Sum_i(Expr_i * exp(-dG_total_i / RT)) + W_on_s
          P_s(t) = (Expr_t * exp(-dG_total_t / RT)) / Z_s

        Returns:
            Tuple of (DataFrame with P_off_target per siRNA, metadata dict).
        """
        if "energy" not in df.columns:
            logger.warning("No 'energy' column found. Cannot calculate probabilities.")
            return df, {}

        if "sirna_id" not in df.columns:
            logger.info("No sirna_id column; falling back to single partition function")
            return self.calculate_probabilities(df)

        # Check if single siRNA (optimization)
        # Only fallback if no custom parameters are active (since standard calc doesn't support them)
        unique_sirnas = df["sirna_id"].unique()
        if (
            len(unique_sirnas) == 1
            and not alpha_gamma_pairs
            and not theta_values
            and not on_target_map
        ):
            logger.info("Single siRNA detected; using standard calculation")
            return self.calculate_probabilities(df)

        logger.info(
            f"Processing {len(unique_sirnas)} siRNAs with per-siRNA partition functions"
        )

        if self.accessibility_service:
            df = self._annotate_opening_energy(df)

        # 2. Calculate dG_total (Base Case: alpha=1.0, gamma=1.0)
        if "opening_energy" in df.columns:
            df = df.with_columns(
                (pl.col("energy") + pl.col("opening_energy")).alias("dG_total")
            )
        else:
            df = df.with_columns(
                pl.lit(0.0).alias("opening_energy"),
                pl.col("energy").alias("dG_total"),
            )

        # 3. Expression value (default 1.0)
        if "exp_value" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("exp_value"))

        # Extract on-target rows if on_target_map is provided
        on_target_data: dict[str, dict] = {}
        if on_target_map and "transcript_id" in df.columns:
            for sirna_id, target_tid in on_target_map.items():
                # Find matching rows for this siRNA and transcript
                matches = df.filter(
                    (pl.col("sirna_id") == sirna_id)
                    & (pl.col("transcript_id") == target_tid)
                )
                if matches.height > 0:
                    # Use best (lowest energy) match as on-target
                    best = matches.sort("energy").head(1).to_dicts()[0]
                    on_target_data[sirna_id] = {
                        "transcript_id": target_tid,
                        "energy": best["energy"],
                        "opening_energy": best.get("opening_energy", 0.0),
                        "dG_total": best.get("dG_total", best["energy"]),
                    }
                    logger.debug(
                        f"On-target for {sirna_id}: {target_tid} (dG={best['energy']:.2f})"
                    )
                else:
                    logger.warning(
                        f"On-target transcript '{target_tid}' not found for siRNA '{sirna_id}'"
                    )

            # Remove on-target rows from off-target set to avoid double-counting
            if on_target_data:
                exclude_conditions = [
                    (pl.col("sirna_id") == sid)
                    & (pl.col("transcript_id") == data["transcript_id"])
                    for sid, data in on_target_data.items()
                ]
                exclude_expr = exclude_conditions[0]
                for cond in exclude_conditions[1:]:
                    exclude_expr = exclude_expr | cond
                df = df.filter(~exclude_expr)
                logger.info(
                    f"Excluded {len(on_target_data)} on-target entries from off-target set"
                )

        # Calculate E_min per siRNA (required for alpha/gamma clamping)
        # We assume E_min is the minimum hybridization energy observed for that siRNA
        min_energies = df.group_by("sirna_id").agg(
            pl.col("energy").min().alias("E_min")
        )
        df = df.join(min_energies, on="sirna_id", how="left")

        # Define list of calculations to perform
        # format: (name_suffix, energy_expression_col)
        calc_configs = [("", "dG_total")]

        # Prepare expressions for additional parameters
        if alpha_gamma_pairs:
            for alpha, gamma in alpha_gamma_pairs:
                # Skip default base case if present
                if alpha == 1.0 and gamma == 1.0:
                    continue

                suffix = f":alpha={alpha},gamma={gamma}"
                col_name = f"dG_total{suffix}"

                # Logic: if energy < alpha * E_min -> use gamma * E_min + open
                # otherwise -> use energy + open
                df = df.with_columns(
                    pl.when(pl.col("energy") < alpha * pl.col("E_min"))
                    .then(gamma * pl.col("E_min") + pl.col("opening_energy"))
                    .otherwise(pl.col("energy") + pl.col("opening_energy"))
                    .alias(col_name)
                )
                calc_configs.append((suffix, col_name))

        if theta_values:
            for theta in theta_values:
                suffix = f":theta={theta}"
                col_name = f"dG_total{suffix}"

                # Logic: ((theta * (energy + 10)) - 10) + open
                df = df.with_columns(
                    (
                        ((theta * (pl.col("energy") + 10.0)) - 10.0)
                        + pl.col("opening_energy")
                    ).alias(col_name)
                )
                calc_configs.append((suffix, col_name))

        # 4. & 5. & 6. Calculate W, Z, and P for all configurations
        # We build a list of aggregations to do them all in one group_by pass for efficiency

        z_aggs = []
        for suffix, dG_col in calc_configs:
            # W = Expr * exp(-dG / RT)
            weight_col = f"boltzmann_weight{suffix}"
            df = df.with_columns(
                (pl.col("exp_value") * ((-pl.col(dG_col) / RT).exp())).alias(weight_col)
            )
            # Z aggregation
            z_col = f"Z_sirna{suffix}"
            z_aggs.append(pl.col(weight_col).sum().alias(z_col))

        # Compute all Zs (parallelized reduction)
        z_df = df.group_by("sirna_id").agg(z_aggs)

        # Add on-target weights to Z for each siRNA
        if on_target_data:
            import numpy as np

            on_target_weights = {}
            for sirna_id, data in on_target_data.items():
                dG_on = data["dG_total"]
                w_on = on_target_expression * np.exp(-dG_on / RT)
                on_target_weights[sirna_id] = w_on

            # Update Z_sirna to include on-target weight
            z_df = z_df.with_columns(
                [
                    (
                        pl.col("Z_sirna")
                        + pl.col("sirna_id").replace_strict(
                            on_target_weights, default=0.0
                        )
                    ).alias("Z_sirna")
                ]
            )
            logger.info(
                f"Added on-target weights for {len(on_target_weights)} siRNAs to partition functions"
            )

        # Join Zs back
        df = df.join(z_df, on="sirna_id", how="left")

        # Compute Probabilities for all configs
        probs_exprs = []
        for suffix, _ in calc_configs:
            weight_col = f"boltzmann_weight{suffix}"
            z_col = f"Z_sirna{suffix}"
            p_col = f"P_off_target{suffix}"

            probs_exprs.append(
                pl.when(pl.col(z_col) > 0)
                .then(pl.col(weight_col) / pl.col(z_col))
                .otherwise(0.0)
                .alias(p_col)
            )

        df = df.with_columns(probs_exprs)

        # Cleanup intermediate columns (optional, but good for memory)
        # We verify by checking if dG_total is base dict
        # We keep the base P_off_target and others.

        # Collect metadata (base Z)
        # For metadata return, we use the base config Z
        z_stats = z_df.select(["sirna_id", "Z_sirna"]).to_dicts()
        metadata = {
            "n_sirnas": len(unique_sirnas),
            "z_per_sirna": {row["sirna_id"]: float(row["Z_sirna"]) for row in z_stats},
            "z_total": float(df["boltzmann_weight"].sum()),
            "on_target_count": len(on_target_data) if on_target_data else 0,
        }

        return df, metadata

    def calculate_legacy_format(
        self,
        df: pl.DataFrame,
        sirna_id: str = "siRNA",
        on_target_path: Optional[Path] = None,
        query_path: Optional[Path] = None,
        on_target_expression: float = 1000.0,
        on_target_accessibility_path: Optional[Path] = None,
        on_target_risearch_path: Optional[Path] = None,
        alpha_gamma_pairs: list[tuple[float, float]] = [(1.0, 1.0), (0.8, 0.8)],
        verbose: bool = False,
    ) -> str:
        """
        Calculate probabilities in legacy output format.

        Groups by transcript and applies multiple alpha/gamma scaling factors.
        Returns formatted string matching old pipeline gw.results format.
        """
        if "energy" not in df.columns:
            return "# Error: No energy column found\n"

        # Annotate with opening energy if service available
        if self.accessibility_service:
            df = self._annotate_opening_energy(df)

        # Ensure opening_energy exists
        if "opening_energy" not in df.columns:
            df = df.with_columns(pl.lit(0.0).alias("opening_energy"))

        # Ensure exp_value exists
        if "exp_value" not in df.columns:
            df = df.with_columns(pl.lit(1.0).alias("exp_value"))

        # Calculate on-target dG components
        dG_hyb_on = 0.0
        dG_open_on = 0.0
        if on_target_path and query_path:
            dG_hyb_on, dG_open_on = self._calculate_on_target_components(
                on_target_path,
                query_path,
                on_target_accessibility_path,
                on_target_risearch_path,
            )

        results_by_ag = {}
        on_target_stats = {}

        # Get minimum energy for clamping logic (old pipeline feature)
        # Old pipeline uses the minimum HYBRIDIZATION energy (noacc) for clamping condition
        # and includes both off-targets and on-target in this minimum finding.

        current_min = dG_hyb_on if (on_target_path and query_path) else 0.0

        if df.height > 0:
            off_min = df["energy"].min()
            if off_min < current_min:
                current_min = off_min

        min_energy = current_min
        # min_energy_noacc is effectively the same as min_energy in this context
        # (both refer to hybridization energy minimum)

        # Filter out on-target row if it exists (added by calculate_probabilities)
        df = df.filter(pl.col("transcript_id") != "onTarget")

        for alpha, gamma in alpha_gamma_pairs:
            # --- WITH ACCESSIBILITY ---
            # OLD PIPELINE FORMULA (for tuple alpha,gamma):
            #   assigned_energy = energy + open_eng
            #   if energy < alpha * minimum: assigned_energy = gamma * minimum + open_eng
            # The alpha/gamma are for CLAMPING, not multiplication!

            df_scaled = df.with_columns(
                [
                    # Apply clamping: if energy < alpha * min_energy, use gamma * min_energy
                    pl.when(pl.col("energy") < alpha * min_energy)
                    .then(gamma * min_energy + pl.col("opening_energy"))
                    .otherwise(pl.col("energy") + pl.col("opening_energy"))
                    .alias("dG_scaled"),
                    # NoAcc version: same clamping but without opening_energy
                    pl.when(pl.col("energy") < alpha * min_energy)
                    .then(gamma * min_energy)
                    .otherwise(pl.col("energy"))
                    .alias("dG_scaled_noacc"),
                ]
            )

            # Boltzmann weight per row (WITH ACC)
            df_scaled = df_scaled.with_columns(
                [
                    (pl.col("exp_value") * ((-pl.col("dG_scaled") / RT).exp())).alias(
                        "W"
                    ),
                    (
                        pl.col("exp_value") * ((-pl.col("dG_scaled_noacc") / RT).exp())
                    ).alias("W_noacc"),
                ]
            )

            # Aggregate by transcript (sum weights)
            df_agg = df_scaled.group_by(["chrom", "transcript_id"]).agg(
                [
                    pl.col("W").sum().alias("W_transcript"),
                    pl.col("W_noacc").sum().alias("W_transcript_noacc"),
                ]
            )

            # Z_off = sum of all off-target weights
            z_off = df_agg["W_transcript"].sum()
            z_off_noacc = df_agg["W_transcript_noacc"].sum()

            # On-target weight (apply same clamping logic)
            if dG_hyb_on < alpha * min_energy:
                dG_total_on = gamma * min_energy + dG_open_on
                dG_on_noacc = gamma * min_energy
            else:
                dG_total_on = dG_hyb_on + dG_open_on
                dG_on_noacc = dG_hyb_on

            w_on = on_target_expression * np.exp(-dG_total_on / RT)
            w_on_noacc = on_target_expression * np.exp(-dG_on_noacc / RT)

            # Total partition function
            z_total = z_off + w_on
            z_total_noacc = z_off_noacc + w_on_noacc

            # P(on) and P(off)
            p_on = w_on / z_total if z_total > 0 else 0.0
            p_off_total = z_off / z_total if z_total > 0 else 0.0
            p_on_noacc = w_on_noacc / z_total_noacc if z_total_noacc > 0 else 0.0
            p_off_total_noacc = (
                z_off_noacc / z_total_noacc if z_total_noacc > 0 else 0.0
            )

            on_target_stats[(alpha, gamma)] = {
                "Pon": p_on,
                "Pon_noacc": p_on_noacc,
                "Poff": p_off_total,
                "Poff_noacc": p_off_total_noacc,
                "Z": w_on,
                "Z_noacc": w_on_noacc,
                "Zoff": z_off,
                "Zoff_noacc": z_off_noacc,
            }

            # Per-transcript probabilities (With Acc and NoAcc)
            df_agg = df_agg.with_columns(
                [
                    (pl.col("W_transcript") / z_total).alias("P_transcript")
                    if z_total > 0
                    else pl.lit(0.0).alias("P_transcript"),
                    (pl.col("W_transcript_noacc") / z_total_noacc).alias(
                        "P_transcript_noacc"
                    )
                    if z_total_noacc > 0
                    else pl.lit(0.0).alias("P_transcript_noacc"),
                ]
            )

            results_by_ag[(alpha, gamma)] = df_agg

        # Build output string
        lines = []

        # On-target header
        lines.append(f"# On-target info for {sirna_id} #")
        for (alpha, gamma), stats in on_target_stats.items():
            line = f"# For alpha={alpha} and gamma={gamma}; "
            line += f"Pon: {stats['Pon']:.12g}; Pon_noacc: {stats['Pon_noacc']:.12g}; "
            line += (
                f"Poff: {stats['Poff']:.12g}; Poff_noacc: {stats['Poff_noacc']:.12g}; "
            )
            line += f"Z: {stats['Z']:.12g}; Z_noacc: {stats['Z_noacc']:.12g}; "
            line += (
                f"Zoff: {stats['Zoff']:.12g}; Zoff_noacc: {stats['Zoff_noacc']:.12g}"
            )
            lines.append(line)
        lines.append("## End of on-target info ##")

        # Off-target section
        if verbose:
            lines.append("# Off-target info#")

            # Column headers
            header_parts = ["# ID1", "ID2"]
            for alpha, gamma in alpha_gamma_pairs:
                header_parts.append(f"alpha:{alpha},gamma:{gamma}")
                header_parts.append(f"alpha:{alpha},gamma:{gamma}(NoAcc)")
            lines.append("\t".join(header_parts))

            # Merge all (alpha, gamma) results
            # Calculate base set from first pair
            first_ag = alpha_gamma_pairs[0]
            merged = (
                results_by_ag[first_ag]
                .select(
                    ["chrom", "transcript_id", "P_transcript", "P_transcript_noacc"]
                )
                .rename(
                    {
                        "P_transcript": f"P_{first_ag[0]}_{first_ag[1]}",
                        "P_transcript_noacc": f"P_{first_ag[0]}_{first_ag[1]}_noacc",
                    }
                )
            )

            for ag in alpha_gamma_pairs[1:]:
                other = (
                    results_by_ag[ag]
                    .select(
                        ["chrom", "transcript_id", "P_transcript", "P_transcript_noacc"]
                    )
                    .rename(
                        {
                            "P_transcript": f"P_{ag[0]}_{ag[1]}",
                            "P_transcript_noacc": f"P_{ag[0]}_{ag[1]}_noacc",
                        }
                    )
                )
                merged = merged.join(other, on=["chrom", "transcript_id"], how="outer")

            # Sort by chrom, transcript_id
            merged = merged.sort(["chrom", "transcript_id"])

            # Output rows
            for row in merged.iter_rows(named=True):
                parts = [row["chrom"], row["transcript_id"]]
                for ag in alpha_gamma_pairs:
                    # With Acc
                    col = f"P_{ag[0]}_{ag[1]}"
                    val = row.get(col, 0.0) or 0.0
                    parts.append(f"{val:.12g}")

                    # No Acc
                    col_noacc = f"P_{ag[0]}_{ag[1]}_noacc"
                    val_noacc = row.get(col_noacc, 0.0) or 0.0
                    parts.append(f"{val_noacc:.12g}")

                lines.append("\t".join(parts))

            lines.append("## End of off-target info")

        lines.append("")

        return "\n".join(lines)

    def _annotate_opening_energy(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Lookup opening energy for each row using the AccessibilityService.
        """
        if df.height == 0:
            return df

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
        # Exclude 'onTarget' as it is handled separately
        targets = (
            df.filter(pl.col("chrom") != "onTarget")
            .select(["chrom", "start", "end", "strand"])
            .unique()
        )

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
            start_in = row["start"]
            end_in = row["end"]
            strand = row["strand"]

            # Convert 1-based RIsearch coords to 0-based for internal query
            start0 = start_in - 1

            try:
                # Query region
                opening_energies = self.accessibility_service.query(
                    chrom, start0, end_in, strand
                )

                # Interaction length 'u'
                interaction_len = end_in - start0

                # Check for 2D (matrix) profile
                if opening_energies.ndim == 2:
                    # Matrix: rows = positions, cols = u-lengths
                    # We want the column corresponding to our interaction length 'u'
                    # col_idx = u - 1
                    matrix_width = opening_energies.shape[1]
                    col_idx = min(interaction_len, matrix_width) - 1
                    if col_idx < 0:
                        col_idx = 0

                    # Pick 3' end value
                    if len(opening_energies) > 0:
                        row_idx = -1 if strand == "+" else 0
                        val = opening_energies[row_idx, col_idx]
                        mean_e_open = int(round(float(val) * 10.0)) / 10.0
                    else:
                        mean_e_open = 10.0
                else:
                    # 1D Legacy/Simplified (Assumes fixed large u, e.g. 30)
                    if len(opening_energies) > 0:
                        val = (
                            opening_energies[-1]
                            if strand == "+"
                            else opening_energies[0]
                        )
                        mean_e_open = int(round(float(val) * 10.0)) / 10.0
                    else:
                        mean_e_open = 10.0

                results.append(
                    {
                        "chrom": chrom,
                        "start": start_in,
                        "end": end_in,
                        "strand": strand,
                        "opening_energy": mean_e_open,
                    }
                )
            except AccessibilityError as e:
                # Log first few errors per chromosome to avoid spam
                result_key = f"acc_error_{chrom}"
                if not hasattr(self, "_logged_errors"):
                    self._logged_errors = set()

                if result_key not in self._logged_errors:
                    logger.warning(
                        f"Accessibility lookup failed for {chrom}:{start_in}-{end_in} ({strand}): {e}. Defaulting to 10.0"
                    )
                    self._logged_errors.add(result_key)

                results.append(
                    {
                        "chrom": chrom,
                        "start": start_in,
                        "end": end_in,
                        "strand": strand,
                        "opening_energy": 10.0,
                    }
                )

        # Join back
        df_acc = pl.DataFrame(results)

        # Join on keys
        return df.join(df_acc, on=["chrom", "start", "end", "strand"], how="left")

    def _calculate_on_target_dg(
        self,
        on_target_path: Path,
        query_path: Path,
        accessibility_path: Optional[Path] = None,
        risearch_path: Optional[Path] = None,
    ) -> float:
        """
        Calculate dG_total for the on-target sequence.
        """
        dG_hyb, dG_open = self._calculate_on_target_components(
            on_target_path, query_path, accessibility_path, risearch_path
        )
        return dG_hyb + dG_open

    def _calculate_on_target_components(
        self,
        on_target_path: Path,
        query_path: Path,
        accessibility_path: Optional[Path] = None,
        risearch_path: Optional[Path] = None,
    ) -> tuple[float, float]:
        """
        Calculate (dG_hyb, dG_open) for On-Target.

        If accessibility_path is not provided, computes accessibility on-the-fly
        using ViennaRNA's RNA.pfl_fold_up.
        """
        # 1. Hybridization Energy + Coords
        if risearch_path:
            dG_hyb, t_start, t_end, strand = self._parse_risearch_file(risearch_path)
        else:
            dG_hyb, t_start, t_end, strand = self._run_risearch_binary(
                query_path, on_target_path
            )

        # 2. Accessibility Energy
        dG_open = 0.0

        if accessibility_path and accessibility_path.exists() and t_end > 0:
            # Use pre-computed accessibility file
            try:
                # Manual parsing of the accessibility file (Text or Binary)
                # Since we don't have a service instance pointing here necessarily.
                # Use AccessibilityService static-like logic?
                # For now, implement simplified text parser if file is text.
                # Assuming text format from tests/output/old_accessibility.

                # We need to reuse the sophisticated logic from AccessibilityService for matrix/coords.
                # Best way: Check if we have a service. If not, separate logic.

                # Use a temp service to parse/query
                temp_service = GenomeAccessibilityService(accessibility_path.parent)
                # Force load this specific file into cache with a known key
                # Query using the file stem as "chrom"?
                # If accessibility_path is "path/to/onTarget_openen", stem is "onTarget_openen".
                # The user "onTarget" ID might be "onTarget".
                # RIsearch target ID is "onTarget" (from FASTA).

                # Hack: Determine "chrom" name expected by query.
                # The _run_risearch_binary returns target ID e.g. "onTarget".

                # If we use temp_service.query("onTarget", ...):
                # It will look for various files.
                # If we rename/symlink? No.

                # Direct parse:
                if "openen" in accessibility_path.name:
                    profile = temp_service._parse_openen_text(accessibility_path)
                elif accessibility_path.suffix == ".npy":
                    profile = np.load(accessibility_path, mmap_mode="r")
                elif accessibility_path.suffix == ".bin":
                    # Legacy binary format: flattened 2D matrix (Nx30)
                    raw_data = np.memmap(accessibility_path, dtype=np.uint8, mode="r")
                    # Reshape to 2D: each row has 30 columns (u-values from 1 to 30)
                    if len(raw_data) % 30 == 0:
                        profile = raw_data.reshape(-1, 30)
                    else:
                        logger.warning(
                            "Binary file size not divisible by 30, using as 1D"
                        )
                        profile = raw_data
                else:
                    profile = np.array([])

                # Now replicate the _annotate_opening_energy logic
                # Profile is likely 2D (Nx30) or 1D.
                # Coords: t_start (1-based), t_end (1-based).
                start0 = t_start - 1
                interaction_len = t_end - start0

                if profile.size > 0:
                    if profile.ndim == 2:
                        matrix_width = profile.shape[1]
                        col_idx = min(interaction_len, matrix_width) - 1
                        if col_idx < 0:
                            col_idx = 0

                        # Bounds check
                        if start0 >= 0 and t_end <= profile.shape[0]:
                            # Endpoint index
                            # If +, end of site (t_end-1).
                            # If -, start of site (start0).
                            target_idx = (t_end - 1) if strand == "+" else start0

                            if 0 <= target_idx < profile.shape[0]:
                                val = profile[target_idx, col_idx]
                                # Binary files store uint8 * 10, so divide by 10
                                if profile.dtype == np.uint8:
                                    dG_open = int(round(float(val))) / 10.0
                                else:
                                    dG_open = int(round(float(val) * 10.0)) / 10.0
                            else:
                                dG_open = 10.0
                        else:
                            dG_open = 10.0
                    else:
                        # 1D
                        target_idx = (t_end - 1) if strand == "+" else start0
                        if 0 <= target_idx < profile.shape[0]:
                            val = profile[target_idx]
                            # If legacy bin (uint8)
                            if profile.dtype == np.uint8:
                                val = float(val) / 10.0

                            dG_open = int(round(float(val) * 10.0)) / 10.0
                        else:
                            dG_open = 10.0

            except Exception as e:
                logger.error(
                    f"Failed to calculate on-target accessibility from file: {e}"
                )
                dG_open = 0.0  # Fallback

        elif t_end > 0:
            # Compute accessibility on-the-fly from the on-target FASTA
            try:
                from RIsearch_pipeline.services.helpers import read_fasta

                # Read sequence from on-target FASTA
                sequence = None
                for seq_id, seq in read_fasta(on_target_path):
                    sequence = seq
                    break  # Use first sequence

                if sequence:
                    logger.info(
                        f"Computing on-target accessibility on-the-fly (len={len(sequence)})"
                    )

                    # Create temp service for computation
                    import tempfile

                    with tempfile.TemporaryDirectory() as tmpdir:
                        temp_service = GenomeAccessibilityService(Path(tmpdir))
                        profile = temp_service.compute_sequence_accessibility(
                            sequence, use_cli=self.use_rnaplfold_cli
                        )

                        # Extract opening energy at binding site
                        start0 = t_start - 1
                        interaction_len = t_end - start0

                        if profile.size > 0 and profile.ndim == 2:
                            col_idx = min(interaction_len, profile.shape[1]) - 1
                            if col_idx < 0:
                                col_idx = 0

                            # Use appropriate endpoint based on strand
                            target_idx = (t_end - 1) if strand == "+" else start0

                            if 0 <= target_idx < profile.shape[0]:
                                dG_open = float(profile[target_idx, col_idx])
                                logger.info(
                                    f"On-target opening energy: {dG_open:.2f} kcal/mol"
                                )
                            else:
                                logger.warning(
                                    f"Binding site index {target_idx} out of bounds (len={profile.shape[0]})"
                                )
                                dG_open = 10.0
                        else:
                            logger.warning("Empty or 1D profile from ViennaRNA")
                            dG_open = 10.0
                else:
                    logger.warning(
                        f"No sequence found in on-target FASTA: {on_target_path}"
                    )

            except ImportError:
                logger.warning(
                    "ViennaRNA not available for on-the-fly accessibility computation"
                )
            except Exception as e:
                logger.error(f"Failed to compute on-target accessibility: {e}")
                dG_open = 0.0

        return dG_hyb, dG_open

    def _run_risearch_binary(
        self, query_path: Path, target_path: Path
    ) -> tuple[float, int, int, str]:
        """
        Run RIsearch binary to get hybridization energy and coordinates.

        Delegates to RIsearchService for actual execution.

        Returns:
            Tuple of (energy, start, end, strand) for the best hit.
            Returns (0.0, 0, 0, "+") if no hits found.
        """
        from RIsearch_pipeline.services.risearch_service import RIsearchService

        service = RIsearchService()
        return service.search_single_sirna(query_path, target_path)

    def _parse_risearch_file(self, path: Path) -> tuple[float, int, int, str]:
        """
        Parse a pre-computed RIsearch output file.
        Returns (best_energy, start, end, strand).
        """
        if not path.exists():
            logger.error(f"Provided RIsearch file not found: {path}")
            return 0.0, 0, 0, "+"

        energies = []
        try:
            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    # Output format (8 cols): QueryID, QStart, QEnd, TargetID, TStart, TEnd, Strand, Energy
                    parts = line.split()
                    if len(parts) >= 8:
                        try:
                            energy = float(parts[7])
                            t_start = int(parts[4])
                            t_end = int(parts[5])
                            strand = parts[6]
                            energies.append((energy, t_start, t_end, strand))
                        except ValueError:
                            pass
        except Exception as e:
            logger.error(f"Error reading RIsearch file {path}: {e}")
            return 0.0, 0, 0, "+"

        if energies:
            # Find min energy
            return min(energies, key=lambda x: x[0])
        else:
            logger.warning(f"No valid predictions found in {path}")
            return 0.0, 0, 0, "+"
