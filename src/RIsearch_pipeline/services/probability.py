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
        precomputed_accessibility: Optional[pl.DataFrame] = None,
    ):
        self.accessibility_service = accessibility_service
        self.precomputed_accessibility = precomputed_accessibility

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
            "z_total": z_total,
            "z_off_target": z_partial,
            "w_on_target": w_on,
            "dG_on_target": dG_total_on,
            "has_on_target": bool(on_target_path and query_path),
        }

        # Add P_on to metadata if applicable
        if on_target_path and query_path:
            metadata["p_on_target"] = w_on / z_total if z_total > 0 else 0.0

        return df, metadata

    def calculate_probabilities_per_sirna(
        self,
        df: pl.DataFrame,
        alpha_gamma_pairs: Optional[list[tuple[float, float]]] = None,
        theta_values: Optional[list[float]] = None,
        on_target_map: Optional[dict[str, str]] = None,
        on_target_expression: float = 1000.0,
        on_target_fasta_data: Optional[dict[str, dict]] = None,
    ) -> Tuple[pl.DataFrame, Dict]:
        """
        Calculate P(OT) per siRNA using per-siRNA partition functions.

        Two on-target modes:

        1. **In-predictions mode** (on_target_map provided, no fasta_data):
           Finds on-target binding sites within the predictions by gene/transcript
           ID matching. Rows stay in the DataFrame and use their BED expression.
           Matches old pipeline behavior.

        2. **FASTA mode** (on_target_fasta_data provided):
           On-target energy is computed separately (from FASTA + RIsearch/ViennaRNA).
           A fixed on_target_expression is used. The on-target weight is added
           to Z as a separate term.

        Args:
            df: DataFrame with predictions.
            alpha_gamma_pairs: List of (alpha, gamma) tuples for clamping.
            theta_values: List of theta scaling values.
            on_target_map: Mapping of siRNA_id -> on_target_gene_or_transcript_id.
            on_target_expression: Expression for FASTA mode on-targets (default 1000.0).
            on_target_fasta_data: Dict of {siRNA_id: {energy, opening_energy}} from
                separate FASTA computation. When provided, activates FASTA mode.

        Returns:
            Tuple of (DataFrame with P_off_target per siRNA, metadata dict).

        Complexity: O(N) time where N = number of predictions.
        """
        if "energy" not in df.columns:
            logger.warning("No 'energy' column found. Cannot calculate probabilities.")
            return df, {}

        if "sirna_id" not in df.columns:
            logger.info("No sirna_id column; falling back to single partition function")
            return self.calculate_probabilities(df)

        # Check if single siRNA (optimization)
        unique_sirnas = df["sirna_id"].unique()
        if (
            len(unique_sirnas) == 1
            and not alpha_gamma_pairs
            and not theta_values
            and not on_target_map
            and not on_target_fasta_data
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

        # --- Determine on-target mode ---
        # Mode A: In-predictions (old pipeline behavior) — tag rows, keep in DF
        # Mode B: FASTA mode — remove rows and re-add with fixed expression
        use_fasta_mode = on_target_fasta_data is not None
        in_predictions_mode = (
            on_target_map is not None
            and not use_fasta_mode
            and "transcript_id" in df.columns
        )

        # Extracted on-target data for FASTA mode only
        on_target_data_fasta: dict[str, dict] = {}

        if in_predictions_mode:
            # --- Mode A: Tag on-target rows in place ---
            # Build a condition that matches on-target rows by gene_id or transcript_id
            # The old pipeline checks: gid == on_id OR tid == on_id
            #
            # Filter map to siRNAs present in this DataFrame before building
            # expressions — avoids an O(N_corpus)-deep OR chain when called
            # per-siRNA from a worker (each df has exactly one unique siRNA_id).
            sirna_set = set(unique_sirnas.to_list())
            local_map = {k: v for k, v in on_target_map.items() if k in sirna_set}

            tag_conditions = []
            for sirna_id, target_id in local_map.items():
                cond = pl.col("sirna_id") == sirna_id
                id_match = pl.col("transcript_id") == target_id
                if "gene_id" in df.columns:
                    id_match = id_match | (pl.col("gene_id") == target_id)
                tag_conditions.append(cond & id_match)

            if tag_conditions:
                combined = tag_conditions[0]
                for c in tag_conditions[1:]:
                    combined = combined | c
                df = df.with_columns(combined.alias("is_on_target"))
            else:
                df = df.with_columns(pl.lit(False).alias("is_on_target"))

            n_on = df.filter(pl.col("is_on_target")).height
            logger.info(
                f"In-predictions on-target mode: tagged {n_on} rows as on-target "
                f"(from {len(local_map)} mappings)"
            )

        elif use_fasta_mode and on_target_map and "transcript_id" in df.columns:
            # --- Mode B: FASTA mode — extract and remove on-target rows ---
            on_target_data_fasta = on_target_fasta_data

            # Remove on-target rows from off-target set to avoid double-counting
            exclude_conditions = [
                (pl.col("sirna_id") == sid)
                & (
                    (pl.col("transcript_id") == on_target_map[sid])
                    | (
                        pl.col("gene_id") == on_target_map[sid]
                        if "gene_id" in df.columns
                        else pl.lit(False)
                    )
                )
                for sid in on_target_fasta_data
                if sid in on_target_map
            ]
            if exclude_conditions:
                exclude_expr = exclude_conditions[0]
                for cond in exclude_conditions[1:]:
                    exclude_expr = exclude_expr | cond
                df = df.filter(~exclude_expr)
                logger.info(
                    f"FASTA mode: excluded on-target rows for {len(on_target_fasta_data)} siRNAs"
                )

            df = df.with_columns(pl.lit(False).alias("is_on_target"))
        else:
            df = df.with_columns(pl.lit(False).alias("is_on_target"))

        # Calculate E_min per siRNA (required for alpha/gamma clamping).
        # Replicate old-pipeline behaviour: anchor clamping to the minimum
        # energy from the RAW RIsearch file (all genome hits, including those
        # that don't overlap any annotated transcript).  The parser attaches
        # this as `raw_e_min` when loading per-siRNA directory files.
        # Fall back to the post-intersection minimum when `raw_e_min` is absent
        # (e.g. when loading a single merged file).
        if "raw_e_min" in df.columns:
            min_energies = df.group_by("sirna_id").agg(
                pl.col("raw_e_min").min().alias("E_min")
            )
        else:
            min_energies = df.group_by("sirna_id").agg(
                pl.col("energy").min().alias("E_min")
            )
        df = df.join(min_energies, on="sirna_id", how="left")

        # Define list of calculations to perform
        # format: (name_suffix, energy_expression_col, noacc_energy_expression_col)
        calc_configs = [("", "dG_total", "energy")]

        # Collect all dG variant expressions so they can be applied in a single
        # with_columns() call (one pass over the data regardless of parameter count).
        all_dG_exprs = []

        if alpha_gamma_pairs:
            for alpha, gamma in alpha_gamma_pairs:
                # Skip default base case if present
                if alpha == 1.0 and gamma == 1.0:
                    continue

                suffix = f":alpha={alpha},gamma={gamma}"
                col_name = f"dG_total{suffix}"
                col_name_noacc = f"energy{suffix}"

                # Logic: if energy < alpha * E_min -> use gamma * E_min + open
                # otherwise -> use energy + open
                all_dG_exprs.extend([
                    pl.when(pl.col("energy") < alpha * pl.col("E_min"))
                    .then(gamma * pl.col("E_min") + pl.col("opening_energy"))
                    .otherwise(pl.col("energy") + pl.col("opening_energy"))
                    .alias(col_name),
                    pl.when(pl.col("energy") < alpha * pl.col("E_min"))
                    .then(gamma * pl.col("E_min"))
                    .otherwise(pl.col("energy"))
                    .alias(col_name_noacc),
                ])
                calc_configs.append((suffix, col_name, col_name_noacc))

        if theta_values:
            for theta in theta_values:
                suffix = f":theta={theta}"
                col_name = f"dG_total{suffix}"
                col_name_noacc = f"energy{suffix}"

                # Logic: ((theta * (energy + 10)) - 10) + open
                all_dG_exprs.extend([
                    (
                        ((theta * (pl.col("energy") + 10.0)) - 10.0)
                        + pl.col("opening_energy")
                    ).alias(col_name),
                    ((theta * (pl.col("energy") + 10.0)) - 10.0).alias(col_name_noacc),
                ])
                calc_configs.append((suffix, col_name, col_name_noacc))

        # Single pass for all dG variant columns
        if all_dG_exprs:
            df = df.with_columns(all_dG_exprs)

        # 4. & 5. & 6. Calculate W, Z, and P for all configurations
        # Build all weight expressions and aggregations up-front, then apply in
        # two single passes (weights) and one group_by (Z + on-target combined).

        all_weight_exprs = []
        z_aggs = []
        on_aggs = []  # populated below when in_predictions_mode
        z_col_names = []

        for suffix, dG_col, dG_noacc_col in calc_configs:
            # W = Expr * exp(-dG / RT)
            weight_col = f"boltzmann_weight{suffix}"
            weight_noacc_col = f"boltzmann_weight_noacc{suffix}"
            all_weight_exprs.extend([
                (pl.col("exp_value") * ((-pl.col(dG_col) / RT).exp())).alias(weight_col),
                (pl.col("exp_value") * ((-pl.col(dG_noacc_col) / RT).exp())).alias(weight_noacc_col),
            ])
            # Z aggregation
            z_col = f"Z_sirna{suffix}"
            z_noacc_col = f"Z_sirna_noacc{suffix}"
            z_aggs.extend([
                pl.col(weight_col).sum().alias(z_col),
                pl.col(weight_noacc_col).sum().alias(z_noacc_col),
            ])
            z_col_names.extend([z_col, z_noacc_col])

            if in_predictions_mode:
                w_on_col = f"W_on{suffix}"
                w_on_noacc_col = f"W_on_noacc{suffix}"
                on_aggs.extend([
                    pl.col(weight_col).filter(pl.col("is_on_target")).sum().alias(w_on_col),
                    pl.col(weight_noacc_col).filter(pl.col("is_on_target")).sum().alias(w_on_noacc_col),
                ])

        # Single pass for all Boltzmann weight columns
        df = df.with_columns(all_weight_exprs)

        # Single group_by for Z sums and (when applicable) on-target weight sums
        agg_df = df.group_by("sirna_id").agg(z_aggs + on_aggs)
        z_df = agg_df.select(["sirna_id"] + z_col_names)

        # --- Compute on-target weights ---
        on_target_weight_metadata = {}

        if in_predictions_mode:
            # Mode A: On-target rows are IN the DataFrame.
            # Z already includes them from the group_by sum above.
            # Extract W_on columns from the combined aggregation result.
            on_df = agg_df.drop(z_col_names)

            for row in on_df.to_dicts():
                sid = row["sirna_id"]
                on_target_weight_metadata[sid] = {
                    k: v for k, v in row.items() if k != "sirna_id"
                }

            logger.info(
                f"In-predictions mode: computed on-target weights for "
                f"{on_df.height} siRNAs"
            )

        elif on_target_data_fasta:
            # Mode B: FASTA mode — add separate W_on to Z
            import numpy as np

            min_e_dict = {
                row["sirna_id"]: row["E_min"] for row in min_energies.to_dicts()
            }

            z_modifications = []

            for suffix, dG_col, dG_noacc_col in calc_configs:
                on_target_w = {}
                on_target_w_noacc = {}

                for sirna_id, data in on_target_data_fasta.items():
                    orig_energy = data["energy"]
                    orig_open = data.get("opening_energy", 0.0)
                    e_min = min_e_dict.get(sirna_id, orig_energy)
                    if orig_energy < e_min:
                        e_min = orig_energy

                    if suffix.startswith(":alpha="):
                        parts = suffix.split(",")
                        alpha = float(parts[0].split("=")[1])
                        gamma = float(parts[1].split("=")[1])
                        if orig_energy < alpha * e_min:
                            assigned_e = gamma * e_min + orig_open
                            assigned_e_noacc = gamma * e_min
                        else:
                            assigned_e = orig_energy + orig_open
                            assigned_e_noacc = orig_energy
                    elif suffix.startswith(":theta="):
                        theta = float(suffix.split("=")[1])
                        assigned_e = ((theta * (orig_energy + 10.0)) - 10.0) + orig_open
                        assigned_e_noacc = (theta * (orig_energy + 10.0)) - 10.0
                    else:
                        assigned_e = orig_energy + orig_open
                        assigned_e_noacc = orig_energy

                    w_on = on_target_expression * np.exp(-assigned_e / RT)
                    w_on_noacc = on_target_expression * np.exp(-assigned_e_noacc / RT)

                    on_target_w[sirna_id] = w_on
                    on_target_w_noacc[sirna_id] = w_on_noacc

                    if sirna_id not in on_target_weight_metadata:
                        on_target_weight_metadata[sirna_id] = {}
                    on_target_weight_metadata[sirna_id][f"W_on{suffix}"] = w_on
                    on_target_weight_metadata[sirna_id][f"W_on_noacc{suffix}"] = (
                        w_on_noacc
                    )

                z_col = f"Z_sirna{suffix}"
                z_noacc_col = f"Z_sirna_noacc{suffix}"

                z_modifications.extend(
                    [
                        (
                            pl.col(z_col)
                            + pl.col("sirna_id").replace_strict(
                                on_target_w, default=0.0
                            )
                        ).alias(z_col),
                        (
                            pl.col(z_noacc_col)
                            + pl.col("sirna_id").replace_strict(
                                on_target_w_noacc, default=0.0
                            )
                        ).alias(z_noacc_col),
                    ]
                )

            z_df = z_df.with_columns(z_modifications)

            logger.info(
                f"FASTA mode: added on-target weights for "
                f"{len(on_target_data_fasta)} siRNAs to partition functions"
            )

        # Join Zs back
        df = df.join(z_df, on="sirna_id", how="left")

        # Compute Probabilities for all configs
        probs_exprs = []
        for suffix, _, _ in calc_configs:
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

        # Collect metadata (FULL Z DATAFRAME)
        z_stats = z_df.to_dicts()
        z_stats_dict = {row["sirna_id"]: row for row in z_stats}

        metadata = {
            "n_sirnas": len(unique_sirnas),
            "z_per_sirna": z_stats_dict,
            "on_target_weights": on_target_weight_metadata,
            "on_target_count": len(on_target_weight_metadata),
        }

        return df, metadata

    def format_legacy_summary(
        self,
        sirna_id: str,
        z_per_config: dict[str, float],
        on_target_weights: dict[str, float],
        alpha_gamma_pairs: list[tuple[float, float]],
        theta_values: list[float],
    ) -> str:
        """Format the aggregated metadata into the exact legacy .off header block."""
        lines = []
        lines.append(f"# On-target info for {sirna_id} #")

        suffixes = [("", 1.0, 1.0)]
        for theta in theta_values or []:
            suffixes.append((f":theta={theta}", theta, None))
        for alpha, gamma in alpha_gamma_pairs or []:
            if alpha == 1.0 and gamma == 1.0:
                continue
            suffixes.append((f":alpha={alpha},gamma={gamma}", alpha, gamma))

        for suffix_info in suffixes:
            suffix = suffix_info[0]

            if suffix == "":
                preamble = "# For alpha=1.0 and gamma=1.0;"
            elif suffix.startswith(":theta="):
                preamble = f"# For theta={suffix_info[1]};"
            else:
                preamble = f"# For alpha={suffix_info[1]} and gamma={suffix_info[2]};"

            z = z_per_config.get(f"Z_sirna{suffix}", 0.0)
            z_noacc = z_per_config.get(f"Z_sirna_noacc{suffix}", 0.0)
            w_on = on_target_weights.get(f"W_on{suffix}", 0.0)
            w_on_noacc = on_target_weights.get(f"W_on_noacc{suffix}", 0.0)

            z_off = z - w_on
            z_off_noacc = z_noacc - w_on_noacc

            p_on = w_on / z if z > 0 else 0.0
            p_on_noacc = w_on_noacc / z_noacc if z_noacc > 0 else 0.0

            lines.append(
                f"{preamble} P: {p_on!s}; P_noacc: {p_on_noacc!s}; "
                f"Z: {z!s}; Z_noacc: {z_noacc!s}; "
                f"Zoff: {z_off!s}; Zoff_noacc: {z_off_noacc!s}"
            )

        lines.append("## End of on-target info ##\n")
        return "\n".join(lines)

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
        Lookup opening energy for each row.

        Two paths:
        1. Precomputed Parquet (fast): single Polars join on (chrom, start, end, strand).
        2. Profile-based (LRU cache): per-row query_single() with group-by-chromosome locality.

        Complexity:
            Path 1: O(N log N) for hash join.
            Path 2: O(N) with O(1) per lookup via mmap.
        """
        if df.height == 0:
            return df

        if not all(c in df.columns for c in ["chrom", "start", "end", "strand"]):
            logger.warning(
                "Missing coordinate or strand columns. Cannot lookup accessibility."
            )
            return df

        # --- Path 1: Precomputed Parquet join ---
        if self.precomputed_accessibility is not None:
            acc_df = self.precomputed_accessibility
            join_keys = ["chrom", "start", "end", "strand"]
            result = df.join(
                acc_df.select(join_keys + ["opening_energy"]),
                on=join_keys,
                how="left",
            )
            # Fill nulls (sites not in precomputed file) with default penalty
            if "opening_energy" in result.columns:
                result = result.with_columns(pl.col("opening_energy").fill_null(10.0))
            else:
                result = result.with_columns(pl.lit(10.0).alias("opening_energy"))
            logger.info(
                f"Annotated {result.height} rows with precomputed accessibility"
            )
            return result

        # --- Path 2: Profile-based lookup (vectorized numpy gather per chrom/strand) ---
        if not self.accessibility_service:
            return df

        return self.accessibility_service.annotate_opening_energy_vectorized(df)

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
            # Use pre-computed accessibility Parquet file
            try:
                import polars as pl

                df_acc = (
                    pl.read_parquet(accessibility_path)
                    .filter(pl.col("strand") == strand)
                    .sort("position")
                )
                u_cols = sorted(
                    [c for c in df_acc.columns if c.startswith("u") and c[1:].isdigit()],
                    key=lambda x: int(x[1:]),
                )
                max_pos = int(df_acc["position"].max())
                n_u = len(u_cols)
                profile = np.full((max_pos, n_u), 25.5, dtype=np.float32)
                pos_0 = df_acc["position"].to_numpy() - 1
                for col_idx, col_name in enumerate(u_cols):
                    profile[pos_0, col_idx] = df_acc[col_name].to_numpy()

                start0 = t_start - 1
                interaction_len = t_end - start0
                matrix_width = profile.shape[1]
                col_idx = min(interaction_len, matrix_width) - 1
                if col_idx < 0:
                    col_idx = 0
                target_idx = (t_end - 1) if strand == "+" else start0
                if 0 <= target_idx < profile.shape[0]:
                    dG_open = int(round(float(profile[target_idx, col_idx]) * 10.0)) / 10.0
                else:
                    dG_open = 10.0

            except Exception as e:
                logger.error(
                    f"Failed to read on-target accessibility parquet: {e}"
                )
                dG_open = 0.0

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
                            sequence
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
