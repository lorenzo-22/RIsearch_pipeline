#!/usr/bin/env python3
"""
Profiling script for parallelization comparison.
Run with: uv run scalene tests/off-targets/profile_benchmark.py
"""

import time
import polars as pl
from concurrent.futures import ThreadPoolExecutor
from RIsearch_pipeline.services.risearch_parser import RIsearchParser
from RIsearch_pipeline.services.probability import RT


def polars_native(df: pl.DataFrame) -> pl.DataFrame:
    """Polars Native approach."""
    # Step 1: Add dG_total
    df = df.with_columns(pl.col("energy").alias("dG_total"))

    # Step 2: Boltzmann weights
    df = df.with_columns(((-pl.col("dG_total") / RT).exp()).alias("boltzmann_weight"))

    # Step 3: Group by siRNA, compute Z
    z_per_sirna = df.group_by("sirna_id").agg(
        pl.col("boltzmann_weight").sum().alias("Z_sirna")
    )

    # Step 4: Join Z back
    df = df.join(z_per_sirna, on="sirna_id", how="left")

    # Step 5: Calculate P_off_target
    df = df.with_columns(
        (pl.col("boltzmann_weight") / pl.col("Z_sirna")).alias("P_off_target")
    )
    return df


def threads_approach(df: pl.DataFrame, workers: int = 4) -> pl.DataFrame:
    """Threads approach."""

    def process_group(group_df: pl.DataFrame) -> pl.DataFrame:
        group_df = group_df.with_columns(pl.col("energy").alias("dG_total"))
        group_df = group_df.with_columns(
            ((-pl.col("dG_total") / RT).exp()).alias("boltzmann_weight")
        )
        z = group_df["boltzmann_weight"].sum()
        return group_df.with_columns(
            (pl.col("boltzmann_weight") / z).alias("P_off_target")
        )

    # Step 1: Partition
    groups = df.partition_by("sirna_id", maintain_order=True)

    # Step 2: Execute in parallel
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(process_group, groups))

    # Step 3: Concatenate
    return pl.concat(results)


if __name__ == "__main__":
    # Load data
    parser = RIsearchParser()
    df = parser.load("tests/off-targets/10k_dataset/risearch_10k.out")
    print(f"Dataset: {len(df):,} predictions, {df['sirna_id'].n_unique()} siRNAs")

    # Profile Polars Native
    print("\n--- Polars Native ---")
    t0 = time.perf_counter()
    result1 = polars_native(df.clone())
    print(f"Total: {(time.perf_counter() - t0) * 1000:.1f} ms")

    # Profile Threads
    print("\n--- Threads (4 workers) ---")
    t0 = time.perf_counter()
    result2 = threads_approach(df.clone(), workers=4)
    print(f"Total: {(time.perf_counter() - t0) * 1000:.1f} ms")
