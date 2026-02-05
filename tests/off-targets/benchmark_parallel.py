#!/usr/bin/env python3
"""
Benchmark parallel processing strategies for multi-siRNA off-target analysis.

Usage:
    uv run python tests/off-targets/benchmark_parallel.py
    uv run python tests/off-targets/benchmark_parallel.py --config tests/off-targets/benchmark_config.toml
"""

import time
import tomllib
from pathlib import Path
from typing import Optional
import polars as pl
from rich.console import Console
from rich.table import Table
from rich.progress import track
import typer

console = Console()


def load_config(config_path: Path) -> dict:
    """Load benchmark configuration from TOML file."""
    with open(config_path, "rb") as f:
        return tomllib.load(f)


def benchmark_polars_native_detailed(df: pl.DataFrame) -> dict[str, float]:
    """Benchmark native Polars with step-by-step timing."""
    from RIsearch_pipeline.services.probability import RT

    timings = {}

    # Step 1: Add dG_total column
    t0 = time.perf_counter()
    df = df.with_columns(pl.col("energy").alias("dG_total"))
    timings["1_add_dG_total"] = time.perf_counter() - t0

    # Step 2: Calculate Boltzmann weights
    t0 = time.perf_counter()
    df = df.with_columns(((-pl.col("dG_total") / RT).exp()).alias("boltzmann_weight"))
    timings["2_boltzmann_weight"] = time.perf_counter() - t0

    # Step 3: Group by siRNA and compute Z per siRNA
    t0 = time.perf_counter()
    z_per_sirna = df.group_by("sirna_id").agg(
        pl.col("boltzmann_weight").sum().alias("Z_sirna")
    )
    timings["3_group_by_Z"] = time.perf_counter() - t0

    # Step 4: Join Z back to main DataFrame
    t0 = time.perf_counter()
    df = df.join(z_per_sirna, on="sirna_id", how="left")
    timings["4_join_Z"] = time.perf_counter() - t0

    # Step 5: Calculate P_off_target
    t0 = time.perf_counter()
    df = df.with_columns(
        (pl.col("boltzmann_weight") / pl.col("Z_sirna")).alias("P_off_target")
    )
    timings["5_calc_P_off"] = time.perf_counter() - t0

    timings["total"] = sum(v for k, v in timings.items() if k != "total")
    return timings


def benchmark_threads_detailed(df: pl.DataFrame, workers: int) -> dict[str, float]:
    """Benchmark ThreadPoolExecutor with step-by-step timing."""
    from concurrent.futures import ThreadPoolExecutor
    from RIsearch_pipeline.services.probability import RT

    timings = {}

    def process_sirna_group(group_df: pl.DataFrame) -> pl.DataFrame:
        """Process a single siRNA's predictions."""
        group_df = group_df.with_columns(pl.col("energy").alias("dG_total"))
        group_df = group_df.with_columns(
            ((-pl.col("dG_total") / RT).exp()).alias("boltzmann_weight")
        )
        z = group_df["boltzmann_weight"].sum()
        return group_df.with_columns(
            (pl.col("boltzmann_weight") / z).alias("P_off_target")
        )

    # Step 1: Partition by siRNA
    t0 = time.perf_counter()
    groups = df.partition_by("sirna_id", maintain_order=True)
    timings["1_partition_by"] = time.perf_counter() - t0

    # Step 2: Create ThreadPoolExecutor
    t0 = time.perf_counter()
    executor = ThreadPoolExecutor(max_workers=workers)
    timings["2_create_executor"] = time.perf_counter() - t0

    # Step 3: Submit and execute tasks
    t0 = time.perf_counter()
    results = list(executor.map(process_sirna_group, groups))
    timings["3_execute_tasks"] = time.perf_counter() - t0

    # Step 4: Shutdown executor
    t0 = time.perf_counter()
    executor.shutdown(wait=True)
    timings["4_shutdown"] = time.perf_counter() - t0

    # Step 5: Concatenate results
    t0 = time.perf_counter()
    result_df = pl.concat(results)
    timings["5_concat_results"] = time.perf_counter() - t0

    timings["total"] = sum(v for k, v in timings.items() if k != "total")
    return timings


def benchmark_polars_native(df: pl.DataFrame, iterations: int) -> list[float]:
    """Benchmark native Polars group_by parallelism."""
    from RIsearch_pipeline.services.probability import ProbabilityService

    prob_service = ProbabilityService(None)
    times = []

    for _ in range(iterations):
        start = time.perf_counter()
        prob_service.calculate_probabilities_per_sirna(df.clone())
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return times


def benchmark_threads(df: pl.DataFrame, workers: int, iterations: int) -> list[float]:
    """Benchmark ThreadPoolExecutor approach."""
    from concurrent.futures import ThreadPoolExecutor
    from RIsearch_pipeline.services.probability import RT

    def process_sirna_group(group_df: pl.DataFrame) -> pl.DataFrame:
        """Process a single siRNA's predictions."""
        group_df = group_df.with_columns(pl.col("energy").alias("dG_total"))
        group_df = group_df.with_columns(
            ((-pl.col("dG_total") / RT).exp()).alias("boltzmann_weight")
        )
        z = group_df["boltzmann_weight"].sum()
        return group_df.with_columns(
            (pl.col("boltzmann_weight") / z).alias("P_off_target")
        )

    times = []
    for _ in range(iterations):
        groups = df.partition_by("sirna_id", maintain_order=True)
        start = time.perf_counter()
        with ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(process_sirna_group, groups))
        pl.concat(results)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return times


def benchmark_processes(df: pl.DataFrame, workers: int, iterations: int) -> list[float]:
    """Benchmark ProcessPoolExecutor approach."""
    from concurrent.futures import ProcessPoolExecutor
    import pickle

    times = []

    for _ in range(iterations):
        groups = df.partition_by("sirna_id", maintain_order=True)

        start = time.perf_counter()
        with ProcessPoolExecutor(max_workers=workers) as executor:
            serialized = [pickle.dumps(g) for g in groups]
            [pickle.loads(s) for s in serialized]
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return times


def run_benchmark(
    config_path: Path = typer.Option(
        Path("tests/off-targets/benchmark_config.toml"),
        "--config",
        "-c",
        help="Path to benchmark configuration file",
    ),
    quick: bool = typer.Option(
        False,
        "--quick",
        "-q",
        help="Run quick benchmark with fewer iterations",
    ),
):
    """Run parallel processing benchmarks."""
    console.print("[bold cyan]Multi-siRNA Parallelization Benchmark[/bold cyan]\n")

    # Load config
    config = load_config(config_path)
    iterations = 1 if quick else config["benchmark"]["iterations"]

    # Load test data
    from RIsearch_pipeline.services.risearch_parser import RIsearchParser

    parser = RIsearchParser()

    # Check if we have multi-siRNA test data
    risearch_path = Path(config["data"]["risearch_file"])
    if not risearch_path.exists():
        # Generate test data from existing multi-siRNA results
        sirna_fasta = Path(config["data"]["sirna_fasta"])
        target_fasta = Path(config["data"]["target_fasta"])

        if sirna_fasta.exists() and target_fasta.exists():
            console.print("[yellow]Generating test RIsearch output...[/yellow]")
            from RIsearch_pipeline.services.risearch_service import RIsearchService

            service = RIsearchService()
            service.run_search(
                query_path=sirna_fasta,
                index_path=service.index_target(target_fasta),
                output_path=risearch_path,
            )
        else:
            console.print("[red]Test data not found. Using single siRNA test.[/red]")
            risearch_path = Path("tests/off-targets/single/input/risearch_siRNAID.out")

    df = parser.load(risearch_path)
    n_sirnas = df["sirna_id"].n_unique()
    n_predictions = len(df)

    console.print(
        f"[green]✓[/green] Loaded {n_predictions:,} predictions for {n_sirnas} siRNAs\n"
    )

    # Run benchmarks for each scale
    sirna_counts = config["data"]["sirna_counts"]
    worker_counts = config["parallel_modes"]["worker_counts"]

    for target_count in sirna_counts:
        console.print(
            f"\n[bold yellow]--- Benchmarking Scale: {target_count} siRNAs ---[/bold yellow]"
        )

        # Scale data if needed
        current_count = df["sirna_id"].n_unique()
        if target_count > current_count:
            # Replicate data to match target count
            multiplier = (target_count // current_count) + 1

            # Create duplicates with new IDs
            dfs = []
            for i in range(multiplier):
                suffix = f"_copy{i}"
                df_copy = df.clone()
                df_copy = df_copy.with_columns(
                    (pl.col("sirna_id") + suffix).alias("sirna_id")
                )
                dfs.append(df_copy)

            scaled_df = pl.concat(dfs)

            # Trim to exact count
            unique_ids = scaled_df["sirna_id"].unique().head(target_count)
            scaled_df = scaled_df.filter(pl.col("sirna_id").is_in(unique_ids))
        else:
            # Downsample if needed
            unique_ids = df["sirna_id"].unique().head(target_count)
            scaled_df = df.filter(pl.col("sirna_id").is_in(unique_ids))

        n_predictions = len(scaled_df)
        console.print(
            f"Dataset: {n_predictions:,} predictions for {target_count} siRNAs"
        )

        # Results table for this scale
        results = Table(title=f"Results ({target_count} siRNAs)")
        results.add_column("Mode", style="cyan")
        results.add_column("Workers", justify="right")
        results.add_column("Mean Time (s)", justify="right")
        results.add_column("Std Dev", justify="right")
        results.add_column("Speedup", justify="right", style="green")

        baseline_time = None

        # Benchmark Polars native
        console.print("[bold]Testing: Polars Native[/bold]")
        times = benchmark_polars_native(scaled_df, iterations)
        mean_time = sum(times) / len(times)
        std_dev = (sum((t - mean_time) ** 2 for t in times) / len(times)) ** 0.5
        baseline_time = mean_time
        results.add_row(
            "Polars Native", "-", f"{mean_time:.4f}", f"{std_dev:.4f}", "1.00x"
        )

        # Benchmark Threads
        for workers in worker_counts:
            console.print(f"Testing: Threads (workers={workers})")
            times = benchmark_threads(scaled_df, workers, iterations)
            mean_time = sum(times) / len(times)
            std_dev = (sum((t - mean_time) ** 2 for t in times) / len(times)) ** 0.5
            speedup = baseline_time / mean_time if mean_time > 0 else 0
            results.add_row(
                "Threads",
                str(workers),
                f"{mean_time:.4f}",
                f"{std_dev:.4f}",
                f"{speedup:.2f}x",
            )

        # Print table for this scale
        console.print(results)


if __name__ == "__main__":
    typer.run(run_benchmark)
