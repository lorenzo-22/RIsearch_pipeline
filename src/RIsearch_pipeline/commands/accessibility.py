from loguru import logger
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.progress import (
    Progress, SpinnerColumn, TextColumn, BarColumn,
    TaskProgressColumn, TimeElapsedColumn, TimeRemainingColumn,
)

from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService


def _make_progress(console):
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=False,
    )


def run(
    genome: Optional[Path] = typer.Option(None, "--fasta", "-f", help="Path to genome or transcriptome FASTA file"),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output directory — one {chrom}.accessibility.parquet per chromosome is written here",
    ),
    window_size: int = typer.Option(80, "--window", "-W", help="Window size (W)"),
    max_span: int = typer.Option(40, "--span", "-L", help="Max base pair span (L)"),
    unpaired_prob: int = typer.Option(
        30, "--unpaired", "-u", help="Unpaired probability length (u)"
    ),
    temperature: float = typer.Option(
        37.0, "--temperature", "-T", help="Folding temperature in °C (default 37.0)."
    ),
    workers: int = typer.Option(
        1,
        "--workers",
        "-j",
        help="Number of parallel workers (one per chromosome).",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Show detailed progress information."
    ),
):
    """
    Pre-compute per-chromosome accessibility profiles.

    Folds every chromosome in the FASTA with ViennaRNA RNAplfold and writes
    one {chrom}.accessibility.parquet per chromosome to the output directory.
    Each Parquet has columns [position (int32), strand (utf8), u1..u{u} (float32)].

    Profiles can later be consumed by the off-targets command via
    --accessibility-dir.
    """
    console = Console(stderr=True)

    if genome is None:
        console.print("[red]Error: --fasta is required[/red]")
        raise typer.Exit(code=1)

    output.mkdir(parents=True, exist_ok=True)
    service = GenomeAccessibilityService(output)

    console.print("\n[bold]Full-genome accessibility computation[/bold]")
    console.print(f"  Genome:       {genome}")
    console.print(f"  Output dir:   {output}")
    console.print(
        f"  Parameters:   W={window_size}, L={max_span}, u={unpaired_prob}, T={temperature}°C"
    )
    console.print(f"  Workers:      {workers}\n")

    try:
        with _make_progress(console) as prog:
            results = service.compute_genome_accessibility(
                genome,
                window_size=window_size,
                max_span=max_span,
                unpaired_prob=unpaired_prob,
                workers=workers,
                progress=prog,
                temperature=temperature,
            )

        console.print(
            f"\n[green bold]✓[/green bold] Processed {len(results)} chromosome(s)"
        )
        for chrom, path in results.items():
            console.print(f"  {chrom}: [cyan]{path}[/cyan]")

    except Exception as e:
        logger.exception("Failed to compute accessibility")
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(code=1)
