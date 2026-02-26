from loguru import logger
from pathlib import Path
from typing import Optional
import typer
from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService


def run(
    genome: Path = typer.Option(..., "--fasta", "-f", help="Path to genome FASTA file"),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output directory (full-genome mode) or .parquet file (binding-site mode)",
    ),
    risearch_dir: Optional[Path] = typer.Option(
        None,
        "--risearch-dir",
        "-r",
        help="Directory of per-siRNA RIsearch output files. "
        "When provided, computes accessibility only for binding site regions "
        "and saves to Parquet instead of full-genome profiles.",
    ),
    window_size: int = typer.Option(80, "--window", "-W", help="Window size (W)"),
    max_span: int = typer.Option(40, "--span", "-L", help="Max base pair span (L)"),
    unpaired_prob: int = typer.Option(
        30, "--unpaired", "-u", help="Unpaired probability length (u)"
    ),
    workers: int = typer.Option(
        1,
        "--workers",
        "-j",
        help="Number of parallel workers for island folding (binding-site mode only).",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Show detailed progress information."
    ),
):
    """
    Pre-compute accessibility profiles for a genome.

    Two modes:
    - Full-genome: computes profiles for every position (no --risearch-dir).
    - Binding-site: computes only for predicted binding sites (--risearch-dir).
      This is ~25x faster and produces a compact Parquet file.
    """
    from rich.console import Console
    from rich.progress import (
        Progress,
        SpinnerColumn,
        TextColumn,
        BarColumn,
        TaskProgressColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )

    console = Console(stderr=True)

    if risearch_dir is not None:
        # --- Binding-site mode ---
        if not risearch_dir.is_dir():
            console.print(f"[red]Error: {risearch_dir} is not a directory[/red]")
            raise typer.Exit(code=1)

        out_path = (
            output if output.suffix == ".parquet" else output / "accessibility.parquet"
        )
        service = GenomeAccessibilityService(out_path.parent)

        console.print(f"\n[bold]Binding-site accessibility computation[/bold]")
        console.print(f"  Genome:       {genome}")
        console.print(f"  RIsearch dir: {risearch_dir}")
        console.print(f"  Output:       {out_path}")
        console.print(
            f"  Parameters:   W={window_size}, L={max_span}, u={unpaired_prob}"
        )
        console.print(f"  Workers:      {workers}\n")

        try:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
                console=console,
                transient=False,
            ) as progress:
                result_path = service.compute_binding_site_accessibility(
                    genome_path=genome,
                    risearch_dir=risearch_dir,
                    output_path=out_path,
                    window_size=window_size,
                    max_span=max_span,
                    unpaired_prob=unpaired_prob,
                    workers=workers,
                    progress=progress,
                    verbose=verbose,
                )

            console.print(
                f"\n[green bold]✓[/green bold] Saved to [cyan]{result_path}[/cyan]"
            )

        except Exception as e:
            logger.exception("Failed to compute binding-site accessibility")
            console.print(f"[red]Error: {e}[/red]")
            raise typer.Exit(code=1)
    else:
        # --- Full-genome mode (existing behavior) ---
        service = GenomeAccessibilityService(output)
        try:
            results = service.compute_genome_accessibility(
                genome,
                window_size=window_size,
                max_span=max_span,
                unpaired_prob=unpaired_prob,
            )
            typer.echo(f"Successfully processed {len(results)} sequences.")
            for chrom, path in results.items():
                typer.echo(f"  {chrom}: {path}")
        except Exception as e:
            logger.exception("Failed to compute accessibility")
            typer.echo(f"Error: {e}", err=True)
            raise typer.Exit(code=1)
