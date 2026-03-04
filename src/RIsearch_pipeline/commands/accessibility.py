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
        help="Directory of per-siRNA RIsearch output files.",
    ),
    risearch_file: Optional[Path] = typer.Option(
        None,
        "--risearch-file",
        "-R",
        help="Single merged RIsearch TSV file (alternative to --risearch-dir).",
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
        help="Number of parallel workers for chromosome/island folding.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Show detailed progress information."
    ),
):
    """
    Pre-compute accessibility profiles for a genome.

    Two modes:
    - Full-genome: computes profiles for every position.
    - Binding-site: computes only for predicted binding sites.
      Accepts --risearch-dir (per-siRNA files) or --risearch-file (single merged TSV).
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

    if risearch_dir is not None or risearch_file is not None:
        # --- Binding-site mode ---
        if risearch_dir and not risearch_dir.is_dir():
            console.print(f"[red]Error: {risearch_dir} is not a directory[/red]")
            raise typer.Exit(code=1)
        if risearch_file and not risearch_file.is_file():
            console.print(f"[red]Error: {risearch_file} is not a file[/red]")
            raise typer.Exit(code=1)

        out_path = (
            output if output.suffix == ".parquet" else output / "accessibility.parquet"
        )
        service = GenomeAccessibilityService(out_path.parent)

        source = risearch_dir or risearch_file
        console.print("\n[bold]Binding-site accessibility computation[/bold]")
        console.print(f"  Genome:       {genome}")
        console.print(f"  Input:        {source}")
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
            ) as prog:
                result_path = service.compute_binding_site_accessibility(
                    genome_path=genome,
                    output_path=out_path,
                    risearch_dir=risearch_dir,
                    risearch_file=risearch_file,
                    window_size=window_size,
                    max_span=max_span,
                    unpaired_prob=unpaired_prob,
                    workers=workers,
                    progress=prog,
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
        # --- Full-genome mode ---
        output.mkdir(parents=True, exist_ok=True)
        service = GenomeAccessibilityService(output)

        console.print("\n[bold]Full-genome accessibility computation[/bold]")
        console.print(f"  Genome:       {genome}")
        console.print(f"  Output dir:   {output}")
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
            ) as prog:
                results = service.compute_genome_accessibility(
                    genome,
                    window_size=window_size,
                    max_span=max_span,
                    unpaired_prob=unpaired_prob,
                    workers=workers,
                    progress=prog,
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
