from loguru import logger
from pathlib import Path
import typer
from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService


def run(
    genome: Path = typer.Option(..., "--fasta", "-f", help="Path to genome FASTA file"),
    output_dir: Path = typer.Option(
        ..., "--output", "-o", help="Directory to save accessibility profiles"
    ),
    window_size: int = typer.Option(80, "--window", "-W", help="Window size (W)"),
    max_span: int = typer.Option(40, "--span", "-L", help="Max base pair span (L)"),
    unpaired_prob: int = typer.Option(
        30, "--unpaired", "-u", help="Unpaired probability length (u)"
    ),
):
    """
    Pre-compute accessibility profiles for a whole genome.

    Uses ViennaRNA to calculate opening energies/unpaired probabilities and stores them
    in an efficient binary format for fast random access.
    """
    service = GenomeAccessibilityService(output_dir)
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
