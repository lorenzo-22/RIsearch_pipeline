"""CLI commands for running RIsearch index and search operations."""

from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from RIsearch_pipeline.services.risearch_service import RIsearchService


def index(
    target: Path = typer.Argument(..., help="Target FASTA file to index."),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output index path. Defaults to <target>.idx next to the target file.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose logging."
    ),
) -> None:
    """Build a RIsearch index from a target FASTA file."""
    if verbose:
        import sys
        logger.remove()
        logger.add(
            sys.stderr,
            level="DEBUG",
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>",
        )

    svc = RIsearchService()
    index_path = svc.index_target(target, output)
    typer.echo(f"Index written to: {index_path}")


def search(
    query: Path = typer.Argument(..., help="Query siRNA FASTA file."),
    index: Path = typer.Argument(..., help="Pre-built RIsearch index (.idx)."),
    target: Optional[Path] = typer.Option(
        None,
        "--target",
        "-t",
        help="Target FASTA used to build the index (required if index was built externally).",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output TSV file. Defaults to stdout.",
    ),
    seed_length: int = typer.Option(6, "--seed", "-s", help="Seed length."),
    max_extension: int = typer.Option(
        20, "--max-extension", "-e", help="Max extension length on each side."
    ),
    energy_threshold: float = typer.Option(
        -10.0,
        "--energy",
        "-E",
        help="Energy threshold in kcal/mol (only hits below this value are kept).",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose logging."
    ),
) -> None:
    """Run a RIsearch search and output hits as TSV."""
    if verbose:
        import sys
        logger.remove()
        logger.add(
            sys.stderr,
            level="DEBUG",
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>",
        )

    svc = RIsearchService()
    df = svc.run_search(
        query_path=query,
        index_path=index,
        target_fasta=target,
        seed_length=seed_length,
        max_extension=max_extension,
        energy_threshold=energy_threshold,
    )

    if output is not None:
        df.write_csv(output, separator="\t")
        typer.echo(f"Results written to: {output} ({df.height} hits)")
    else:
        typer.echo(df.write_csv(separator="\t"), nl=False)
