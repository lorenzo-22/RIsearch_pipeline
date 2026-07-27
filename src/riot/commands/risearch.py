"""CLI commands for running RIsearch index and search operations."""

from pathlib import Path
from typing import Annotated, Optional

import polars as pl
import typer

from riot._logging import setup_logging
from riot.core import risearch as core


def index(
    target: Annotated[Path, typer.Argument(help="Target FASTA file to index.")],
    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Output index path. Defaults to <target>.idx next to the target file.",
        ),
    ] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging.")
    ] = False,
) -> Path:
    """Build a RIsearch index from a target FASTA file."""
    setup_logging(verbose)
    index_path = core.build_index(target, output)
    typer.echo(f"Index written to: {index_path}")
    return index_path


def search(
    query: Annotated[Path, typer.Argument(help="Query siRNA FASTA file.")],
    index: Annotated[Path, typer.Argument(help="Pre-built RIsearch index (.idx).")],
    target: Annotated[
        Optional[Path],
        typer.Option(
            "--target",
            "-t",
            help="Target FASTA used to build the index (required if index was built externally).",
        ),
    ] = None,
    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Output TSV file. Defaults to stdout.",
        ),
    ] = None,
    seed_length: Annotated[int, typer.Option("--seed", "-s", help="Seed length.")] = 6,
    max_extension: Annotated[
        int,
        typer.Option(
            "--max-extension", "-e", help="Max extension length on each side."
        ),
    ] = 20,
    energy_threshold: Annotated[
        float,
        typer.Option(
            "--energy",
            "-E",
            help="Energy threshold in kcal/mol (only hits below this value are kept).",
        ),
    ] = -10.0,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging.")
    ] = False,
) -> pl.DataFrame:
    """Run a RIsearch search and output hits as TSV."""
    setup_logging(verbose)
    df = core.run_search(
        query=query,
        index=index,
        target=target,
        seed_length=seed_length,
        max_extension=max_extension,
        energy_threshold=energy_threshold,
    )

    output = Path(output) if isinstance(output, str) else output
    if output is not None:
        df.write_csv(output, separator="\t")
        typer.echo(f"Results written to: {output} ({df.height} hits)")
    else:
        typer.echo(df.write_csv(separator="\t"), nl=False)

    return df
