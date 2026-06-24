"""CLI commands for running RIsearch index and search operations."""

from pathlib import Path
from typing import Annotated, Optional

import polars as pl
import typer

from riot._logging import setup_logging
from riot.services.risearch_service import RIsearchService


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
    # Accept str paths (Python API) as well as Path objects.
    target = Path(target) if isinstance(target, str) else target
    output = Path(output) if isinstance(output, str) else output
    svc = RIsearchService()
    index_path = svc.index_target(target, output)
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
    seed_length: Annotated[
        int, typer.Option("--seed", "-s", help="Seed length.")
    ] = 6,
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
    # Accept str paths (Python API) as well as Path objects.
    query = Path(query) if isinstance(query, str) else query
    index = Path(index) if isinstance(index, str) else index
    target = Path(target) if isinstance(target, str) else target
    output = Path(output) if isinstance(output, str) else output
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

    return df
