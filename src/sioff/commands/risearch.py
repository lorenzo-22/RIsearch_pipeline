"""CLI commands for running RIsearch index and search operations."""

from pathlib import Path
from typing import Annotated, Optional

import polars as pl
import typer

from sioff._logging import setup_logging
from sioff.core import risearch as core
from sioff.services.risearch_service import RIsearchError


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
    seed: Annotated[
        str,
        typer.Option(
            "--seed",
            "-s",
            help=(
                "Seed spec, RIsearch2 syntax: 'l' (length only), 'n:m' or 'n:m/l' "
                "(seed must lie in 1-based query positions n..m and be l nt). "
                "E.g. '7', '2:8/7'."
            ),
        ),
    ] = "6",
    no_gu_seed: Annotated[
        bool,
        typer.Option(
            "--no-gu-seed",
            help=(
                "Forbid G:U wobble pairs inside the seed (RIsearch2 --noGUseed). "
                "Affects seed location only, not the energy model, so this is a "
                "candidate-generation knob. Removes ~85% of hits on fly 3'UTRs."
            ),
        ),
    ] = False,
    matrix: Annotated[
        str,
        typer.Option(
            "--matrix",
            "-z",
            help=(
                "Nearest-neighbour energy parameter set: 't04' (Turner 2004, "
                "RNA-RNA), 't99' (Turner 1999, RNA-RNA), 'slh04' "
                "(SantaLucia-Hicks 2004, DNA-DNA), 's95-rna-dna' or "
                "'s95-dna-rna' (Sugimoto 1995, RNA/DNA hybrids)."
            ),
        ),
    ] = "t04",
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
    # Seed geometry and matrix are validated before the search starts, so a bad
    # spec is a one-line CLI error rather than a traceback out of the bindings.
    try:
        seed_start, seed_end, seed_length = core.parse_seed_spec(seed)
        df = core.run_search(
            query=query,
            index=index,
            target=target,
            seed_length=seed_length,
            seed_start=seed_start,
            seed_end=seed_end,
            seed_wobble=not no_gu_seed,
            matrix=matrix,
            max_extension=max_extension,
            energy_threshold=energy_threshold,
        )
    except RIsearchError as e:
        typer.secho(f"Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    output = Path(output) if isinstance(output, str) else output
    if output is not None:
        df.write_csv(output, separator="\t")
        typer.echo(f"Results written to: {output} ({df.height} hits)")
    else:
        typer.echo(df.write_csv(separator="\t"), nl=False)

    return df
