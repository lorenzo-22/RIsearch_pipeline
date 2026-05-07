"""Typer CLI for the siRNA off-target discovery pipeline."""

from pathlib import Path
from typing import Optional

import typer

from RIsearch_pipeline._logging import setup_logging
from RIsearch_pipeline.commands import accessibility, off_targets, risearch


app = typer.Typer(
    name="risearch-pipeline",
    help="siRNA off-target discovery pipeline — analyze RIsearch2 predictions.",
    add_completion=False,
)
app.command(name="accessibility")(accessibility.run)
app.command(name="off-targets")(off_targets.run)
app.command(name="index")(risearch.index)
app.command(name="search")(risearch.search)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to YAML config file. Runs the command specified in the config.",
        exists=True,
        readable=True,
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose logging (DEBUG level).",
    ),
) -> None:
    """
    siRNA off-target discovery pipeline.

    Use 'accessibility' to pre-compute profiles or 'off-targets' to analyze predictions.
    Alternatively, use --config to run from a YAML configuration file.
    """
    setup_logging(verbose)

    if config is not None:
        from RIsearch_pipeline.config import load_config, config_to_kwargs

        try:
            cfg = load_config(config)
            kwargs = config_to_kwargs(cfg, cfg.command)

            final_verbose = verbose or cfg.verbose
            setup_logging(final_verbose)
            kwargs["verbose"] = final_verbose

            if cfg.command == "off-targets":
                off_targets.run(**kwargs)
            elif cfg.command == "accessibility":
                accessibility.run(**kwargs)
        except Exception as e:
            typer.echo(f"Error loading config: {e}", err=True)
            raise typer.Exit(code=1)

        raise typer.Exit()

    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


def entrypoint() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    entrypoint()
