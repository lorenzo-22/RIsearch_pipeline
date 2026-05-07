"""Typer CLI for the siRNA off-target discovery pipeline."""

import sys
from pathlib import Path
from typing import Optional

import typer
from loguru import logger


from RIsearch_pipeline.commands import accessibility, off_targets, risearch


def setup_logging(verbose: bool) -> None:
    """Configure loguru based on verbosity."""
    logger.remove()  # Remove default handler
    if verbose:
        logger.add(
            sys.stderr,
            level="DEBUG",
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan> - <level>{message}</level>",
        )
    else:
        logger.add(sys.stderr, level="WARNING", format="<level>{message}</level>")


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

            # Pass verbose from CLI or config
            # Pass verbose from CLI or config (config overrides default CLI if not explicitly set by user, but explicit CLI beats config? usually global CLI flag beats config)
            # Actually, let's say CLI flag overrides config if set to True.
            # But if config has verbose=True, we should honor it even if CLI is False (default).
            final_verbose = verbose or cfg.verbose
            if final_verbose:
                setup_logging(True)
                kwargs["verbose"] = True
            else:
                kwargs["verbose"] = False

            if cfg.command == "off-targets":
                off_targets.run(**kwargs)
            elif cfg.command == "accessibility":
                accessibility.run(**kwargs)
        except Exception as e:
            typer.echo(f"Error loading config: {e}", err=True)
            raise typer.Exit(code=1)

        raise typer.Exit()

    # If no config and no subcommand, show help
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


def entrypoint() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    entrypoint()
