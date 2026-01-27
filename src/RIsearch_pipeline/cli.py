"""Typer CLI for the siRNA off-target discovery pipeline."""

from pathlib import Path
from typing import Optional

import typer


from RIsearch_pipeline.commands import accessibility, off_targets

app = typer.Typer(
    name="risearch-pipeline",
    help="siRNA off-target discovery pipeline — analyze RIsearch2 predictions.",
    add_completion=False,
)
app.command(name="accessibility")(accessibility.run)
app.command(name="off-targets")(off_targets.run)


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
) -> None:
    """
    siRNA off-target discovery pipeline.

    Use 'accessibility' to pre-compute profiles or 'off-targets' to analyze predictions.
    Alternatively, use --config to run from a YAML configuration file.
    """
    if config is not None:
        from RIsearch_pipeline.config import load_config, config_to_kwargs

        try:
            cfg = load_config(config)
            kwargs = config_to_kwargs(cfg, cfg.command)

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
