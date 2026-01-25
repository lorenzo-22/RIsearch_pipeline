"""Typer CLI for the siRNA off-target discovery pipeline."""

import typer


from RIsearch_pipeline.commands import accessibility, off_targets

app = typer.Typer(
    name="risearch-pipeline",
    help="siRNA off-target discovery pipeline — analyze RIsearch2 predictions.",
    add_completion=False,
)
app.command(name="accessibility")(accessibility.run)
app.command(name="off-targets")(off_targets.run)


@app.callback()
def main(ctx: typer.Context) -> None:
    """
    siRNA off-target discovery pipeline.

    Use 'accessibility' to pre-compute profiles or 'off-targets' to analyze predictions.
    """
    pass


def entrypoint() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    entrypoint()
