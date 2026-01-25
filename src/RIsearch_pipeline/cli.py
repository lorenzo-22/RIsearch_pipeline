"""Typer CLI for the siRNA off-target discovery pipeline."""

from pathlib import Path

import typer

from RIsearch_pipeline.services.risearch_parser import RIsearchParser

app = typer.Typer(
    name="risearch-pipeline",
    help="siRNA off-target discovery pipeline — analyze RIsearch2 predictions.",
    add_completion=False,
)


@app.callback(invoke_without_command=True)
def analyze(
    ctx: typer.Context,
    risearch_file: Path = typer.Option(
        ...,
        "-r",
        "--risearch-file",
        help="Path to RIsearch2 output TSV file.",
        exists=True,
        readable=True,
    ),
    gtf_file: Path = typer.Option(
        None,
        "-t",
        "--transcriptome",
        help="Path to transcriptome annotation file (.gtf).",
        exists=True,
        readable=True,
    ),
    feature_type: str = typer.Option(
        "exon",
        "--feature",
        help="Feature type to select from GTF (default: exon).",
    ),
    expression_metric: str = typer.Option(
        "RPKM",
        "--expression-metric",
        help="Attribute to use for expression score (default: RPKM).",
    ),
) -> None:
    """Analyze siRNA off-target predictions."""
    if ctx.invoked_subcommand is not None:
        return

    risearch_parser = RIsearchParser()

    try:
        # Load RIsearch2 predictions
        if not isinstance(risearch_file, Path):
            risearch_file = Path(risearch_file)

        df_risearch = risearch_parser.load(risearch_file)
        summary_ris = risearch_parser.summary(df_risearch)

        typer.echo(
            f"✓ Loaded {summary_ris['row_count']} predictions from {risearch_file.name}"
        )
        typer.echo(
            f"  Energy range: {summary_ris['energy_min']:.2f} to {summary_ris['energy_max']:.2f} kcal/mol"
        )

        # Load Transcriptome if provided
        if gtf_file:
            from RIsearch_pipeline.services.transcriptome_parser import (
                TranscriptomeParser,
            )
            from RIsearch_pipeline.services.intersection_service import (
                IntersectionService,
            )

            if not isinstance(gtf_file, Path):
                gtf_file = Path(gtf_file)

            trans_parser = TranscriptomeParser()
            df_trans = trans_parser.load_gtf(
                gtf_file, feature=feature_type, score_col=expression_metric
            )
            summary_trans = trans_parser.summary(df_trans)

            typer.echo(
                f"✓ Loaded {summary_trans['row_count']} {feature_type}s from {gtf_file.name}"
            )

            # Perform intersection
            typer.echo("  Intersecting predictions with transcriptome...")
            intersector = IntersectionService()
            df_intersect = intersector.intersect(df_risearch, df_trans)

            typer.echo(
                f"✓ Found {df_intersect.height} intersecting off-target candidates"
            )

            if df_intersect.height > 0:
                typer.echo("\nFirst 5 intersecting candidates:")
                # Select useful columns to show
                display_cols = [
                    "sirna_id",
                    "chrom",
                    "start",
                    "end",
                    "gene_id",
                    "transcript_id",
                    "energy",
                    "exp_value",
                ]
                typer.echo(df_intersect.select(display_cols).head(5))
        else:
            typer.echo("\nFirst 5 rows (Predictions):")
            typer.echo(df_risearch.head(5))

    except FileNotFoundError as e:
        typer.echo(f"✗ Error: {e}", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"✗ Failed to parse file: {e}", err=True)
        raise typer.Exit(code=1)


def entrypoint() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    entrypoint()
