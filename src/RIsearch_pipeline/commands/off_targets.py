from pathlib import Path
import typer
from RIsearch_pipeline.services.risearch_parser import RIsearchParser


def run(
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
    accessibility_dir: Path = typer.Option(
        None,
        "-a",
        "--accessibility-dir",
        help="Directory containing pre-computed accessibility profiles to calculate P(OT).",
        exists=True,
        file_okay=False,
    ),
    output_file: Path = typer.Option(
        None,
        "-o",
        "--output",
        help="Path to save the analysis results (TSV).",
    ),
    genome_file: Path = typer.Option(
        None,
        "-f",
        "--fasta",
        help="Path to genome FASTA file (computes accessibility on-the-fly).",
        exists=True,
        readable=True,
    ),
    window_size: int = typer.Option(80, "--window", "-W", help="Window size (W)"),
    max_span: int = typer.Option(40, "--span", "-L", help="Max base pair span (L)"),
    unpaired_prob: int = typer.Option(
        30, "--unpaired", "-u", help="Unpaired probability length (u)"
    ),
) -> None:
    """
    Analyze siRNA off-target predictions.

    Integrates RIsearch2 predictions with transcriptome data and optionally
    calculates off-target probabilities.

    Accessibility can be provided via pre-computed directory (-a) OR
    computed on-the-fly from a genome FASTA (-f).
    """
    risearch_parser = RIsearchParser()

    try:
        # Load RIsearch2 predictions
        if not isinstance(risearch_file, Path):
            risearch_file = Path(risearch_file)

        df = risearch_parser.load(risearch_file)
        summary_ris = risearch_parser.summary(df)

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
            df = intersector.intersect(df, df_trans)

            typer.echo(f"✓ Found {df.height} intersecting off-target candidates")

        # Accessibility and Probability Calculation
        from RIsearch_pipeline.services.probability import ProbabilityService
        from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService
        import tempfile

        # Keep reference to temp dir object so it persists until function exit
        temp_dir_obj = None

        if accessibility_dir:
            typer.echo(
                f"  Calculating probabilities using accessibility profiles from {accessibility_dir}..."
            )
            acc_service = GenomeAccessibilityService(accessibility_dir)
            prob_service = ProbabilityService(acc_service)

        elif genome_file:
            typer.echo(f"  Computing accessibility on-the-fly from {genome_file}...")
            typer.echo("  (This may take some time for large genomes)")

            # Create temp dir
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="risearch_accessibility_")
            temp_dir = Path(temp_dir_obj.name)

            acc_service = GenomeAccessibilityService(temp_dir)
            acc_service.compute_genome_accessibility(
                genome_file,
                window_size=window_size,
                max_span=max_span,
                unpaired_prob=unpaired_prob,
            )
            prob_service = ProbabilityService(acc_service)

        else:
            typer.echo(
                "  Warning: No accessibility data provided (-a or -f). P(OT) will be based only on hybridization energy."
            )
            prob_service = ProbabilityService(None)

        # Calculate P(OT)
        df = prob_service.calculate_probabilities(df)

        # Display Results
        if df.height > 0:
            if "P_off_target" in df.columns:
                # Sort by Probability desc
                df = df.sort("P_off_target", descending=True)
                typer.echo("\nTop 10 candidates (Sorted by P(OT)):")
            else:
                typer.echo("\nFirst 5 rows (Predictions):")

            # Select useful columns to show
            potential_cols = [
                "sirna_id",
                "chrom",
                "start",
                "end",
                "gene_id",
                "transcript_id",
                "energy",
                "opening_energy",
                "dG_total",
                "P_off_target",
                "exp_value",
            ]
            # specific intersection selection
            display_cols = [c for c in potential_cols if c in df.columns]

            typer.echo(df.select(display_cols).head(10))

        # Save to Output File
        if output_file:
            df.write_csv(output_file, separator="\t")
            typer.echo(f"\n✓ Results saved to {output_file}")

    except FileNotFoundError as e:
        typer.echo(f"✗ Error: {e}", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.echo(f"✗ Failed to parse file: {e}", err=True)
        raise typer.Exit(code=1)
