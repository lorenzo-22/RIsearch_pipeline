from pathlib import Path
import typer
import polars as pl
from typing import Optional
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
        help="Path to transcriptome annotation file (.gtf) or .bed file.",
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
    on_target_file: Path = typer.Option(
        None,
        "-on",
        "--on-target",
        help="Path to On-Target sequence FASTA (for Partition Function).",
        exists=True,
        readable=True,
    ),
    on_target_risearch_file: Optional[Path] = typer.Option(
        None,
        "--on-target-risearch-file",
        "-on-ris",
        help="Path to pre-computed RIsearch output file for On-Target (skips on-the-fly calculation).",
        exists=True,
        readable=True,
    ),
    query_file: Path = typer.Option(
        None,
        "-q",
        "--query",
        help="Path to siRNA query FASTA (required if --on-target is used).",
        exists=True,
        readable=True,
    ),
    on_target_expression: float = typer.Option(
        1000.0,
        "--on-target-expression",
        "-oexp",
        help="Expression level for On-Target (default: 1000.0).",
    ),
    on_target_accessibility: Optional[Path] = typer.Option(
        None,
        "--on-target-accessibility",
        help="Path to accessibility file for On-Target (text or binary).",
    ),
    legacy_format: bool = typer.Option(
        False,
        "--legacy-format",
        help="Output in legacy format (gw.results style, aggregated by transcript).",
    ),
    detailed_report: bool = typer.Option(
        False,
        "--detailed-report",
        help="Report off-target probabilities for individual transcripts (Legacy Format).",
    ),
    sense_only: bool = typer.Option(
        False,
        "--sense-only",
        help="Limit RIsearch2 predictions to sense strand (+), ignore antisense.",
    ),
    predictions_type: str = typer.Option(
        "gw",
        "--type",
        help="Type of RIsearch2 predictions. gw for genome-wide and tw for transcriptome-wide.",
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

        # Apply --sense-only filter (Sense strand only)
        if sense_only:
            df = df.filter(pl.col("strand") == "+")

        typer.echo("--- RIsearch Predictions (First 5 rows) ---")
        typer.echo(df.head(5))
        typer.echo("------------------------------------------")

        summary_ris = risearch_parser.summary(df)

        typer.echo(
            f"✓ Loaded {summary_ris['row_count']} predictions from {risearch_file.name}"
        )
        typer.echo(
            f"  Energy range: {summary_ris['energy_min']:.2f} to {summary_ris['energy_max']:.2f} kcal/mol"
        )
        if sense_only:
            typer.echo("  (Filtered to sense strand only)")

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
            # load_gtf handles BED files too now (by extension)
            df_trans = trans_parser.load_gtf(
                gtf_file, feature=feature_type, score_col=expression_metric
            )
            typer.echo("--- Transcriptome Data (First 5 rows) ---")
            typer.echo(df_trans.head(5))
            typer.echo("----------------------------------------")
            summary_trans = trans_parser.summary(df_trans)

            typer.echo(
                f"✓ Loaded {summary_trans['row_count']} features from {gtf_file.name}"
            )

            # Perform intersection
            typer.echo(
                f"  Intersecting predictions with transcriptome (mode={predictions_type})..."
            )
            intersector = IntersectionService()
            df = intersector.intersect(df, df_trans, mode=predictions_type)

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
        df = prob_service.calculate_probabilities(
            df,
            on_target_path=on_target_file,
            query_path=query_file,
            on_target_expression=on_target_expression,
            on_target_accessibility_path=on_target_accessibility,
            on_target_risearch_path=on_target_risearch_file,
        )

        # Legacy Format Output
        if legacy_format:
            # Extract siRNA ID from query or use default
            sirna_id = "siRNA"
            if query_file:
                sirna_id = query_file.stem

            legacy_output = prob_service.calculate_legacy_format(
                df,
                sirna_id=sirna_id,
                on_target_path=on_target_file,
                query_path=query_file,
                on_target_expression=on_target_expression,
                on_target_accessibility_path=on_target_accessibility,
                on_target_risearch_path=on_target_risearch_file,
                verbose=detailed_report,
            )

            if output_file:
                legacy_path = output_file.with_suffix(".results")
                with open(legacy_path, "w") as f:
                    f.write(legacy_output)
                typer.echo(f"✓ Legacy format results saved to {legacy_path}")
            else:
                typer.echo(legacy_output)

            # Also save TSV if requested (Standard flow continues below? No, return in orig)
            if output_file:
                df.write_csv(output_file, separator="\t")
                typer.echo(f"✓ Detailed results saved to {output_file}")
            return

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
