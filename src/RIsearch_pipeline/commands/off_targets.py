from pathlib import Path
import typer
import polars as pl
from typing import Optional
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TimeElapsedColumn,
)

from RIsearch_pipeline.services.risearch_parser import RIsearchParser

console = Console()


def run(
    risearch_file: Optional[Path] = typer.Option(
        None,
        "-r",
        "--risearch-file",
        help="Path to pre-computed RIsearch output TSV file.",
        exists=True,
        readable=True,
    ),
    sirna_fasta: Optional[Path] = typer.Option(
        None,
        "-s",
        "--sirna-fasta",
        help="Path to siRNA FASTA file (one or more sequences). Runs RIsearch internally.",
        exists=True,
        readable=True,
    ),
    target_fasta: Optional[Path] = typer.Option(
        None,
        "--target-fasta",
        "--genome",
        help="Path to target FASTA (genome or transcriptome) for RIsearch.",
        exists=True,
        readable=True,
    ),
    target_index: Optional[Path] = typer.Option(
        None,
        "-idx",
        "--target-index",
        help="Pre-built RIsearch index (optional, speeds up repeated runs).",
        exists=True,
        readable=True,
    ),
    workers: Optional[int] = typer.Option(
        None,
        "-j",
        "--workers",
        help="Number of parallel workers for multi-siRNA processing (default: CPU count).",
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
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Show detailed output (tables, stats).",
    ),
) -> None:
    """
    Analyze siRNA off-target predictions.

    Integrates RIsearch2 predictions with transcriptome data and optionally
    calculates off-target probabilities.
    """
    risearch_parser = RIsearchParser()

    console.print(Panel("RIsearch Pipeline", style="bold cyan"))

    try:
        # Validate input options
        if risearch_file is None and sirna_fasta is None:
            console.print(
                "[bold red]Error:[/bold red] Must provide either --risearch-file OR --sirna-fasta"
            )
            raise typer.Exit(code=1)

        if sirna_fasta is not None and target_fasta is None:
            console.print(
                "[bold red]Error:[/bold red] --sirna-fasta requires --target-fasta"
            )
            raise typer.Exit(code=1)

        # Mode 1: Integrated RIsearch execution
        if sirna_fasta is not None:
            from RIsearch_pipeline.services.risearch_service import RIsearchService
            import tempfile

            risearch_service = RIsearchService()

            # Validate siRNA FASTA (checks for duplicates)
            try:
                sirna_ids = risearch_service.validate_sirna_fasta(sirna_fasta)
                console.print(
                    f"[green]✓[/green] Validated [bold]{len(sirna_ids)}[/bold] siRNA(s) from {sirna_fasta.name}"
                )
            except ValueError as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
                raise typer.Exit(code=1)

            # Index target (or reuse existing index)
            with console.status("[bold green]Indexing target...") as status:
                if target_index is not None:
                    index_path = target_index
                    console.print(
                        f"[green]✓[/green] Using pre-built index: {index_path.name}"
                    )
                else:
                    index_path = risearch_service.index_target(target_fasta)
                    console.print(f"[green]✓[/green] Created index: {index_path.name}")

            # Run RIsearch
            with tempfile.TemporaryDirectory(prefix="risearch_") as tmpdir:
                output_path = Path(tmpdir) / "predictions.out"

                with console.status("[bold green]Running RIsearch...") as status:
                    risearch_service.run_search(
                        query_path=sirna_fasta,
                        index_path=index_path,
                        output_path=output_path,
                    )

                # Load predictions
                df = risearch_parser.load(output_path)

        # Mode 2: Pre-computed RIsearch file
        else:
            if not isinstance(risearch_file, Path):
                risearch_file = Path(risearch_file)
            df = risearch_parser.load(risearch_file)

        # Apply --sense-only filter (Sense strand only)
        if sense_only:
            df = df.filter(pl.col("strand") == "+")

        # Verbose: Raw DF
        if verbose:
            console.print("[dim]--- RIsearch Predictions (First 5 rows) ---[/dim]")
            console.print(df.head(5))

        summary_ris = risearch_parser.summary(df)
        source_name = sirna_fasta.name if sirna_fasta else risearch_file.name
        console.print(
            f"[green]✓[/green] Loaded [bold]{summary_ris['row_count']}[/bold] predictions from {source_name}"
        )
        console.print(
            f"  └─ Energy range: {summary_ris['energy_min']:.2f} to {summary_ris['energy_max']:.2f} kcal/mol"
        )

        if sense_only:
            console.print("  └─ (Filtered to sense strand only)")

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

            if verbose:
                console.print(f"[dim]--- Transcriptome Data (First 5 rows) ---[/dim]")
                console.print(df_trans.head(5))

            summary_trans = trans_parser.summary(df_trans)
            console.print(
                f"[green]✓[/green] Loaded [bold]{summary_trans['row_count']}[/bold] features from {gtf_file.name}"
            )

            # Perform intersection
            with console.status(
                f"[bold green]Intersecting predictions (mode={predictions_type})..."
            ) as status:
                intersector = IntersectionService()
                df = intersector.intersect(df, df_trans, mode=predictions_type)

            console.print(
                f"[green]✓[/green] Found [bold]{df.height}[/bold] intersecting off-target candidates"
            )

        # Accessibility and Probability Calculation
        from RIsearch_pipeline.services.probability import ProbabilityService
        from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService
        import tempfile

        # Keep reference to temp dir object so it persists until function exit
        temp_dir_obj = None

        if accessibility_dir:
            console.print(
                f"  [dim]Calculating probabilities using profiles from {accessibility_dir}...[/dim]"
            )
            acc_service = GenomeAccessibilityService(accessibility_dir)
            prob_service = ProbabilityService(acc_service)

        elif genome_file:
            console.print(
                f"  [dim]Computing accessibility on-the-fly from {genome_file} (this may take time)...[/dim]"
            )

            # Create temp dir
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="risearch_accessibility_")
            temp_dir = Path(temp_dir_obj.name)

            acc_service = GenomeAccessibilityService(temp_dir)

            # Use progress bar for accessibility
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeElapsedColumn(),
                console=console,
            ) as progress:
                task = progress.add_task(
                    "Calculating Accessibility...", total=None
                )  # Total unknown initially

                def progress_callback(advance=0, description=""):
                    progress.update(task, advance=advance, description=description)

                acc_service.compute_genome_accessibility(
                    genome_file,
                    window_size=window_size,
                    max_span=max_span,
                    unpaired_prob=unpaired_prob,
                    progress_callback=progress_callback,
                )

            prob_service = ProbabilityService(acc_service)

        else:
            console.print(
                "[yellow]Warning:[/yellow] No accessibility data provided. P(OT) based on energy only."
            )
            prob_service = ProbabilityService(None)

        # Calculate P(OT)
        # Calculate P(OT)
        console.print("[bold green]Calculating probabilities...[/bold green]")

        # 1. Annotate Opening Energies
        if accessibility_dir or genome_file:
            console.print(
                "  └─ Annotating opening energies from accessibility profiles..."
            )
        else:
            console.print(
                "  └─ [yellow]Skipping accessibility annotation (energy only)[/yellow]"
            )

        df, meta = prob_service.calculate_probabilities(
            df,
            on_target_path=on_target_file,
            query_path=query_file,
            on_target_expression=on_target_expression,
            on_target_accessibility_path=on_target_accessibility,
            on_target_risearch_path=on_target_risearch_file,
        )

        # 2. Boltzmann Weights
        console.print("  └─ Computing Boltzmann weights (α=1.0, γ=1.0)...")

        # 3. Partition Function Info
        z_fmt = f"{meta['z_total']:.2e}"
        z_off = f"{meta['z_off_target']:.2e}"
        z_on = f"{meta['w_on_target']:.2e}"
        console.print(
            f"  └─ Partition Function Z = [bold]{z_fmt}[/bold] (Off-Target={z_off}, On-Target={z_on})"
        )

        # 4. On-Target Details
        if meta["has_on_target"]:
            dg_on = f"{meta['dG_on_target']:.2f}"
            # Use high precision to show deviation from 1.0
            p_on = f"{meta.get('p_on_target', 0.0):.10f}"
            console.print(f"  └─ On-Target: ΔG_total={dg_on} kcal/mol, P(on)={p_on}")

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
                console.print(
                    f"[green]✓[/green] Legacy format results saved to {legacy_path}"
                )
            else:
                console.print(legacy_output)

            # Also save TSV if requested (Standard flow continues below? No, return in orig)
            if output_file:
                detailed_path = output_file.parent / "detailed_results.tsv"
                df.write_csv(detailed_path, separator="\t")
                console.print(
                    f"[green]✓[/green] Detailed results saved to {detailed_path}"
                )
            return

        # Display Results Table
        if df.height > 0:
            if "P_off_target" in df.columns:
                df = df.sort("P_off_target", descending=True)
                title = "Top 10 Candidates (By P_off_target)"
            else:
                title = "Predictions Preview"

            table = Table(title=title)

            # Select useful columns to show
            potential_cols = [
                "sirna_id",
                "chrom",
                "gene_id",
                "transcript_id",
                "energy",
                "dG_total",
                "P_off_target",
                "exp_value",
            ]
            display_cols = [c for c in potential_cols if c in df.columns]

            for col in display_cols:
                table.add_column(col)

            for row in df.select(display_cols).head(10).iter_rows():
                # Format floats nicely
                formatted_row = []
                for val in row:
                    if isinstance(val, float):
                        formatted_row.append(f"{val:.4g}")
                    else:
                        formatted_row.append(str(val))
                table.add_row(*formatted_row)

            console.print(table)

        # Save to Output File
        if output_file:
            df.write_csv(output_file, separator="\t")
            console.print(
                f"\n[green]✓[/green] Results saved to [bold]{output_file}[/bold]"
            )

    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]Failed:[/bold red] {e}")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1)
