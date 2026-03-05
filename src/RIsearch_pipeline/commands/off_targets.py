from pathlib import Path
import tempfile
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


def _downcast_schema(df: pl.DataFrame) -> pl.DataFrame:
    """Downcast columns to minimize Arrow IPC file size.

    - Coordinates: start/end → UInt32 (covers positions up to 4.3B)
    - Strings: sirna_id, transcript_id, chrom, strand, gene_id → Categorical
    - Floats: energy, opening_energy, dG_total, P_off_target, exp_value → Float32

    Complexity: O(n) time, O(1) extra space.
    """
    casts = {}
    for col in ["start", "end"]:
        if col in df.columns:
            casts[col] = pl.UInt32
    for col in ["sirna_id", "transcript_id", "chrom", "strand", "gene_id"]:
        if col in df.columns:
            casts[col] = pl.Categorical
    for col in [
        "energy",
        "opening_energy",
        "dG_total",
        "P_off_target",
        "exp_value",
    ]:
        if col in df.columns:
            casts[col] = pl.Float32
    if casts:
        df = df.cast(casts)
    return df


def run(
    risearch_file: Optional[Path] = typer.Option(
        None,
        "-r",
        "--risearch-file",
        help="Path to pre-computed RIsearch output file or directory of files.",
        exists=True,
        readable=True,
        dir_okay=True,
    ),
    input_dir: Optional[Path] = typer.Option(
        None,
        "-d",
        "--input-dir",
        help="Directory containing per-siRNA RIsearch output files (*.gz or *.tsv).",
        exists=True,
        file_okay=False,
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
    transcriptome_format: str = typer.Option(
        "auto",
        "--transcriptome-format",
        help="Transcriptome file format: auto, bed6, bed7, or gtf (default: auto-detect).",
    ),
    accessibility_dir: Path = typer.Option(
        None,
        "-a",
        "--accessibility-dir",
        help="Directory containing pre-computed accessibility profiles to calculate P(OT).",
        exists=True,
        file_okay=False,
    ),
    accessibility_file: Path = typer.Option(
        None,
        "--accessibility-file",
        help="Parquet file with pre-computed binding-site accessibility (from 'accessibility --risearch-dir').",
        exists=True,
        dir_okay=False,
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
    on_target_ids_file: Optional[Path] = typer.Option(
        None,
        "-oi",
        "--on-target-ids",
        help="TSV file mapping siRNA IDs to on-target transcript IDs (sirna_id \\t transcript_id).",
        exists=True,
        readable=True,
    ),
    alpha: str = typer.Option(
        "1.0",
        "--alpha",
        help="Alpha clamping parameter(s). Separate multiple values with ';' (e.g., '0.8;1.0').",
    ),
    gamma: str = typer.Option(
        "1.0",
        "--gamma",
        help="Gamma clamping parameter(s). Separate multiple values with ';' (e.g., '0.8;1.0').",
    ),
    theta: str = typer.Option(
        "",
        "--theta",
        help="Theta scaling parameter(s). Separate multiple values with ';' (e.g., '0.5;0.7').",
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
    use_rnaplfold_cli: bool = typer.Option(
        False,
        "--use-rnaplfold-cli",
        help="Use RNAplfold binary instead of ViennaRNA Python bindings for accessibility.",
    ),
    chunk_mode: bool = typer.Option(
        False,
        "--chunk-mode",
        help="Process siRNAs in batches for large files (reduces memory usage).",
    ),
    batch_size: int = typer.Option(
        50,
        "--batch-size",
        "-b",
        help="Number of siRNAs to process per batch in chunk mode (default: 50).",
    ),
    scratch_dir: Optional[Path] = typer.Option(
        None,
        "--scratch-dir",
        help="Directory for intermediate Arrow IPC files. Defaults to system temp.",
    ),
    output_format: str = typer.Option(
        "tsv",
        "--output-format",
        help="Final output format: 'tsv', 'csv', or 'parquet'.",
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
        # Handle the case where the user passes a directory to --risearch-file
        if risearch_file is not None and risearch_file.is_dir():
            input_dir = risearch_file
            risearch_file = None

        if risearch_file is None and sirna_fasta is None and input_dir is None:
            console.print(
                "[bold red]Error:[/bold red] Must provide either --risearch-file, --input-dir, OR --sirna-fasta"
            )
            raise typer.Exit(code=1)

        # Integrated RIsearch execution requires target_fasta
        is_running_risearch = (
            (sirna_fasta is not None)
            and (risearch_file is None)
            and (input_dir is None)
        )
        if is_running_risearch and target_fasta is None:
            console.print(
                "[bold red]Error:[/bold red] --sirna-fasta requires --target-fasta when running RIsearch dynamically"
            )
            raise typer.Exit(code=1)

        # Mode 1: Integrated RIsearch execution
        if is_running_risearch:
            from RIsearch_pipeline.services.risearch_service import RIsearchService

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

                with console.status("[bold green]Running RIsearch..."):
                    risearch_service.run_search(
                        query_path=sirna_fasta,
                        index_path=index_path,
                        output_path=output_path,
                    )

                # Load predictions
                df = risearch_parser.load(output_path)

        # Mode 3: Directory of per-siRNA files
        elif input_dir is not None:
            if not isinstance(input_dir, Path):
                input_dir = Path(input_dir)

            console.print(
                f"[bold cyan]Directory mode[/bold cyan] - loading files from {input_dir}"
            )

            # List all files
            all_files = risearch_parser.list_directory_files(input_dir)
            if not all_files:
                console.print(
                    f"[bold red]Error:[/bold red] No RIsearch files found in {input_dir}"
                )
                raise typer.Exit(code=1)

            console.print(
                f"[green]✓[/green] Found [bold]{len(all_files)}[/bold] files to process"
            )

            # Prepare services
            from RIsearch_pipeline.services.probability import ProbabilityService
            from RIsearch_pipeline.services.accessibility import (
                GenomeAccessibilityService,
            )

            acc_service = None
            precomputed_acc = None
            if accessibility_file:
                precomputed_acc = pl.read_parquet(accessibility_file)
                console.print(
                    f"  [dim]Loaded {precomputed_acc.height} precomputed accessibility values[/dim]"
                )
            elif accessibility_dir:
                acc_service = GenomeAccessibilityService(
                    accessibility_dir, max_cached=4
                )
            prob_service = ProbabilityService(
                acc_service,
                use_rnaplfold_cli=use_rnaplfold_cli,
                precomputed_accessibility=precomputed_acc,
            )

            # Prepare transcriptome if provided
            df_trans = None
            intersector = None
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
                    gtf_file,
                    feature=feature_type,
                    score_col=expression_metric,
                    format=transcriptome_format,
                )
                intersector = IntersectionService()
                console.print(
                    f"[green]✓[/green] Loaded transcriptome with {df_trans.height} features"
                )

            # Parse on-target map
            on_target_map = {}
            if on_target_ids_file is not None:
                with open(on_target_ids_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith("#"):
                            parts = line.split("\t")
                            if len(parts) >= 2:
                                on_target_map[parts[0]] = parts[1]
                console.print(
                    f"  [dim]Loaded {len(on_target_map)} on-target mappings[/dim]"
                )

            # Parse parameter lists
            alpha_vals = [float(x) for x in alpha.split(";") if x.strip()]
            gamma_vals = [float(x) for x in gamma.split(";") if x.strip()]
            theta_vals = [float(x) for x in theta.split(";") if x.strip()]

            alpha_gamma_pairs = [(1.0, 1.0)]
            for a in alpha_vals:
                for g in gamma_vals:
                    if a == 1.0 and g == 1.0:
                        continue
                    if a <= g:
                        alpha_gamma_pairs.append((a, g))
            alpha_gamma_pairs = list(dict.fromkeys(alpha_gamma_pairs))

            # Process in batches - write Arrow IPC intermediates to scratch disk
            num_batches = (len(all_files) + batch_size - 1) // batch_size
            total_rows = 0

            console.print(
                f"  └─ Processing in {num_batches} batches of up to {batch_size} files"
            )

            # Resolve scratch directory (physical disk, not /dev/shm)
            scratch_base = str(scratch_dir) if scratch_dir else None

            # Storage for global metadata
            global_z_stats = {}
            global_on_w_stats = {}

            with tempfile.TemporaryDirectory(
                dir=scratch_base, prefix="risearch_ipc_"
            ) as tmp_scratch:
                scratch_path = Path(tmp_scratch)
                ipc_idx = 0

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TextColumn("[progress.percentage]{task.completed}/{task.total}"),
                    TimeElapsedColumn(),
                    console=console,
                ) as progress:
                    task = progress.add_task("Processing batches...", total=num_batches)

                    for batch_idx in range(num_batches):
                        start_idx = batch_idx * batch_size
                        end_idx = min(start_idx + batch_size, len(all_files))
                        batch_files = all_files[start_idx:end_idx]

                        # Load entire batch - Polars parallelizes via Rayon
                        df_chunk = risearch_parser.load_directory_batch(
                            input_dir, batch_files
                        )

                        if sense_only:
                            df_chunk = df_chunk.filter(pl.col("strand") == "+")

                        if df_chunk.height == 0:
                            progress.update(
                                task,
                                advance=1,
                                description=f"Batch {batch_idx + 1}: empty",
                            )
                            continue

                        # Accumulate intersection for the full batch to compute probabilities correctly
                        if df_trans is not None and intersector is not None:
                            intersected_chunks = []
                            for intersect_batch in intersector.intersect_streaming(
                                df_chunk, df_trans, mode=predictions_type
                            ):
                                if intersect_batch.height > 0:
                                    intersected_chunks.append(intersect_batch)

                            if not intersected_chunks:
                                progress.update(
                                    task,
                                    advance=1,
                                    description=f"Batch {batch_idx + 1}/{num_batches}",
                                )
                                continue

                            df_chunk = pl.concat(intersected_chunks)
                            del intersected_chunks

                        df_chunk, batch_metadata = (
                            prob_service.calculate_probabilities_per_sirna(
                                df_chunk,
                                alpha_gamma_pairs=alpha_gamma_pairs,
                                theta_values=theta_vals,
                                on_target_expression=on_target_expression,
                                on_target_map=on_target_map,
                            )
                        )

                        # Accumulate Z metadata
                        for sid, z_dict in batch_metadata["z_per_sirna"].items():
                            global_z_stats[sid] = z_dict
                        for sid, w_dict in batch_metadata["on_target_weights"].items():
                            global_on_w_stats[sid] = w_dict

                        # Drop heavy intermediate columns before serialization
                        drop_cols = [
                            c
                            for c in df_chunk.columns
                            if c.startswith("boltzmann_weight")
                            or c.startswith("Z_sirna")
                            or c == "E_min"
                        ]
                        if drop_cols:
                            df_chunk = df_chunk.drop(drop_cols)

                        # Downcast and write Arrow IPC to scratch disk
                        df_chunk = _downcast_schema(df_chunk)
                        ipc_file = scratch_path / f"batch_{ipc_idx:06d}.arrow"
                        df_chunk.write_ipc(ipc_file)
                        total_rows += df_chunk.height
                        ipc_idx += 1

                        # Free memory
                        del df_chunk

                        progress.update(
                            task,
                            advance=1,
                            description=f"Batch {batch_idx + 1}/{num_batches}",
                        )

                # Final merge: lazily scan all IPC files and write output
                if ipc_idx > 0 and output_file:
                    console.print(
                        f"  └─ Merging {ipc_idx} intermediate files ({total_rows} rows)..."
                    )
                    ipc_pattern = sorted(scratch_path.glob("batch_*.arrow"))
                    lf = pl.scan_ipc(ipc_pattern, memory_map=True)

                    if output_format == "parquet":
                        out_path = output_file.with_suffix(".parquet")
                        lf.sink_parquet(out_path)
                    else:
                        sep = "\t" if output_format == "tsv" else ","
                        suffix = ".tsv" if output_format == "tsv" else ".csv"
                        out_path = output_file.with_suffix(suffix)
                        final_df = lf.collect(streaming=True)
                        final_df.write_csv(out_path, separator=sep)
                        del final_df

                    console.print(
                        f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions → [bold]{out_path}[/bold]"
                    )

                    # Write .summary files next to the output_file's directory
                    out_dir = output_file.parent
                    out_dir.mkdir(parents=True, exist_ok=True)
                    for sid in global_z_stats:
                        summary_text = prob_service.format_legacy_summary(
                            sid,
                            global_z_stats[sid],
                            global_on_w_stats.get(sid, {}),
                            alpha_gamma_pairs,
                            theta_vals,
                        )
                        summary_file = out_dir / f"{sid}.summary"
                        summary_file.write_text(summary_text)

                    console.print(
                        f"[green]✓[/green] Wrote {len(global_z_stats)} .summary files to [bold]{out_dir}[/bold]"
                    )
                elif ipc_idx > 0:
                    console.print(
                        f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions (no output file specified)"
                    )
                elif total_rows == 0:
                    console.print("[yellow]Warning:[/yellow] No predictions to process")

            return  # Exit after directory mode processing

        # Mode 2: Pre-computed RIsearch file
        elif risearch_file is not None:
            if not isinstance(risearch_file, Path):
                risearch_file = Path(risearch_file)

            # CHUNK MODE: Process siRNAs one at a time for large files
            if chunk_mode:
                console.print(
                    "[bold cyan]Chunk mode enabled[/bold cyan] - processing siRNAs individually"
                )

                # Step 1: Stream-scan to get unique siRNA IDs (low memory)
                with console.status("[bold green]Scanning for siRNA IDs..."):
                    sirna_ids = risearch_parser.scan_sirna_ids(risearch_file)

                console.print(
                    f"[green]✓[/green] Found [bold]{len(sirna_ids)}[/bold] unique siRNAs to process"
                )

                # Prepare services
                from RIsearch_pipeline.services.probability import ProbabilityService
                from RIsearch_pipeline.services.accessibility import (
                    GenomeAccessibilityService,
                )

                acc_service = None
                precomputed_acc = None
                if accessibility_file:
                    precomputed_acc = pl.read_parquet(accessibility_file)
                elif accessibility_dir:
                    acc_service = GenomeAccessibilityService(
                        accessibility_dir, max_cached=4
                    )
                prob_service = ProbabilityService(
                    acc_service,
                    use_rnaplfold_cli=use_rnaplfold_cli,
                    precomputed_accessibility=precomputed_acc,
                )

                # Prepare transcriptome if provided
                df_trans = None
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
                        gtf_file,
                        feature=feature_type,
                        score_col=expression_metric,
                        format=transcriptome_format,
                    )
                    intersector = IntersectionService()
                    console.print(
                        f"[green]✓[/green] Loaded transcriptome with {df_trans.height} features"
                    )

                # Parse on-target IDs map
                on_target_map = {}
                if on_target_ids_file is not None:
                    with open(on_target_ids_file, "r") as f:
                        for line in f:
                            line = line.strip()
                            if line and not line.startswith("#"):
                                parts = line.split("\t")
                                if len(parts) >= 2:
                                    on_target_map[parts[0]] = parts[1]

                # Parse parameter lists
                alpha_vals = [float(x) for x in alpha.split(";") if x.strip()]
                gamma_vals = [float(x) for x in gamma.split(";") if x.strip()]
                theta_vals = [float(x) for x in theta.split(";") if x.strip()]

                alpha_gamma_pairs = [(1.0, 1.0)]
                for a in alpha_vals:
                    for g in gamma_vals:
                        if a == 1.0 and g == 1.0:
                            continue
                        if a <= g:
                            alpha_gamma_pairs.append((a, g))
                alpha_gamma_pairs = list(dict.fromkeys(alpha_gamma_pairs))

                # Storage for global metadata
                global_z_stats = {}
                global_on_w_stats = {}

                # Step 2: Process siRNAs in batches
                all_results = []
                num_batches = (
                    len(sirna_ids) + batch_size - 1
                ) // batch_size  # Ceiling division

                console.print(
                    f"  └─ Processing {len(sirna_ids)} siRNAs in {num_batches} batches of up to {batch_size}"
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TextColumn("[progress.percentage]{task.completed}/{task.total}"),
                    TimeElapsedColumn(),
                    console=console,
                ) as progress:
                    task = progress.add_task("Processing batches...", total=num_batches)

                    for batch_idx in range(num_batches):
                        # Get siRNAs for this batch
                        start_idx = batch_idx * batch_size
                        end_idx = min(start_idx + batch_size, len(sirna_ids))
                        batch_sirnas = sirna_ids[start_idx:end_idx]

                        # Load entire batch at once (enables Polars/Rayon parallelization)
                        df_chunk = risearch_parser.load_by_sirna_batch(
                            risearch_file, batch_sirnas
                        )

                        if sense_only:
                            df_chunk = df_chunk.filter(pl.col("strand") == "+")

                        if df_chunk.height == 0:
                            progress.update(
                                task,
                                advance=1,
                                description=f"Batch {batch_idx + 1}: empty",
                            )
                            continue

                        # Intersect with transcriptome - Polars parallelizes this
                        if df_trans is not None:
                            df_chunk = intersector.intersect(
                                df_chunk, df_trans, mode=predictions_type
                            )

                        # Build on-target map for all siRNAs in this batch
                        batch_map = None
                        if on_target_map:
                            batch_map = {
                                sid: on_target_map.get(sid)
                                for sid in batch_sirnas
                                if sid in on_target_map
                            }

                        # Calculate probabilities - Polars parallelizes across all rows
                        df_chunk, batch_metadata = (
                            prob_service.calculate_probabilities_per_sirna(
                                df_chunk,
                                alpha_gamma_pairs=alpha_gamma_pairs,
                                theta_values=theta_vals,
                                on_target_map=batch_map if batch_map else None,
                                on_target_expression=on_target_expression,
                            )
                        )

                        # Accumulate metadata
                        for sid, z_dict in batch_metadata["z_per_sirna"].items():
                            global_z_stats[sid] = z_dict
                        for sid, w_dict in batch_metadata["on_target_weights"].items():
                            global_on_w_stats[sid] = w_dict

                        all_results.append(df_chunk)
                        progress.update(
                            task,
                            advance=1,
                            description=f"Batch {batch_idx + 1}/{num_batches}",
                        )

                # Step 3: Combine all results
                if all_results:
                    df = pl.concat(all_results, how="diagonal")
                    console.print(
                        f"[green]✓[/green] Processed [bold]{df.height}[/bold] total predictions"
                    )
                else:
                    console.print("[yellow]Warning:[/yellow] No predictions to process")
                    return

                # Save output
                if output_file:
                    df.write_csv(output_file, separator="\t")
                    console.print(
                        f"[green]✓[/green] Results saved to [bold]{output_file}[/bold]"
                    )

                    # Write .summary files
                    out_dir = output_file.parent
                    out_dir.mkdir(parents=True, exist_ok=True)
                    for sid in global_z_stats:
                        summary_text = prob_service.format_legacy_summary(
                            sid,
                            global_z_stats[sid],
                            global_on_w_stats.get(sid, {}),
                            alpha_gamma_pairs,
                            theta_vals,
                        )
                        summary_file = out_dir / f"{sid}.summary"
                        summary_file.write_text(summary_text)

                    console.print(
                        f"[green]✓[/green] Wrote {len(global_z_stats)} .summary files to [bold]{out_dir}[/bold]"
                    )
                else:
                    console.print(df.head(10))

                return  # Exit after chunk mode processing

            # Standard mode: Load entire file
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
                gtf_file,
                feature=feature_type,
                score_col=expression_metric,
                format=transcriptome_format,
            )

            if verbose:
                console.print("[dim]--- Transcriptome Data (First 5 rows) ---[/dim]")
                console.print(df_trans.head(5))

            summary_trans = trans_parser.summary(df_trans)
            console.print(
                f"[green]✓[/green] Loaded [bold]{summary_trans['row_count']}[/bold] features from {gtf_file.name}"
            )

            # Perform intersection
            with console.status(
                f"[bold green]Intersecting predictions (mode={predictions_type})..."
            ):
                intersector = IntersectionService()
                df = intersector.intersect(df, df_trans, mode=predictions_type)

            console.print(
                f"[green]✓[/green] Found [bold]{df.height}[/bold] intersecting off-target candidates"
            )

        # Accessibility and Probability Calculation
        from RIsearch_pipeline.services.probability import ProbabilityService
        from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService

        # Keep reference to temp dir object so it persists until function exit
        temp_dir_obj = None

        if accessibility_file:
            precomputed_acc = pl.read_parquet(accessibility_file)
            console.print(
                f"  [dim]Loaded {precomputed_acc.height} precomputed accessibility values from {accessibility_file}...[/dim]"
            )
            prob_service = ProbabilityService(
                None,
                use_rnaplfold_cli=use_rnaplfold_cli,
                precomputed_accessibility=precomputed_acc,
            )

        elif accessibility_dir:
            console.print(
                f"  [dim]Calculating probabilities using profiles from {accessibility_dir}...[/dim]"
            )
            acc_service = GenomeAccessibilityService(accessibility_dir, max_cached=4)
            prob_service = ProbabilityService(
                acc_service, use_rnaplfold_cli=use_rnaplfold_cli
            )

        elif genome_file:
            console.print(
                f"  [dim]Computing accessibility on-the-fly from {genome_file} (this may take time)...[/dim]"
            )

            # Create temp dir
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="risearch_accessibility_")
            temp_dir = Path(temp_dir_obj.name)

            acc_service = GenomeAccessibilityService(temp_dir, max_cached=4)

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

            prob_service = ProbabilityService(
                acc_service, use_rnaplfold_cli=use_rnaplfold_cli
            )

        else:
            console.print(
                "[yellow]Warning:[/yellow] No accessibility data provided. P(OT) based on energy only."
            )
            prob_service = ProbabilityService(None, use_rnaplfold_cli=use_rnaplfold_cli)

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

        # Parse parameter lists
        alpha_vals = [float(x) for x in alpha.split(";") if x.strip()]
        gamma_vals = [float(x) for x in gamma.split(";") if x.strip()]
        theta_vals = [float(x) for x in theta.split(";") if x.strip()]

        # Construct alpha-gamma pairs (Cartesian product with a <= g constraint)
        alpha_gamma_pairs = []
        # Always include baseline (1.0, 1.0) first
        alpha_gamma_pairs.append((1.0, 1.0))

        for a in alpha_vals:
            for g in gamma_vals:
                if a == 1.0 and g == 1.0:
                    continue
                if a <= g:
                    alpha_gamma_pairs.append((a, g))

        # Remove duplicates while preserving order
        alpha_gamma_pairs = list(dict.fromkeys(alpha_gamma_pairs))

        # Detect multi-siRNA mode OR if custom parameters are used (forcing vectorized path)
        unique_sirnas = df["sirna_id"].unique() if "sirna_id" in df.columns else []
        has_custom_params = len(alpha_gamma_pairs) > 1 or len(theta_vals) > 0

        is_multi_sirna = len(unique_sirnas) > 1 or has_custom_params

        # Parse on-target IDs file if provided
        on_target_map: dict[str, str] = {}
        if on_target_ids_file is not None:
            try:
                mapping_df = pl.read_csv(
                    on_target_ids_file, separator="\t", has_header=False
                )
                on_target_map = dict(
                    zip(
                        mapping_df.get_column("column_1"),
                        mapping_df.get_column("column_2"),
                    )
                )
                console.print(
                    f"[green]✓[/green] Loaded {len(on_target_map)} on-target ID mappings from {on_target_ids_file.name}"
                )
            except Exception as e:
                console.print(
                    f"[bold red]Error parsing on-target mapping file:[/bold red] {e}"
                )
                raise typer.Exit(code=1)

        if is_multi_sirna:
            mode_msg = (
                "Multi-siRNA"
                if len(unique_sirnas) > 1
                else "Parameterized Single-siRNA"
            )
            console.print(f"  └─ {mode_msg} detection: {len(unique_sirnas)} siRNAs")

            df, meta = prob_service.calculate_probabilities_per_sirna(
                df,
                alpha_gamma_pairs=alpha_gamma_pairs,
                theta_values=theta_vals,
                on_target_map=on_target_map if on_target_map else None,
                on_target_expression=on_target_expression,
            )

            # Write .summary files if output_file is provided
            if output_file:
                out_dir = output_file.parent
                out_dir.mkdir(parents=True, exist_ok=True)
                global_z_stats = meta["z_per_sirna"]
                global_on_w_stats = meta["on_target_weights"]
                for sid in global_z_stats:
                    summary_text = prob_service.format_legacy_summary(
                        sid,
                        global_z_stats[sid],
                        global_on_w_stats.get(sid, {}),
                        alpha_gamma_pairs,
                        theta_vals,
                    )
                    summary_file = out_dir / f"{sid}.summary"
                    summary_file.write_text(summary_text)

            # Display per-siRNA Z summary
            console.print("  └─ Computing per-siRNA Boltzmann weights...")
            if has_custom_params:
                console.print(
                    f"  └─ Applied {len(alpha_gamma_pairs) - 1} alpha/gamma pairs and {len(theta_vals)} theta values"
                )

            # Compute total Z from per-siRNA data
            total_z = sum(
                z_dict.get("Z_sirna", 0.0) for z_dict in meta["z_per_sirna"].values()
            )
            console.print(
                f"  └─ Total Z across all siRNAs (base): [bold]{total_z:.2e}[/bold]"
            )
        else:
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

        # 4. On-Target Details (single siRNA mode only)
        if meta.get("has_on_target", False):
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
