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

from loguru import logger

from RIsearch_pipeline.services.risearch_parser import RIsearchParser
from RIsearch_pipeline.services.profiling import PipelineProfiler

console = Console()

# ---------------------------------------------------------------------------
# Per-worker state for directory-mode parallel processing (spawn-safe).
# _init_worker() runs once per spawned worker process and populates these
# module-level globals.  spawn (not fork) is used to avoid Rayon/Polars
# thread-pool deadlocks that occur when forking a process that already has
# background threads active.
# ---------------------------------------------------------------------------
_WORKER_DF_TRANS: Optional[pl.DataFrame] = None
_WORKER_INTERSECTOR: Optional[object] = None
_WORKER_PROB_SERVICE: Optional[object] = None


def _init_worker(
    gtf_path: str,
    feature: str,
    score_col: str,
    transcriptome_format: str,
    accessibility_dir: str,
    accessibility_file: str,
    use_rnaplfold_cli: bool,
    polars_max_threads: int,
) -> None:
    """Initialise per-worker state.  Called once per spawned worker process."""
    import os
    os.environ["POLARS_MAX_THREADS"] = str(polars_max_threads)

    global _WORKER_DF_TRANS, _WORKER_INTERSECTOR, _WORKER_PROB_SERVICE

    from RIsearch_pipeline.services.transcriptome_parser import TranscriptomeParser
    from RIsearch_pipeline.services.intersection_service import IntersectionService
    from RIsearch_pipeline.services.probability import ProbabilityService

    if gtf_path:
        _WORKER_DF_TRANS = TranscriptomeParser().load_gtf(
            Path(gtf_path),
            feature=feature,
            score_col=score_col,
            format=transcriptome_format,
        )
        _WORKER_INTERSECTOR = IntersectionService()

    acc_service = None
    precomputed_acc = None
    if accessibility_file:
        precomputed_acc = pl.read_parquet(accessibility_file)
    elif accessibility_dir:
        from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService
        acc_service = GenomeAccessibilityService(Path(accessibility_dir), max_cached=4)

    _WORKER_PROB_SERVICE = ProbabilityService(
        acc_service,
        use_rnaplfold_cli=use_rnaplfold_cli,
        precomputed_accessibility=precomputed_acc,
    )


def _process_one_sirna(
    file_path_str: str,
    on_target_map: dict,
    self_hyb_emin: dict,
    alpha_gamma_pairs: list,
    theta_vals: list,
    on_target_expression: float,
    sense_only: bool,
    predictions_type: str,
) -> tuple:
    """Process a single siRNA file through the full pipeline.

    Runs in a spawned worker process.  Reads worker state from module-level
    globals populated by _init_worker().

    Returns:
        (PyArrow Table, metadata dict) or (None, {}) if no predictions remain.
    """
    import sys
    _m = sys.modules[__name__]

    df_trans = _m._WORKER_DF_TRANS
    intersector = _m._WORKER_INTERSECTOR
    prob_service = _m._WORKER_PROB_SERVICE

    from RIsearch_pipeline.services.risearch_parser import RIsearchParser

    parser = RIsearchParser()
    df = parser.load_single_file(Path(file_path_str))

    if sense_only:
        df = df.filter(pl.col("strand") == "+")
    if df.height == 0:
        return None, {}

    # Apply self-hybridisation E_min override (one value per siRNA file)
    if self_hyb_emin and "raw_e_min" in df.columns:
        sirna_id = str(df["sirna_id"][0])
        if sirna_id in self_hyb_emin:
            df = df.with_columns(
                pl.lit(self_hyb_emin[sirna_id]).cast(pl.Float32).alias("raw_e_min")
            )

    # Intersect with transcriptome (sequential — each worker owns its CPU)
    if df_trans is not None and intersector is not None:
        df = intersector.intersect(df, df_trans, mode=predictions_type, workers=1)
        if df.height == 0:
            return None, {}
        df = df.unique(
            subset=["sirna_id", "chrom", "start", "end", "strand", "energy", "transcript_id"],
            keep="first",
        )

    # Probabilities
    df, meta = prob_service.calculate_probabilities_per_sirna(
        df,
        alpha_gamma_pairs=alpha_gamma_pairs,
        theta_values=theta_vals,
        on_target_map=on_target_map,
        on_target_expression=on_target_expression,
    )

    # Drop heavy intermediate columns before returning
    drop_cols = [
        c for c in df.columns
        if c.startswith("boltzmann_weight") or c.startswith("Z_sirna") or c == "E_min"
    ]
    if drop_cols:
        df = df.drop(drop_cols)

    return df.to_arrow(), meta


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
        help="Path to siRNA FASTA file (one or more sequences). Runs RIsearch internally when used alone; also used in directory mode to compute self-hybridisation E_min.",
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
        help="Number of parallel threads for intersection (default: CPU count).",
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
    profile: bool = typer.Option(
        False,
        "--profile",
        help="Print a per-stage timing and memory profile table at the end of the run.",
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

    profiler = PipelineProfiler(enabled=profile)

    import os
    n_workers: int = workers if workers is not None else os.cpu_count() or 1

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
            with console.status("[bold green]Running RIsearch..."):
                df = risearch_service.run_search(
                    query_path=sirna_fasta,
                    index_path=index_path,
                    target_fasta=target_fasta,
                )

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
                with profiler.stage("Load transcriptome") as _s:
                    df_trans = trans_parser.load_gtf(
                        gtf_file,
                        feature=feature_type,
                        score_col=expression_metric,
                        format=transcriptome_format,
                    )
                    _s.rows_out = df_trans.height
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

            # Compute self-hybridisation E_min for each siRNA if FASTA provided.
            # This replicates the old pipeline's behaviour where E_min is the minimum
            # energy of the siRNA hybridised against itself (not the genome-wide min).
            self_hyb_emin: dict[str, float] = {}
            if sirna_fasta is not None:
                from RIsearch_pipeline.services.risearch_service import RIsearchService
                _ris_svc = RIsearchService()
                console.print(
                    f"  └─ Computing self-hybridisation E_min from [bold]{sirna_fasta.name}[/bold]..."
                )
                self_hyb_emin = _ris_svc.self_hybridization_emin_batch(sirna_fasta)
                console.print(
                    f"  [dim]Self-hyb E_min computed for {len(self_hyb_emin)} siRNA(s)[/dim]"
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

            total_rows = 0
            n_proc = min(n_workers, len(all_files))

            console.print(
                f"  └─ Processing {len(all_files)} siRNAs with up to {n_proc} parallel workers"
            )

            # Resolve output path up-front so each batch can stream directly
            out_path = None
            _csv_first_batch = True
            _pq_writer = None
            if output_file:
                if output_format == "parquet":
                    out_path = output_file.with_suffix(".parquet")
                elif output_format == "tsv":
                    out_path = output_file.with_suffix(".tsv")
                else:
                    out_path = output_file.with_suffix(".csv")
                out_path.parent.mkdir(parents=True, exist_ok=True)

            # Storage for global metadata
            global_z_stats = {}
            global_on_w_stats = {}

            try:
                from concurrent.futures import ProcessPoolExecutor, as_completed
                import multiprocessing as _mp

                # Polars threads per worker: 2 keeps N_workers × 2 ≤ available cores.
                polars_threads_per_worker = max(1, n_workers // n_proc)

                # Build initializer args so each spawned worker loads the
                # transcriptome independently (avoids fork+Rayon deadlocks).
                _init_args = (
                    str(gtf_file) if gtf_file else "",
                    feature_type,
                    expression_metric,
                    transcriptome_format,
                    str(accessibility_dir) if accessibility_dir else "",
                    str(accessibility_file) if accessibility_file else "",
                    use_rnaplfold_cli,
                    polars_threads_per_worker,
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TextColumn("[progress.percentage]{task.completed}/{task.total}"),
                    TimeElapsedColumn(),
                    console=console,
                ) as progress:
                    task = progress.add_task(
                        f"Processing siRNAs...", total=len(all_files)
                    )

                    with ProcessPoolExecutor(
                        max_workers=n_proc,
                        mp_context=_mp.get_context("spawn"),
                        initializer=_init_worker,
                        initargs=_init_args,
                    ) as pool:
                        futures = {
                            pool.submit(
                                _process_one_sirna,
                                str(f),
                                on_target_map,
                                self_hyb_emin,
                                alpha_gamma_pairs,
                                theta_vals,
                                on_target_expression,
                                sense_only,
                                predictions_type,
                            ): f
                            for f in all_files
                        }

                        completed = 0
                        for future in as_completed(futures):
                            completed += 1
                            progress.update(
                                task,
                                advance=1,
                                description=f"Batch {completed}/{len(all_files)}",
                            )

                            try:
                                arrow_table, batch_metadata = future.result()
                            except Exception as exc:
                                f = futures[future]
                                logger.error(f"Worker failed for {f.name}: {exc}")
                                continue

                            if arrow_table is None:
                                continue

                            df_chunk = pl.from_arrow(arrow_table)
                            total_rows += df_chunk.height

                            # Accumulate Z metadata
                            for sid, z_dict in batch_metadata.get("z_per_sirna", {}).items():
                                global_z_stats[sid] = z_dict
                            for sid, w_dict in batch_metadata.get("on_target_weights", {}).items():
                                global_on_w_stats[sid] = w_dict

                            if out_path is not None:
                                if output_format == "parquet":
                                    arrow_tbl = _downcast_schema(df_chunk).to_arrow()
                                    if _pq_writer is None:
                                        import pyarrow.parquet as pq
                                        _pq_writer = pq.ParquetWriter(out_path, arrow_tbl.schema)
                                    _pq_writer.write_table(arrow_tbl)
                                else:
                                    sep = "\t" if output_format == "tsv" else ","
                                    with open(out_path, "w" if _csv_first_batch else "a") as _f:
                                        df_chunk.write_csv(_f, separator=sep, include_header=_csv_first_batch)
                                    _csv_first_batch = False

                            del df_chunk

                # Report results
                if total_rows > 0 and out_path is not None:
                    console.print(
                        f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions → [bold]{out_path}[/bold]"
                    )

                    # Write .summary files (Opt 8: generate all text first, then
                    # flush each file in one write_text call per siRNA)
                    if global_z_stats:
                        out_dir = output_file.parent
                        out_dir.mkdir(parents=True, exist_ok=True)
                        summary_texts = {
                            sid: prob_service.format_legacy_summary(
                                sid,
                                global_z_stats[sid],
                                global_on_w_stats.get(sid, {}),
                                alpha_gamma_pairs,
                                theta_vals,
                            )
                            for sid in global_z_stats
                        }
                        for sid, text in summary_texts.items():
                            (out_dir / f"{sid}.summary").write_text(text)
                        console.print(
                            f"[green]✓[/green] Wrote {len(global_z_stats)} .summary files to [bold]{out_dir}[/bold]"
                        )
                elif total_rows > 0:
                    console.print(
                        f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions (no output file specified)"
                    )
                elif total_rows == 0:
                    console.print("[yellow]Warning:[/yellow] No predictions to process")

            finally:
                # Ensure Parquet writer is always closed (flushes footer)
                if _pq_writer is not None:
                    _pq_writer.close()

            profiler.print_summary(console)
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
                    with profiler.stage("Load transcriptome") as _s:
                        df_trans = trans_parser.load_gtf(
                            gtf_file,
                            feature=feature_type,
                            score_col=expression_metric,
                            format=transcriptome_format,
                        )
                        _s.rows_out = df_trans.height
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
                        with profiler.stage(f"[chunk] Load batch {batch_idx + 1}") as _s:
                            df_chunk = risearch_parser.load_by_sirna_batch(
                                risearch_file, batch_sirnas
                            )
                            if sense_only:
                                df_chunk = df_chunk.filter(pl.col("strand") == "+")
                            _s.rows_out = df_chunk.height

                        if df_chunk.height == 0:
                            progress.update(
                                task,
                                advance=1,
                                description=f"Batch {batch_idx + 1}: empty",
                            )
                            continue

                        # Intersect with transcriptome - Polars parallelizes this
                        if df_trans is not None:
                            with profiler.stage(f"[chunk] Intersect batch {batch_idx + 1}", rows_in=df_chunk.height) as _s:
                                df_chunk = intersector.intersect(
                                    df_chunk, df_trans, mode=predictions_type, workers=n_workers
                                )
                                _s.rows_out = df_chunk.height

                        # Build on-target map for all siRNAs in this batch
                        batch_map = None
                        if on_target_map:
                            batch_map = {
                                sid: on_target_map.get(sid)
                                for sid in batch_sirnas
                                if sid in on_target_map
                            }

                        # Calculate probabilities - Polars parallelizes across all rows
                        with profiler.stage(f"[chunk] Probabilities batch {batch_idx + 1}", rows_in=df_chunk.height) as _s:
                            df_chunk, batch_metadata = (
                                prob_service.calculate_probabilities_per_sirna(
                                    df_chunk,
                                    alpha_gamma_pairs=alpha_gamma_pairs,
                                    theta_values=theta_vals,
                                    on_target_map=batch_map if batch_map else None,
                                    on_target_expression=on_target_expression,
                                )
                            )
                            _s.rows_out = df_chunk.height

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

                profiler.print_summary(console)
                return  # Exit after chunk mode processing

            # Standard mode: Load entire file
            with profiler.stage("Load predictions") as _s:
                df = risearch_parser.load(risearch_file)
                _s.rows_out = df.height

        logger.info(f"Loaded {df.height:,} predictions for {df['sirna_id'].n_unique()} siRNAs")

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
            with profiler.stage("Load transcriptome") as _s:
                df_trans = trans_parser.load_gtf(
                    gtf_file,
                    feature=feature_type,
                    score_col=expression_metric,
                    format=transcriptome_format,
                )
                _s.rows_out = df_trans.height

            if verbose:
                console.print("[dim]--- Transcriptome Data (First 5 rows) ---[/dim]")
                console.print(df_trans.head(5))

            summary_trans = trans_parser.summary(df_trans)
            console.print(
                f"[green]✓[/green] Loaded [bold]{summary_trans['row_count']}[/bold] features from {gtf_file.name}"
            )
            logger.info(f"Transcriptome: {df_trans.height:,} rows, {df_trans['gene_id'].n_unique()} genes")

            # Perform intersection
            with profiler.stage("Intersection", rows_in=df.height) as _s:
                with console.status(
                    f"[bold green]Intersecting predictions (mode={predictions_type})..."
                ):
                    intersector = IntersectionService()
                    df = intersector.intersect(df, df_trans, mode=predictions_type, workers=n_workers)
                    _s.rows_out = df.height

            console.print(
                f"[green]✓[/green] Found [bold]{df.height}[/bold] intersecting off-target candidates"
            )
            logger.info(f"After intersection: {df.height:,} rows")

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

            with profiler.stage("Probabilities (per-siRNA)", rows_in=df.height) as _s:
                df, meta = prob_service.calculate_probabilities_per_sirna(
                    df,
                    alpha_gamma_pairs=alpha_gamma_pairs,
                    theta_values=theta_vals,
                    on_target_map=on_target_map if on_target_map else None,
                    on_target_expression=on_target_expression,
                )
                _s.rows_out = df.height
            logger.info(f"After probabilities: {df.height:,} rows")

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
            with profiler.stage("Probabilities (single-siRNA)", rows_in=df.height) as _s:
                df, meta = prob_service.calculate_probabilities(
                    df,
                    on_target_path=on_target_file,
                    query_path=query_file,
                    on_target_expression=on_target_expression,
                    on_target_accessibility_path=on_target_accessibility,
                    on_target_risearch_path=on_target_risearch_file,
                )
                _s.rows_out = df.height

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
            logger.info(f"Writing {df.height:,} rows to {output_file}")
            with profiler.stage("Write output", rows_in=df.height):
                df.write_csv(output_file, separator="\t")
            console.print(
                f"\n[green]✓[/green] Results saved to [bold]{output_file}[/bold]"
            )

        profiler.print_summary(console)

    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]Failed:[/bold red] {e}")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1)
