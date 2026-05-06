from pathlib import Path
import tempfile
import typer
import polars as pl
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

# Per-worker state for directory-mode parallel processing.
# Populated by _init_worker() once per spawned process. spawn (not fork)
# avoids Rayon/Polars thread-pool deadlocks in the parent process.
_WORKER_DF_TRANS: pl.DataFrame | None = None
_WORKER_INTERSECTOR: object | None = None
_WORKER_PROB_SERVICE: object | None = None
_WORKER_ON_TARGET_MAP: dict = {}
_WORKER_SELF_HYB_EMIN: dict = {}


def _parse_alpha_gamma_theta(
    alpha: str, gamma: str, theta: str
) -> tuple[list[tuple[float, float]], list[float]]:
    """Parse semicolon-delimited alpha/gamma/theta strings into typed lists.

    Always includes the (1.0, 1.0) baseline pair first.
    Only adds (a, g) pairs where a <= g to preserve clamping semantics.
    """
    alpha_vals = [float(x) for x in alpha.split(";") if x.strip()]
    gamma_vals = [float(x) for x in gamma.split(";") if x.strip()]
    theta_vals = [float(x) for x in theta.split(";") if x.strip()]

    alpha_gamma_pairs: list[tuple[float, float]] = [(1.0, 1.0)]
    for a in alpha_vals:
        for g in gamma_vals:
            if not (a == 1.0 and g == 1.0) and a <= g:
                alpha_gamma_pairs.append((a, g))
    alpha_gamma_pairs = list(dict.fromkeys(alpha_gamma_pairs))

    return alpha_gamma_pairs, theta_vals


def _load_on_target_map(path: Path) -> dict[str, str]:
    """Parse a TSV mapping file (sirna_id \\t on_target_id) into a dict."""
    mapping: dict[str, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
    return mapping


def _init_worker(
    ipc_path: str,
    accessibility_dir: str,
    accessibility_file: str,
    use_rnaplfold_cli: bool,
    polars_max_threads: int,
    on_target_map: dict,
    self_hyb_emin: dict,
) -> None:
    """Initialise per-worker state. Called once per spawned worker process.

    Transcriptome loaded from Arrow IPC (memory-mapped) so all workers share
    the same physical pages via the OS page cache.
    """
    import os
    os.environ["POLARS_MAX_THREADS"] = str(polars_max_threads)

    global _WORKER_DF_TRANS, _WORKER_INTERSECTOR, _WORKER_PROB_SERVICE
    global _WORKER_ON_TARGET_MAP, _WORKER_SELF_HYB_EMIN

    from RIsearch_pipeline.services.intersection_service import IntersectionService
    from RIsearch_pipeline.services.probability import ProbabilityService

    if ipc_path:
        _WORKER_DF_TRANS = pl.read_ipc(Path(ipc_path), memory_map=True)
        _WORKER_INTERSECTOR = IntersectionService()
        _WORKER_INTERSECTOR.preload_transcriptome(_WORKER_DF_TRANS)

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
    _WORKER_ON_TARGET_MAP = on_target_map
    _WORKER_SELF_HYB_EMIN = self_hyb_emin


def _process_one_sirna(
    file_path_str: str,
    alpha_gamma_pairs: list,
    theta_vals: list,
    on_target_expression: float,
    sense_only: bool,
    predictions_type: str,
) -> tuple:
    """Process a single siRNA file through the full pipeline.

    Runs in a spawned worker process. Returns (PyArrow Table, metadata dict)
    or (None, {}) if no predictions remain after filtering.
    """
    import sys
    import time
    _m = sys.modules[__name__]

    df_trans = _m._WORKER_DF_TRANS
    intersector = _m._WORKER_INTERSECTOR
    prob_service = _m._WORKER_PROB_SERVICE
    on_target_map = _m._WORKER_ON_TARGET_MAP
    self_hyb_emin = _m._WORKER_SELF_HYB_EMIN

    from RIsearch_pipeline.services.risearch_parser import RIsearchParser

    _t_start = time.perf_counter()
    parser = RIsearchParser()
    df = parser.load_single_file(Path(file_path_str))
    _t_load = time.perf_counter()

    if sense_only:
        df = df.filter(pl.col("strand") == "+")
    if df.height == 0:
        return None, {}

    if self_hyb_emin and "raw_e_min" in df.columns:
        sirna_id = str(df["sirna_id"][0])
        if sirna_id in self_hyb_emin:
            df = df.with_columns(
                pl.lit(self_hyb_emin[sirna_id]).cast(pl.Float32).alias("raw_e_min")
            )

    _intersect_subtimings: dict = {}
    if df_trans is not None and intersector is not None:
        df = intersector.intersect(df, df_trans, mode=predictions_type, workers=1, _timings=_intersect_subtimings)
        if df.height == 0:
            return None, {}
    _t_intersect = time.perf_counter()

    df, meta = prob_service.calculate_probabilities_per_sirna(
        df,
        alpha_gamma_pairs=alpha_gamma_pairs,
        theta_values=theta_vals,
        on_target_map=on_target_map,
        on_target_expression=on_target_expression,
    )
    _t_prob = time.perf_counter()

    drop_cols = [
        c for c in df.columns
        if c.startswith("boltzmann_weight") or c.startswith("Z_sirna") or c == "E_min"
    ]
    if drop_cols:
        df = df.drop(drop_cols)

    arrow = df.to_arrow()
    _t_end = time.perf_counter()

    # Return freed pages to OS to prevent RSS growth across sequential siRNAs
    # in the same worker. gc + malloc_trim release pages Python's allocator retains.
    import gc
    gc.collect()
    try:
        import ctypes
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass

    meta["_timings"] = {
        "load_s": _t_load - _t_start,
        "intersect_s": _t_intersect - _t_load,
        "prob_s": _t_prob - _t_intersect,
        "serialize_s": _t_end - _t_prob,
        "total_s": _t_end - _t_start,
        **_intersect_subtimings,
    }

    return arrow, meta


def _downcast_schema(df: pl.DataFrame) -> pl.DataFrame:
    """Downcast columns to minimize Arrow IPC file size.

    Coordinates → UInt32, strings → Categorical, floats → Float32.
    """
    casts = {}
    for col in ["start", "end"]:
        if col in df.columns:
            casts[col] = pl.UInt32
    for col in ["sirna_id", "transcript_id", "chrom", "strand", "gene_id"]:
        if col in df.columns:
            casts[col] = pl.Categorical
    for col in ["energy", "opening_energy", "dG_total", "P_off_target", "exp_value"]:
        if col in df.columns:
            casts[col] = pl.Float32
    if casts:
        df = df.cast(casts)
    return df


def run(
    risearch_file: Path | None = typer.Option(
        None, "-r", "--risearch-file",
        help="Path to pre-computed RIsearch output file or directory of files.",
        exists=True, readable=True, dir_okay=True,
    ),
    input_dir: Path | None = typer.Option(
        None, "-d", "--input-dir",
        help="Directory containing per-siRNA RIsearch output files (*.gz or *.tsv).",
        exists=True, file_okay=False,
    ),
    sirna_fasta: Path | None = typer.Option(
        None, "-s", "--sirna-fasta",
        help="Path to siRNA FASTA file. Runs RIsearch internally when used alone; also used in directory mode to compute self-hybridisation E_min.",
        exists=True, readable=True,
    ),
    target_fasta: Path | None = typer.Option(
        None, "--target-fasta", "--genome",
        help="Path to target FASTA (genome or transcriptome) for RIsearch.",
        exists=True, readable=True,
    ),
    target_index: Path | None = typer.Option(
        None, "-idx", "--target-index",
        help="Pre-built RIsearch index (optional, speeds up repeated runs).",
        exists=True, readable=True,
    ),
    workers: int | None = typer.Option(
        None, "-j", "--workers",
        help="Number of parallel threads for intersection (default: CPU count).",
    ),
    gtf_file: Path = typer.Option(
        None, "-t", "--transcriptome",
        help="Path to transcriptome annotation file (.gtf) or .bed file.",
        exists=True, readable=True,
    ),
    feature_type: str = typer.Option("exon", "--feature",
        help="Feature type to select from GTF (default: exon)."),
    expression_metric: str = typer.Option("RPKM", "--expression-metric",
        help="Attribute to use for expression score (default: RPKM)."),
    transcriptome_format: str = typer.Option("auto", "--transcriptome-format",
        help="Transcriptome file format: auto, bed6, bed7, or gtf (default: auto-detect)."),
    accessibility_dir: Path = typer.Option(
        None, "-a", "--accessibility-dir",
        help="Directory containing pre-computed accessibility profiles.",
        exists=True, file_okay=False,
    ),
    accessibility_file: Path = typer.Option(
        None, "--accessibility-file",
        help="Parquet file with pre-computed binding-site accessibility.",
        exists=True, dir_okay=False,
    ),
    output_file: Path = typer.Option(None, "-o", "--output",
        help="Path to save the analysis results (TSV)."),
    genome_file: Path = typer.Option(
        None, "-f", "--fasta",
        help="Path to genome FASTA file (computes accessibility on-the-fly).",
        exists=True, readable=True,
    ),
    window_size: int = typer.Option(80, "--window", "-W", help="Window size (W)"),
    max_span: int = typer.Option(40, "--span", "-L", help="Max base pair span (L)"),
    unpaired_prob: int = typer.Option(30, "--unpaired", "-u", help="Unpaired probability length (u)"),
    on_target_file: Path = typer.Option(
        None, "-on", "--on-target",
        help="Path to On-Target sequence FASTA (for Partition Function).",
        exists=True, readable=True,
    ),
    on_target_risearch_file: Path | None = typer.Option(
        None, "--on-target-risearch-file", "-on-ris",
        help="Path to pre-computed RIsearch output file for On-Target.",
        exists=True, readable=True,
    ),
    query_file: Path = typer.Option(
        None, "-q", "--query",
        help="Path to siRNA query FASTA (required if --on-target is used).",
        exists=True, readable=True,
    ),
    on_target_expression: float = typer.Option(
        1000.0, "--on-target-expression", "-oexp",
        help="Expression level for On-Target (default: 1000.0)."),
    on_target_accessibility: Path | None = typer.Option(
        None, "--on-target-accessibility",
        help="Path to accessibility file for On-Target (text or binary)."),
    on_target_ids_file: Path | None = typer.Option(
        None, "-oi", "--on-target-ids",
        help="TSV file mapping siRNA IDs to on-target transcript IDs (sirna_id \\t transcript_id).",
        exists=True, readable=True,
    ),
    alpha: str = typer.Option("1.0", "--alpha",
        help="Alpha clamping parameter(s). Separate multiple values with ';' (e.g., '0.8;1.0')."),
    gamma: str = typer.Option("1.0", "--gamma",
        help="Gamma clamping parameter(s). Separate multiple values with ';' (e.g., '0.8;1.0')."),
    theta: str = typer.Option("", "--theta",
        help="Theta scaling parameter(s). Separate multiple values with ';' (e.g., '0.5;0.7')."),
    legacy_format: bool = typer.Option(False, "--legacy-format",
        help="Output in legacy format (gw.results style, aggregated by transcript)."),
    detailed_report: bool = typer.Option(False, "--detailed-report",
        help="Report off-target probabilities for individual transcripts (Legacy Format)."),
    sense_only: bool = typer.Option(False, "--sense-only",
        help="Limit RIsearch2 predictions to sense strand (+), ignore antisense."),
    predictions_type: str = typer.Option("gw", "--type",
        help="Type of RIsearch2 predictions. gw for genome-wide and tw for transcriptome-wide."),
    verbose: bool = typer.Option(False, "--verbose", "-v",
        help="Show detailed output (tables, stats)."),
    profile: bool = typer.Option(False, "--profile",
        help="Print a per-stage timing and memory profile table at the end of the run."),
    use_rnaplfold_cli: bool = typer.Option(False, "--use-rnaplfold-cli",
        help="Use RNAplfold binary instead of ViennaRNA Python bindings for accessibility."),
    chunk_mode: bool = typer.Option(False, "--chunk-mode",
        help="Process siRNAs in batches for large files (reduces memory usage)."),
    batch_size: int = typer.Option(50, "--batch-size", "-b",
        help="Number of siRNAs to process per batch in chunk mode (default: 50)."),
    scratch_dir: Path | None = typer.Option(None, "--scratch-dir",
        help="Directory for intermediate Arrow IPC files. Defaults to system temp."),
    output_format: str = typer.Option("tsv", "--output-format",
        help="Final output format: 'tsv', 'csv', or 'parquet'."),
    summary_only: bool = typer.Option(False, "--summary-only",
        help="Write only .summary files per siRNA; skip the per-prediction output file."),
) -> None:
    """Analyze siRNA off-target predictions.

    Integrates RIsearch2 predictions with transcriptome data and optionally
    calculates off-target probabilities.
    """
    import os
    risearch_parser = RIsearchParser()
    console.print(Panel("RIsearch Pipeline", style="bold cyan"))
    profiler = PipelineProfiler(enabled=profile)
    n_workers: int = workers if workers is not None else os.cpu_count() or 1

    try:
        if risearch_file is not None and risearch_file.is_dir():
            input_dir = risearch_file
            risearch_file = None

        if risearch_file is None and sirna_fasta is None and input_dir is None:
            console.print("[bold red]Error:[/bold red] Must provide either --risearch-file, --input-dir, OR --sirna-fasta")
            raise typer.Exit(code=1)

        is_running_risearch = sirna_fasta is not None and risearch_file is None and input_dir is None
        if is_running_risearch and target_fasta is None:
            console.print("[bold red]Error:[/bold red] --sirna-fasta requires --target-fasta when running RIsearch dynamically")
            raise typer.Exit(code=1)

        if is_running_risearch:
            from RIsearch_pipeline.services.risearch_service import RIsearchService
            risearch_service = RIsearchService()

            try:
                sirna_ids = risearch_service.validate_sirna_fasta(sirna_fasta)
                console.print(f"[green]✓[/green] Validated [bold]{len(sirna_ids)}[/bold] siRNA(s) from {sirna_fasta.name}")
            except ValueError as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
                raise typer.Exit(code=1)

            with console.status("[bold green]Indexing target..."):
                if target_index is not None:
                    index_path = target_index
                    console.print(f"[green]✓[/green] Using pre-built index: {index_path.name}")
                else:
                    index_path = risearch_service.index_target(target_fasta)
                    console.print(f"[green]✓[/green] Created index: {index_path.name}")

            with console.status("[bold green]Running RIsearch..."):
                df = risearch_service.run_search(
                    query_path=sirna_fasta, index_path=index_path, target_fasta=target_fasta,
                )

        elif input_dir is not None:
            console.print(f"[bold cyan]Directory mode[/bold cyan] - loading files from {input_dir}")

            all_files = risearch_parser.list_directory_files(input_dir)
            if not all_files:
                console.print(f"[bold red]Error:[/bold red] No RIsearch files found in {input_dir}")
                raise typer.Exit(code=1)

            console.print(f"[green]✓[/green] Found [bold]{len(all_files)}[/bold] files to process")

            from RIsearch_pipeline.services.probability import ProbabilityService
            from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService

            acc_service = None
            precomputed_acc = None
            if accessibility_file:
                precomputed_acc = pl.read_parquet(accessibility_file)
                console.print(f"  [dim]Loaded {precomputed_acc.height} precomputed accessibility values[/dim]")
            elif accessibility_dir:
                acc_service = GenomeAccessibilityService(accessibility_dir, max_cached=4)
            prob_service = ProbabilityService(acc_service, use_rnaplfold_cli=use_rnaplfold_cli,
                                              precomputed_accessibility=precomputed_acc)

            df_trans = None
            _ipc_path = None
            _ipc_tmp = None
            if gtf_file:
                from RIsearch_pipeline.services.transcriptome_parser import TranscriptomeParser
                from RIsearch_pipeline.services.intersection_service import IntersectionService

                trans_parser = TranscriptomeParser()
                with profiler.stage("Load transcriptome") as _s:
                    df_trans = trans_parser.load_gtf(gtf_file, feature=feature_type,
                                                     score_col=expression_metric, format=transcriptome_format)
                    _s.rows_out = df_trans.height
                console.print(f"[green]✓[/green] Loaded transcriptome with {df_trans.height} features")

                _ipc_tmp = tempfile.mkdtemp(prefix="risearch_ipc_")
                _ipc_path = Path(_ipc_tmp) / "transcriptome.arrow"
                df_trans.write_ipc(_ipc_path)

            on_target_map: dict[str, str] = {}
            if on_target_ids_file is not None:
                on_target_map = _load_on_target_map(on_target_ids_file)
                console.print(f"  [dim]Loaded {len(on_target_map)} on-target mappings[/dim]")

            self_hyb_emin: dict[str, float] = {}
            if sirna_fasta is not None:
                from RIsearch_pipeline.services.risearch_service import RIsearchService
                console.print(f"  └─ Computing self-hybridisation E_min from [bold]{sirna_fasta.name}[/bold]...")
                self_hyb_emin = RIsearchService().self_hybridization_emin_batch(sirna_fasta)
                console.print(f"  [dim]Self-hyb E_min computed for {len(self_hyb_emin)} siRNA(s)[/dim]")

            alpha_gamma_pairs, theta_vals = _parse_alpha_gamma_theta(alpha, gamma, theta)

            total_rows = 0
            n_proc = min(n_workers, len(all_files))
            console.print(f"  └─ Processing {len(all_files)} siRNAs with up to {n_proc} parallel workers")

            out_path = None
            _csv_first_batch = True
            _pq_writer = None
            if output_file and not summary_only:
                if output_format == "parquet":
                    out_path = output_file.with_suffix(".parquet")
                elif output_format == "tsv":
                    out_path = output_file.with_suffix(".tsv")
                else:
                    out_path = output_file.with_suffix(".csv")
                out_path.parent.mkdir(parents=True, exist_ok=True)

            summary_out_dir = None
            if output_file is not None:
                if output_file.suffix == "" or str(output_file).endswith("/") or output_file.is_dir():
                    summary_out_dir = output_file
                else:
                    summary_out_dir = output_file.parent
                summary_out_dir.mkdir(parents=True, exist_ok=True)
            summary_count = 0

            try:
                from concurrent.futures import ProcessPoolExecutor, as_completed
                import multiprocessing as _mp

                polars_threads_per_worker = max(1, n_workers // n_proc)
                _ctx = _mp.get_context("spawn")

                _init_args = (
                    str(_ipc_path) if gtf_file else "",
                    str(accessibility_dir) if accessibility_dir else "",
                    str(accessibility_file) if accessibility_file else "",
                    use_rnaplfold_cli,
                    polars_threads_per_worker,
                    on_target_map,
                    self_hyb_emin,
                )

                with Progress(
                    SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                    BarColumn(), TextColumn("[progress.percentage]{task.completed}/{task.total}"),
                    TimeElapsedColumn(), console=console,
                ) as progress:
                    task = progress.add_task("Processing siRNAs...", total=len(all_files))

                    with ProcessPoolExecutor(
                        max_workers=n_proc, mp_context=_ctx,
                        initializer=_init_worker, initargs=_init_args,
                    ) as pool:
                        futures = {
                            pool.submit(
                                _process_one_sirna, str(f), alpha_gamma_pairs,
                                theta_vals, on_target_expression, sense_only, predictions_type,
                            ): f
                            for f in all_files
                        }

                        import time as _time
                        _drain_stage_totals: dict[str, float] = {}
                        _drain_n = 0
                        _drain_wait_total = 0.0
                        completed = 0

                        for future in as_completed(futures):
                            _t_drain_start = _time.perf_counter()
                            completed += 1
                            progress.update(task, advance=1, description=f"Batch {completed}/{len(all_files)}")

                            # Pop immediately so Future._result (arrow_table + metadata)
                            # can be GC'd as soon as the loop variable is overwritten —
                            # without this, all N results accumulate in memory.
                            f = futures.pop(future)
                            try:
                                _t_wait = _time.perf_counter()
                                arrow_table, batch_metadata = future.result()
                                _drain_wait_total += _time.perf_counter() - _t_wait
                            except Exception as exc:
                                logger.error(f"Worker failed for {f.name}: {exc}")
                                continue

                            if arrow_table is None:
                                continue

                            _t_deserialize = _time.perf_counter()
                            df_chunk = pl.from_arrow(arrow_table)
                            del arrow_table
                            total_rows += df_chunk.height
                            _dt_deserialize = _time.perf_counter() - _t_deserialize

                            _t_summary = _time.perf_counter()
                            if summary_out_dir is not None:
                                z_per_sirna = batch_metadata.get("z_per_sirna", {})
                                w_per_sirna = batch_metadata.get("on_target_weights", {})
                                for sid, z_dict in z_per_sirna.items():
                                    text = prob_service.format_legacy_summary(
                                        sid, z_dict, w_per_sirna.get(sid, {}),
                                        alpha_gamma_pairs, theta_vals,
                                    )
                                    (summary_out_dir / f"{sid}.summary").write_text(text)
                                    summary_count += 1
                            _dt_summary = _time.perf_counter() - _t_summary

                            _t_write = _time.perf_counter()
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
                            _dt_write = _time.perf_counter() - _t_write

                            del df_chunk

                            worker_timings = batch_metadata.get("_timings", {})
                            for k, v in worker_timings.items():
                                _drain_stage_totals[k] = _drain_stage_totals.get(k, 0.0) + v
                            _drain_stage_totals["drain_deserialize_s"] = _drain_stage_totals.get("drain_deserialize_s", 0.0) + _dt_deserialize
                            _drain_stage_totals["drain_summary_s"] = _drain_stage_totals.get("drain_summary_s", 0.0) + _dt_summary
                            _drain_stage_totals["drain_write_s"] = _drain_stage_totals.get("drain_write_s", 0.0) + _dt_write
                            _drain_n += 1

                    if _drain_n > 0:
                        perf_lines = [f"[perf] per-siRNA stage breakdown (N={_drain_n}, times are sums over all siRNAs):"]
                        for k, v in sorted(_drain_stage_totals.items()):
                            perf_lines.append(f"  {k}: {v:.3f}s  (avg {v/_drain_n*1000:.1f}ms/siRNA)")
                        perf_lines.append(f"  drain_wait_for_future_s: {_drain_wait_total:.3f}s  (avg {_drain_wait_total/_drain_n*1000:.1f}ms/siRNA)")
                        console.print("\n".join(perf_lines))

                if total_rows > 0 and out_path is not None:
                    console.print(f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions → [bold]{out_path}[/bold]")
                elif total_rows > 0:
                    console.print(f"[green]✓[/green] Processed [bold]{total_rows}[/bold] predictions (no output file specified)")
                else:
                    console.print("[yellow]Warning:[/yellow] No predictions to process")

                if summary_count > 0:
                    console.print(f"[green]✓[/green] Wrote {summary_count} .summary files to [bold]{summary_out_dir}[/bold]")

            finally:
                if _pq_writer is not None:
                    _pq_writer.close()
                import shutil
                if _ipc_tmp:
                    shutil.rmtree(_ipc_tmp, ignore_errors=True)

            profiler.print_summary(console)
            return

        elif risearch_file is not None:
            if chunk_mode:
                console.print("[bold cyan]Chunk mode enabled[/bold cyan] - processing siRNAs individually")

                with console.status("[bold green]Scanning for siRNA IDs..."):
                    sirna_ids = risearch_parser.scan_sirna_ids(risearch_file)
                console.print(f"[green]✓[/green] Found [bold]{len(sirna_ids)}[/bold] unique siRNAs to process")

                from RIsearch_pipeline.services.probability import ProbabilityService
                from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService

                acc_service = None
                precomputed_acc = None
                if accessibility_file:
                    precomputed_acc = pl.read_parquet(accessibility_file)
                elif accessibility_dir:
                    acc_service = GenomeAccessibilityService(accessibility_dir, max_cached=4)
                prob_service = ProbabilityService(acc_service, use_rnaplfold_cli=use_rnaplfold_cli,
                                                  precomputed_accessibility=precomputed_acc)

                df_trans = None
                intersector = None
                if gtf_file:
                    from RIsearch_pipeline.services.transcriptome_parser import TranscriptomeParser
                    from RIsearch_pipeline.services.intersection_service import IntersectionService

                    trans_parser = TranscriptomeParser()
                    with profiler.stage("Load transcriptome") as _s:
                        df_trans = trans_parser.load_gtf(gtf_file, feature=feature_type,
                                                         score_col=expression_metric, format=transcriptome_format)
                        _s.rows_out = df_trans.height
                    intersector = IntersectionService()
                    console.print(f"[green]✓[/green] Loaded transcriptome with {df_trans.height} features")

                on_target_map: dict[str, str] = {}
                if on_target_ids_file is not None:
                    on_target_map = _load_on_target_map(on_target_ids_file)

                alpha_gamma_pairs, theta_vals = _parse_alpha_gamma_theta(alpha, gamma, theta)

                global_z_stats: dict = {}
                global_on_w_stats: dict = {}
                all_results = []
                num_batches = (len(sirna_ids) + batch_size - 1) // batch_size

                console.print(f"  └─ Processing {len(sirna_ids)} siRNAs in {num_batches} batches of up to {batch_size}")

                with Progress(
                    SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                    BarColumn(), TextColumn("[progress.percentage]{task.completed}/{task.total}"),
                    TimeElapsedColumn(), console=console,
                ) as progress:
                    task = progress.add_task("Processing batches...", total=num_batches)

                    for batch_idx in range(num_batches):
                        start_idx = batch_idx * batch_size
                        batch_sirnas = sirna_ids[start_idx:start_idx + batch_size]

                        with profiler.stage(f"[chunk] Load batch {batch_idx + 1}") as _s:
                            df_chunk = risearch_parser.load_by_sirna_batch(risearch_file, batch_sirnas)
                            if sense_only:
                                df_chunk = df_chunk.filter(pl.col("strand") == "+")
                            _s.rows_out = df_chunk.height

                        if df_chunk.height == 0:
                            progress.update(task, advance=1, description=f"Batch {batch_idx + 1}: empty")
                            continue

                        if df_trans is not None:
                            with profiler.stage(f"[chunk] Intersect batch {batch_idx + 1}", rows_in=df_chunk.height) as _s:
                                df_chunk = intersector.intersect(df_chunk, df_trans, mode=predictions_type, workers=n_workers)
                                _s.rows_out = df_chunk.height

                        batch_map = {sid: on_target_map[sid] for sid in batch_sirnas if sid in on_target_map} or None

                        with profiler.stage(f"[chunk] Probabilities batch {batch_idx + 1}", rows_in=df_chunk.height) as _s:
                            df_chunk, batch_metadata = prob_service.calculate_probabilities_per_sirna(
                                df_chunk, alpha_gamma_pairs=alpha_gamma_pairs, theta_values=theta_vals,
                                on_target_map=batch_map, on_target_expression=on_target_expression,
                            )
                            _s.rows_out = df_chunk.height

                        for sid, z_dict in batch_metadata["z_per_sirna"].items():
                            global_z_stats[sid] = z_dict
                        for sid, w_dict in batch_metadata["on_target_weights"].items():
                            global_on_w_stats[sid] = w_dict

                        all_results.append(df_chunk)
                        progress.update(task, advance=1, description=f"Batch {batch_idx + 1}/{num_batches}")

                if not all_results:
                    console.print("[yellow]Warning:[/yellow] No predictions to process")
                    return

                df = pl.concat(all_results, how="diagonal")
                console.print(f"[green]✓[/green] Processed [bold]{df.height}[/bold] total predictions")

                if output_file:
                    df.write_csv(output_file, separator="\t")
                    console.print(f"[green]✓[/green] Results saved to [bold]{output_file}[/bold]")

                    out_dir = output_file.parent
                    out_dir.mkdir(parents=True, exist_ok=True)
                    for sid in global_z_stats:
                        summary_text = prob_service.format_legacy_summary(
                            sid, global_z_stats[sid], global_on_w_stats.get(sid, {}),
                            alpha_gamma_pairs, theta_vals,
                        )
                        (out_dir / f"{sid}.summary").write_text(summary_text)
                    console.print(f"[green]✓[/green] Wrote {len(global_z_stats)} .summary files to [bold]{out_dir}[/bold]")
                else:
                    console.print(df.head(10))

                profiler.print_summary(console)
                return

            with profiler.stage("Load predictions") as _s:
                df = risearch_parser.load(risearch_file)
                _s.rows_out = df.height

        logger.info(f"Loaded {df.height:,} predictions for {df['sirna_id'].n_unique()} siRNAs")

        if sense_only:
            df = df.filter(pl.col("strand") == "+")

        if verbose:
            console.print("[dim]--- RIsearch Predictions (First 5 rows) ---[/dim]")
            console.print(df.head(5))

        summary_ris = risearch_parser.summary(df)
        source_name = sirna_fasta.name if sirna_fasta else risearch_file.name
        console.print(f"[green]✓[/green] Loaded [bold]{summary_ris['row_count']}[/bold] predictions from {source_name}")
        console.print(f"  └─ Energy range: {summary_ris['energy_min']:.2f} to {summary_ris['energy_max']:.2f} kcal/mol")
        if sense_only:
            console.print("  └─ (Filtered to sense strand only)")

        if gtf_file:
            from RIsearch_pipeline.services.transcriptome_parser import TranscriptomeParser
            from RIsearch_pipeline.services.intersection_service import IntersectionService

            trans_parser = TranscriptomeParser()
            with profiler.stage("Load transcriptome") as _s:
                df_trans = trans_parser.load_gtf(gtf_file, feature=feature_type,
                                                 score_col=expression_metric, format=transcriptome_format)
                _s.rows_out = df_trans.height

            if verbose:
                console.print("[dim]--- Transcriptome Data (First 5 rows) ---[/dim]")
                console.print(df_trans.head(5))

            summary_trans = trans_parser.summary(df_trans)
            console.print(f"[green]✓[/green] Loaded [bold]{summary_trans['row_count']}[/bold] features from {gtf_file.name}")
            logger.info(f"Transcriptome: {df_trans.height:,} rows, {df_trans['gene_id'].n_unique()} genes")

            with profiler.stage("Intersection", rows_in=df.height) as _s:
                with console.status(f"[bold green]Intersecting predictions (mode={predictions_type})..."):
                    df = IntersectionService().intersect(df, df_trans, mode=predictions_type, workers=n_workers)
                    _s.rows_out = df.height

            console.print(f"[green]✓[/green] Found [bold]{df.height}[/bold] intersecting off-target candidates")
            logger.info(f"After intersection: {df.height:,} rows")

        from RIsearch_pipeline.services.probability import ProbabilityService
        from RIsearch_pipeline.services.accessibility import GenomeAccessibilityService

        temp_dir_obj = None

        if accessibility_file:
            precomputed_acc = pl.read_parquet(accessibility_file)
            console.print(f"  [dim]Loaded {precomputed_acc.height} precomputed accessibility values from {accessibility_file}...[/dim]")
            prob_service = ProbabilityService(None, use_rnaplfold_cli=use_rnaplfold_cli, precomputed_accessibility=precomputed_acc)

        elif accessibility_dir:
            console.print(f"  [dim]Calculating probabilities using profiles from {accessibility_dir}...[/dim]")
            prob_service = ProbabilityService(
                GenomeAccessibilityService(accessibility_dir, max_cached=4),
                use_rnaplfold_cli=use_rnaplfold_cli,
            )

        elif genome_file:
            console.print(f"  [dim]Computing accessibility on-the-fly from {genome_file} (this may take time)...[/dim]")
            temp_dir_obj = tempfile.TemporaryDirectory(prefix="risearch_accessibility_")
            acc_service = GenomeAccessibilityService(Path(temp_dir_obj.name), max_cached=4)

            with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                          BarColumn(), TimeElapsedColumn(), console=console) as progress:
                task = progress.add_task("Calculating Accessibility...", total=None)

                def progress_callback(advance=0, description=""):
                    progress.update(task, advance=advance, description=description)

                acc_service.compute_genome_accessibility(
                    genome_file, window_size=window_size, max_span=max_span,
                    unpaired_prob=unpaired_prob, progress_callback=progress_callback,
                )

            prob_service = ProbabilityService(acc_service, use_rnaplfold_cli=use_rnaplfold_cli)

        else:
            console.print("[yellow]Warning:[/yellow] No accessibility data provided. P(OT) based on energy only.")
            prob_service = ProbabilityService(None, use_rnaplfold_cli=use_rnaplfold_cli)

        console.print("[bold green]Calculating probabilities...[/bold green]")
        if accessibility_dir or genome_file:
            console.print("  └─ Annotating opening energies from accessibility profiles...")
        else:
            console.print("  └─ [yellow]Skipping accessibility annotation (energy only)[/yellow]")

        alpha_gamma_pairs, theta_vals = _parse_alpha_gamma_theta(alpha, gamma, theta)

        unique_sirnas = df["sirna_id"].unique() if "sirna_id" in df.columns else []
        is_multi_sirna = len(unique_sirnas) > 1 or len(alpha_gamma_pairs) > 1 or len(theta_vals) > 0

        on_target_map: dict[str, str] = {}
        if on_target_ids_file is not None:
            try:
                mapping_df = pl.read_csv(on_target_ids_file, separator="\t", has_header=False)
                on_target_map = dict(zip(mapping_df["column_1"], mapping_df["column_2"]))
                console.print(f"[green]✓[/green] Loaded {len(on_target_map)} on-target ID mappings from {on_target_ids_file.name}")
            except Exception as e:
                console.print(f"[bold red]Error parsing on-target mapping file:[/bold red] {e}")
                raise typer.Exit(code=1)

        if is_multi_sirna:
            mode_msg = "Multi-siRNA" if len(unique_sirnas) > 1 else "Parameterized Single-siRNA"
            console.print(f"  └─ {mode_msg} detection: {len(unique_sirnas)} siRNAs")

            with profiler.stage("Probabilities (per-siRNA)", rows_in=df.height) as _s:
                df, meta = prob_service.calculate_probabilities_per_sirna(
                    df, alpha_gamma_pairs=alpha_gamma_pairs, theta_values=theta_vals,
                    on_target_map=on_target_map if on_target_map else None,
                    on_target_expression=on_target_expression,
                )
                _s.rows_out = df.height
            logger.info(f"After probabilities: {df.height:,} rows")

            if output_file:
                out_dir = output_file.parent
                out_dir.mkdir(parents=True, exist_ok=True)
                global_z_stats = meta["z_per_sirna"]
                global_on_w_stats = meta["on_target_weights"]
                for sid in global_z_stats:
                    summary_text = prob_service.format_legacy_summary(
                        sid, global_z_stats[sid], global_on_w_stats.get(sid, {}),
                        alpha_gamma_pairs, theta_vals,
                    )
                    (out_dir / f"{sid}.summary").write_text(summary_text)

            console.print("  └─ Computing per-siRNA Boltzmann weights...")
            if len(alpha_gamma_pairs) > 1 or theta_vals:
                console.print(f"  └─ Applied {len(alpha_gamma_pairs) - 1} alpha/gamma pairs and {len(theta_vals)} theta values")

            total_z = sum(z_dict.get("Z_sirna", 0.0) for z_dict in meta["z_per_sirna"].values())
            console.print(f"  └─ Total Z across all siRNAs (base): [bold]{total_z:.2e}[/bold]")
        else:
            with profiler.stage("Probabilities (single-siRNA)", rows_in=df.height) as _s:
                df, meta = prob_service.calculate_probabilities(
                    df, on_target_path=on_target_file, query_path=query_file,
                    on_target_expression=on_target_expression,
                    on_target_accessibility_path=on_target_accessibility,
                    on_target_risearch_path=on_target_risearch_file,
                )
                _s.rows_out = df.height

            console.print("  └─ Computing Boltzmann weights (α=1.0, γ=1.0)...")
            z_fmt = f"{meta['z_total']:.2e}"
            console.print(
                f"  └─ Partition Function Z = [bold]{z_fmt}[/bold] "
                f"(Off-Target={meta['z_off_target']:.2e}, On-Target={meta['w_on_target']:.2e})"
            )

        if meta.get("has_on_target", False):
            console.print(
                f"  └─ On-Target: ΔG_total={meta['dG_on_target']:.2f} kcal/mol, "
                f"P(on)={meta.get('p_on_target', 0.0):.10f}"
            )

        if legacy_format:
            sirna_id = query_file.stem if query_file else "siRNA"
            legacy_output = prob_service.calculate_legacy_format(
                df, sirna_id=sirna_id, on_target_path=on_target_file, query_path=query_file,
                on_target_expression=on_target_expression,
                on_target_accessibility_path=on_target_accessibility,
                on_target_risearch_path=on_target_risearch_file, verbose=detailed_report,
            )

            if output_file:
                legacy_path = output_file.with_suffix(".results")
                with open(legacy_path, "w") as f:
                    f.write(legacy_output)
                console.print(f"[green]✓[/green] Legacy format results saved to {legacy_path}")
                detailed_path = output_file.parent / "detailed_results.tsv"
                df.write_csv(detailed_path, separator="\t")
                console.print(f"[green]✓[/green] Detailed results saved to {detailed_path}")
            else:
                console.print(legacy_output)
            return

        if df.height > 0:
            if "P_off_target" in df.columns:
                df = df.sort("P_off_target", descending=True)
                title = "Top 10 Candidates (By P_off_target)"
            else:
                title = "Predictions Preview"

            table = Table(title=title)
            potential_cols = ["sirna_id", "chrom", "gene_id", "transcript_id",
                              "energy", "dG_total", "P_off_target", "exp_value"]
            display_cols = [c for c in potential_cols if c in df.columns]
            for col in display_cols:
                table.add_column(col)
            for row in df.select(display_cols).head(10).iter_rows():
                table.add_row(*[f"{v:.4g}" if isinstance(v, float) else str(v) for v in row])
            console.print(table)

        if output_file:
            logger.info(f"Writing {df.height:,} rows to {output_file}")
            with profiler.stage("Write output", rows_in=df.height):
                df.write_csv(output_file, separator="\t")
            console.print(f"\n[green]✓[/green] Results saved to [bold]{output_file}[/bold]")

        profiler.print_summary(console)

    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]Failed:[/bold red] {e}")
        if verbose:
            console.print_exception()
        raise typer.Exit(code=1)
