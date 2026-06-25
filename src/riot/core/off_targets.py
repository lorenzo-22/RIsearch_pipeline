"""Pure core for off-target analysis.

Two entry points, both free of CLI concerns (no stdout, no Typer/Click
exceptions, no file writes):

- :func:`compute_off_targets_single` — single-file or inline-RIsearch input;
  returns ``(pl.DataFrame, meta)``.
- :func:`compute_off_targets_directory` — a directory of per-siRNA prediction
  files; a generator yielding ``(pl.DataFrame, meta)`` per siRNA. The pool and
  the temp Arrow-IPC transcriptome live inside the generator's ``try/finally``,
  so callers should consume it fully (or close it) for deterministic cleanup.

``meta`` carries the partition-function state downstream formatters need:
``z_per_sirna``/``on_target_weights`` (multi-siRNA), or ``z_total`` etc.
(single-siRNA). When ``legacy_format`` is requested the single core also stashes
the rendered legacy report under ``meta["legacy_text"]`` (it reuses the already
built accessibility service, avoiding a second on-the-fly fold).
"""

import ctypes
import gc
import multiprocessing
import os
import shutil
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Generator, Optional

import polars as pl
from loguru import logger

from riot.services.accessibility import GenomeAccessibilityService
from riot.services.annotation_parser import AnnotationParser
from riot.services.intersection_service import IntersectionService
from riot.services.probability import ProbabilityService
from riot.services.profiling import PipelineProfiler
from riot.services.risearch_parser import RIsearchParser
from riot.services.risearch_service import RIsearchService

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
_WORKER_ON_TARGET_MAP: dict = {}
_WORKER_SELF_HYB_EMIN: dict = {}


def _init_worker(
    ipc_path: str,
    accessibility_dir: str,
    polars_max_threads: int,
    on_target_map: dict,
    self_hyb_emin: dict,
    temperature: float = 37.0,
) -> None:
    """Initialise per-worker state.  Called once per spawned worker process.

    The transcriptome is loaded from a pre-written Arrow IPC file (memory-mapped)
    rather than being re-parsed from the original BED/GTF.  All workers share the
    same physical memory pages via the OS page cache.
    """
    os.environ["POLARS_MAX_THREADS"] = str(polars_max_threads)

    global _WORKER_DF_TRANS, _WORKER_INTERSECTOR, _WORKER_PROB_SERVICE
    global _WORKER_ON_TARGET_MAP, _WORKER_SELF_HYB_EMIN

    if ipc_path:
        _WORKER_DF_TRANS = pl.read_ipc(Path(ipc_path), memory_map=True)
        _WORKER_INTERSECTOR = IntersectionService()
        _WORKER_INTERSECTOR.preload_transcriptome(_WORKER_DF_TRANS)

    acc_service = None
    if accessibility_dir:
        acc_service = GenomeAccessibilityService(Path(accessibility_dir), max_cached=4)

    _WORKER_PROB_SERVICE = ProbabilityService(acc_service, temperature=temperature)

    # Store per-run constants so they don't need to be pickled per task.
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

    Runs in a spawned worker process.  Reads worker state from module-level
    globals populated by _init_worker().

    Returns:
        (PyArrow Table, metadata dict) or (None, {}) if no predictions remain.
    """
    _m = sys.modules[__name__]

    df_trans = _m._WORKER_DF_TRANS
    intersector = _m._WORKER_INTERSECTOR
    prob_service = _m._WORKER_PROB_SERVICE
    on_target_map = _m._WORKER_ON_TARGET_MAP
    self_hyb_emin = _m._WORKER_SELF_HYB_EMIN

    _t_start = time.perf_counter()

    parser = RIsearchParser()
    df = parser.load_single_file(Path(file_path_str))
    _t_load = time.perf_counter()

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

    # Intersect with transcriptome (sequential — each worker owns its CPU).
    # intersect() already deduplicates per (chrom, strand) pair internally.
    _intersect_subtimings: dict = {}
    if df_trans is not None and intersector is not None:
        df = intersector.intersect(
            df,
            df_trans,
            mode=predictions_type,
            workers=1,
            _timings=_intersect_subtimings,
        )
        if df.height == 0:
            return None, {}
    _t_intersect = time.perf_counter()

    # Probabilities
    df, meta = prob_service.calculate_probabilities_per_sirna(
        df,
        alpha_gamma_pairs=alpha_gamma_pairs,
        theta_values=theta_vals,
        on_target_map=on_target_map,
        on_target_expression=on_target_expression,
    )
    _t_prob = time.perf_counter()

    # Drop heavy intermediate columns before returning
    drop_cols = [
        c
        for c in df.columns
        if c.startswith("boltzmann_weight") or c.startswith("Z_sirna") or c == "E_min"
    ]
    if drop_cols:
        df = df.drop(drop_cols)

    arrow = df.to_arrow()
    _t_end = time.perf_counter()

    # Return freed pages to OS to prevent RSS growth across sequential siRNAs
    # in the same worker. Python's allocator retains pages by default; gc +
    # malloc_trim release them, keeping per-siRNA memory cost constant.
    gc.collect()
    try:
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


def _build_alpha_gamma_pairs(alpha: str, gamma: str) -> list[tuple[float, float]]:
    """Parse semicolon-separated alpha/gamma strings into (alpha, gamma) pairs.

    Always includes the baseline (1.0, 1.0). Enforces alpha ≤ gamma. Deduplicates.
    """
    alpha_vals = [float(x) for x in alpha.split(";") if x.strip()]
    gamma_vals = [float(x) for x in gamma.split(";") if x.strip()]
    pairs: list[tuple[float, float]] = [(1.0, 1.0)]
    for a in alpha_vals:
        for g in gamma_vals:
            if a == 1.0 and g == 1.0:
                continue
            if a <= g:
                pairs.append((a, g))
    return list(dict.fromkeys(pairs))


def _parse_theta(theta: str) -> list[float]:
    return [float(x) for x in theta.split(";") if x.strip()]


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


# ---------------------------------------------------------------------------
# Single-file / inline-RIsearch core
# ---------------------------------------------------------------------------
def compute_off_targets_single(
    *,
    risearch_file: Optional[Path] = None,
    sirna_fasta: Optional[Path] = None,
    target_fasta: Optional[Path] = None,
    target_index: Optional[Path] = None,
    gtf_file: Optional[Path] = None,
    feature_type: str = "exon",
    expression_metric: str = "RPKM",
    transcriptome_format: str = "auto",
    accessibility_dir: Optional[Path] = None,
    genome_file: Optional[Path] = None,
    window_size: int = 80,
    max_span: int = 40,
    unpaired_prob: int = 30,
    temperature: float = 37.0,
    on_target_file: Optional[Path] = None,
    on_target_risearch_file: Optional[Path] = None,
    query_file: Optional[Path] = None,
    on_target_expression: float = 1000.0,
    on_target_accessibility: Optional[Path] = None,
    on_target_ids_file: Optional[Path] = None,
    alpha: str = "1.0",
    gamma: str = "1.0",
    theta: str = "",
    sense_only: bool = False,
    predictions_type: str = "gw",
    legacy_format: bool = False,
    detailed_report: bool = False,
    n_workers: int = 1,
    profiler: Optional[PipelineProfiler] = None,
    accessibility_progress_callback=None,
) -> tuple[pl.DataFrame, dict]:
    """Compute off-target predictions for a single predictions file (or inline RIsearch).

    Returns ``(df, meta)``. Writes no files and prints nothing. Raises ``ValueError``
    on bad inputs.
    """
    profiler = profiler if profiler is not None else PipelineProfiler(enabled=False)
    risearch_parser = RIsearchParser()

    is_running_risearch = (sirna_fasta is not None) and (risearch_file is None)
    if risearch_file is None and sirna_fasta is None:
        raise ValueError("Must provide either risearch_file (a file) or sirna_fasta")
    if is_running_risearch and target_fasta is None:
        raise ValueError(
            "sirna_fasta requires target_fasta when running RIsearch dynamically"
        )

    # --- Acquire predictions ---
    if is_running_risearch:
        # Mode 1: integrated RIsearch execution
        risearch_service = RIsearchService()
        risearch_service.validate_sirna_fasta(sirna_fasta)  # raises ValueError on dups
        index_path = (
            target_index
            if target_index is not None
            else risearch_service.index_target(target_fasta)
        )
        df = risearch_service.run_search(
            query_path=sirna_fasta,
            index_path=index_path,
            target_fasta=target_fasta,
        )
    else:
        # Mode 2: pre-computed RIsearch file
        with profiler.stage("Load predictions") as _s:
            df = risearch_parser.load(Path(risearch_file))
            _s.rows_out = df.height

    logger.info(
        f"Loaded {df.height:,} predictions for {df['sirna_id'].n_unique()} siRNAs"
    )

    if sense_only:
        df = df.filter(pl.col("strand") == "+")

    # Diagnostics the CLI prints (loaded count / energy range / intersection
    # counts). Collected here because the wrapper no longer holds these frames.
    _summary_ris = risearch_parser.summary(df)
    report: dict = {
        "n_loaded": _summary_ris["row_count"],
        "energy_min": _summary_ris["energy_min"],
        "energy_max": _summary_ris["energy_max"],
        "sense_only": sense_only,
        "n_features": None,
        "n_intersected": None,
        # Raw predictions head (6 narrow columns) for the CLI's --verbose preview;
        # the final frame has too many columns to render without truncation.
        "preview": df.head(5),
    }

    # --- Transcriptome + intersection ---
    if gtf_file:
        trans_parser = AnnotationParser()
        with profiler.stage("Load transcriptome") as _s:
            df_trans = trans_parser.load_gtf(
                Path(gtf_file),
                feature=feature_type,
                score_col=expression_metric,
                format=transcriptome_format,
            )
            _s.rows_out = df_trans.height
        logger.info(
            f"Transcriptome: {df_trans.height:,} rows, {df_trans['gene_id'].n_unique()} genes"
        )
        with profiler.stage("Intersection", rows_in=df.height) as _s:
            intersector = IntersectionService()
            df = intersector.intersect(
                df, df_trans, mode=predictions_type, workers=n_workers
            )
            _s.rows_out = df.height
        report["n_features"] = df_trans.height
        report["n_intersected"] = df.height
        logger.info(f"After intersection: {df.height:,} rows")

    alpha_gamma_pairs = _build_alpha_gamma_pairs(alpha, gamma)
    theta_vals = _parse_theta(theta)

    unique_sirnas = df["sirna_id"].unique() if "sirna_id" in df.columns else []
    has_custom_params = len(alpha_gamma_pairs) > 1 or len(theta_vals) > 0
    is_multi_sirna = len(unique_sirnas) > 1 or has_custom_params

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
        except Exception as e:
            raise ValueError(f"Error parsing on-target mapping file: {e}") from e

    def _finish(prob_service: ProbabilityService, frame: pl.DataFrame):
        if is_multi_sirna:
            with profiler.stage(
                "Probabilities (per-siRNA)", rows_in=frame.height
            ) as _s:
                frame, meta = prob_service.calculate_probabilities_per_sirna(
                    frame,
                    alpha_gamma_pairs=alpha_gamma_pairs,
                    theta_values=theta_vals,
                    on_target_map=on_target_map if on_target_map else None,
                    on_target_expression=on_target_expression,
                )
                _s.rows_out = frame.height
        else:
            with profiler.stage(
                "Probabilities (single-siRNA)", rows_in=frame.height
            ) as _s:
                frame, meta = prob_service.calculate_probabilities(
                    frame,
                    on_target_path=on_target_file,
                    query_path=query_file,
                    on_target_expression=on_target_expression,
                    on_target_accessibility_path=on_target_accessibility,
                    on_target_risearch_path=on_target_risearch_file,
                )
                _s.rows_out = frame.height
        logger.info(f"After probabilities: {frame.height:,} rows")

        # Legacy report must be rendered here (reuses prob_service's accessibility
        # service, avoiding a second on-the-fly fold). Replicates the CLI's call,
        # which uses calculate_legacy_format's default alpha/gamma pairs.
        if legacy_format:
            sirna_id = query_file.stem if query_file else "siRNA"
            meta["legacy_text"] = prob_service.calculate_legacy_format(
                frame,
                sirna_id=sirna_id,
                on_target_path=on_target_file,
                query_path=query_file,
                on_target_expression=on_target_expression,
                on_target_accessibility_path=on_target_accessibility,
                on_target_risearch_path=on_target_risearch_file,
                verbose=detailed_report,
            )
        meta["_report"] = report
        return frame, meta

    # --- Accessibility service selection, then probabilities ---
    if accessibility_dir:
        acc_service = GenomeAccessibilityService(Path(accessibility_dir), max_cached=4)
        prob_service = ProbabilityService(acc_service, temperature=temperature)
        return _finish(prob_service, df)
    elif genome_file:
        # Compute accessibility on-the-fly into a temp dir; the dir must outlive the
        # probability + legacy computation (both read the profiles), so do them
        # inside the context manager before it is torn down.
        with tempfile.TemporaryDirectory(prefix="risearch_accessibility_") as temp_dir:
            acc_service = GenomeAccessibilityService(Path(temp_dir), max_cached=4)
            acc_service.compute_genome_accessibility(
                Path(genome_file),
                window_size=window_size,
                max_span=max_span,
                unpaired_prob=unpaired_prob,
                progress_callback=accessibility_progress_callback,
                temperature=temperature,
            )
            prob_service = ProbabilityService(acc_service, temperature=temperature)
            return _finish(prob_service, df)
    else:
        return _finish(ProbabilityService(None, temperature=temperature), df)


# ---------------------------------------------------------------------------
# Directory core (generator, one DataFrame per siRNA)
# ---------------------------------------------------------------------------
def compute_off_targets_directory(
    *,
    input_dir: Path,
    sirna_fasta: Optional[Path] = None,
    gtf_file: Optional[Path] = None,
    feature_type: str = "exon",
    expression_metric: str = "RPKM",
    transcriptome_format: str = "auto",
    accessibility_dir: Optional[Path] = None,
    temperature: float = 37.0,
    on_target_ids_file: Optional[Path] = None,
    on_target_expression: float = 1000.0,
    alpha: str = "1.0",
    gamma: str = "1.0",
    theta: str = "",
    sense_only: bool = False,
    predictions_type: str = "gw",
    n_workers: int = 1,
) -> Generator[tuple[pl.DataFrame, dict], None, None]:
    """Yield ``(df, meta)`` for each siRNA file in *input_dir*.

    A generator: the process pool and the temp Arrow-IPC transcriptome are held
    open across yields and cleaned up in a ``finally``. Consume it fully (or call
    ``.close()`` / wrap in ``contextlib.closing``) for deterministic teardown.
    Raises ``FileNotFoundError`` if the directory contains no prediction files.
    """
    input_dir = Path(input_dir)
    risearch_parser = RIsearchParser()
    all_files = risearch_parser.list_directory_files(input_dir)
    if not all_files:
        raise FileNotFoundError(f"No RIsearch files found in {input_dir}")

    # Serialize transcriptome to a temp Arrow IPC file so workers can memory-map
    # it instead of re-parsing the BED/GTF. Workers share pages via the OS cache.
    _ipc_tmp: Optional[str] = None
    _ipc_path: Optional[Path] = None
    if gtf_file:
        trans_parser = AnnotationParser()
        df_trans = trans_parser.load_gtf(
            Path(gtf_file),
            feature=feature_type,
            score_col=expression_metric,
            format=transcriptome_format,
        )
        _ipc_tmp = tempfile.mkdtemp(prefix="risearch_ipc_")
        _ipc_path = Path(_ipc_tmp) / "transcriptome.arrow"
        df_trans.write_ipc(_ipc_path)

    on_target_map: dict[str, str] = {}
    if on_target_ids_file is not None:
        with open(on_target_ids_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        on_target_map[parts[0]] = parts[1]

    # Self-hybridisation E_min per siRNA (matches the legacy E_min semantics).
    self_hyb_emin: dict[str, float] = {}
    if sirna_fasta is not None:
        self_hyb_emin = RIsearchService().self_hybridization_emin_batch(sirna_fasta)

    alpha_gamma_pairs = _build_alpha_gamma_pairs(alpha, gamma)
    theta_vals = _parse_theta(theta)

    n_proc = min(n_workers, len(all_files))
    # Polars threads per worker: keeps N_workers × threads ≤ available cores.
    polars_threads_per_worker = max(1, n_workers // n_proc)

    # spawn: each worker starts a fresh interpreter, avoiding fork+Rayon deadlocks.
    _ctx = multiprocessing.get_context("spawn")
    _init_args = (
        str(_ipc_path) if gtf_file else "",
        str(accessibility_dir) if accessibility_dir else "",
        polars_threads_per_worker,
        on_target_map,
        self_hyb_emin,
        temperature,
    )

    try:
        with ProcessPoolExecutor(
            max_workers=n_proc,
            mp_context=_ctx,
            initializer=_init_worker,
            initargs=_init_args,
        ) as pool:
            futures = {
                pool.submit(
                    _process_one_sirna,
                    str(f),
                    alpha_gamma_pairs,
                    theta_vals,
                    on_target_expression,
                    sense_only,
                    predictions_type,
                ): f
                for f in all_files
            }

            for future in as_completed(futures):
                # Pop the future immediately so its held result (arrow table +
                # metadata) can be GC'd as soon as the loop variable rotates.
                f = futures.pop(future)
                try:
                    arrow_table, batch_metadata = future.result()
                except Exception as exc:
                    logger.error(f"Worker failed for {f.name}: {exc}")
                    continue

                if arrow_table is None:
                    continue

                df_chunk = pl.from_arrow(arrow_table)
                del arrow_table  # Arrow copy no longer needed; Polars owns the data
                yield df_chunk, batch_metadata
    finally:
        # Clean up the temp Arrow IPC dir used to share the transcriptome.
        if _ipc_tmp:
            shutil.rmtree(_ipc_tmp, ignore_errors=True)
