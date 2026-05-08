from loguru import logger
from pathlib import Path
from collections import OrderedDict
from typing import Dict
import numpy as np

# Try to import ViennaRNA
try:
    import RNA

    HAS_VIENNA_BINDINGS = True
except ImportError:
    HAS_VIENNA_BINDINGS = False


class AccessibilityError(Exception):
    """Base exception for accessibility service."""

    pass


_COMPLEMENT_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def _fold_full_chromosome(
    sequence: str,
    chrom: str,
    window_size: int,
    max_span: int,
    unpaired_prob: int,
    output_path: str,
    chunk_size: int = 1_000_000,
) -> str:
    """Fold both strands of a chromosome in chunks and stream to Parquet.

    Top-level (not a method) so ProcessPoolExecutor can pickle it; only the
    output path string is pickled back, avoiding O(N×u) OOM on large chromosomes.
    """
    import RNA as _RNA
    import math
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq

    RT = 0.616  # kcal/mol at 37°C
    seq_len = len(sequence)

    schema = pa.schema(
        [("position", pa.int32()), ("strand", pa.utf8())]
        + [(f"u{i}", pa.float32()) for i in range(1, unpaired_prob + 1)]
    )

    def _fold_strand(writer: pq.ParquetWriter, seq: str, strand_char: str) -> None:
        rna_seq = seq.upper().replace("T", "U")
        pos = 0

        while pos < seq_len:
            ext_start = max(0, pos - window_size)
            ext_end = min(seq_len, pos + chunk_size + window_size)
            chunk_seq = rna_seq[ext_start:ext_end]
            chunk_len = len(chunk_seq)

            data_start_in_chunk = pos - ext_start
            data_end_in_chunk = min(seq_len, pos + chunk_size) - ext_start

            w = min(window_size, chunk_len)
            l_adj = min(max_span, chunk_len)

            try:
                probs_matrix = _RNA.pfl_fold_up(chunk_seq, unpaired_prob, w, l_adj)
                probs_ok = True
            except Exception:
                probs_ok = False
                probs_matrix = None

            n_pos = data_end_in_chunk - data_start_in_chunk
            positions = np.arange(pos + 1, pos + n_pos + 1, dtype=np.int32)
            u_matrix = np.full((n_pos, unpaired_prob), 25.5, dtype=np.float32)

            if probs_ok:
                for i in range(n_pos):
                    chunk_pos_1based = data_start_in_chunk + i + 1
                    if chunk_pos_1based >= len(probs_matrix):
                        continue
                    row_probs = probs_matrix[chunk_pos_1based]
                    for u in range(1, unpaired_prob + 1):
                        if u >= len(row_probs):
                            break
                        p = row_probs[u]
                        if p is not None and p > 0:
                            ev = -RT * math.log(p)
                            u_matrix[i, u - 1] = round(ev * 10.0) / 10.0
            else:
                u_matrix[:] = 10.0  # folding error

            arrays = [
                pa.array(positions),
                pa.array([strand_char] * n_pos, type=pa.utf8()),
            ] + [
                pa.array(u_matrix[:, u - 1])
                for u in range(1, unpaired_prob + 1)
            ]
            writer.write_batch(pa.RecordBatch.from_arrays(arrays, schema=schema))

            pos += chunk_size

    with pq.ParquetWriter(output_path, schema) as writer:
        _fold_strand(writer, sequence, "+")
        _fold_strand(writer, _reverse_complement(sequence), "-")

    return output_path


def _fold_island(
    rna_seq: str,
    binding_sites: list[tuple[int, int]],
    ctx_start: int,
    seq_len: int,
    strand: str,
    rev_start: int,
    window_size: int,
    max_span: int,
    unpaired_prob: int,
) -> dict[tuple[int, int], float]:
    """Fold one island and return {(orig_s, orig_e): opening_energy}. Top-level for pickling."""
    import RNA as _RNA

    RT = 0.616  # kcal/mol at 37°C
    sub_len = len(rna_seq)
    w = min(window_size, sub_len)
    l_adj = min(max_span, sub_len)

    result: dict[tuple[int, int], float] = {}

    try:
        probs_matrix = _RNA.pfl_fold_up(rna_seq, unpaired_prob, w, l_adj)
    except Exception:
        # Default penalty for all sites in this island
        for orig_s, orig_e in binding_sites:
            result[(orig_s, orig_e)] = 10.0
        return result

    for orig_s, orig_e in binding_sites:
        interaction_len = orig_e - orig_s

        if strand == "+":
            pos_in_sub = (orig_e - 1) - ctx_start
        else:
            pos_in_sub = (seq_len - orig_s - 1) - rev_start

        idx_1based = pos_in_sub + 1

        u_col = min(interaction_len, unpaired_prob)
        if u_col < 1:
            u_col = 1

        energy_val = 10.0
        if 0 < idx_1based < len(probs_matrix) and u_col < len(probs_matrix[idx_1based]):
            p = probs_matrix[idx_1based][u_col]
            if p is not None and p > 0:
                import math

                energy_val = -RT * math.log(p)
                energy_val = int(round(energy_val * 10.0)) / 10.0

        result[(orig_s, orig_e)] = energy_val

    return result


def _fold_island_batch(
    tasks: list[tuple],
) -> dict[str, dict[tuple[int, int], float]]:
    """Fold a batch of islands; reduces ProcessPoolExecutor IPC overhead. Top-level for pickling."""
    merged: dict[str, dict[tuple[int, int], float]] = {}
    for strand, fold_args in tasks:
        result = _fold_island(*fold_args)
        if strand not in merged:
            merged[strand] = {}
        merged[strand].update(result)
    return merged


class GenomeAccessibilityService:
    """Pre-compute and query genomic accessibility profiles (memory-mapped numpy arrays)."""

    def __init__(self, data_dir: Path, max_cached: int = 4):
        self.data_dir = data_dir
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._max_cached = max_cached
        self._profiles: OrderedDict[str, np.ndarray] = OrderedDict()
        self._profile_flags: Dict[str, bool] = {}

    def get_profile_path(self, chrom: str, strand: str = "+") -> Path:
        """Get path for accessibility profile. Strand can be '+' or '-'."""
        suffix = "plus" if strand == "+" else "minus"
        return self.data_dir / f"{chrom}_{suffix}.access.npy"

    def compute_genome_accessibility(
        self,
        genome_path: Path,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
        workers: int = 1,
        progress=None,
        progress_callback=None,
    ) -> Dict[str, Path]:
        """Fold all chromosomes; return {chrom: parquet_path}."""
        from concurrent.futures import ProcessPoolExecutor, as_completed

        if not HAS_VIENNA_BINDINGS:
            raise AccessibilityError(
                "ViennaRNA Python bindings ('import RNA') not found. "
                "Please install ViennaRNA or use the CLI fallback."
            )

        from RIsearch_pipeline.services.helpers import read_fasta

        chromosomes = list(read_fasta(genome_path))
        # longest-first scheduling minimises tail latency
        chromosomes.sort(key=lambda cs: len(cs[1]), reverse=True)

        logger.info(
            f"Computing full-genome accessibility for {len(chromosomes)} "
            f"chromosome(s) using {workers} worker(s)"
        )

        total_tasks = len(chromosomes)
        prog_task = None
        if progress is not None:
            prog_task = progress.add_task("Folding chromosomes...", total=total_tasks)

        results: Dict[str, Path] = {}
        futures = {}

        with ProcessPoolExecutor(max_workers=workers) as pool:
            for chrom, sequence in chromosomes:
                out_path = str(self.data_dir / f"{chrom}.accessibility.parquet")
                fut = pool.submit(
                    _fold_full_chromosome,
                    sequence,
                    chrom,
                    window_size,
                    max_span,
                    unpaired_prob,
                    out_path,
                )
                futures[fut] = chrom

            for fut in as_completed(futures):
                chrom = futures[fut]
                try:
                    path = fut.result()
                except Exception as e:
                    logger.error(f"Error folding {chrom}: {e}")
                    raise

                results[chrom] = Path(path)
                logger.info(f"Completed {chrom} → {path}")

                if progress is not None and prog_task is not None:
                    progress.update(
                        prog_task,
                        advance=1,
                        description=f"Folded {chrom}",
                    )
                elif progress_callback:
                    progress_callback(advance=1, description=f"Folded {chrom}")

        return results

    def compute_binding_site_accessibility(
        self,
        genome_path: Path,
        output_path: Path,
        risearch_dir: Path = None,
        risearch_file: Path = None,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
        workers: int = 1,
        progress=None,
        verbose: bool = False,
    ) -> Path:
        """Fold only binding-site islands; write per-site opening energies to Parquet."""
        if not HAS_VIENNA_BINDINGS:
            raise AccessibilityError(
                "ViennaRNA Python bindings ('import RNA') not found. "
                "Please install ViennaRNA or use the CLI fallback."
            )

        if not risearch_dir and not risearch_file:
            raise AccessibilityError("Provide either risearch_dir or risearch_file")

        import polars as pl
        import pyarrow as pa
        import pyarrow.parquet as pq
        from RIsearch_pipeline.services.helpers import read_fasta, reverse_complement, merge_intervals

        if progress:
            scan_task = progress.add_task("Scanning RIsearch files...", total=None)

        if risearch_file:
            import gzip

            opener = gzip.open if str(risearch_file).endswith(".gz") else open
            with opener(risearch_file, "rt") as fh:
                first_line = fh.readline().strip()
            n_cols = len(first_line.split("\t"))

            logger.info(f"Scanning {risearch_file.name} ({n_cols}-column format)...")

            if n_cols <= 4:
                chrom_col, start_col, end_col, strand_col = "column_1", "column_2", "column_3", "column_4"
            else:
                chrom_col, start_col, end_col, strand_col = "column_4", "column_5", "column_6", "column_7"

            base_lf = pl.scan_csv(risearch_file, separator="\t", has_header=False).select([
                pl.col(chrom_col).alias("chrom"),
                pl.col(start_col).cast(pl.Int32).alias("start"),
                pl.col(end_col).cast(pl.Int32).alias("end"),
                pl.col(strand_col).alias("strand"),
            ])
            source_desc = risearch_file.name
        else:
            from RIsearch_pipeline.services.risearch_parser import RIsearchParser

            parser = RIsearchParser()
            all_files = parser.list_directory_files(risearch_dir)
            if not all_files:
                raise AccessibilityError(f"No RIsearch files found in {risearch_dir}")

            logger.info(f"Scanning {len(all_files)} RIsearch files for binding site coordinates...")

            lazy_frames = []
            skipped = 0
            for f in all_files:
                if f.stat().st_size == 0:
                    skipped += 1
                    continue
                lazy_frames.append(
                    pl.scan_csv(f, separator="\t", has_header=False).select([
                        pl.col("column_4").alias("chrom"),
                        pl.col("column_5").cast(pl.Int32).alias("start"),
                        pl.col("column_6").cast(pl.Int32).alias("end"),
                        pl.col("column_7").alias("strand"),
                    ])
                )

            if skipped:
                logger.info(f"Skipped {skipped} empty files")
            if not lazy_frames:
                raise AccessibilityError("All RIsearch files are empty")

            base_lf = pl.concat(lazy_frames)
            source_desc = f"{len(all_files)} files"

        chroms_in_predictions = set(
            base_lf.select("chrom").unique().collect()["chrom"].to_list()
        )
        n_chroms = len(chroms_in_predictions)
        logger.info(f"Found binding sites on {n_chroms} chromosomes from {source_desc}")

        if progress:
            progress.update(
                scan_task,
                completed=True,
                description=f"Scanned {source_desc} → {n_chroms} chromosomes",
            )
            progress.remove_task(scan_task)

        arrow_schema = pa.schema(
            [
                ("chrom", pa.string()),
                ("start", pa.int32()),
                ("end", pa.int32()),
                ("strand", pa.string()),
                ("opening_energy", pa.float64()),
            ]
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(str(output_path), arrow_schema)
        total_written = 0

        if progress:
            chrom_task = progress.add_task(
                f"Processing chromosomes (0/{n_chroms})",
                total=n_chroms,
            )
            island_task = None
        chroms_done = 0
        total_islands = 0
        for chrom, chrom_seq in read_fasta(genome_path):
            if chrom not in chroms_in_predictions:
                continue

            chrom_sites = base_lf.filter(pl.col("chrom") == chrom).unique().collect()
            seq_len = len(chrom_seq)

            logger.info(
                f"Processing {chrom} ({seq_len:,} bp, "
                f"{chrom_sites.height} unique binding sites, "
                f"{chrom_sites.estimated_size('mb'):.1f} MB)..."
            )
            if progress:
                progress.update(
                    chrom_task,
                    description=f"Processing {chrom} ({chrom_sites.height} sites)",
                )

            all_island_tasks = []
            for strand in ["+", "-"]:
                strand_sites = chrom_sites.filter(pl.col("strand") == strand)
                if strand_sites.height == 0:
                    continue

                starts = strand_sites["start"].to_numpy()
                ends = strand_sites["end"].to_numpy()
                intervals = [(int(s) - 1, int(e)) for s, e in zip(starts, ends)]

                islands = merge_intervals(intervals, padding=window_size)
                logger.info(f"  {strand} strand: {len(intervals)} sites → {len(islands)} islands")

                strand_seq = reverse_complement(chrom_seq) if strand == "-" else chrom_seq

                for island_start, island_end in islands:
                    ctx_start = max(0, island_start - window_size)
                    ctx_end = min(seq_len, island_end + window_size)

                    if strand == "-":
                        rev_start = seq_len - ctx_end
                        rev_end = seq_len - ctx_start
                        island_subseq = strand_seq[rev_start:rev_end]
                    else:
                        rev_start = 0
                        island_subseq = strand_seq[ctx_start:ctx_end]

                    rna_seq = island_subseq.replace("T", "U").replace("t", "u")

                    sites_in_island = [
                        (s, e)
                        for s, e in intervals
                        if island_start <= s and e <= island_end
                    ]

                    all_island_tasks.append(
                        (
                            strand,
                            starts,
                            ends,
                            strand_sites,
                            (
                                rna_seq,
                                sites_in_island,
                                ctx_start,
                                seq_len,
                                strand,
                                rev_start,
                                window_size,
                                max_span,
                                unpaired_prob,
                            ),
                        )
                    )

            n_tasks = len(all_island_tasks)
            logger.info(f"  {chrom}: {n_tasks} total islands across both strands")
            total_islands += n_tasks

            if n_tasks == 0:
                chroms_done += 1
                if progress:
                    progress.update(
                        chrom_task,
                        advance=1,
                        description=f"Processing chromosomes ({chroms_done}/{n_chroms})",
                    )
                continue

            if progress:
                island_task = progress.add_task(
                    f"  {chrom}: folding {n_tasks} islands",
                    total=n_tasks,
                )

            pos_energy_by_strand: dict[str, dict[tuple[int, int], float]] = {"+": {}, "-": {}}

            if workers > 1 and n_tasks > 1:
                from concurrent.futures import ProcessPoolExecutor, as_completed
                from itertools import islice as _islice

                packed = [(strand, fold_args) for strand, _, _, _, fold_args in all_island_tasks]

                # workers×4 chunks: amortise IPC overhead while keeping cores balanced
                chunk_size = max(1, (n_tasks + workers * 4 - 1) // (workers * 4))
                chunks = []
                it = iter(packed)
                while batch_chunk := list(_islice(it, chunk_size)):
                    chunks.append(batch_chunk)

                with ProcessPoolExecutor(max_workers=workers) as pool:
                    futures = {
                        pool.submit(_fold_island_batch, chunk): len(chunk)
                        for chunk in chunks
                    }

                    for future in as_completed(futures):
                        batch_results = future.result()
                        n_done = futures[future]
                        for strand, energy_map in batch_results.items():
                            pos_energy_by_strand[strand].update(energy_map)
                        if progress:
                            for _ in range(n_done):
                                progress.advance(island_task)
            else:
                for idx, (strand, _, _, _, fold_args) in enumerate(all_island_tasks):
                    pos_energy_by_strand[strand].update(_fold_island(*fold_args))
                    if progress:
                        progress.advance(island_task)

            if progress and island_task is not None:
                progress.remove_task(island_task)
                island_task = None

            for strand in ["+", "-"]:
                strand_sites = chrom_sites.filter(pl.col("strand") == strand)
                if strand_sites.height == 0:
                    continue
                pos_energy = pos_energy_by_strand[strand]
                starts = strand_sites["start"].to_numpy()
                ends = strand_sites["end"].to_numpy()

                chrom_results = []
                for row_idx in range(strand_sites.height):
                    s_1based = int(starts[row_idx])
                    e_1based = int(ends[row_idx])
                    s0 = s_1based - 1
                    e0 = e_1based
                    oe = pos_energy.get((s0, e0), 10.0)
                    chrom_results.append(
                        {
                            "chrom": chrom,
                            "start": s_1based,
                            "end": e_1based,
                            "strand": strand,
                            "opening_energy": oe,
                        }
                    )

                # Write this chrom/strand batch as a row group immediately
                if chrom_results:
                    batch_df = pl.DataFrame(chrom_results)
                    writer.write_table(batch_df.to_arrow().cast(arrow_schema))
                    total_written += len(chrom_results)
                    logger.info(f"  Wrote {len(chrom_results)} results for {chrom} {strand} (total: {total_written})")

            chroms_done += 1
            if progress:
                progress.update(
                    chrom_task,
                    advance=1,
                    description=f"Processing chromosomes ({chroms_done}/{n_chroms})",
                )

        logger.info(f"Processed {total_islands} total islands across all chromosomes")
        writer.close()
        logger.info(f"Saved {total_written} binding site accessibility values to {output_path}")
        return output_path

    def compute_sequence_accessibility(
        self,
        sequence: str,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
        use_cli: bool = False,
    ) -> np.ndarray:
        """Return 2D opening-energy array [seq_len, unpaired_prob] for a single sequence."""
        seq_len = len(sequence)
        if seq_len == 0:
            return np.array([], dtype=np.float32)

        if use_cli or not HAS_VIENNA_BINDINGS:
            if use_cli:
                logger.info("Using RNAplfold CLI (--use-rnaplfold-cli flag)")
            else:
                logger.info("ViennaRNA Python not available, falling back to RNAplfold CLI")
            return self._run_rnaplfold_cli(sequence, window_size, max_span, unpaired_prob)

        w = min(window_size, seq_len)
        max_span_adj = min(max_span, seq_len)
        RT = 0.616  # kcal/mol at 37°C
        profile = np.full((seq_len, unpaired_prob), 25.5, dtype=np.float32)

        try:
            probs_matrix = RNA.pfl_fold_up(sequence, unpaired_prob, w, max_span_adj)
            for i in range(1, seq_len + 1):
                if i < len(probs_matrix):
                    for u in range(1, unpaired_prob + 1):
                        if u < len(probs_matrix[i]):
                            p = probs_matrix[i][u]
                            if p is not None and p > 0:
                                profile[i - 1, u - 1] = -RT * np.log(p)
                            else:
                                profile[i - 1, u - 1] = 25.5
            logger.debug(f"Computed accessibility for sequence of length {seq_len}")
        except Exception as e:
            logger.error(f"ViennaRNA error computing accessibility: {e}")
            raise AccessibilityError(f"Failed to compute accessibility: {e}") from e

        return profile

    def _run_rnaplfold_cli(
        self,
        sequence: str,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
    ) -> np.ndarray:
        import subprocess
        import tempfile
        import shutil

        seq_len = len(sequence)
        seq_id = "rnaplfold_seq"
        w = min(window_size, seq_len)
        max_span_adj = min(max_span, seq_len)

        if not shutil.which("RNAplfold"):
            raise AccessibilityError(
                "RNAplfold binary not found in PATH. "
                "Please install ViennaRNA or add RNAplfold to PATH."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "input.fa"
            fasta_path.write_text(f">{seq_id}\n{sequence}\n")

            cmd = f"RNAplfold -W {w} -L {max_span_adj} -u {unpaired_prob} -O"
            logger.debug(f"Running: {cmd}")

            try:
                result = subprocess.run(
                    cmd,
                    shell=True,
                    cwd=tmpdir,
                    stdin=open(fasta_path),
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                if result.returncode != 0:
                    logger.error(f"RNAplfold stderr: {result.stderr}")
                    raise AccessibilityError(f"RNAplfold failed with exit code {result.returncode}")
            except subprocess.TimeoutExpired:
                raise AccessibilityError("RNAplfold timed out after 5 minutes")
            except FileNotFoundError:
                raise AccessibilityError("RNAplfold binary not found")

            openen_path = Path(tmpdir) / f"{seq_id}_openen"
            if not openen_path.exists():
                candidates = list(Path(tmpdir).glob("*_openen"))
                if candidates:
                    openen_path = candidates[0]
                else:
                    raise AccessibilityError(
                        f"RNAplfold output file not found. Files in tmpdir: {list(Path(tmpdir).iterdir())}"
                    )

            profile = self._parse_openen_text(openen_path)

        logger.info(f"Computed accessibility via RNAplfold CLI (len={seq_len})")
        return profile

    def _ensure_profile(self, chrom: str, strand: str) -> str:
        """Load profile for (chrom, strand) into LRU cache; return profile_key."""
        profile_key = f"{chrom}_{strand}"

        if profile_key in self._profiles:
            self._profiles.move_to_end(profile_key)
            return profile_key

        path = self._find_profile(chrom, strand)
        if not path:
            raise AccessibilityError(
                f"Profile for {chrom} {strand} not found in {self.data_dir}. "
                "Expected .access.npy, .access.bin, or legacy open.acc.bin files."
            )

        if path.suffix == ".npy":
            self._profiles[profile_key] = np.load(path, mmap_mode="r")
        elif path.suffix == ".bin":
            raw_data = np.memmap(path, dtype=np.uint8, mode="r")
            if len(raw_data) % 30 == 0:
                self._profiles[profile_key] = raw_data.reshape(-1, 30)
            else:
                logger.warning(f"Binary file {path} size not divisible by 30, using as 1D")
                self._profiles[profile_key] = raw_data
            self._profile_flags[f"{profile_key}_is_legacy_bin"] = True
        elif "openen" in path.name:
            logger.info(f"Parsing legacy text file: {path}")
            self._profiles[profile_key] = self._parse_openen_text(path)

        logger.info(f"Loaded accessibility profile for {chrom} {strand} from {path}")

        while len(self._profiles) > self._max_cached:
            evicted_key, _ = self._profiles.popitem(last=False)
            self._profile_flags.pop(f"{evicted_key}_is_legacy_bin", None)
            logger.debug(f"Evicted profile cache entry: {evicted_key}")

        return profile_key

    def query(self, chrom: str, start: int, end: int, strand: str = "+") -> np.ndarray:
        """Return opening energies for a region (0-based, half-open [start, end))."""
        profile_key = self._ensure_profile(chrom, strand)
        profile = self._profiles[profile_key]
        is_legacy_bin = self._profile_flags.get(f"{profile_key}_is_legacy_bin", False)
        seq_len = len(profile)

        if start < 0 or end > seq_len:
            raise AccessibilityError(
                f"Query {start}-{end} out of bounds for {chrom} {strand} (len {seq_len}, path={self._find_profile(chrom, strand)})"
            )

        data = profile[start:end]
        if is_legacy_bin:
            return data.astype(np.float32) / 10.0
        return data

    def query_single(
        self, chrom: str, start: int, end: int, strand: str = "+"
    ) -> float:
        """Return one opening energy for a 1-based RIsearch prediction (quantized to 0.1 kcal/mol)."""
        profile_key = self._ensure_profile(chrom, strand)
        profile = self._profiles[profile_key]
        is_legacy_bin = self._profile_flags.get(f"{profile_key}_is_legacy_bin", False)
        seq_len = len(profile)

        start0 = start - 1
        end0 = end  # half-open

        if start0 < 0 or end0 > seq_len:
            return 10.0

        interaction_len = end0 - start0

        if profile.ndim == 2:
            matrix_width = profile.shape[1]
            col_idx = max(min(interaction_len, matrix_width) - 1, 0)
            row_idx = end0 - 1 if strand == "+" else start0
            raw_val = float(profile[row_idx, col_idx])
            if is_legacy_bin:
                raw_val /= 10.0
        else:
            row_idx = end0 - 1 if strand == "+" else start0
            raw_val = float(profile[row_idx])
            if is_legacy_bin:
                raw_val /= 10.0

        return int(round(raw_val * 10.0)) / 10.0

    def annotate_opening_energy_vectorized(self, df: "pl.DataFrame") -> "pl.DataFrame":
        """Add `opening_energy` column by vectorized numpy gather per (chrom, strand) group."""
        import polars as pl
        from collections import defaultdict

        n_rows = df.height
        opening_energies = np.full(n_rows, 10.0, dtype=np.float32)

        chroms = df["chrom"].to_numpy()
        starts = df["start"].to_numpy().astype(np.int64)
        ends = df["end"].to_numpy().astype(np.int64)
        strands = df["strand"].to_numpy()

        groups: dict[tuple, list] = defaultdict(list)
        for i in range(n_rows):
            c = str(chroms[i])
            if c == "onTarget":
                continue
            groups[(c, str(strands[i]))].append(i)

        for (chrom, strand), row_indices in groups.items():
            try:
                profile_key = self._ensure_profile(chrom, strand)
            except Exception:
                continue

            profile = self._profiles[profile_key]
            is_legacy_bin = self._profile_flags.get(f"{profile_key}_is_legacy_bin", False)
            seq_len = len(profile)

            idx = np.array(row_indices, dtype=np.int64)
            s = starts[idx]
            e = ends[idx]
            start0 = s - 1
            end0 = e
            interaction_len = end0 - start0
            in_bounds = (start0 >= 0) & (end0 <= seq_len)

            if profile.ndim == 2:
                matrix_width = profile.shape[1]
                col_idx = np.maximum(np.minimum(interaction_len, matrix_width) - 1, 0)
                row_idx = np.where(strand == "+", end0 - 1, start0)
                row_idx = np.clip(row_idx, 0, seq_len - 1)
                raw_vals = profile[row_idx, col_idx].astype(np.float32)
            else:
                row_idx = np.where(strand == "+", end0 - 1, start0)
                row_idx = np.clip(row_idx, 0, seq_len - 1)
                raw_vals = profile[row_idx].astype(np.float32)

            if is_legacy_bin:
                raw_vals = raw_vals / 10.0

            raw_vals = (np.round(raw_vals * 10.0) / 10.0).astype(np.float32)
            opening_energies[idx[in_bounds]] = raw_vals[in_bounds]

        return df.with_columns(pl.Series("opening_energy", opening_energies, dtype=pl.Float32))

    def _parse_openen_text(self, path: Path) -> np.ndarray:
        """
        Parse RNAplfold -O text output.
        Returns 2D array [seq_len, stride] (usually stride=30).
        """
        vals = []
        max_idx = 0
        stride = 0

        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if (
                    not line
                    or line.startswith("#")
                    or line.lower().startswith("position")
                ):
                    continue

                parts = line.split()
                try:
                    pos = int(parts[0])
                    max_idx = max(max_idx, pos)

                    # Values are parts[1:]
                    row_vals = []
                    for s in parts[1:]:
                        if s == "NA" or s == "nan":
                            row_vals.append(25.5)
                        else:
                            row_vals.append(float(s))

                    if stride == 0:
                        stride = len(row_vals)

                    vals.append((pos, row_vals))
                except (ValueError, IndexError):
                    continue

        if max_idx == 0:
            return np.array([], dtype=np.float32)

        # Create 2D array
        arr = np.full((max_idx, stride), 25.5, dtype=np.float32)
        for pos, r_vals in vals:
            if 0 <= pos - 1 < max_idx:
                # Truncate or pad if length mismatch (though unlikely if stride constant)
                # Just take min len
                n = min(len(r_vals), stride)
                arr[pos - 1, :n] = r_vals[:n]

        return arr

    def _find_profile(self, chrom: str, strand: str) -> Path | None:
        """Find profile file: .npy > .bin > legacy bin > legacy text."""
        suffix = "plus" if strand == "+" else "minus"

        p = self.data_dir / f"{chrom}_{suffix}.access.npy"
        if p.exists():
            return p

        p = self.data_dir / f"{chrom}_{suffix}.access.bin"
        if p.exists():
            return p

        candidates = list(self.data_dir.glob(f"{chrom}*{suffix}*.bin"))
        if candidates:
            return candidates[0]

        if strand == "+":
            for name in (f"{chrom}.open.acc.bin", f"{chrom}_open.acc.bin"):
                p = self.data_dir / name
                if p.exists():
                    return p
        else:
            for name in (f"{chrom}.rev.open.acc.bin", f"{chrom}_rev.open.acc.bin"):
                p = self.data_dir / name
                if p.exists():
                    return p

        # Use _ or . separator to avoid prefix collision (transcript_3 vs transcript_35)
        candidates_txt = list(self.data_dir.glob(f"{chrom}_*openen")) or \
                         list(self.data_dir.glob(f"{chrom}.*openen"))
        if candidates_txt:
            return candidates_txt[0]

        return None

    def precompute_parquet_from_profiles(
        self,
        risearch_dir: Path,
        output_path: Path,
    ) -> Path:
        """Look up opening energies for all binding sites in risearch_dir; write to Parquet."""
        import polars as pl

        from RIsearch_pipeline.services.risearch_parser import RIsearchParser

        parser = RIsearchParser()
        all_files = parser.list_directory_files(risearch_dir)
        if not all_files:
            raise AccessibilityError(f"No RIsearch files found in {risearch_dir}")

        logger.info(f"Scanning {len(all_files)} files for unique binding sites...")
        lazy_frames = [
            pl.scan_csv(f, separator="\t", has_header=False).select([
                pl.col("column_4").alias("chrom"),
                pl.col("column_5").cast(pl.Int32).alias("start"),
                pl.col("column_6").cast(pl.Int32).alias("end"),
                pl.col("column_7").alias("strand"),
            ])
            for f in all_files
            if f.stat().st_size > 0
        ]
        sites = (
            pl.concat(lazy_frames)
            .unique(subset=["chrom", "start", "end", "strand"])
            .collect(engine="streaming")
        )
        logger.info(f"Found {sites.height:,} unique binding sites across {sites['chrom'].n_unique()} chromosomes")

        results: list[pl.DataFrame] = []

        for (chrom, strand), group in sites.partition_by(["chrom", "strand"], as_dict=True).items():
            path = self._find_profile(chrom, strand)
            if path is None:
                logger.warning(f"No profile found for {chrom} {strand}, using default penalty 10.0")
                results.append(group.with_columns(pl.lit(10.0).cast(pl.Float32).alias("opening_energy")))
                continue

            is_legacy_bin = False
            if path.suffix == ".npy":
                profile = np.load(path, mmap_mode="r")
            elif path.suffix == ".bin":
                raw = np.memmap(path, dtype=np.uint8, mode="r")
                profile = raw.reshape(-1, 30) if len(raw) % 30 == 0 else raw
                is_legacy_bin = True
            else:
                profile = self._parse_openen_text(path)

            seq_len = len(profile)
            starts = group["start"].to_numpy().astype(np.int64)
            ends = group["end"].to_numpy().astype(np.int64)
            start0 = starts - 1
            end0 = ends
            interaction_len = end0 - start0

            if profile.ndim == 2:
                matrix_width = profile.shape[1]
                col_idx = np.maximum(np.minimum(interaction_len, matrix_width) - 1, 0)
                row_idx = np.where(strand == "+", end0 - 1, start0)
                row_idx = np.clip(row_idx, 0, seq_len - 1)
                raw_vals = profile[row_idx, col_idx].astype(np.float32)
            else:
                row_idx = np.where(strand == "+", end0 - 1, start0)
                row_idx = np.clip(row_idx, 0, seq_len - 1)
                raw_vals = profile[row_idx].astype(np.float32)

            if is_legacy_bin:
                raw_vals = raw_vals / 10.0

            raw_vals = (np.round(raw_vals * 10.0) / 10.0).astype(np.float32)
            results.append(group.with_columns(pl.Series("opening_energy", raw_vals, dtype=pl.Float32)))

        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_df = pl.concat(results, how="diagonal")
        out_df.write_parquet(output_path)
        logger.info(f"Wrote {out_df.height:,} rows to {output_path}")
        return output_path
