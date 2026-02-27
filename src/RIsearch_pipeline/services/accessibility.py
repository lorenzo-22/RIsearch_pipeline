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
    """
    Worker function: fold a single island and extract opening energies.

    Must be a top-level function (not a method) for ProcessPoolExecutor
    pickling. Receives only small, pickle-friendly arguments.

    Args:
        rna_seq: Island subsequence (already T→U converted).
        binding_sites: List of (orig_s_0based, orig_e_0based) sites in this island.
        ctx_start: 0-based start of the context window (for + strand positioning).
        seq_len: Full chromosome length (for - strand positioning).
        strand: '+' or '-'.
        rev_start: Reversed start position (only used for - strand).
        window_size: -W parameter.
        max_span: -L parameter.
        unpaired_prob: -u parameter.

    Returns:
        Dict mapping (orig_s, orig_e) → opening_energy (float).

    Complexity: O(L × W) for pfl_fold_up where L = len(rna_seq).
    """
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


class GenomeAccessibilityService:
    """
    Service to pre-compute and query genomic accessibility profiles.

    Stores profiles as memory-mapped numpy arrays (float16) for efficient random access.
    """

    def __init__(self, data_dir: Path, max_cached: int = 4):
        self.data_dir = data_dir
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._max_cached = max_cached
        # OrderedDict for LRU eviction: profile_key -> numpy array
        self._profiles: OrderedDict[str, np.ndarray] = OrderedDict()
        # Metadata flags stored separately (not counted toward LRU slots)
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
        progress_callback=None,
    ) -> Dict[str, Path]:
        """
        Compute accessibility for all sequences in the genome FASTA.
        Computes profiles for both forward (+) and reverse (-) strands.

        Args:
            genome_path: Path to genome FASTA file.
            window_size: -W parameter (default 80)
            max_span: -L parameter (default 40)
            unpaired_prob: -u parameter (default 30)

        Returns:
            Dictionary mapping chromosome names to their profile file paths.
        """
        if not HAS_VIENNA_BINDINGS:
            raise AccessibilityError(
                "ViennaRNA Python bindings ('import RNA') not found. "
                "Please install ViennaRNA or use the CLI fallback (not yet implemented)."
            )

        from RIsearch_pipeline.services.helpers import read_fasta, reverse_complement

        results = {}

        logger.info(f"Computing accessibility for genome: {genome_path}")

        for chrom, sequence in read_fasta(genome_path):
            seq_len = len(sequence)
            logger.info(f"Processing {chrom} (length {seq_len})...")

            # Process both strands
            for strand in ["+", "-"]:
                logger.info(f"  Computing {strand} strand...")

                # Use reverse complement for minus strand
                seq_to_process = (
                    sequence if strand == "+" else reverse_complement(sequence)
                )

                # Create profile array for opening energies
                profile = np.zeros(seq_len, dtype=np.float32)
                RT = 0.616  # kcal/mol at 37°C

                try:
                    # Use RNA.pfl_fold_up - direct equivalent of RNAplfold -u
                    # Returns 2D array: result[i][u] = P(segment of size u starting at position i is unpaired)
                    # Note: result is 1-based indexing
                    probs_matrix = RNA.pfl_fold_up(
                        seq_to_process, unpaired_prob, window_size, max_span
                    )

                    # Extract probabilities for our target unpaired length
                    # probs_matrix[i][u] is probability for segment of length u at position i
                    for i in range(1, seq_len + 1):
                        if i < len(probs_matrix) and unpaired_prob < len(
                            probs_matrix[i]
                        ):
                            p = probs_matrix[i][unpaired_prob]
                            if p is not None and p > 0:
                                # Convert probability to opening energy: E = -RT * ln(P)
                                profile[i - 1] = -RT * np.log(p)
                            else:
                                profile[i - 1] = 25.5  # Max storable value
                        else:
                            profile[i - 1] = 0.0

                except Exception as e:
                    logger.error(f"Error calling ViennaRNA on {chrom} {strand}: {e}")
                    raise

                # For minus strand, reverse the profile like old pipeline does
                if strand == "-":
                    profile = profile[::-1]

                # Save to disk
                out_path = self.get_profile_path(chrom, strand)
                np.save(out_path, profile)

                # Save readable TSV for FULL matrix validation (all u lengths)
                # This matches raw RNAplfold output format
                tsv_path = out_path.with_suffix(".tsv")
                with open(tsv_path, "w") as f:
                    # Header
                    header = "position" + "".join(
                        [f"\tu{u}" for u in range(1, unpaired_prob + 1)]
                    )
                    f.write(header + "\n")

                    for i in range(1, seq_len + 1):
                        row_vals = [str(i)]
                        for u in range(1, unpaired_prob + 1):
                            val = 0.0  # Default prob
                            if i < len(probs_matrix) and u < len(probs_matrix[i]):
                                p = probs_matrix[i][u]
                                if p is not None and p > 0:
                                    val = p  # Write raw probability to TSV
                            row_vals.append(f"{val:.6f}")
                        f.write("\t".join(row_vals) + "\n")

                results[f"{chrom}_{strand}"] = out_path

            if progress_callback:
                progress_callback(advance=1, description=f"Processing {chrom}")

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
        """
        Compute accessibility only for regions with predicted binding sites.

        Scans all per-siRNA RIsearch files in risearch_dir to extract unique
        binding site coordinates, groups by chromosome, merges nearby sites
        into islands (with W-bp flanking), folds only those islands using
        ViennaRNA, and saves per-site opening energies to Parquet.

        Args:
            genome_path: Path to genome FASTA file.
            output_path: Path for the output .parquet file.
            risearch_dir: Directory of per-siRNA RIsearch output files.
            risearch_file: Single merged RIsearch TSV file (alternative to risearch_dir).
            window_size: -W parameter (default 80).
            max_span: -L parameter (default 40).
            unpaired_prob: -u parameter (default 30).
            workers: Number of parallel workers (default 1 = serial).
            progress: Optional Rich Progress object for live progress bars.
            verbose: Log per-island folding details.

        Returns:
            Path to the output Parquet file.

        Complexity:
            O(S) lazy scan where S = total predictions across all files.
            O(I × W) folding time where I = number of merged islands.
            O(C) memory where C = results for one chromosome (streamed to disk).
        """
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
        from RIsearch_pipeline.services.helpers import (
            read_fasta,
            reverse_complement,
            merge_intervals,
        )

        # --- Step 1: Build a LazyFrame and extract chromosome names ---
        if progress:
            scan_task = progress.add_task("Scanning RIsearch files...", total=None)

        if risearch_file:
            # Single file — auto-detect column layout
            import gzip

            opener = gzip.open if str(risearch_file).endswith(".gz") else open
            with opener(risearch_file, "rt") as fh:
                first_line = fh.readline().strip()
            n_cols = len(first_line.split("\t"))

            logger.info(f"Scanning {risearch_file.name} ({n_cols}-column format)...")

            if n_cols <= 4:
                col_map = {
                    "column_1": "chrom",
                    "column_2": "start",
                    "column_3": "end",
                    "column_4": "strand",
                }
            else:
                col_map = {
                    "column_4": "chrom",
                    "column_5": "start",
                    "column_6": "end",
                    "column_7": "strand",
                }

            base_lf = pl.scan_csv(
                risearch_file, separator="\t", has_header=False
            ).select(
                [
                    pl.col(list(col_map.keys())[0]).alias("chrom"),
                    pl.col(list(col_map.keys())[1]).cast(pl.Int32).alias("start"),
                    pl.col(list(col_map.keys())[2]).cast(pl.Int32).alias("end"),
                    pl.col(list(col_map.keys())[3]).alias("strand"),
                ]
            )
            source_desc = risearch_file.name
        else:
            # Directory of per-siRNA files
            from RIsearch_pipeline.services.risearch_parser import RIsearchParser

            parser = RIsearchParser()
            all_files = parser.list_directory_files(risearch_dir)
            if not all_files:
                raise AccessibilityError(f"No RIsearch files found in {risearch_dir}")

            logger.info(
                f"Scanning {len(all_files)} RIsearch files for binding "
                "site coordinates..."
            )

            lazy_frames = []
            skipped = 0
            for f in all_files:
                if f.stat().st_size == 0:
                    skipped += 1
                    continue
                lf = pl.scan_csv(f, separator="\t", has_header=False).select(
                    [
                        pl.col("column_4").alias("chrom"),
                        pl.col("column_5").cast(pl.Int32).alias("start"),
                        pl.col("column_6").cast(pl.Int32).alias("end"),
                        pl.col("column_7").alias("strand"),
                    ]
                )
                lazy_frames.append(lf)

            if skipped:
                logger.info(f"Skipped {skipped} empty files")
            if not lazy_frames:
                raise AccessibilityError("All RIsearch files are empty")

            base_lf = pl.concat(lazy_frames)
            source_desc = f"{len(all_files)} files"

        # Extract only the unique chromosome names (tiny: ~24 strings)
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

        # --- Step 2: Stream genome FASTA — process one chromosome at a time ---
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

            # Lazy scan + filter for this chromosome only — collect just this chrom's sites
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

                # Convert 1-based RIsearch coords to 0-based intervals
                starts = strand_sites["start"].to_numpy()
                ends = strand_sites["end"].to_numpy()
                intervals = [(int(s) - 1, int(e)) for s, e in zip(starts, ends)]

                # Merge nearby intervals with W-bp padding
                islands = merge_intervals(intervals, padding=window_size)
                logger.info(
                    f"  {strand} strand: {len(intervals)} sites → "
                    f"{len(islands)} islands"
                )

                # Prepare sequence for this strand
                if strand == "-":
                    strand_seq = reverse_complement(chrom_seq)
                else:
                    strand_seq = chrom_seq

                # Prepare island folding tasks
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

            # Log total island count for this chromosome
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

            # --- Fold all islands (both strands, parallel or serial) ---
            if progress:
                island_task = progress.add_task(
                    f"  {chrom}: folding {n_tasks} islands",
                    total=n_tasks,
                )

            # Collect pos_energy per strand
            pos_energy_by_strand: dict[str, dict[tuple[int, int], float]] = {
                "+": {},
                "-": {},
            }

            if workers > 1 and n_tasks > 1:
                from concurrent.futures import (
                    ProcessPoolExecutor,
                    as_completed,
                )

                with ProcessPoolExecutor(max_workers=workers) as pool:
                    futures = {}
                    for idx, (strand, _, _, _, fold_args) in enumerate(
                        all_island_tasks
                    ):
                        fut = pool.submit(_fold_island, *fold_args)
                        futures[fut] = (idx, strand)

                    for future in as_completed(futures):
                        idx, strand = futures[future]
                        pos_energy_by_strand[strand].update(future.result())
                        if progress:
                            progress.advance(island_task)
            else:
                for idx, (strand, _, _, _, fold_args) in enumerate(all_island_tasks):
                    pos_energy_by_strand[strand].update(_fold_island(*fold_args))
                    if progress:
                        progress.advance(island_task)

            if progress and island_task is not None:
                progress.remove_task(island_task)
                island_task = None

            # --- Write results per strand ---
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
                    logger.info(
                        f"  Wrote {len(chrom_results)} results for "
                        f"{chrom} {strand} (total: {total_written})"
                    )

            # chrom_seq goes out of scope here — GC can reclaim
            chroms_done += 1
            if progress:
                progress.update(
                    chrom_task,
                    advance=1,
                    description=f"Processing chromosomes ({chroms_done}/{n_chroms})",
                )

        # --- Step 3: Finalize Parquet ---
        logger.info(f"Processed {total_islands} total islands across all chromosomes")
        writer.close()
        logger.info(
            f"Saved {total_written} binding site accessibility values to {output_path}"
        )

        return output_path

    def compute_sequence_accessibility(
        self,
        sequence: str,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
        use_cli: bool = False,
    ) -> np.ndarray:
        """
        Compute accessibility for a single sequence (e.g., on-target).

        Uses ViennaRNA's RNA.pfl_fold_up to compute unpaired probabilities,
        then converts to opening energies. Falls back to RNAplfold CLI if
        Python bindings are unavailable or use_cli=True.

        Args:
            sequence: RNA/DNA sequence string.
            window_size: -W parameter (default 80).
            max_span: -L parameter (default 40).
            unpaired_prob: -u parameter (default 30).
            use_cli: Force using RNAplfold binary instead of Python bindings.

        Returns:
            2D numpy array [seq_len, unpaired_prob] of opening energies.
            Use result[pos, u-1] to get opening energy for length u at position pos.
        """
        seq_len = len(sequence)
        if seq_len == 0:
            return np.array([], dtype=np.float32)

        # Decide whether to use CLI or Python bindings
        if use_cli or not HAS_VIENNA_BINDINGS:
            if use_cli:
                logger.info("Using RNAplfold CLI (--use-rnaplfold-cli flag)")
            else:
                logger.info(
                    "ViennaRNA Python not available, falling back to RNAplfold CLI"
                )
            return self._run_rnaplfold_cli(
                sequence, window_size, max_span, unpaired_prob
            )

        # Use ViennaRNA Python bindings
        # Adjust window/span if sequence is shorter
        w = min(window_size, seq_len)
        max_span_adj = min(max_span, seq_len)

        RT = 0.616  # kcal/mol at 37°C

        # Create 2D profile array [seq_len, unpaired_prob]
        profile = np.full((seq_len, unpaired_prob), 25.5, dtype=np.float32)

        try:
            # RNA.pfl_fold_up returns 2D: result[i][u] = P(segment of size u at position i is unpaired)
            # 1-based indexing
            probs_matrix = RNA.pfl_fold_up(sequence, unpaired_prob, w, max_span_adj)

            for i in range(1, seq_len + 1):
                if i < len(probs_matrix):
                    for u in range(1, unpaired_prob + 1):
                        if u < len(probs_matrix[i]):
                            p = probs_matrix[i][u]
                            if p is not None and p > 0:
                                # Convert probability to opening energy: E = -RT * ln(P)
                                profile[i - 1, u - 1] = -RT * np.log(p)
                            else:
                                profile[i - 1, u - 1] = (
                                    25.5  # High penalty for inaccessible
                                )

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
        """
        Run RNAplfold binary to compute accessibility.

        Command: RNAplfold -W {w} -L {l} -u {u} -O < input.fa
        Produces: {seq_id}_openen file with opening energies.

        Args:
            sequence: RNA/DNA sequence string.
            window_size: -W parameter.
            max_span: -L parameter.
            unpaired_prob: -u parameter.

        Returns:
            2D numpy array [seq_len, unpaired_prob] of opening energies.
        """
        import subprocess
        import tempfile
        import shutil

        seq_len = len(sequence)
        seq_id = "rnaplfold_seq"

        # Adjust parameters for short sequences
        w = min(window_size, seq_len)
        max_span_adj = min(max_span, seq_len)

        # Check if RNAplfold is available
        if not shutil.which("RNAplfold"):
            raise AccessibilityError(
                "RNAplfold binary not found in PATH. "
                "Please install ViennaRNA or add RNAplfold to PATH."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Write sequence to temp FASTA
            fasta_path = Path(tmpdir) / "input.fa"
            with open(fasta_path, "w") as f:
                f.write(f">{seq_id}\n{sequence}\n")

            # Run RNAplfold
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
                    timeout=300,  # 5 minute timeout
                )

                if result.returncode != 0:
                    logger.error(f"RNAplfold stderr: {result.stderr}")
                    raise AccessibilityError(
                        f"RNAplfold failed with exit code {result.returncode}"
                    )

            except subprocess.TimeoutExpired:
                raise AccessibilityError("RNAplfold timed out after 5 minutes")
            except FileNotFoundError:
                raise AccessibilityError("RNAplfold binary not found")

            # Parse output file
            openen_path = Path(tmpdir) / f"{seq_id}_openen"
            if not openen_path.exists():
                # Try alternate naming
                candidates = list(Path(tmpdir).glob("*_openen"))
                if candidates:
                    openen_path = candidates[0]
                else:
                    raise AccessibilityError(
                        f"RNAplfold output file not found. Files in tmpdir: {list(Path(tmpdir).iterdir())}"
                    )

            logger.debug(f"Parsing RNAplfold output: {openen_path}")
            profile = self._parse_openen_text(openen_path)

            # Clean up dp.ps file (RNAplfold generates this)
            dp_ps = Path(tmpdir) / f"{seq_id}_dp.ps"
            if dp_ps.exists():
                dp_ps.unlink()

        logger.info(f"Computed accessibility via RNAplfold CLI (len={seq_len})")
        return profile

    def _ensure_profile(self, chrom: str, strand: str) -> str:
        """
        Ensure the profile for (chrom, strand) is loaded, applying LRU eviction.

        Returns the profile_key string.

        Complexity: O(1) amortized.  Eviction frees the oldest cached
        profile when the cache exceeds max_cached slots, bounding
        resident memory to max_cached chromosome-sized arrays.
        """
        profile_key = f"{chrom}_{strand}"

        if profile_key in self._profiles:
            # Move to end (most-recently-used)
            self._profiles.move_to_end(profile_key)
            return profile_key

        # --- Load from disk ---
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
                logger.warning(
                    f"Binary file {path} size not divisible by 30, using as 1D"
                )
                self._profiles[profile_key] = raw_data
            self._profile_flags[f"{profile_key}_is_legacy_bin"] = True
        elif "openen" in path.name:
            logger.info(f"Parsing legacy text file: {path}")
            self._profiles[profile_key] = self._parse_openen_text(path)

        logger.info(f"Loaded accessibility profile for {chrom} {strand} from {path}")

        # --- LRU eviction ---
        while len(self._profiles) > self._max_cached:
            evicted_key, _ = self._profiles.popitem(last=False)
            # Clean up associated flags
            self._profile_flags.pop(f"{evicted_key}_is_legacy_bin", None)
            logger.debug(f"Evicted profile cache entry: {evicted_key}")

        return profile_key

    def query(self, chrom: str, start: int, end: int, strand: str = "+") -> np.ndarray:
        """
        Query accessibility for a region (0-based, half-open [start, end)).

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (0-based, exclusive)
            strand: Strand ('+' or '-')

        Returns:
            Numpy array of opening energies for the region.
        """
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
        """
        Fast-path: return a single opening energy value for a prediction.

        Mirrors the old pipeline's get_opening_energy() logic:
        - For 2D profiles (matrix [pos, u]): picks the value at the 3'-end
          position using the interaction length as the u-column index.
        - For 1D profiles: picks the 3'-end value directly.
        - Quantizes to 0.1 resolution (match legacy behavior).

        Args:
            chrom: Chromosome name.
            start: 1-based start position (RIsearch convention).
            end:   1-based end position (inclusive, RIsearch convention).
            strand: '+' or '-'.

        Returns:
            Opening energy as float, quantized to 0.1 kcal/mol.

        Complexity: O(1) time, O(1) space per call (profile already mmap'd).
        """
        profile_key = self._ensure_profile(chrom, strand)
        profile = self._profiles[profile_key]
        is_legacy_bin = self._profile_flags.get(f"{profile_key}_is_legacy_bin", False)
        seq_len = len(profile)

        # Convert 1-based RIsearch coords to 0-based
        start0 = start - 1
        end0 = end  # half-open

        if start0 < 0 or end0 > seq_len:
            return 10.0  # Default penalty for out-of-bounds

        interaction_len = end0 - start0

        if profile.ndim == 2:
            matrix_width = profile.shape[1]
            col_idx = min(interaction_len, matrix_width) - 1
            if col_idx < 0:
                col_idx = 0
            # 3' end for +, 5' end for -
            row_idx = end0 - 1 if strand == "+" else start0
            raw_val = float(profile[row_idx, col_idx])
            if is_legacy_bin:
                raw_val /= 10.0
        else:
            # 1D profile
            row_idx = end0 - 1 if strand == "+" else start0
            raw_val = float(profile[row_idx])
            if is_legacy_bin:
                raw_val /= 10.0

        # Quantize to 0.1 resolution (legacy compatibility)
        return int(round(raw_val * 10.0)) / 10.0

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
        """Find profile file, prioritizing .npy, .bin, then text."""
        suffix = "plus" if strand == "+" else "minus"

        # 1. Standard NPY
        p1 = self.data_dir / f"{chrom}_{suffix}.access.npy"
        if p1.exists():
            return p1

        # 2. Standard BIN
        p2 = self.data_dir / f"{chrom}_{suffix}.access.bin"
        if p2.exists():
            return p2

        # 3. Legacy BIN (glob)
        # e.g. chr1_plus.bin, chr1.open.acc.bin, chr1_rev.open.acc.bin
        candidates = list(self.data_dir.glob(f"{chrom}*{suffix}*.bin"))
        if candidates:
            return candidates[0]

        # 3b. Legacy open/rev naming (open.acc.bin / rev.open.acc.bin)
        if strand == "+":
            for candidate in self.data_dir.glob(f"{chrom}*.open.acc.bin"):
                if "rev.open.acc.bin" not in candidate.name:
                    return candidate
        else:
            candidates_rev = list(self.data_dir.glob(f"{chrom}*rev.open.acc.bin"))
            if candidates_rev:
                return candidates_rev[0]

        # 4. Legacy Text (openen)
        # e.g. chr1_0_75631_openen
        # Use _ or . as separator to avoid prefix matching (e.g. transcript_3 vs transcript_35)
        candidates_txt = list(self.data_dir.glob(f"{chrom}_*openen"))
        if not candidates_txt:
            candidates_txt = list(self.data_dir.glob(f"{chrom}.*openen"))

        if candidates_txt:
            # Pick the shortest or first?
            return candidates_txt[0]

        return None
