from loguru import logger
from pathlib import Path
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


class GenomeAccessibilityService:
    """
    Service to pre-compute and query genomic accessibility profiles.

    Stores profiles as memory-mapped numpy arrays (float16) for efficient random access.
    """

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._profiles: Dict[str, np.memmap] = {}

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

    def query(self, chrom: str, start: int, end: int, strand: str = "+") -> np.ndarray:
        """
        Query accessibility for a region (0-based, half-open [start, end)).

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (0-based, exclusive)
            strand: Strand ('+' or '-')

        Returns:
            Numpy array of probabilities.
        """
        profile_key = f"{chrom}_{strand}"

        if profile_key not in self._profiles:
            # Try to resolve path (support .npy and legacy .bin)
            path = self._find_profile(chrom, strand)

            if not path:
                raise AccessibilityError(
                    f"Profile for {chrom} {strand} not found in {self.data_dir}. "
                    "Expected .access.npy, .access.bin, or legacy open.acc.bin files."
                )

            # Load based on extension
            if path.suffix == ".npy":
                self._profiles[profile_key] = np.load(path, mmap_mode="r")
            elif path.suffix == ".bin":
                # Legacy format: uint8, flattened 2D matrix (Nx30)
                # val = round(energy * 10)
                raw_data = np.memmap(path, dtype=np.uint8, mode="r")
                # Reshape to 2D if divisible by 30
                if len(raw_data) % 30 == 0:
                    self._profiles[profile_key] = raw_data.reshape(-1, 30)
                else:
                    logger.warning(
                        f"Binary file {path} size not divisible by 30, using as 1D"
                    )
                    self._profiles[profile_key] = raw_data
                self._profiles[f"{profile_key}_is_legacy_bin"] = True
            elif "openen" in path.name:
                # Text format: RNAplfold output
                # We need to parse this into a float32 array
                logger.info(f"Parsing legacy text file: {path}")
                self._profiles[profile_key] = self._parse_openen_text(path)

            logger.info(
                f"Loaded accessibility profile for {chrom} {strand} from {path}"
            )

        profile = self._profiles[profile_key]
        is_legacy_bin = self._profiles.get(f"{profile_key}_is_legacy_bin", False)
        seq_len = len(profile)

        # Check bounds
        if start < 0 or end > seq_len:
            raise AccessibilityError(
                f"Query {start}-{end} out of bounds for {chrom} {strand} (len {seq_len}, path={self._find_profile(chrom, strand)})"
            )

        # Retrieve and decode if needed
        data = profile[start:end]
        if is_legacy_bin:
            return data.astype(np.float32) / 10.0
        return data

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
