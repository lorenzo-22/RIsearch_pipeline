import logging
from pathlib import Path
from typing import Dict
import numpy as np

# Try to import ViennaRNA
try:
    import RNA

    HAS_VIENNA_BINDINGS = True
except ImportError:
    HAS_VIENNA_BINDINGS = False

logger = logging.getLogger(__name__)


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
                results[f"{chrom}_{strand}"] = out_path

        return results

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
            # Try to load
            path = self.get_profile_path(chrom, strand)
            if not path.exists():
                raise AccessibilityError(
                    f"Profile for {chrom} {strand} not found at {path}"
                )

            # Memory map
            self._profiles[profile_key] = np.load(path, mmap_mode="r")

        profile = self._profiles[profile_key]
        seq_len = len(profile)

        # Check bounds
        if start < 0 or end > seq_len:
            raise AccessibilityError(
                f"Query {start}-{end} out of bounds for {chrom} {strand} (len {seq_len})"
            )

        # Profiles are already oriented correctly (minus strand reversed during compute)
        return profile[start:end]
