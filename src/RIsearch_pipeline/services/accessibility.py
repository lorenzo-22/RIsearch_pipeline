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

    def get_profile_path(self, chrom: str) -> Path:
        return self.data_dir / f"{chrom}.access.npy"

    def compute_genome_accessibility(
        self,
        genome_path: Path,
        window_size: int = 80,
        max_span: int = 40,
        unpaired_prob: int = 30,
    ) -> Dict[str, Path]:
        """
        Compute accessibility for all sequences in the genome FASTA.

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

        # We need a FASTA parser. Using a simple generator to avoid heavy dependencies if possible,
        # but pysam is better. We'll use a simple parser for now to minimize dependnecies.
        from RIsearch_pipeline.services.helpers import read_fasta

        results = {}

        logger.info(f"Computing accessibility for genome: {genome_path}")

        for chrom, sequence in read_fasta(genome_path):
            seq_len = len(sequence)
            logger.info(f"Processing {chrom} (length {seq_len})...")

            # Create profile array (initialized to 0.0)
            # Use float32 for compact storage but decent precision
            profile = np.zeros(seq_len, dtype=np.float32)

            try:
                # Configure Model
                md = RNA.md()
                md.window_size = window_size
                md.max_bp_span = max_span

                fc = RNA.fold_compound(sequence, md)

                # Callback to capture probabilities
                # Signature: (probs, size, i, max_bg, data)
                # probs: list of probabilities. probs[k] is P(unpaired) for length k.
                # i: current position (1-based, ending position of the window)
                def prob_callback(probs, size, i, *args):
                    if i <= seq_len:
                        # Capture the probability for the specific length 'unpaired_prob'
                        # i is 1-based, so use i-1 for 0-based index
                        if (
                            unpaired_prob < len(probs)
                            and probs[unpaired_prob] is not None
                        ):
                            profile[i - 1] = probs[unpaired_prob]

                # Run sliding window
                fc.probs_window(unpaired_prob, RNA.PROBS_WINDOW_UP, prob_callback)

            except Exception as e:
                logger.error(f"Error calling ViennaRNA on {chrom}: {e}")
                raise

            # Save to disk
            out_path = self.get_profile_path(chrom)
            np.save(out_path, profile)
            results[chrom] = out_path

        return results

    def query(self, chrom: str, start: int, end: int) -> np.ndarray:
        """
        Query accessibility for a region (0-based, half-open [start, end)).

        Returns:
            Numpy array of probabilities.
        """
        if chrom not in self._profiles:
            # Try to load
            path = self.get_profile_path(chrom)
            if not path.exists():
                raise AccessibilityError(f"Profile for {chrom} not found at {path}")

            # Memory map
            self._profiles[chrom] = np.load(path, mmap_mode="r")

        profile = self._profiles[chrom]

        # Check bounds
        if start < 0 or end > len(profile):
            raise AccessibilityError(
                f"Query {start}-{end} out of bounds for {chrom} (len {len(profile)})"
            )

        return profile[start:end]
