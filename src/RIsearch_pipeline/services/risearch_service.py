"""Service for executing RIsearch binary and validating siRNA inputs.

Provides a clean interface to run RIsearch searches, index targets, and validate
siRNA FASTA files before processing. Designed for both single-siRNA and batch
multi-siRNA workflows.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

from loguru import logger


class RIsearchError(Exception):
    """Raised when RIsearch execution fails."""

    pass


class RIsearchService:
    """Wrapper for RIsearch binary execution.

    Encapsulates index creation, search execution, and input validation.
    Supports both single and multi-siRNA FASTA inputs.

    Attributes:
        binary_path: Path to the RIsearch binary.
    """

    def __init__(self, binary_path: Optional[Path] = None):
        """Initialize with optional custom binary path.

        Args:
            binary_path: Path to RIsearch binary. If None, uses default location
                         at ~/.cargo/bin/RIsearch.
        """
        if binary_path is None:
            binary_path = Path.home() / ".cargo" / "bin" / "RIsearch"

        self.binary_path = binary_path

        if not self.binary_path.exists():
            logger.warning(f"RIsearch binary not found at {self.binary_path}")

    def validate_sirna_fasta(self, path: Path) -> List[str]:
        """Validate siRNA FASTA and return list of sequence IDs.

        Checks for:
        - File existence
        - Valid FASTA format
        - Duplicate sequence IDs (raises error if found)

        Args:
            path: Path to siRNA FASTA file.

        Returns:
            List of unique sequence IDs in order of appearance.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If duplicate IDs are found.
        """
        from Bio import SeqIO

        if not path.exists():
            raise FileNotFoundError(f"siRNA FASTA file not found: {path}")

        ids: List[str] = []
        seen: set[str] = set()

        for record in SeqIO.parse(path, "fasta"):
            seq_id = record.id
            if seq_id in seen:
                raise ValueError(
                    f"Duplicate siRNA ID found: '{seq_id}'. "
                    f"Each siRNA must have a unique identifier."
                )
            seen.add(seq_id)
            ids.append(seq_id)

        logger.info(f"Validated {len(ids)} siRNA(s) from {path.name}")
        return ids

    def index_target(
        self,
        target_path: Path,
        index_path: Optional[Path] = None,
    ) -> Path:
        """Create or retrieve RIsearch index for a target FASTA.

        If index_path is None, creates a temporary index file.
        If index already exists and is newer than target, reuses it.

        Args:
            target_path: Path to target FASTA (genome/transcriptome).
            index_path: Optional path for index output. If None, uses temp file.

        Returns:
            Path to the created/existing index.

        Raises:
            RIsearchError: If indexing fails.
            FileNotFoundError: If target FASTA doesn't exist.
        """
        if not target_path.exists():
            raise FileNotFoundError(f"Target FASTA not found: {target_path}")

        # Generate default index path if not provided
        if index_path is None:
            # Use target name with .idx suffix in same directory
            index_path = target_path.with_suffix(".idx")

        # Check if index exists and is up-to-date
        if index_path.exists():
            if index_path.stat().st_mtime > target_path.stat().st_mtime:
                logger.info(f"Reusing existing index: {index_path}")
                return index_path
            else:
                logger.info(f"Index outdated, rebuilding: {index_path}")

        # Build index
        cmd = [str(self.binary_path), "index", str(target_path), str(index_path)]
        logger.debug(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )
            logger.info(
                f"Created index: {index_path} ({index_path.stat().st_size} bytes)"
            )
            return index_path

        except subprocess.CalledProcessError as e:
            raise RIsearchError(f"RIsearch index failed: {e.stderr}") from e

    def run_search(
        self,
        query_path: Path,
        index_path: Path,
        output_path: Optional[Path] = None,
        seed_length: str = "6",
        extension_length: int = 20,
        energy_threshold: float = -10.0,
        threads: Optional[int] = None,
    ) -> Path:
        """Run RIsearch search and return path to output.

        Args:
            query_path: Path to siRNA FASTA (single or multiple sequences).
            index_path: Path to pre-built RIsearch index.
            output_path: Path for output file. If None, uses stdout to temp file.
            seed_length: Seed length or range (e.g., "6" or "6:12").
            extension_length: Maximum extension length on each side.
            energy_threshold: Energy cutoff in kcal/mol (default -10.0).
            threads: Number of threads. If None, uses RIsearch default.

        Returns:
            Path to the output file containing predictions.

        Raises:
            RIsearchError: If search fails.
            FileNotFoundError: If query or index don't exist.
        """
        if not query_path.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_path}")
        if not index_path.exists():
            raise FileNotFoundError(f"Index not found: {index_path}")

        # Use temp file if no output path specified
        if output_path is None:
            output_fd, output_path_str = tempfile.mkstemp(
                suffix=".out", prefix="risearch_"
            )
            os.close(output_fd)
            output_path = Path(output_path_str)

        # Build command
        cmd = [
            str(self.binary_path),
            "search",
            "-q",
            str(query_path),
            "-i",
            str(index_path),
            "-o",
            str(output_path),
            "-s",
            seed_length,
            "-l",
            str(extension_length),
            "-e",
            str(energy_threshold),
        ]

        if threads is not None:
            cmd.extend(["-t", str(threads)])

        logger.debug(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )
            logger.info(f"Search complete: {output_path}")
            return output_path

        except subprocess.CalledProcessError as e:
            raise RIsearchError(f"RIsearch search failed: {e.stderr}") from e

    def search_single_sirna(
        self,
        query_path: Path,
        target_path: Path,
    ) -> Tuple[float, int, int, str]:
        """Run RIsearch for a single siRNA and return best hit.

        Convenience method for on-target calculations. Creates temporary
        index, runs search, and returns the best (lowest energy) result.

        Args:
            query_path: Path to single-siRNA FASTA.
            target_path: Path to target sequence FASTA.

        Returns:
            Tuple of (energy, start, end, strand) for the best hit.
            Returns (0.0, 0, 0, "+") if no hits found.
        """
        with tempfile.TemporaryDirectory(prefix="risearch_") as tmpdir:
            tmpdir_path = Path(tmpdir)
            index_path = tmpdir_path / "target.idx"
            output_path = tmpdir_path / "results.out"

            try:
                # Index and search
                self.index_target(target_path, index_path)
                self.run_search(query_path, index_path, output_path)

                # Parse results
                return self._parse_best_hit(output_path)

            except RIsearchError as e:
                logger.error(f"RIsearch failed: {e}")
                return 0.0, 0, 0, "+"

    def _parse_best_hit(self, output_path: Path) -> Tuple[float, int, int, str]:
        """Parse RIsearch output and return best (lowest energy) hit.

        Args:
            output_path: Path to RIsearch output file.

        Returns:
            Tuple of (energy, start, end, strand) for best hit.
            Returns (0.0, 0, 0, "+") if no valid hits.
        """
        best_hit: Optional[Tuple[float, int, int, str]] = None

        with open(output_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                # Format: QueryID, QStart, QEnd, TargetID, TStart, TEnd, Strand, Energy
                if len(parts) >= 8:
                    try:
                        energy = float(parts[7])
                        t_start = int(parts[4])
                        t_end = int(parts[5])
                        strand = parts[6]

                        if best_hit is None or energy < best_hit[0]:
                            best_hit = (energy, t_start, t_end, strand)
                    except ValueError:
                        continue

        if best_hit:
            return best_hit

        logger.warning(f"No valid hits found in {output_path}")
        return 0.0, 0, 0, "+"
