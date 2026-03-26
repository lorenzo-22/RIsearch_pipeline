"""Service for executing RIsearch and validating siRNA inputs.

Provides a clean interface to run RIsearch searches, index targets, and validate
siRNA FASTA files before processing. Uses the risearch Python bindings (PyO3/maturin)
for in-process execution — no subprocess or intermediate TSV files.
"""

import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
from loguru import logger


class RIsearchError(Exception):
    """Raised when RIsearch execution fails."""

    pass


class RIsearchService:
    """Wrapper for the risearch Python bindings.

    Encapsulates index creation, search execution, and input validation.
    Supports both single and multi-siRNA FASTA inputs.

    An internal registry maps each index path to its source target FASTA so
    that ``run_search()`` can resolve integer target indices to sequence names.
    """

    def __init__(self) -> None:
        # Maps str(index_path) → target_fasta path for name resolution
        self._target_registry: Dict[str, Path] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

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

        If index_path is None, creates an index next to the target file.
        If index already exists and is newer than target, reuses it.
        Registers the target_path so run_search() can resolve target names.

        Args:
            target_path: Path to target FASTA (genome/transcriptome).
            index_path: Optional path for index output.

        Returns:
            Path to the created/existing index.

        Raises:
            RIsearchError: If indexing fails.
            FileNotFoundError: If target FASTA doesn't exist.
        """
        import risearch

        if not target_path.exists():
            raise FileNotFoundError(f"Target FASTA not found: {target_path}")

        if index_path is None:
            index_path = target_path.with_suffix(".idx")

        # Always register so run_search() can resolve target names
        self._target_registry[str(index_path)] = target_path

        if index_path.exists():
            if index_path.stat().st_mtime > target_path.stat().st_mtime:
                logger.info(f"Reusing existing index: {index_path}")
                return index_path
            else:
                logger.info(f"Index outdated, rebuilding: {index_path}")

        try:
            risearch.index(target_path, index_path)
            logger.info(
                f"Created index: {index_path} ({index_path.stat().st_size} bytes)"
            )
            return index_path
        except Exception as e:
            raise RIsearchError(f"RIsearch index failed: {e}") from e

    def run_search(
        self,
        query_path: Path,
        index_path: Path,
        target_fasta: Optional[Path] = None,
        seed_length: int = 6,
        max_extension: int = 20,
        energy_threshold: float = -10.0,
    ) -> pl.DataFrame:
        """Run RIsearch search and return results as a DataFrame.

        Args:
            query_path: Path to siRNA FASTA (single or multiple sequences).
            index_path: Path to pre-built RIsearch index.
            target_fasta: Path to the target FASTA used to build the index.
                          Required when the index was not created in this
                          session via index_target(). Needed for name resolution.
            seed_length: Seed length integer (default 6).
            max_extension: Max extension length on each side (default 20).
            energy_threshold: Energy cutoff in kcal/mol (default -10.0).

        Returns:
            Polars DataFrame with columns: sirna_id, chrom, start, end, strand, energy.

        Raises:
            RIsearchError: If search fails or target names cannot be resolved.
            FileNotFoundError: If query or index don't exist.
        """
        import risearch

        if not query_path.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_path}")
        if not index_path.exists():
            raise FileNotFoundError(f"Index not found: {index_path}")

        # Resolve target FASTA for name lookup
        resolved_target = target_fasta or self._target_registry.get(str(index_path))
        if resolved_target is None:
            raise RIsearchError(
                "Cannot resolve target sequence names: pass target_fasta= or "
                "build the index via index_target() first."
            )

        try:
            store = risearch.TargetStore.open(index_path)
            raw = risearch.search(
                query_path,
                store,
                seed_length=seed_length,
                max_extension=max_extension,
                energy_threshold=energy_threshold,
            )
        except Exception as e:
            raise RIsearchError(f"RIsearch search failed: {e}") from e

        if raw.is_empty():
            logger.info("Search complete: 0 hits")
            return pl.DataFrame(
                schema={
                    "sirna_id": pl.Utf8,
                    "chrom": pl.Utf8,
                    "start": pl.Int32,
                    "end": pl.Int32,
                    "strand": pl.Utf8,
                    "energy": pl.Float32,
                }
            )

        # Build name lookup arrays
        query_names = _fasta_names(query_path)
        target_names = _fasta_names(resolved_target)

        sirna_ids = [query_names[i] for i in raw["query_idx"].to_list()]
        chroms = [target_names[i] for i in raw["target_idx"].to_list()]

        df = (
            raw.with_columns([
                pl.Series("sirna_id", sirna_ids),
                pl.Series("chrom", chroms),
            ])
            .rename({"t_start": "start", "t_end": "end"})
            .select(["sirna_id", "chrom", "start", "end", "strand", "energy"])
            .cast({"start": pl.Int32, "end": pl.Int32, "energy": pl.Float32})
        )

        logger.info(f"Search complete: {df.height} hits")
        return df

    def search_single_sirna(
        self,
        query_path: Path,
        target_path: Path,
    ) -> Tuple[float, int, int, str]:
        """Run RIsearch for a single siRNA and return best hit.

        Convenience method for on-target calculations. Creates a temporary
        index, runs search, and returns the best (lowest energy) result.

        Args:
            query_path: Path to single-siRNA FASTA.
            target_path: Path to target sequence FASTA.

        Returns:
            Tuple of (energy, start, end, strand) for the best hit.
            Returns (0.0, 0, 0, "+") if no hits found.
        """
        with tempfile.TemporaryDirectory(prefix="risearch_") as tmpdir:
            index_path = Path(tmpdir) / "target.idx"

            try:
                self.index_target(target_path, index_path)
                df = self.run_search(query_path, index_path)

                if df.is_empty():
                    logger.warning("No valid hits found")
                    return 0.0, 0, 0, "+"

                best = df.sort("energy").row(0, named=True)
                return (
                    float(best["energy"]),
                    int(best["start"]),
                    int(best["end"]),
                    str(best["strand"]),
                )

            except RIsearchError as e:
                logger.error(f"RIsearch failed: {e}")
                return 0.0, 0, 0, "+"


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _fasta_names(path: Path) -> List[str]:
    """Return sequence IDs from a FASTA file in order of appearance."""
    from Bio import SeqIO
    return [r.id for r in SeqIO.parse(path, "fasta")]
