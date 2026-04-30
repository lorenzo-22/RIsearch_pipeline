import tempfile
from pathlib import Path
from typing import Optional

import polars as pl
from loguru import logger


class RIsearchError(Exception):
    pass


class RIsearchService:
    """Wrapper for the risearch PyO3 bindings (in-process, no subprocess)."""

    def __init__(self) -> None:
        self._target_registry: dict[str, Path] = {}

    def validate_sirna_fasta(self, path: Path) -> list[str]:
        """Return unique siRNA IDs; raise ValueError on duplicates."""
        from Bio import SeqIO

        if not path.exists():
            raise FileNotFoundError(f"siRNA FASTA file not found: {path}")

        ids: list[str] = []
        seen: set[str] = set()
        for record in SeqIO.parse(path, "fasta"):
            if record.id in seen:
                raise ValueError(f"Duplicate siRNA ID found: '{record.id}'.")
            seen.add(record.id)
            ids.append(record.id)

        logger.info(f"Validated {len(ids)} siRNA(s) from {path.name}")
        return ids

    def index_target(self, target_path: Path, index_path: Optional[Path] = None) -> Path:
        """Create or reuse RIsearch index. Register target for name resolution."""
        import risearch

        if not target_path.exists():
            raise FileNotFoundError(f"Target FASTA not found: {target_path}")

        if index_path is None:
            index_path = target_path.with_suffix(".idx")

        self._target_registry[str(index_path)] = target_path

        if index_path.exists() and index_path.stat().st_mtime > target_path.stat().st_mtime:
            logger.info(f"Reusing existing index: {index_path}")
            return index_path

        try:
            risearch.index(target_path, index_path)
            logger.info(f"Created index: {index_path} ({index_path.stat().st_size} bytes)")
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
        """Run RIsearch and return results as a DataFrame.

        target_fasta required when index was not built in this session via index_target().
        """
        import risearch

        if not query_path.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_path}")
        if not index_path.exists():
            raise FileNotFoundError(f"Index not found: {index_path}")

        resolved_target = target_fasta or self._target_registry.get(str(index_path))
        if resolved_target is None:
            raise RIsearchError(
                "Cannot resolve target sequence names: pass target_fasta= or "
                "build the index via index_target() first."
            )

        try:
            store = risearch.TargetStore.open(index_path)
            raw = risearch.search(
                query_path, store,
                seed_length=seed_length,
                max_extension=max_extension,
                energy_threshold=energy_threshold,
            )
        except Exception as e:
            raise RIsearchError(f"RIsearch search failed: {e}") from e

        if raw.is_empty():
            logger.info("Search complete: 0 hits")
            return pl.DataFrame(schema={
                "sirna_id": pl.Utf8, "chrom": pl.Utf8,
                "start": pl.Int32, "end": pl.Int32,
                "strand": pl.Utf8, "energy": pl.Float32,
            })

        query_names = _fasta_names(query_path)
        target_names = _fasta_names(resolved_target)

        return (
            raw.with_columns([
                pl.Series("sirna_id", [query_names[i] for i in raw["query_idx"].to_list()]),
                pl.Series("chrom", [target_names[i] for i in raw["target_idx"].to_list()]),
            ])
            .rename({"t_start": "start", "t_end": "end"})
            .select(["sirna_id", "chrom", "start", "end", "strand", "energy"])
            .cast({"start": pl.Int32, "end": pl.Int32, "energy": pl.Float32})
        )

    def self_hybridization_emin(self, sequence: str, sirna_id: str = "query") -> float:
        """Compute E_min via self-hybridization (seed = len-1, threshold = 0.0).

        Replicates old pipeline: ``risearch2.x -q siRNA.fa -i siRNA.pksuf -s len-1 -e 0``
        """
        seq_dna = sequence.upper().replace("U", "T")
        with tempfile.TemporaryDirectory(prefix="risearch_self_") as tmpdir:
            fasta_path = Path(tmpdir) / "sirna.fa"
            index_path = Path(tmpdir) / "sirna.idx"
            fasta_path.write_text(f">{sirna_id}\n{seq_dna}\n")
            self.index_target(fasta_path, index_path)
            df = self.run_search(
                fasta_path, index_path, target_fasta=fasta_path,
                seed_length=len(seq_dna) - 1, energy_threshold=0.0,
            )
            if df.is_empty():
                logger.warning(f"No self-hybridisation hits for {sirna_id}, using E_min=0.0")
                return 0.0
            emin = float(df["energy"].min())
            logger.debug(f"Self-hyb E_min {sirna_id}: {emin:.4f} kcal/mol")
            return emin

    def self_hybridization_emin_batch(self, fasta_path: Path) -> dict[str, float]:
        """Compute self-hybridisation E_min for all siRNAs in a FASTA file."""
        from Bio import SeqIO

        if not fasta_path.exists():
            raise FileNotFoundError(f"siRNA FASTA not found: {fasta_path}")

        records = list(SeqIO.parse(fasta_path, "fasta"))
        logger.info(f"Computing self-hybridisation E_min for {len(records)} siRNA(s)...")
        result = {r.id: self.self_hybridization_emin(str(r.seq), r.id) for r in records}
        logger.info(f"Self-hybridisation complete: {len(result)} siRNA(s) processed")
        return result

    def search_single_sirna(self, query_path: Path, target_path: Path) -> tuple[float, int, int, str]:
        """Run RIsearch for a single siRNA; return best (energy, start, end, strand)."""
        with tempfile.TemporaryDirectory(prefix="risearch_") as tmpdir:
            index_path = Path(tmpdir) / "target.idx"
            try:
                self.index_target(target_path, index_path)
                df = self.run_search(query_path, index_path)
                if df.is_empty():
                    logger.warning("No valid hits found")
                    return 0.0, 0, 0, "+"
                best = df.sort("energy").row(0, named=True)
                return float(best["energy"]), int(best["start"]), int(best["end"]), str(best["strand"])
            except RIsearchError as e:
                logger.error(f"RIsearch failed: {e}")
                return 0.0, 0, 0, "+"


def _fasta_names(path: Path) -> list[str]:
    from Bio import SeqIO
    return [r.id for r in SeqIO.parse(path, "fasta")]
