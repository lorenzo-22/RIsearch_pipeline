"""Pure core for RIsearch index/search.

Wraps :class:`RIsearchService` without any CLI side-effects (no stdout, no file
writes beyond the index binary itself). Returns in-memory objects.
"""

from pathlib import Path
from typing import Optional

import polars as pl

from riot.services.risearch_service import RIsearchService


def build_index(target: Path, output: Optional[Path] = None) -> Path:
    """Build (or reuse) a RIsearch index for a target FASTA.

    The index is a binary on-disk artifact consumed by :func:`run_search`, so this
    returns its :class:`Path` rather than in-memory data. Reuses an existing index
    if it is newer than the target.
    """
    target = Path(target)
    output = Path(output) if output is not None else None
    return RIsearchService().index_target(target, output)


def run_search(
    query: Path,
    index: Path,
    target: Optional[Path] = None,
    seed_length: int = 6,
    max_extension: int = 20,
    energy_threshold: float = -10.0,
) -> pl.DataFrame:
    """Run a RIsearch search and return hits as a Polars DataFrame.

    Columns: ``sirna_id, chrom, start, end, strand, energy``. ``target`` is required
    when the index was built outside this process (needed to resolve target names).
    """
    return RIsearchService().run_search(
        query_path=Path(query),
        index_path=Path(index),
        target_fasta=Path(target) if target is not None else None,
        seed_length=seed_length,
        max_extension=max_extension,
        energy_threshold=energy_threshold,
    )
