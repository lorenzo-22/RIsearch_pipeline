"""Pure core for accessibility computation — returns in-memory DataFrames.

Folds each sequence in the FASTA with ViennaRNA and returns one Polars DataFrame
per chromosome (schema ``[position, strand, u1..u{unpaired_prob}]`` — the same
schema the CLI streams to ``{chrom}.accessibility.parquet``), holding nothing on
disk. For genome-scale inputs prefer the CLI's streaming path; the in-memory dict
keeps every profile resident at once.
"""

from pathlib import Path

import numpy as np
import polars as pl

from riot.services.accessibility import _reverse_complement, fold_sequence
from riot.services.helpers import read_fasta


def _profile_to_df(
    profile: np.ndarray, strand: str, unpaired_prob: int
) -> pl.DataFrame:
    """Turn a [seq_len, u] opening-energy array into the parquet-equivalent schema."""
    n_pos = 0 if profile.size == 0 else profile.shape[0]
    columns: dict[str, pl.Series] = {
        "position": pl.Series("position", np.arange(1, n_pos + 1, dtype=np.int32)),
        "strand": pl.Series("strand", [strand] * n_pos, dtype=pl.Utf8),
    }
    for u in range(1, unpaired_prob + 1):
        col = (
            profile[:, u - 1].astype(np.float32)
            if n_pos
            else np.array([], dtype=np.float32)
        )
        columns[f"u{u}"] = pl.Series(f"u{u}", col, dtype=pl.Float32)
    return pl.DataFrame(columns)


def compute_accessibility(
    genome: Path,
    window_size: int = 80,
    max_span: int = 40,
    unpaired_prob: int = 30,
    temperature: float = 37.0,
) -> dict[str, pl.DataFrame]:
    """Compute per-chromosome accessibility profiles in memory (no files written).

    Returns a dict mapping chromosome name to a Polars DataFrame with columns
    ``[position, strand, u1..u{unpaired_prob}]``, both strands stacked (``+`` then
    ``-``), matching the on-disk parquet schema.
    """
    genome = Path(genome)
    if not genome.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome}")

    profiles: dict[str, pl.DataFrame] = {}
    for chrom, sequence in read_fasta(genome):
        plus = fold_sequence(
            sequence,
            window_size=window_size,
            max_span=max_span,
            unpaired_prob=unpaired_prob,
            temperature=temperature,
        )
        minus = fold_sequence(
            _reverse_complement(sequence),
            window_size=window_size,
            max_span=max_span,
            unpaired_prob=unpaired_prob,
            temperature=temperature,
        )
        profiles[chrom] = pl.concat(
            [
                _profile_to_df(plus, "+", unpaired_prob),
                _profile_to_df(minus, "-", unpaired_prob),
            ]
        )
    return profiles
