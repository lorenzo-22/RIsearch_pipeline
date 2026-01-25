"""Polars schema definitions for RIsearch2 data structures."""

import polars as pl

# RIsearch2 output TSV schema (tab-separated, no header)
# Columns: siRNA_id, si_start, si_end, chrom, start, end, strand, energy, structure, query_seq, left_flank, right_flank
RISEARCH_COLUMNS = [
    "sirna_id",
    "chrom",
    "start",
    "end",
    "strand",
    "energy",
]

RISEARCH_SCHEMA = {
    "sirna_id": pl.Utf8,
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "strand": pl.Utf8,
    "energy": pl.Float32,
}

TRANSCRIPTOME_SCHEMA = {
    "chrom": pl.Utf8,
    "source": pl.Utf8,
    "feature": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "score_raw": pl.Utf8,  # GTF column 6 is often ".", so read as string first
    "strand": pl.Utf8,
    "frame": pl.Utf8,
    "attributes": pl.Utf8,
}
