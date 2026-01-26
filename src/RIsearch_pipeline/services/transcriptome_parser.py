"""Service for parsing transcriptome annotation files (GTF/GFF/BED)."""

from pathlib import Path
from typing import Optional

import polars as pl

from RIsearch_pipeline.models import TRANSCRIPTOME_SCHEMA


class TranscriptomeParser:
    """Parser for transcriptome files to extract gene/transcript locations and expression."""

    def load_gtf(
        self, path: Path, feature: str = "exon", score_col: str = "RPKM"
    ) -> pl.DataFrame:
        """Load GTF/GFF or BED file and extract fields.

        Args:
            path: Path to .gtf, .gff, or .bed file
            feature: Feature type to filter by (GTF only, e.g. 'exon')
            score_col: Attribute for score (GTF) or use score column (BED)

        Returns:
            DataFrame with columns: chrom, start, end, strand, gene_id, transcript_id, exp_value
        """
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Transcriptome file not found: {path}")

        # Detect format by suffix
        if path.suffix.lower() in [".bed", ".bed.gz"]:
            return self._load_bed(path)
        else:
            return self._load_gtf_impl(path, feature, score_col)

    def _load_bed(self, path: Path) -> pl.DataFrame:
        """Load BED file (7 columns expected)."""
        # User format: chrom, start, end, transcript_id, score_placeholder, strand, exp_value
        headers = [
            "chrom",
            "start",
            "end",
            "transcript_id",
            "score_placeholder",
            "strand",
            "exp_value",
        ]

        # Read 7 columns
        df = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            columns=range(7),
            new_columns=headers,
            truncate_ragged_lines=True,
        )

        # Select and cast
        df = df.select(
            [
                pl.col("chrom"),
                pl.col("start"),
                pl.col("end"),
                pl.col("strand"),
                pl.col("transcript_id").alias(
                    "gene_id"
                ),  # Use transcript_id as gene_id for BED
                pl.col("transcript_id"),
                pl.col("exp_value").cast(pl.Float32, strict=False).fill_null(0.0),
            ]
        )

        return df

    def _load_gtf_impl(self, path: Path, feature: str, score_col: str) -> pl.DataFrame:
        """Internal GTF loading logic."""
        # Load raw GTF (9 columns)
        # We rename them immediately to match schema names
        df = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=list(TRANSCRIPTOME_SCHEMA.keys()),
            schema_overrides=TRANSCRIPTOME_SCHEMA,
            truncate_ragged_lines=True,
        )

        # Filter by feature type if needed
        if feature:
            df = df.filter(pl.col("feature") == feature)

        # Regex patterns to extract attributes
        # Looking for: key "value"; or key value;
        # Example: gene_id "gene_1"; transcript_id "transcript_21"; RPKM "1000"

        # Helper string expression for extracting value by key
        def extract_attr(key: str) -> pl.Expr:
            # Matches: key followed by space, then quote?, capture group, quote?, then semicolon
            # We handle both quoted and unquoted values, ending with ; or end of string
            return pl.col("attributes").str.extract(rf'{key}\s+"?([^";]+)"?', 1)

        # Extract required columns
        df = df.select(
            [
                pl.col("chrom"),
                pl.col("start"),
                pl.col("end"),
                pl.col("strand"),
                extract_attr("gene_id").alias("gene_id"),
                extract_attr("transcript_id").alias("transcript_id"),
                extract_attr(score_col)
                .cast(pl.Float32, strict=False)
                .fill_null(0.0)
                .alias("exp_value"),
            ]
        )

        return df

    def summary(self, df: pl.DataFrame) -> dict:
        """Generate summary statistics for loaded transcriptome."""
        return {
            "row_count": df.height,
            "genes": df["gene_id"].n_unique(),
            "transcripts": df["transcript_id"].n_unique(),
            "chromosomes": df["chrom"].unique().to_list(),
            "score_mean": df["exp_value"].mean(),
        }
