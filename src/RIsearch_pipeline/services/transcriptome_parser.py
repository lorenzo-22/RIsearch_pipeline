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
        """Load GTF file and extract fields.

        Args:
            path: Path to .gtf or .gff file
            feature: Feature type to filter by (e.g. 'exon', 'transcript')
            score_col: Name of the attribute to use as expression score (e.g. 'RPKM', 'FPKM')

        Returns:
            DataFrame with columns: chrom, start, end, strand, gene_id, transcript_id, score
        """
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Transcriptome file not found: {path}")

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
