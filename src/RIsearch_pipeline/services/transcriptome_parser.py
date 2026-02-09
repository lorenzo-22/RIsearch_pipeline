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
        """Load BED file (supports BED6 and BED7 formats).

        BED6: chrom, start, end, name, score, strand
        BED7: chrom, start, end, name, score, strand, exp_value
        """
        # First, detect number of columns by reading first line
        import gzip

        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f:
            first_line = f.readline().strip()
            num_cols = len(first_line.split("\t"))

        if num_cols >= 7:
            # BED7+: chrom, start, end, transcript_id, score_placeholder, strand, exp_value
            headers = [
                "chrom",
                "start",
                "end",
                "transcript_id",
                "score_placeholder",
                "strand",
                "exp_value",
            ]
            df = pl.read_csv(
                path,
                separator="\t",
                has_header=False,
                columns=range(7),
                new_columns=headers,
                truncate_ragged_lines=True,
            )
            df = df.select(
                [
                    pl.col("chrom"),
                    pl.col("start"),
                    pl.col("end"),
                    pl.col("strand"),
                    pl.col("transcript_id").alias("gene_id"),
                    pl.col("transcript_id"),
                    pl.col("exp_value").cast(pl.Float32, strict=False).fill_null(0.0),
                ]
            )
        else:
            # BED6: chrom, start, end, name, score, strand
            headers = ["chrom", "start", "end", "transcript_id", "score", "strand"]
            df = pl.read_csv(
                path,
                separator="\t",
                has_header=False,
                columns=range(6),
                new_columns=headers,
                truncate_ragged_lines=True,
                schema_overrides={
                    h: pl.Utf8 for h in headers
                },  # All strings to avoid type inference issues
            )
            # Use score column as exp_value for BED6
            df = df.select(
                [
                    pl.col("chrom"),
                    pl.col("start").cast(pl.Int64),
                    pl.col("end").cast(pl.Int64),
                    pl.col("strand"),
                    pl.col("transcript_id").alias("gene_id"),
                    pl.col("transcript_id"),
                    pl.col("score")
                    .cast(pl.Float32, strict=False)
                    .fill_null(0.0)
                    .alias("exp_value"),
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
