"""Service for parsing transcriptome annotation files (GTF/GFF/BED)."""

from pathlib import Path

import polars as pl

from RIsearch_pipeline.models import TRANSCRIPTOME_SCHEMA


class TranscriptomeParser:
    """Parser for transcriptome files to extract gene/transcript locations and expression."""

    def load_gtf(
        self,
        path: Path,
        feature: str = "exon",
        score_col: str = "RPKM",
        format: str = "auto",
    ) -> pl.DataFrame:
        """Load GTF/GFF or BED file and extract fields.

        Returns DataFrame with: chrom, start, end, strand, gene_id, transcript_id, exp_value
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Transcriptome file not found: {path}")

        fmt = format.lower()
        if fmt == "auto":
            if path.suffix.lower() in [".bed", ".bed.gz"]:
                return self._load_bed(path)
            return self._load_gtf_impl(path, feature, score_col)
        elif fmt == "bed6":
            return self._load_bed(path, force_bed6=True)
        elif fmt == "bed7":
            return self._load_bed(path, force_bed7=True)
        else:
            return self._load_gtf_impl(path, feature, score_col)

    def _load_bed(
        self, path: Path, force_bed6: bool = False, force_bed7: bool = False
    ) -> pl.DataFrame:
        import gzip

        if force_bed6:
            num_cols = 6
        elif force_bed7:
            num_cols = 7
        else:
            opener = gzip.open if str(path).endswith(".gz") else open
            try:
                with opener(path, "rt") as f:
                    num_cols = len(f.readline().strip().split("\t"))
            except Exception:
                num_cols = 6

        if num_cols >= 7:
            headers = ["chrom", "start", "end", "transcript_id", "score_placeholder", "strand", "exp_value"]
            df = pl.read_csv(
                path, separator="\t", has_header=False,
                columns=range(7), new_columns=headers, truncate_ragged_lines=True,
            )
            return df.select([
                pl.col("chrom"),
                pl.col("start"),
                pl.col("end"),
                pl.col("strand"),
                pl.col("transcript_id").alias("gene_id"),
                pl.col("transcript_id"),
                pl.col("exp_value").cast(pl.Float32, strict=False).fill_null(0.0),
            ])
        else:
            headers = ["chrom", "start", "end", "transcript_id", "score", "strand"]
            df = pl.read_csv(
                path, separator="\t", has_header=False,
                columns=range(6), new_columns=headers, truncate_ragged_lines=True,
                schema_overrides={h: pl.Utf8 for h in headers},
            )
            return df.select([
                pl.col("chrom"),
                pl.col("start").cast(pl.Int64),
                pl.col("end").cast(pl.Int64),
                pl.col("strand"),
                pl.col("transcript_id").alias("gene_id"),
                pl.col("transcript_id"),
                pl.col("score").cast(pl.Float32, strict=False).fill_null(0.0).alias("exp_value"),
            ])

    def _load_gtf_impl(self, path: Path, feature: str, score_col: str) -> pl.DataFrame:
        df = pl.read_csv(
            path, separator="\t", has_header=False,
            comment_prefix="#", new_columns=list(TRANSCRIPTOME_SCHEMA.keys()),
            schema_overrides=TRANSCRIPTOME_SCHEMA, truncate_ragged_lines=True,
        )

        if feature:
            df = df.filter(pl.col("feature") == feature)

        def extract_attr(key: str) -> pl.Expr:
            return pl.col("attributes").str.extract(rf'{key}\s+"?([^";]+)"?', 1)

        return df.select([
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col("strand"),
            extract_attr("gene_id").alias("gene_id"),
            extract_attr("transcript_id").alias("transcript_id"),
            extract_attr(score_col).cast(pl.Float32, strict=False).fill_null(0.0).alias("exp_value"),
        ])

    def summary(self, df: pl.DataFrame) -> dict:
        return {
            "row_count": df.height,
            "genes": df["gene_id"].n_unique(),
            "transcripts": df["transcript_id"].n_unique(),
            "chromosomes": df["chrom"].unique().to_list(),
            "score_mean": df["exp_value"].mean(),
        }
