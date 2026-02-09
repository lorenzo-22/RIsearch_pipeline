"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl
import numpy as np
from bisect import bisect_right
from typing import Iterator, Optional


class IntersectionService:
    """Service to intersect off-target predictions with genomic features."""

    def intersect(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
        mode: str = "gw",
    ) -> pl.DataFrame:
        """Filter predictions to those contained within transcriptome features.

        Args:
            risearch_df: DataFrame from RIsearchParser.
            transcriptome_df: DataFrame from TranscriptomeParser.
            mode: "gw" (genome-wide) or "tw" (transcriptome-wide).

        Returns:
            DataFrame containing intersected results.
        """
        if mode == "tw":
            return self._transcriptome_wide_join(risearch_df, transcriptome_df)
        else:
            return self._genome_wide_streaming(risearch_df, transcriptome_df)

    def _transcriptome_wide_join(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Transcriptome-Wide Mode - direct ID join."""
        joined = risearch_df.join(
            transcriptome_df,
            left_on="chrom",
            right_on="transcript_id",
            how="inner",
            suffix="_trans",
        )

        if "transcript_id" not in joined.columns:
            joined = joined.with_columns(pl.col("chrom").alias("transcript_id"))

        return joined

    def _genome_wide_streaming(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Genome-wide intersection with streaming to limit memory.

        Processes predictions in small batches and yields results immediately.
        """
        BATCH_SIZE = 50000  # Process 50K predictions at a time

        all_results = []

        # Get unique (chrom, strand) pairs
        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        for pair in chrom_strand_pairs:
            chrom = pair["chrom"]
            strand = pair["strand"]

            # Filter to this chromosome+strand
            preds = risearch_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )
            trans = transcriptome_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )

            if preds.height == 0 or trans.height == 0:
                continue

            # Build transcript index once per chromosome
            trans_sorted = trans.sort("start")
            trans_index = self._build_transcript_index(trans_sorted)

            # Process predictions in batches
            for batch_start in range(0, preds.height, BATCH_SIZE):
                batch = preds.slice(batch_start, BATCH_SIZE)
                batch_result = self._process_batch(batch, trans_index)

                if batch_result is not None and batch_result.height > 0:
                    all_results.append(batch_result)

        if not all_results:
            return self._empty_result_schema(risearch_df)

        return pl.concat(all_results, how="diagonal")

    def _build_transcript_index(self, trans_sorted: pl.DataFrame) -> dict:
        """Pre-build numpy arrays for fast binary search."""
        return {
            "starts": trans_sorted["start"].to_numpy(),
            "ends": trans_sorted["end"].to_numpy(),
            "gene_ids": trans_sorted["gene_id"].to_list(),
            "transcript_ids": trans_sorted["transcript_id"].to_list(),
            "exp_values": trans_sorted["exp_value"].to_numpy(),
        }

    def _process_batch(
        self,
        preds: pl.DataFrame,
        trans_index: dict,
    ) -> Optional[pl.DataFrame]:
        """Process a batch of predictions using binary search."""
        trans_starts = trans_index["starts"]
        trans_ends = trans_index["ends"]
        trans_gene_ids = trans_index["gene_ids"]
        trans_transcript_ids = trans_index["transcript_ids"]
        trans_exp_values = trans_index["exp_values"]

        # Collect matching rows
        result_rows = []
        pred_data = preds.to_dicts()

        for pred in pred_data:
            pred_start = pred["start"]
            pred_end = pred["end"]

            # Find transcripts that could contain this prediction
            right_idx = bisect_right(trans_starts, pred_start)

            # Check all transcripts starting at or before pred_start
            for i in range(right_idx):
                if trans_ends[i] >= pred_end:
                    result_row = pred.copy()
                    result_row["trans_start"] = int(trans_starts[i])
                    result_row["trans_end"] = int(trans_ends[i])
                    result_row["gene_id"] = trans_gene_ids[i]
                    result_row["transcript_id"] = trans_transcript_ids[i]
                    result_row["exp_value"] = float(trans_exp_values[i])
                    result_rows.append(result_row)

        if not result_rows:
            return None

        return pl.DataFrame(result_rows)

    def intersect_streaming(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
        mode: str = "gw",
    ) -> Iterator[pl.DataFrame]:
        """Streaming intersection that yields batches instead of accumulating.

        Use this for large datasets to avoid OOM.
        """
        if mode == "tw":
            yield self._transcriptome_wide_join(risearch_df, transcriptome_df)
            return

        BATCH_SIZE = 50000

        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        for pair in chrom_strand_pairs:
            chrom = pair["chrom"]
            strand = pair["strand"]

            preds = risearch_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )
            trans = transcriptome_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )

            if preds.height == 0 or trans.height == 0:
                continue

            trans_sorted = trans.sort("start")
            trans_index = self._build_transcript_index(trans_sorted)

            for batch_start in range(0, preds.height, BATCH_SIZE):
                batch = preds.slice(batch_start, BATCH_SIZE)
                batch_result = self._process_batch(batch, trans_index)

                if batch_result is not None and batch_result.height > 0:
                    yield batch_result

    def _empty_result_schema(self, risearch_df: pl.DataFrame) -> pl.DataFrame:
        """Return empty DataFrame with expected output schema."""
        return risearch_df.head(0).with_columns(
            [
                pl.lit(None).cast(pl.Int64).alias("trans_start"),
                pl.lit(None).cast(pl.Int64).alias("trans_end"),
                pl.lit(None).cast(pl.Utf8).alias("gene_id"),
                pl.lit(None).cast(pl.Utf8).alias("transcript_id"),
                pl.lit(None).cast(pl.Float32).alias("exp_value"),
            ]
        )
