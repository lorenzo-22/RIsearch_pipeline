"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl
import numpy as np
from bisect import bisect_left, bisect_right


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
            # Transcriptome-Wide Mode - direct ID join
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

        else:
            # Genome-Wide Mode - interval containment via binary search
            return self._interval_containment(risearch_df, transcriptome_df)

    def _interval_containment(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Efficient interval containment using sorted binary search.

        For each prediction, find all transcripts that contain it.
        Uses numpy arrays and binary search for O(n log m) complexity
        instead of O(n * m) cross-join.
        """
        results = []

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

            # Build sorted index on transcript starts and ends
            chunk_result = self._binary_search_containment(preds, trans)

            if chunk_result is not None and chunk_result.height > 0:
                results.append(chunk_result)

        if not results:
            return self._empty_result_schema(risearch_df)

        return pl.concat(results, how="diagonal")

    def _binary_search_containment(
        self,
        preds: pl.DataFrame,
        trans: pl.DataFrame,
    ) -> pl.DataFrame:
        """Find all transcripts containing each prediction using binary search.

        Algorithm:
        1. Sort transcripts by start position
        2. For each prediction, binary search to find transcripts starting <= pred.start
        3. Among those, filter for transcripts ending >= pred.end
        """
        # Convert to numpy for fast iteration
        trans_sorted = trans.sort("start")
        trans_starts = trans_sorted["start"].to_numpy()
        trans_ends = trans_sorted["end"].to_numpy()
        trans_gene_ids = trans_sorted["gene_id"].to_list()
        trans_transcript_ids = trans_sorted["transcript_id"].to_list()
        trans_exp_values = trans_sorted["exp_value"].to_numpy()

        # Collect results
        result_rows = []
        pred_data = preds.to_dicts()

        for pred in pred_data:
            pred_start = pred["start"]
            pred_end = pred["end"]

            # Find transcripts that could contain this prediction
            # Transcript must start at or before prediction start
            right_idx = bisect_right(trans_starts, pred_start)

            # Check all transcripts starting before pred_start
            for i in range(right_idx):
                # Transcript must end at or after prediction end
                if trans_ends[i] >= pred_end:
                    # This transcript contains the prediction
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
