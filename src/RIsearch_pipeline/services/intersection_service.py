"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl


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
            # Transcriptome-Wide Mode
            # RIsearch 'chrom' column is actually Transcript ID.
            # Join on transcript ID directly (no Cartesian explosion).
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
            # Genome-Wide Mode - Per-chromosome range-join
            return self._range_join_per_chrom(risearch_df, transcriptome_df)

    def _range_join_per_chrom(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Efficient range-join by processing per chromosome+strand.

        Strategy:
        1. Process each (chrom, strand) separately to limit memory
        2. For each chunk, use interval overlap detection
        3. A prediction is contained if: pred.start >= trans.start AND pred.end <= trans.end

        Memory complexity: O(max_chunk_size) instead of O(n*m)
        """
        results = []

        # Get unique (chrom, strand) pairs from RIsearch data
        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        # Prepare transcriptome with renamed columns to avoid collision
        trans_renamed = transcriptome_df.rename(
            {
                "start": "trans_start",
                "end": "trans_end",
            }
        )

        for pair in chrom_strand_pairs:
            chrom = pair["chrom"]
            strand = pair["strand"]

            # Filter to this chromosome+strand
            pred_chunk = risearch_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )

            trans_chunk = trans_renamed.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )

            if pred_chunk.height == 0 or trans_chunk.height == 0:
                continue

            # For smaller chunks, cross-join + filter is acceptable
            # This limits peak memory to (chunk_pred * chunk_trans) per chromosome
            if pred_chunk.height * trans_chunk.height < 100_000_000:
                # Small enough for cross-join approach
                chunk_result = self._cross_join_filter(pred_chunk, trans_chunk)
            else:
                # Large chunk - use sorted interval approach
                chunk_result = self._sorted_interval_join(pred_chunk, trans_chunk)

            if chunk_result.height > 0:
                results.append(chunk_result)

        if not results:
            # Return empty DataFrame with expected schema
            return self._empty_result_schema(risearch_df)

        return pl.concat(results, how="diagonal")

    def _cross_join_filter(
        self,
        pred_chunk: pl.DataFrame,
        trans_chunk: pl.DataFrame,
    ) -> pl.DataFrame:
        """Cross-join with containment filter for smaller chunks."""
        # Select only needed columns from transcriptome
        trans_subset = trans_chunk.select(
            ["trans_start", "trans_end", "gene_id", "transcript_id", "exp_value"]
        )

        # Cross join
        crossed = pred_chunk.join(trans_subset, how="cross")

        # Filter for containment
        contained = crossed.filter(
            (pl.col("start") >= pl.col("trans_start"))
            & (pl.col("end") <= pl.col("trans_end"))
        )

        return contained

    def _sorted_interval_join(
        self,
        pred_chunk: pl.DataFrame,
        trans_chunk: pl.DataFrame,
    ) -> pl.DataFrame:
        """Sorted interval join for large chunks using batch processing.

        Process predictions in batches to limit memory.
        """
        BATCH_SIZE = 10000
        results = []

        # Sort transcriptome by start for efficient searching
        trans_sorted = trans_chunk.sort("trans_start")
        trans_subset = trans_sorted.select(
            ["trans_start", "trans_end", "gene_id", "transcript_id", "exp_value"]
        )

        # Process predictions in batches
        for i in range(0, pred_chunk.height, BATCH_SIZE):
            batch = pred_chunk.slice(i, BATCH_SIZE)

            # For each batch, do cross-join + filter
            crossed = batch.join(trans_subset, how="cross")
            contained = crossed.filter(
                (pl.col("start") >= pl.col("trans_start"))
                & (pl.col("end") <= pl.col("trans_end"))
            )

            if contained.height > 0:
                results.append(contained)

        if not results:
            return self._empty_result_schema(pred_chunk)

        return pl.concat(results, how="diagonal")

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
