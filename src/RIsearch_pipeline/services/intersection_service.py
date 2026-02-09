"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl
import gc


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
            return self._range_join_streaming(risearch_df, transcriptome_df)

    def _range_join_streaming(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
        """Streaming range-join that processes and concatenates per-chromosome.

        Uses lazy concatenation to avoid holding all results in memory.
        """
        # Prepare transcriptome with renamed columns
        trans_renamed = transcriptome_df.rename(
            {
                "start": "trans_start",
                "end": "trans_end",
            }
        )

        # Get unique (chrom, strand) pairs
        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        # Process each chromosome in a streaming fashion
        # Use lazy frames for memory efficiency
        lazy_results = []

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

            # Process this chunk
            chunk_result = self._process_chunk(pred_chunk, trans_chunk)

            if chunk_result.height > 0:
                # Convert to lazy for efficient concatenation
                lazy_results.append(chunk_result.lazy())

            # Free memory
            del pred_chunk, trans_chunk, chunk_result
            gc.collect()

        if not lazy_results:
            return self._empty_result_schema(risearch_df)

        # Collect all lazy frames at once - more memory efficient
        return pl.concat(lazy_results).collect()

    def _process_chunk(
        self,
        pred_chunk: pl.DataFrame,
        trans_chunk: pl.DataFrame,
    ) -> pl.DataFrame:
        """Process a single chromosome+strand chunk using batched cross-join."""
        BATCH_SIZE = 5000  # Process predictions in smaller batches
        results = []

        trans_subset = trans_chunk.select(
            ["trans_start", "trans_end", "gene_id", "transcript_id", "exp_value"]
        )

        # If small enough, process directly
        if pred_chunk.height * trans_subset.height < 50_000_000:
            crossed = pred_chunk.join(trans_subset, how="cross")
            return crossed.filter(
                (pl.col("start") >= pl.col("trans_start"))
                & (pl.col("end") <= pl.col("trans_end"))
            )

        # Batch process for large chunks
        for i in range(0, pred_chunk.height, BATCH_SIZE):
            batch = pred_chunk.slice(i, BATCH_SIZE)
            crossed = batch.join(trans_subset, how="cross")
            contained = crossed.filter(
                (pl.col("start") >= pl.col("trans_start"))
                & (pl.col("end") <= pl.col("trans_end"))
            )
            if contained.height > 0:
                results.append(contained)

            # Free intermediate memory
            del batch, crossed, contained

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
