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
        # Rename transcriptome columns to avoid collision and allow comparison
        # We keep chrom and strand for the join key in gw mode

        if mode == "tw":
            # Transcriptome-Wide Mode
            # RIsearch 'chrom' column is actually Transcript ID.
            # We join on Transcript ID.
            # We assume prediction is 'within' the transcript by definition if IDs match.

            # We join on risearch.chrom == transcriptome.transcript_id
            # We can also populate 'transcript_id' in the result to be the same as chrom

            joined = risearch_df.join(
                transcriptome_df,
                left_on="chrom",
                right_on="transcript_id",
                how="inner",
                suffix="_trans",
            )

            # Ensure transcript_id column exists (it might have been consumed or renamed depending on polars version/join)
            # In Polars join, if columns are different names, both are kept.
            # So we will have 'chrom' (from risearch) and 'transcript_id' (from transcriptome, but invalid as join key?? No, left_on/right_on keeps them usually)
            # Actually, let's explicit alias to be safe.

            # But wait, if we join on differnet names, both cols are kept.
            # We want 'transcript_id' to satisfy ProbabilityService which uses it for grouping.
            # And 'chrom' should ideally stay as the ID too.

            if "transcript_id" not in joined.columns:
                joined = joined.with_columns(pl.col("chrom").alias("transcript_id"))

            return joined

        else:
            # Genome-Wide Mode (Default)
            # Join on (chrom, strand)
            joined = risearch_df.join(
                transcriptome_df,
                on=["chrom", "strand"],
                how="inner",
                suffix="_trans",
            )

            # Filter for containment
            # risearch interval: start, end
            # transcriptome interval: start_trans, end_trans

            # Condition: risearch_start >= trans_start AND risearch_end <= trans_end
            intersected = joined.filter(
                (pl.col("start") >= pl.col("start_trans"))
                & (pl.col("end") <= pl.col("end_trans"))
            )

            return intersected
