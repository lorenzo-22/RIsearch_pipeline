"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl


class IntersectionService:
    """Service to intersect off-target predictions with genomic features."""

    def intersect(
        self, risearch_df: pl.DataFrame, transcriptome_df: pl.DataFrame
    ) -> pl.DataFrame:
        """Filter predictions to those contained within transcriptome features.

        Performs an inner join on chromosome and strand, then filters rows where
        the prediction interval [start, end] is fully contained within the
        feature interval [start, end].

        Args:
            risearch_df: DataFrame from RIsearchParser.
            transcriptome_df: DataFrame from TranscriptomeParser.

        Returns:
            DataFrame containing intersected results. Includes columns from both
            inputs, with suffixes if needed.
        """
        # We need to distinguish start/end columns
        # risearch: start, end
        # transcriptome: start, end

        # Rename transcriptome columns to avoid collision and allow comparison
        # We keep chrom and strand for the join key

        # Suffix handling in join:
        # Polars join allows suffix parameter.

        # Perform join
        # This will create a Cartesian product for each (chrom, strand) pair,
        # which is memory intensive if pairs are large, but usually manageable
        # for genome-wide siRNA studies unless thousands of siRNAs vs full genome at once.
        # Given standard typical usage, this is acceptable. Use lazy() for optimization if needed later.

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

        # Clean up:
        # We might want to keep the transcript info
        # Columns: sirna_id, chrom, start, end, strand, energy,
        #          source, feature, start_trans, end_trans, score_raw, frame, attributes,
        #          gene_id, transcript_id, exp_value

        return intersected
