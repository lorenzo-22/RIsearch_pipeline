"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl
import numpy as np
from typing import Iterator, Optional

from loguru import logger

try:
    from ncls import NCLS
except Exception:  # pragma: no cover - optional dependency
    NCLS = None


class IntersectionService:
    """Service to intersect off-target predictions with genomic features."""

    def intersect(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
        mode: str = "gw",
        workers: int = 1,
    ) -> pl.DataFrame:
        """Filter predictions to those overlapping transcriptome features.

        Args:
            risearch_df: DataFrame from RIsearchParser.
            transcriptome_df: DataFrame from TranscriptomeParser.
            mode: "gw" (genome-wide) or "tw" (transcriptome-wide).
            workers: Number of threads for parallel chromosome processing.

        Returns:
            DataFrame containing intersected results.
        """
        if mode == "tw":
            return self._transcriptome_wide_join(risearch_df, transcriptome_df)
        else:
            return self._genome_wide_streaming(risearch_df, transcriptome_df, workers=workers)

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
        workers: int = 1,
    ) -> pl.DataFrame:
        """Genome-wide intersection with streaming to limit memory.

        Processes predictions in small batches per (chrom, strand) pair.
        When workers > 1, pairs are processed in parallel via ThreadPoolExecutor
        (NCLS and Polars/Rust operations release the GIL).
        """
        BATCH_SIZE = 50000  # Process 50K predictions at a time

        # Get unique (chrom, strand) pairs
        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        def _process_pair(pair: dict) -> list[pl.DataFrame]:
            chrom = pair["chrom"]
            strand = pair["strand"]

            preds = risearch_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )
            trans = transcriptome_df.filter(
                (pl.col("chrom") == chrom) & (pl.col("strand") == strand)
            )

            if preds.height == 0 or trans.height == 0:
                return []

            trans_sorted = trans.sort("start")
            trans_index = self._build_transcript_index(trans_sorted)
            logger.debug(f"Intersecting chrom={chrom} strand={strand}: {preds.height} preds × {trans.height} transcripts")

            results = []
            for batch_start in range(0, preds.height, BATCH_SIZE):
                batch = preds.slice(batch_start, BATCH_SIZE)
                batch_result = self._process_batch(batch, trans_index)

                if batch_result is not None and batch_result.height > 0:
                    # Deduplicate each batch before accumulating to keep peak
                    # memory proportional to batch size rather than total rows.
                    batch_result = batch_result.unique(
                        subset=["sirna_id", "chrom", "start", "end", "strand", "energy", "transcript_id"],
                        keep="first",
                    )
                    results.append(batch_result)
            return results

        all_results: list[pl.DataFrame] = []

        if workers > 1 and len(chrom_strand_pairs) > 1:
            from concurrent.futures import ThreadPoolExecutor, as_completed

            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_process_pair, pair): pair for pair in chrom_strand_pairs}
                for future in as_completed(futures):
                    all_results.extend(future.result())
        else:
            for pair in chrom_strand_pairs:
                all_results.extend(_process_pair(pair))

        if not all_results:
            return self._empty_result_schema(risearch_df)

        result = pl.concat(all_results, how="diagonal")
        logger.debug(f"GW intersection: {result.height:,} rows after per-batch dedup")
        return result

    def _build_transcript_index(self, trans_sorted: pl.DataFrame) -> dict:
        """Pre-build numpy arrays and optional interval index for fast lookup."""
        starts = trans_sorted["start"].to_numpy().astype(np.int64, copy=False)
        ends = trans_sorted["end"].to_numpy().astype(np.int64, copy=False)

        index = {
            "starts": starts,
            "ends": ends,
            "gene_ids": trans_sorted["gene_id"].to_list(),
            "transcript_ids": trans_sorted["transcript_id"].to_list(),
            "exp_values": trans_sorted["exp_value"].to_numpy(),
            "trans_df": trans_sorted,
        }

        if NCLS is not None and len(starts) > 0:
            # BED ends are already exclusive (half-open); pass them directly to NCLS.
            interval_ids = np.arange(len(starts), dtype=np.int64)
            index["ncls"] = NCLS(starts, ends, interval_ids)

        return index

    def _process_batch(
        self,
        preds: pl.DataFrame,
        trans_index: dict,
    ) -> Optional[pl.DataFrame]:
        """Process a batch of predictions using the fastest available method."""
        if trans_index.get("ncls") is not None:
            return self._process_batch_ncls(preds, trans_index)

        return self._process_batch_numpy(preds, trans_index)

    def _process_batch_ncls(
        self,
        preds: pl.DataFrame,
        trans_index: dict,
    ) -> Optional[pl.DataFrame]:
        """Process a batch using NCLS interval indexing (fast path)."""
        if preds.height == 0:
            return None

        trans_df = trans_index["trans_df"]
        trans_starts = trans_index["starts"]
        trans_ends = trans_index["ends"]
        ncls = trans_index["ncls"]

        pred_starts = preds["start"].to_numpy().astype(np.int64, copy=False)
        pred_ends = preds["end"].to_numpy().astype(np.int64, copy=False)
        pred_ids = np.arange(preds.height, dtype=np.int64)

        # Query with BED half-open intervals — pred_ends is already the exclusive end.
        # all_overlaps_both returns all pairs where intervals overlap (any overlap),
        # matching bedtools intersect default behaviour.
        query_starts = pred_starts
        query_ends = pred_ends
        pred_idx, trans_idx = ncls.all_overlaps_both(query_starts, query_ends, pred_ids)

        if len(pred_idx) == 0:
            return None

        preds_sel = preds[pred_idx]
        trans_sel = trans_df[trans_idx].select(
            [
                pl.col("start").alias("trans_start"),
                pl.col("end").alias("trans_end"),
                pl.col("gene_id"),
                pl.col("transcript_id"),
                pl.col("exp_value"),
            ]
        )

        return preds_sel.hstack(trans_sel)

    def _process_batch_numpy(
        self,
        preds: pl.DataFrame,
        trans_index: dict,
    ) -> Optional[pl.DataFrame]:
        """Process a batch using vectorized numpy searches (fallback path)."""
        if preds.height == 0:
            return None

        trans_df = trans_index["trans_df"]
        trans_starts = trans_index["starts"]
        trans_ends = trans_index["ends"]

        pred_starts = preds["start"].to_numpy().astype(np.int64, copy=False)
        pred_ends = preds["end"].to_numpy().astype(np.int64, copy=False)

        order = np.argsort(pred_starts, kind="mergesort")
        pred_starts_sorted = pred_starts[order]
        pred_ends_sorted = pred_ends[order]

        pred_match_idx_list = []
        trans_match_idx_list = []

        for t_idx in range(len(trans_starts)):
            t_start = trans_starts[t_idx]
            t_end = trans_ends[t_idx]

            # Any-overlap: pred_start < t_end AND pred_end > t_start.
            # All predictions starting before t_end may overlap; then filter
            # to those whose end is past t_start.
            right = np.searchsorted(pred_starts_sorted, t_end, side="left")

            if right == 0:
                continue

            ends_slice = pred_ends_sorted[:right]
            mask = ends_slice > t_start
            if not np.any(mask):
                continue

            pred_match_idx = order[:right][mask]
            pred_match_idx_list.append(pred_match_idx)
            trans_match_idx_list.append(
                np.full(pred_match_idx.shape[0], t_idx, dtype=np.int64)
            )

        if not pred_match_idx_list:
            return None

        pred_idx = np.concatenate(pred_match_idx_list)
        trans_idx = np.concatenate(trans_match_idx_list)

        preds_sel = preds[pred_idx]
        trans_sel = trans_df[trans_idx].select(
            [
                pl.col("start").alias("trans_start"),
                pl.col("end").alias("trans_end"),
                pl.col("gene_id"),
                pl.col("transcript_id"),
                pl.col("exp_value"),
            ]
        )

        return preds_sel.hstack(trans_sel)

    def intersect_streaming(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
        mode: str = "gw",
        workers: int = 1,
    ) -> Iterator[pl.DataFrame]:
        """Streaming intersection that yields batches instead of accumulating.

        Use this for large datasets to avoid OOM.
        When workers > 1, delegates to _genome_wide_streaming for parallel
        chromosome processing, then yields results in chunks.
        """
        if mode == "tw":
            yield self._transcriptome_wide_join(risearch_df, transcriptome_df)
            return

        if workers > 1:
            # Parallel path: collect all results, then yield in BATCH_SIZE chunks
            result = self._genome_wide_streaming(risearch_df, transcriptome_df, workers=workers)
            BATCH_SIZE = 50000
            for offset in range(0, result.height, BATCH_SIZE):
                yield result.slice(offset, BATCH_SIZE)
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
