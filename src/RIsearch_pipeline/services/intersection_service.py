"""Service for intersecting RIsearch predictions with transcriptome annotations."""

import polars as pl
import numpy as np
from typing import Optional

from loguru import logger

try:
    from ncls import NCLS
except Exception:  # pragma: no cover - optional dependency
    NCLS = None


class IntersectionService:
    """Service to intersect off-target predictions with genomic features."""

    def __init__(self) -> None:
        self._trans_partition: dict[tuple, pl.DataFrame] | None = None
        self._trans_ncls_index: dict[tuple, dict] | None = None

    def preload_transcriptome(self, transcriptome_df: pl.DataFrame) -> None:
        """Pre-partition transcriptome and build NCLS indices for all (chrom,strand) pairs.

        Call once from _init_worker. Subsequent intersect() calls reuse the
        cached data instead of rebuilding per siRNA (~89ms/siRNA saved).
        """
        self._trans_partition = {
            (k, v): sub_df
            for (k, v), sub_df in transcriptome_df.partition_by(
                ["chrom", "strand"], as_dict=True
            ).items()
        }
        self._trans_ncls_index = {
            key: self._build_transcript_index(sub_df.sort("start"))
            for key, sub_df in self._trans_partition.items()
            if sub_df.height > 0
        }

    def intersect(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
        mode: str = "gw",
        workers: int = 1,
        _timings: Optional[dict] = None,
    ) -> pl.DataFrame:
        """Filter predictions to those overlapping transcriptome features.

        Args:
            mode: "gw" (genome-wide) or "tw" (transcriptome-wide).
            workers: Number of threads for parallel chromosome processing.
            _timings: Optional dict filled with per-sub-stage wall-clock seconds.
        """
        if mode == "tw":
            return self._transcriptome_wide_join(risearch_df, transcriptome_df)
        return self._genome_wide_streaming(risearch_df, transcriptome_df, workers=workers, _timings=_timings)

    def _transcriptome_wide_join(
        self,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> pl.DataFrame:
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
        _timings: Optional[dict] = None,
    ) -> pl.DataFrame:
        import time
        BATCH_SIZE = 50000

        chrom_strand_pairs = risearch_df.select(["chrom", "strand"]).unique().to_dicts()

        _t_part_preds = time.perf_counter()
        preds_by_pair: dict[tuple, pl.DataFrame] = {
            (k, v): sub_df
            for (k, v), sub_df in risearch_df.partition_by(
                ["chrom", "strand"], as_dict=True
            ).items()
        }
        _t_part_trans = time.perf_counter()

        trans_by_pair = self._trans_partition if self._trans_partition is not None else {
            (k, v): sub_df
            for (k, v), sub_df in transcriptome_df.partition_by(
                ["chrom", "strand"], as_dict=True
            ).items()
        }
        _t_after_partition = time.perf_counter()

        _ncls_time = [0.0]
        _match_time = [0.0]

        def _process_pair(pair: dict) -> list[pl.DataFrame]:
            chrom, strand = pair["chrom"], pair["strand"]
            preds = preds_by_pair.get((chrom, strand))
            trans = trans_by_pair.get((chrom, strand))

            if preds is None or preds.height == 0 or trans is None or trans.height == 0:
                return []

            _tb = time.perf_counter()
            if self._trans_ncls_index is not None:
                trans_index = self._trans_ncls_index.get((chrom, strand))
                if trans_index is None:
                    return []
            else:
                trans_index = self._build_transcript_index(trans.sort("start"))
            _ncls_time[0] += time.perf_counter() - _tb

            logger.debug(f"Intersecting chrom={chrom} strand={strand}: {preds.height} preds × {trans.height} transcripts")

            results = []
            _tm = time.perf_counter()
            for batch_start in range(0, preds.height, BATCH_SIZE):
                batch_result = self._process_batch(preds.slice(batch_start, BATCH_SIZE), trans_index)
                if batch_result is not None and batch_result.height > 0:
                    results.append(batch_result)
            _match_time[0] += time.perf_counter() - _tm

            if not results:
                return []
            combined = pl.concat(results, how="diagonal")
            return [combined.unique(
                subset=["sirna_id", "chrom", "start", "end", "strand", "energy", "transcript_id"],
                keep="first",
            )]

        all_results: list[pl.DataFrame] = []

        if workers > 1 and len(chrom_strand_pairs) > 1:
            from concurrent.futures import ThreadPoolExecutor, as_completed
            safe_workers = min(workers, 8)
            with ThreadPoolExecutor(max_workers=safe_workers) as pool:
                futures = {pool.submit(_process_pair, pair): pair for pair in chrom_strand_pairs}
                for future in as_completed(futures):
                    all_results.extend(future.result())
        else:
            for pair in chrom_strand_pairs:
                all_results.extend(_process_pair(pair))

        if _timings is not None:
            _timings["intersect_partition_preds_s"] = _t_part_trans - _t_part_preds
            _timings["intersect_partition_trans_s"] = _t_after_partition - _t_part_trans
            _timings["intersect_ncls_build_s"] = _ncls_time[0]
            _timings["intersect_match_s"] = _match_time[0]

        if not all_results:
            return self._empty_result_schema(risearch_df)

        result = pl.concat(all_results, how="diagonal")
        logger.debug(f"GW intersection: {result.height:,} rows after per-batch dedup")
        return result

    def _build_transcript_index(self, trans_sorted: pl.DataFrame) -> dict:
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
            index["ncls"] = NCLS(starts, ends, np.arange(len(starts), dtype=np.int64))

        return index

    def _process_batch(self, preds: pl.DataFrame, trans_index: dict) -> Optional[pl.DataFrame]:
        if trans_index.get("ncls") is not None:
            return self._process_batch_ncls(preds, trans_index)
        return self._process_batch_numpy(preds, trans_index)

    def _process_batch_ncls(self, preds: pl.DataFrame, trans_index: dict) -> Optional[pl.DataFrame]:
        if preds.height == 0:
            return None

        pred_starts = preds["start"].to_numpy().astype(np.int64, copy=False)
        pred_ends = preds["end"].to_numpy().astype(np.int64, copy=False)
        pred_idx, trans_idx = trans_index["ncls"].all_overlaps_both(
            pred_starts, pred_ends, np.arange(preds.height, dtype=np.int64)
        )

        if len(pred_idx) == 0:
            return None

        trans_sel = trans_index["trans_df"][trans_idx].select([
            pl.col("start").alias("trans_start"),
            pl.col("end").alias("trans_end"),
            pl.col("gene_id"),
            pl.col("transcript_id"),
            pl.col("exp_value"),
        ])
        return preds[pred_idx].hstack(trans_sel)

    def _process_batch_numpy(self, preds: pl.DataFrame, trans_index: dict) -> Optional[pl.DataFrame]:
        if preds.height == 0:
            return None

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
            right = np.searchsorted(pred_starts_sorted, trans_ends[t_idx], side="left")
            if right == 0:
                continue
            mask = pred_ends_sorted[:right] > trans_starts[t_idx]
            if not np.any(mask):
                continue
            pred_match_idx = order[:right][mask]
            pred_match_idx_list.append(pred_match_idx)
            trans_match_idx_list.append(np.full(pred_match_idx.shape[0], t_idx, dtype=np.int64))

        if not pred_match_idx_list:
            return None

        pred_idx = np.concatenate(pred_match_idx_list)
        trans_idx = np.concatenate(trans_match_idx_list)

        trans_sel = trans_index["trans_df"][trans_idx].select([
            pl.col("start").alias("trans_start"),
            pl.col("end").alias("trans_end"),
            pl.col("gene_id"),
            pl.col("transcript_id"),
            pl.col("exp_value"),
        ])
        return preds[pred_idx].hstack(trans_sel)

    def _empty_result_schema(self, risearch_df: pl.DataFrame) -> pl.DataFrame:
        return risearch_df.head(0).with_columns([
            pl.lit(None).cast(pl.Int64).alias("trans_start"),
            pl.lit(None).cast(pl.Int64).alias("trans_end"),
            pl.lit(None).cast(pl.Utf8).alias("gene_id"),
            pl.lit(None).cast(pl.Utf8).alias("transcript_id"),
            pl.lit(None).cast(pl.Float32).alias("exp_value"),
        ])
