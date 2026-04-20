#!/usr/bin/env python3
"""
Convert RIsearch2 .out.gz files to per-siRNA Parquet files.

Each .out.gz is a headerless 12-column TSV produced by RIsearch2:
  col 0: sirna_id  col 3: chrom  col 4: start  col 5: end
  col 6: strand    col 7: energy  (remaining cols are discarded)

Each input file gets a corresponding output .parquet with the same name stem,
containing the 6 relevant columns + pre-computed raw_e_min per siRNA.

Workers write directly to disk — no large DataFrames are transferred between
processes, so memory usage stays low regardless of dataset size.

Usage:
    # Convert a directory in-place (output alongside input)
    python3 convert_risearch_to_parquet.py <input_dir>

    # Convert to a separate output directory
    python3 convert_risearch_to_parquet.py <input_dir> --out-dir <output_dir>

    # Control parallelism
    python3 convert_risearch_to_parquet.py <input_dir> --workers 32

    # Convert only a specific subset of siRNA IDs (one ID per line)
    python3 convert_risearch_to_parquet.py <input_dir> --ids-file subset.ids --out-dir <output_dir>
"""

import argparse
import multiprocessing
import os
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


# ── Per-file worker (runs in subprocess — writes to disk, returns nothing big) ─

def _convert_one(src: str, dst: str) -> tuple[str, int]:
    """Convert one .out.gz to .parquet. Returns (dst_path, row_count)."""
    # Must be set before polars is imported so the thread pool is sized correctly
    os.environ["POLARS_MAX_THREADS"] = "1"
    import polars as pl

    df = pl.read_csv(
        src,
        has_header=False,
        separator="\t",
        columns=[0, 3, 4, 5, 6, 7],
        new_columns=["sirna_id", "chrom", "start", "end", "strand", "energy"],
        schema_overrides={
            "sirna_id": pl.Utf8,
            "chrom":    pl.Utf8,
            "start":    pl.Int32,
            "end":      pl.Int32,
            "strand":   pl.Utf8,
            "energy":   pl.Float32,
        },
        truncate_ragged_lines=True,
    ).filter(pl.col("sirna_id").is_not_null())

    if df.height == 0:
        return dst, 0

    emin = df.group_by("sirna_id").agg(pl.col("energy").min().alias("raw_e_min"))
    df = df.join(emin, on="sirna_id", how="left")

    df.write_parquet(dst, compression="zstd", compression_level=3, statistics=True)
    return dst, df.height


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_dir", type=Path,
                        help="Directory containing risearch_*.out.gz files")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory (default: same as input_dir)")
    parser.add_argument("--workers", type=int, default=16,
                        help="Parallel workers (default: 16)")
    parser.add_argument("--ids-file", type=Path, default=None,
                        help="File with siRNA IDs to convert (one per line); "
                             "if omitted all risearch_*.out.gz files are converted")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip files whose .parquet already exists")
    args = parser.parse_args()

    input_dir: Path = args.input_dir
    out_dir: Path   = args.out_dir or input_dir
    n_workers: int  = args.workers

    if not input_dir.is_dir():
        sys.exit(f"Error: {input_dir} is not a directory")

    out_dir.mkdir(parents=True, exist_ok=True)

    if args.ids_file is not None:
        if not args.ids_file.exists():
            sys.exit(f"Error: --ids-file {args.ids_file} not found")
        ids = [line.strip() for line in args.ids_file.read_text().splitlines() if line.strip()]
        files = sorted(
            f for sid in ids
            for f in [input_dir / f"risearch_{sid}.out.gz"]
            if f.exists()
        )
        missing = len(ids) - len(files)
        if missing:
            print(f"Warning: {missing}/{len(ids)} IDs had no matching .out.gz file")
    else:
        files = sorted(input_dir.glob("risearch_*.out.gz"))

    if not files:
        sys.exit(f"Error: no risearch_*.out.gz files found in {input_dir}")

    # Build (src, dst) pairs — optionally skip already-converted files
    pairs = []
    for f in files:
        dst = out_dir / (f.stem.replace(".out", "") + ".parquet")
        if args.skip_existing and dst.exists():
            continue
        pairs.append((str(f), str(dst)))

    skipped = len(files) - len(pairs)
    print(f"Converting {len(pairs)} files → {out_dir}"
          + (f"  ({skipped} skipped)" if skipped else ""))
    print(f"  workers={n_workers}  compression=zstd:3")

    t0 = time.perf_counter()
    done = errors = total_rows = 0

    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
        futs = {pool.submit(_convert_one, src, dst): (src, dst) for src, dst in pairs}
        for fut in as_completed(futs):
            done += 1
            src, dst = futs[fut]
            try:
                _, rows = fut.result()
                total_rows += rows
            except Exception as e:
                errors += 1
                if errors <= 5:
                    print(f"\n  Error: {Path(src).name}: {e}")

            if done % 100 == 0 or done == len(pairs):
                elapsed = time.perf_counter() - t0
                rate = done / elapsed
                print(f"  {done}/{len(pairs)}  ({rate:.0f} files/s)  {errors} errors",
                      end="\r", flush=True)

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.1f}s  "
          f"({len(pairs)} files, {total_rows:,} rows, {errors} errors)")


if __name__ == "__main__":
    main()
