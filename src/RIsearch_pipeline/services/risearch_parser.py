"""Service for parsing RIsearch2 output files."""

from pathlib import Path

import polars as pl

from RIsearch_pipeline.models import RISEARCH_COLUMNS, RISEARCH_SCHEMA


class RIsearchParser:
    """Parser for RIsearch2 prediction output files."""

    def load(self, path: Path) -> pl.DataFrame:
        """Load RIsearch2 TSV output into a validated DataFrame.

        Args:
            path: Path to the RIsearch2 output file (TSV, optionally gzipped).

        Returns:
            Polars DataFrame with typed columns.

        Raises:
            FileNotFoundError: If the file does not exist.
            pl.exceptions.ComputeError: If parsing fails.
        """
        # Ensure path is a Path object (Typer might pass it as something else?)
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        # RIsearch2 outputs are tab-separated with no header
        # We only need specific columns as requested
        df = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            columns=[
                0,
                3,
                4,
                5,
                6,
                7,
            ],  # Select only required columns by index
            new_columns=RISEARCH_COLUMNS,
            schema_overrides=RISEARCH_SCHEMA,
            truncate_ragged_lines=True,
        )

        # Filter out empty rows (trailing newlines)
        df = df.filter(pl.col("sirna_id").is_not_null())

        return df

    def summary(self, df: pl.DataFrame) -> dict:
        """Generate summary statistics for loaded predictions.

        Args:
            df: DataFrame from load().

        Returns:
            Dictionary with row_count, energy_min, energy_max, chromosomes.
        """
        return {
            "row_count": df.height,
            "energy_min": df["energy"].min(),
            "energy_max": df["energy"].max(),
            "chromosomes": df["chrom"].unique().to_list(),
            "strands": df["strand"].unique().to_list(),
        }

    def scan_sirna_ids(self, path: Path) -> list[str]:
        """Extract unique siRNA IDs via streaming scan (low memory).

        Uses Polars lazy evaluation to scan the file without loading
        the entire dataset into memory. Supports gzipped files natively.

        Args:
            path: Path to RIsearch2 output file (TSV, optionally gzipped).

        Returns:
            List of unique siRNA IDs in order of first appearance.
        """
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        # Use scan_csv for lazy evaluation - only materializes unique IDs
        return (
            pl.scan_csv(
                path,
                separator="\t",
                has_header=False,
                schema_overrides={"column_1": pl.Utf8},
            )
            .select(pl.col("column_1").alias("sirna_id"))
            .unique(maintain_order=True)
            .collect(engine="streaming")["sirna_id"]
            .to_list()
        )

    def load_by_sirna(self, path: Path, sirna_id: str) -> pl.DataFrame:
        """Load predictions for a single siRNA via streaming filter.

        Uses Polars lazy evaluation to filter before collecting,
        minimizing memory usage. Supports gzipped files natively.

        Args:
            path: Path to RIsearch2 output file (TSV, optionally gzipped).
            sirna_id: The siRNA ID to filter for.

        Returns:
            Polars DataFrame with predictions for the specified siRNA only.
        """
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        # Lazy scan, filter, then collect - only materializes matching rows
        df = (
            pl.scan_csv(
                path,
                separator="\t",
                has_header=False,
            )
            .filter(pl.col("column_1") == sirna_id)
            .select(
                [
                    pl.col("column_1").alias("sirna_id"),
                    pl.col("column_4").alias("chrom"),
                    pl.col("column_5").cast(pl.Int32).alias("start"),
                    pl.col("column_6").cast(pl.Int32).alias("end"),
                    pl.col("column_7").alias("strand"),
                    pl.col("column_8").cast(pl.Float32).alias("energy"),
                ]
            )
            .collect(engine="streaming")
        )

        return df

    def load_by_sirna_batch(self, path: Path, sirna_ids: list[str]) -> pl.DataFrame:
        """Load predictions for multiple siRNAs via streaming filter.

        Uses Polars lazy evaluation with is_in() filter for efficient batch loading.
        This enables Polars/Rayon parallelization across all rows in the batch.

        Args:
            path: Path to RIsearch2 output file (TSV, optionally gzipped).
            sirna_ids: List of siRNA IDs to filter for.

        Returns:
            Polars DataFrame with predictions for all specified siRNAs.
        """
        if not isinstance(path, Path):
            path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        if not sirna_ids:
            return pl.DataFrame(
                schema={
                    "sirna_id": pl.Utf8,
                    "chrom": pl.Utf8,
                    "start": pl.Int32,
                    "end": pl.Int32,
                    "strand": pl.Utf8,
                    "energy": pl.Float32,
                }
            )

        # Lazy scan, filter by batch, then collect
        df = (
            pl.scan_csv(
                path,
                separator="\t",
                has_header=False,
            )
            .filter(pl.col("column_1").is_in(sirna_ids))
            .select(
                [
                    pl.col("column_1").alias("sirna_id"),
                    pl.col("column_4").alias("chrom"),
                    pl.col("column_5").cast(pl.Int32).alias("start"),
                    pl.col("column_6").cast(pl.Int32).alias("end"),
                    pl.col("column_7").alias("strand"),
                    pl.col("column_8").cast(pl.Float32).alias("energy"),
                ]
            )
            .collect(engine="streaming")
        )

        return df

    def load_directory_batch(
        self, directory: Path, file_paths: list[Path]
    ) -> pl.DataFrame:
        """Load multiple per-siRNA files in a single batch.

        Uses pl.concat with lazy frames to enable Polars/Rayon parallelization
        across all rows from all files in the batch.

        Args:
            directory: Base directory (for error messages).
            file_paths: List of file paths to load.

        Returns:
            Combined Polars DataFrame with predictions from all files.
        """
        if not file_paths:
            return pl.DataFrame(
                schema={
                    "sirna_id": pl.Utf8,
                    "chrom": pl.Utf8,
                    "start": pl.Int32,
                    "end": pl.Int32,
                    "strand": pl.Utf8,
                    "energy": pl.Float32,
                }
            )

        # Create lazy frames for each file - Polars parallelizes this
        # Also compute raw E_min per siRNA from the full file before any
        # intersection filtering.  The old pipeline anchors alpha/gamma clamping
        # to the minimum hybridisation energy across ALL genome hits (including
        # those that don't overlap any annotated transcript), not just the
        # post-intersection subset.  Attaching raw_e_min here lets
        # ProbabilityService replicate that behaviour.
        lazy_frames = []
        raw_emin_frames = []
        for f in file_paths:
            lf = pl.scan_csv(f, separator="\t", has_header=False).select(
                [
                    pl.col("column_1").alias("sirna_id"),
                    pl.col("column_4").alias("chrom"),
                    pl.col("column_5").cast(pl.Int32).alias("start"),
                    pl.col("column_6").cast(pl.Int32).alias("end"),
                    pl.col("column_7").alias("strand"),
                    pl.col("column_8").cast(pl.Float32).alias("energy"),
                ]
            )
            lazy_frames.append(lf)
            # Per-siRNA minimum energy from the raw file
            raw_emin_frames.append(
                pl.scan_csv(f, separator="\t", has_header=False)
                .select(
                    pl.col("column_1").alias("sirna_id"),
                    pl.col("column_8").cast(pl.Float32).alias("energy"),
                )
                .group_by("sirna_id")
                .agg(pl.col("energy").min().alias("raw_e_min"))
            )

        # Concat and collect - Rayon parallelizes across all rows
        df = pl.concat(lazy_frames).collect()
        raw_emin = pl.concat(raw_emin_frames).collect()

        df = df.join(raw_emin, on="sirna_id", how="left")

        return df

    def list_directory_files(self, directory: Path) -> list[Path]:
        """List all RIsearch output files in a directory.

        Args:
            directory: Directory to scan.

        Returns:
            List of file paths (*.gz, *.tsv, *.out).
        """
        if not isinstance(directory, Path):
            directory = Path(directory)

        if not directory.exists():
            raise FileNotFoundError(f"Directory not found: {directory}")

        files = []
        for ext in ["*.gz", "*.tsv", "*.out"]:
            files.extend(directory.glob(ext))

        return sorted(files)
