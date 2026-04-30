"""Service for parsing RIsearch2 output files."""

from pathlib import Path

import polars as pl

from RIsearch_pipeline.models import RISEARCH_COLUMNS, RISEARCH_SCHEMA


class RIsearchParser:
    """Parser for RIsearch2 prediction output files."""

    def load(self, path: Path) -> pl.DataFrame:
        """Load RIsearch2 TSV output into a validated DataFrame."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        df = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            columns=[0, 3, 4, 5, 6, 7],
            new_columns=RISEARCH_COLUMNS,
            schema_overrides=RISEARCH_SCHEMA,
            truncate_ragged_lines=True,
        )
        return df.filter(pl.col("sirna_id").is_not_null())

    def summary(self, df: pl.DataFrame) -> dict:
        """Generate summary statistics for loaded predictions."""
        return {
            "row_count": df.height,
            "energy_min": df["energy"].min(),
            "energy_max": df["energy"].max(),
            "chromosomes": df["chrom"].unique().to_list(),
            "strands": df["strand"].unique().to_list(),
        }

    def scan_sirna_ids(self, path: Path) -> list[str]:
        """Extract unique siRNA IDs via streaming scan (low memory)."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        return (
            pl.scan_csv(path, separator="\t", has_header=False,
                        schema_overrides={"column_1": pl.Utf8})
            .select(pl.col("column_1").alias("sirna_id"))
            .unique(maintain_order=True)
            .collect(engine="streaming")["sirna_id"]
            .to_list()
        )

    def load_by_sirna(self, path: Path, sirna_id: str) -> pl.DataFrame:
        """Load predictions for a single siRNA via streaming filter."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        return (
            pl.scan_csv(path, separator="\t", has_header=False)
            .filter(pl.col("column_1") == sirna_id)
            .select([
                pl.col("column_1").alias("sirna_id"),
                pl.col("column_4").alias("chrom"),
                pl.col("column_5").cast(pl.Int32).alias("start"),
                pl.col("column_6").cast(pl.Int32).alias("end"),
                pl.col("column_7").alias("strand"),
                pl.col("column_8").cast(pl.Float32).alias("energy"),
            ])
            .collect(engine="streaming")
        )

    def load_by_sirna_batch(self, path: Path, sirna_ids: list[str]) -> pl.DataFrame:
        """Load predictions for multiple siRNAs via streaming filter."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"RIsearch2 output file not found: {path}")

        if not sirna_ids:
            return pl.DataFrame(schema={
                "sirna_id": pl.Utf8, "chrom": pl.Utf8,
                "start": pl.Int32, "end": pl.Int32,
                "strand": pl.Utf8, "energy": pl.Float32,
            })

        return (
            pl.scan_csv(path, separator="\t", has_header=False)
            .filter(pl.col("column_1").is_in(sirna_ids))
            .select([
                pl.col("column_1").alias("sirna_id"),
                pl.col("column_4").alias("chrom"),
                pl.col("column_5").cast(pl.Int32).alias("start"),
                pl.col("column_6").cast(pl.Int32).alias("end"),
                pl.col("column_7").alias("strand"),
                pl.col("column_8").cast(pl.Float32).alias("energy"),
            ])
            .collect(engine="streaming")
        )

    def load_single_file(self, file_path: Path) -> pl.DataFrame:
        """Load one per-siRNA RIsearch file.

        Parquet files (from convert_risearch_to_parquet.py) already contain
        named columns and pre-computed raw_e_min; returned without re-parsing.
        For TSV/gz files, raw_e_min is derived in a single scan so the old
        pipeline's alpha/gamma clamping anchor (min energy across all genome
        hits) is preserved through intersection filtering.
        """
        file_path = Path(file_path)

        if file_path.suffix == ".parquet":
            return pl.read_parquet(file_path)

        df = (
            pl.scan_csv(file_path, separator="\t", has_header=False)
            .select([
                pl.col("column_1").alias("sirna_id"),
                pl.col("column_4").alias("chrom"),
                pl.col("column_5").cast(pl.Int32).alias("start"),
                pl.col("column_6").cast(pl.Int32).alias("end"),
                pl.col("column_7").alias("strand"),
                pl.col("column_8").cast(pl.Float32).alias("energy"),
            ])
            .collect()
        )
        raw_emin = df.group_by("sirna_id").agg(pl.col("energy").min().alias("raw_e_min"))
        return df.join(raw_emin, on="sirna_id", how="left")

    def list_directory_files(self, directory: Path) -> list[Path]:
        """List all RIsearch output files in a directory."""
        directory = Path(directory)
        if not directory.exists():
            raise FileNotFoundError(f"Directory not found: {directory}")

        files = []
        for ext in ["*.gz", "*.tsv", "*.out", "*.parquet"]:
            files.extend(directory.glob(ext))
        return sorted(files)
