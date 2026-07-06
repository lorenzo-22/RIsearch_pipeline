"""Service for parsing RIsearch2 output files."""

from pathlib import Path

import polars as pl

from riot.models import RISEARCH_COLUMNS, RISEARCH_SCHEMA

_VALID_SUFFIXES = {".gz", ".tsv", ".out", ".parquet"}

_TSV_SELECT = [
    pl.col("column_1").alias("sirna_id"),
    pl.col("column_4").alias("chrom"),
    pl.col("column_5").cast(pl.Int32).alias("start"),
    pl.col("column_6").cast(pl.Int32).alias("end"),
    pl.col("column_7").alias("strand"),
    pl.col("column_8").cast(pl.Float32).alias("energy"),
]


def _check_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"RIsearch2 output file not found: {path}")


class RIsearchParser:
    """Parser for RIsearch2 prediction output files."""

    def load(self, path: Path) -> pl.DataFrame:
        """Load RIsearch2 TSV output into a validated DataFrame."""
        _check_file(path)
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

    def load_single_file(self, file_path: Path) -> pl.DataFrame:
        """Load one per-siRNA RIsearch file and attach raw_e_min.

        Parquet files from convert_risearch_to_parquet.py already contain
        named columns and pre-computed raw_e_min — returned directly.
        """
        if file_path.suffix == ".parquet":
            return pl.read_parquet(file_path)

        df = (
            pl.scan_csv(file_path, separator="\t", has_header=False)
            .select(_TSV_SELECT)
            .collect()
        )
        raw_emin = df.group_by("sirna_id").agg(pl.col("energy").min().alias("raw_e_min"))
        return df.join(raw_emin, on="sirna_id", how="left")

    def list_directory_files(self, directory: Path) -> list[Path]:
        """List all RIsearch output files in a directory (*.gz, *.tsv, *.out, *.parquet)."""
        if not directory.exists():
            raise FileNotFoundError(f"Directory not found: {directory}")
        return sorted(f for f in directory.iterdir() if f.suffix in _VALID_SUFFIXES)
