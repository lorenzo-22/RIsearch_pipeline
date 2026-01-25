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
