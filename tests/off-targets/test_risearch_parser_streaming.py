"""Tests for RIsearchParser streaming methods (chunk mode support)."""

from pathlib import Path

import polars as pl
import pytest

from RIsearch_pipeline.services.risearch_parser import RIsearchParser


@pytest.fixture
def parser() -> RIsearchParser:
    """Create a parser instance."""
    return RIsearchParser()


@pytest.fixture
def test_tsv_path() -> Path:
    """Path to the test RIsearch2 output file."""
    return Path(__file__).parent / "data" / "risearch_siRNAID.out"


class TestScanSirnaIds:
    """Tests for RIsearchParser.scan_sirna_ids()."""

    def test_returns_list(self, parser: RIsearchParser, test_tsv_path: Path) -> None:
        """scan_sirna_ids returns a list."""
        result = parser.scan_sirna_ids(test_tsv_path)
        assert isinstance(result, list)

    def test_returns_unique_ids(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """scan_sirna_ids returns unique siRNA IDs."""
        result = parser.scan_sirna_ids(test_tsv_path)
        assert len(result) == len(set(result))

    def test_contains_expected_id(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """scan_sirna_ids includes expected siRNA ID from test file."""
        result = parser.scan_sirna_ids(test_tsv_path)
        # Test file uses 'siRNAID' as the siRNA identifier
        assert "siRNAID" in result

    def test_file_not_found_raises(self, parser: RIsearchParser) -> None:
        """scan_sirna_ids raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parser.scan_sirna_ids(Path("/nonexistent/file.tsv"))


class TestLoadBySirna:
    """Tests for RIsearchParser.load_by_sirna()."""

    def test_returns_dataframe(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna returns a Polars DataFrame."""
        df = parser.load_by_sirna(test_tsv_path, "siRNAID")
        assert isinstance(df, pl.DataFrame)

    def test_filters_correctly(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna returns only rows for specified siRNA."""
        df = parser.load_by_sirna(test_tsv_path, "siRNAID")
        assert (df["sirna_id"] == "siRNAID").all()

    def test_correct_columns(self, parser: RIsearchParser, test_tsv_path: Path) -> None:
        """load_by_sirna creates expected column names."""
        df = parser.load_by_sirna(test_tsv_path, "siRNAID")
        expected_cols = ["sirna_id", "chrom", "start", "end", "strand", "energy"]
        assert df.columns == expected_cols

    def test_energy_is_float(self, parser: RIsearchParser, test_tsv_path: Path) -> None:
        """Energy column has Float32 dtype."""
        df = parser.load_by_sirna(test_tsv_path, "siRNAID")
        assert df["energy"].dtype == pl.Float32

    def test_nonexistent_sirna_returns_empty(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna returns empty DataFrame for unknown siRNA."""
        df = parser.load_by_sirna(test_tsv_path, "NONEXISTENT_SIRNA")
        assert df.height == 0

    def test_file_not_found_raises(self, parser: RIsearchParser) -> None:
        """load_by_sirna raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parser.load_by_sirna(Path("/nonexistent/file.tsv"), "siRNAID")


class TestLoadBySirnaBatch:
    """Tests for RIsearchParser.load_by_sirna_batch()."""

    def test_returns_dataframe(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna_batch returns a Polars DataFrame."""
        df = parser.load_by_sirna_batch(test_tsv_path, ["siRNAID"])
        assert isinstance(df, pl.DataFrame)

    def test_filters_correctly(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna_batch returns only rows for specified siRNAs."""
        df = parser.load_by_sirna_batch(test_tsv_path, ["siRNAID"])
        assert (df["sirna_id"] == "siRNAID").all()

    def test_empty_list_returns_empty(
        self, parser: RIsearchParser, test_tsv_path: Path
    ) -> None:
        """load_by_sirna_batch returns empty DataFrame for empty list."""
        df = parser.load_by_sirna_batch(test_tsv_path, [])
        assert df.height == 0

    def test_correct_columns(self, parser: RIsearchParser, test_tsv_path: Path) -> None:
        """load_by_sirna_batch creates expected column names."""
        df = parser.load_by_sirna_batch(test_tsv_path, ["siRNAID"])
        expected_cols = ["sirna_id", "chrom", "start", "end", "strand", "energy"]
        assert df.columns == expected_cols
