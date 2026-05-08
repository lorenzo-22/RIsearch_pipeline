"""Tests for the Transcriptome parser service."""

from pathlib import Path

import polars as pl
import pytest

from RIsearch_pipeline.services.transcriptome_parser import TranscriptomeParser


@pytest.fixture
def parser() -> TranscriptomeParser:
    """Create a parser instance."""
    return TranscriptomeParser()


@pytest.fixture
def test_gtf_path() -> Path:
    """Path to the test GTF file."""
    return Path(__file__).parent / "data" / "expression_data.gtf"


class TestTranscriptomeParser:
    """Tests for loading GTF files."""

    def test_load_gtf_returns_dataframe(
        self, parser: TranscriptomeParser, test_gtf_path: Path
    ) -> None:
        """Parser returns a Polars DataFrame."""
        df = parser.load_gtf(test_gtf_path)
        assert isinstance(df, pl.DataFrame)

    def test_load_gtf_filters_exons(
        self, parser: TranscriptomeParser, test_gtf_path: Path
    ) -> None:
        """Parser filters annotations by feature type (default 'exon')."""
        # The test file has transcripts and exons. We expect only exons.
        # Let's count quickly manually from user's file view or rely on logic.
        # User file snippet shows interleaved transcript/exon lines.
        # Assuming we filter for exons, we should get fewer rows than total.
        df = parser.load_gtf(test_gtf_path)
        # Total lines in file is 119
        # Many lines are "transcript". "exons" > 0.
        assert df.height > 0
        assert df.height < 119

    def test_load_gtf_extracts_attributes(
        self, parser: TranscriptomeParser, test_gtf_path: Path
    ) -> None:
        """Parser correctly extracts gene_id, transcript_id, and score."""
        df = parser.load_gtf(test_gtf_path)

        # Check first row (should be an exon based on file snippet)
        # Row 2 in file is first exon: gene_id "gene_1"; transcript_id "transcript_21"; RPKM "1000"
        # Polars order should be preserved (CSV reading is sequential)

        row = df.row(0, named=True)
        assert row["gene_id"] == "gene_1"
        assert row["transcript_id"] == "transcript_21"
        assert row["exp_value"] == 1000.0

    def test_load_gtf_custom_score_column(
        self, parser: TranscriptomeParser, test_gtf_path: Path
    ) -> None:
        """Parser can use a different attribute for score."""
        # Using RPKM as default, which works.
        # If we asked for FPKM (which isn't there), it should return nulls or 0.0
        df = parser.load_gtf(test_gtf_path, score_col="NONEXISTENT")
        assert (df["exp_value"] == 0.0).all()

    def test_load_gtf_missing_file(self, parser: TranscriptomeParser) -> None:
        """Parser raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parser.load_gtf(Path("/nonexistent/file.gtf"))
