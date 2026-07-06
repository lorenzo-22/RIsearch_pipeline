"""Tests for the Intersection Service."""

import polars as pl
import pytest

from riot.services.intersection_service import IntersectionService


@pytest.fixture
def service() -> IntersectionService:
    """Create service instance."""
    return IntersectionService()


@pytest.fixture
def risearch_df() -> pl.DataFrame:
    """Dummy RIsearch dataframe."""
    return pl.DataFrame(
        {
            "sirna_id": ["si1", "si2", "si3", "si4"],
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [100, 200, 300, 100],
            "end": [120, 220, 320, 120],
            "strand": ["+", "+", "-", "+"],
            "energy": [-15.0, -10.0, -12.0, -10.0],
        }
    )


@pytest.fixture
def transcriptome_df() -> pl.DataFrame:
    """Dummy transcriptome dataframe."""
    return pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [90, 250, 50],
            "end": [150, 350, 150],
            "strand": ["+", "+", "+"],
            "gene_id": ["g1", "g2", "g3"],
            "transcript_id": ["t1", "t2", "t3"],
            "exp_value": [10.0, 20.0, 30.0],
        }
    )


class TestIntersectionService:
    """Tests for intersect logic."""

    def test_intersect_containment(
        self,
        service: IntersectionService,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> None:
        """Test simple containment intersection."""
        # si1 (100-120) is inside t1 (90-150) on chr1 (+) -> MATCH
        # si4 (100-120) is inside t3 (50-150) on chr2 (+) -> MATCH

        result = service.intersect(risearch_df, transcriptome_df)

        # Check si1 match
        si1_match = result.filter(pl.col("sirna_id") == "si1")
        assert si1_match.height == 1
        assert si1_match["gene_id"][0] == "g1"

        # Check si4 match
        si4_match = result.filter(pl.col("sirna_id") == "si4")
        assert si4_match.height == 1
        assert si4_match["gene_id"][0] == "g3"

    def test_intersect_excludes_non_overlap(
        self,
        service: IntersectionService,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> None:
        """Test filtering of non-overlapping regions."""
        # si2 (200-220) is NOT in t1 (90-150) nor t2 (250-350) -> NO MATCH

        result = service.intersect(risearch_df, transcriptome_df)
        assert result.filter(pl.col("sirna_id") == "si2").height == 0

    def test_intersect_strands(
        self,
        service: IntersectionService,
        risearch_df: pl.DataFrame,
        transcriptome_df: pl.DataFrame,
    ) -> None:
        """Test strict strand matching."""
        # si3 (300-320) is on '-', t2 (250-350) is on '+' -> NO MATCH despite coordinate overlap

        result = service.intersect(risearch_df, transcriptome_df)
        assert result.filter(pl.col("sirna_id") == "si3").height == 0
