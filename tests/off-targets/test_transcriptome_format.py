"""Tests for AnnotationParser format flag logic."""

import tempfile
from pathlib import Path

import polars as pl
import pytest

from RIsearch_pipeline.services.annotation_parser import AnnotationParser


@pytest.fixture
def parser():
    return AnnotationParser()


def test_load_gtf_auto_detection(parser):
    """Test auto-detection logic (default)."""
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as f:
        f.write("chr1\t100\t200\tprot1\t0.5\t+\n")
        path = Path(f.name)

    try:
        # Default format="auto"
        df = parser.load_gtf(path)
        assert df.height == 1
        assert df["gene_id"][0] == "prot1"
        assert df["exp_value"][0] == 0.5
    finally:
        path.unlink()


def test_load_gtf_force_bed6(parser):
    """Test forcing BED6 format."""
    # Create a file that might be ambiguous or have extra columns
    with tempfile.NamedTemporaryFile(suffix=".txt", mode="w", delete=False) as f:
        f.write("chr1\t100\t200\tprot1\t0.5\t+\tignored_col\n")
        path = Path(f.name)

    try:
        # Force BED6
        df = parser.load_gtf(path, format="bed6")
        assert df.height == 1
        # BED6 uses 5th column (index 4) as score
        assert df["exp_value"][0] == 0.5
    finally:
        path.unlink()


def test_load_gtf_force_bed7(parser):
    """Test forcing BED7 format."""
    with tempfile.NamedTemporaryFile(suffix=".txt", mode="w", delete=False) as f:
        f.write("chr1\t100\t200\tprot1\t0\t+\t10.5\n")
        path = Path(f.name)

    try:
        # Force BED7
        df = parser.load_gtf(path, format="bed7")
        assert df.height == 1
        # BED7 uses 7th column (index 6) as exp_value
        assert df["exp_value"][0] == 10.5
    finally:
        path.unlink()


def test_load_gtf_explicit_gtf(parser):
    """Test explicit GTF format."""
    with tempfile.NamedTemporaryFile(suffix=".txt", mode="w", delete=False) as f:
        f.write(
            'chr1\t.\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; RPKM "50.0";\n'
        )
        path = Path(f.name)

    try:
        # Explicit GTF, even with .txt extension
        df = parser.load_gtf(path, format="gtf")
        assert df.height == 1
        assert df["gene_id"][0] == "g1"
        assert df["exp_value"][0] == 50.0
    finally:
        path.unlink()
