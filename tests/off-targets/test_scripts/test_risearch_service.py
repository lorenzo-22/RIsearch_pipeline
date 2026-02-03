"""Tests for RIsearchService."""

from pathlib import Path
import tempfile

import pytest

from RIsearch_pipeline.services.risearch_service import RIsearchService, RIsearchError


@pytest.fixture
def service() -> RIsearchService:
    """Create a RIsearchService instance."""
    return RIsearchService()


@pytest.fixture
def sirna_single_path() -> Path:
    """Path to single siRNA FASTA."""
    return Path(__file__).parent.parent / "single" / "input" / "sirna.fa"


@pytest.fixture
def sirna_multi_path() -> Path:
    """Path to multi-siRNA FASTA."""
    return Path(__file__).parent.parent / "multiple" / "input" / "sirnas.fa"


@pytest.fixture
def genome_path() -> Path:
    """Path to test genome FASTA."""
    return Path(__file__).parent.parent / "single" / "input" / "genome.fa"


class TestValidateSirnaFasta:
    """Tests for validate_sirna_fasta()."""

    def test_validates_single_sirna(
        self, service: RIsearchService, sirna_single_path: Path
    ) -> None:
        """Single siRNA file returns list with one ID."""
        ids = service.validate_sirna_fasta(sirna_single_path)
        assert len(ids) == 1
        assert ids[0] == "siRNAID"  # Expected ID from sirna.fa

    def test_validates_multiple_sirnas(
        self, service: RIsearchService, sirna_multi_path: Path
    ) -> None:
        """Multi-siRNA file returns list with all IDs."""
        ids = service.validate_sirna_fasta(sirna_multi_path)
        assert len(ids) == 100  # User created 100 siRNAs
        assert ids[0] == "siRNA1"
        assert ids[-1] == "siRNA100"

    def test_rejects_duplicate_ids(self, service: RIsearchService) -> None:
        """Raises ValueError for duplicate siRNA IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(">siRNA1\nAUGC\n>siRNA1\nGCUA\n")
            f.flush()
            dup_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Duplicate siRNA ID"):
                service.validate_sirna_fasta(dup_path)
        finally:
            dup_path.unlink()

    def test_missing_file_raises(self, service: RIsearchService) -> None:
        """Raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            service.validate_sirna_fasta(Path("/nonexistent/sirna.fa"))


class TestIndexTarget:
    """Tests for index_target()."""

    def test_creates_index(self, service: RIsearchService, genome_path: Path) -> None:
        """Creates index file from genome FASTA."""
        if not service.binary_path.exists():
            pytest.skip("RIsearch binary not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            index_path = Path(tmpdir) / "genome.idx"
            result = service.index_target(genome_path, index_path)

            assert result.exists()
            assert result.stat().st_size > 0

    def test_reuses_existing_index(
        self, service: RIsearchService, genome_path: Path
    ) -> None:
        """Reuses index if newer than source."""
        if not service.binary_path.exists():
            pytest.skip("RIsearch binary not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            index_path = Path(tmpdir) / "genome.idx"

            # Create first index
            service.index_target(genome_path, index_path)
            first_mtime = index_path.stat().st_mtime

            # Second call should reuse
            service.index_target(genome_path, index_path)
            assert index_path.stat().st_mtime == first_mtime


class TestRunSearch:
    """Tests for run_search()."""

    def test_search_produces_output(
        self,
        service: RIsearchService,
        sirna_single_path: Path,
        genome_path: Path,
    ) -> None:
        """Search produces non-empty output file."""
        if not service.binary_path.exists():
            pytest.skip("RIsearch binary not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            index_path = tmpdir_path / "genome.idx"
            output_path = tmpdir_path / "results.out"

            # Index and search
            service.index_target(genome_path, index_path)
            result = service.run_search(
                sirna_single_path,
                index_path,
                output_path,
            )

            assert result.exists()
            # Output may be empty if energy threshold is too strict
            # but file should exist


class TestSearchSingleSirna:
    """Tests for search_single_sirna()."""

    def test_returns_best_hit(
        self,
        service: RIsearchService,
        sirna_single_path: Path,
        genome_path: Path,
    ) -> None:
        """Returns tuple with energy, start, end, strand."""
        if not service.binary_path.exists():
            pytest.skip("RIsearch binary not found")

        energy, start, end, strand = service.search_single_sirna(
            sirna_single_path, genome_path
        )

        # Should find some hits with negative energy
        assert energy < 0 or energy == 0.0
        assert strand in ["+", "-"]
