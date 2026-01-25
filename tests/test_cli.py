"""Tests for the CLI module."""

from pathlib import Path

from typer.testing import CliRunner

from RIsearch_pipeline.cli import app

runner = CliRunner()


class TestCLIMain:
    """Tests for the main command."""

    def test_run_valid_file(self) -> None:
        """CLI loads valid TSV and prints summary."""
        tsv_path = Path(__file__).parent / "input" / "risearch_output.tsv"
        result = runner.invoke(app, ["-r", str(tsv_path)])

        assert result.exit_code == 0
        assert "71 predictions" in result.stdout
        assert "Energy range" in result.stdout

    def test_run_missing_file(self) -> None:
        """CLI exits with error for missing file."""
        result = runner.invoke(app, ["-r", "/nonexistent/file.tsv"])

        assert result.exit_code != 0

    def test_run_shows_chromosomes(self) -> None:
        """CLI output includes chromosome information."""
        tsv_path = Path(__file__).parent / "input" / "risearch_output.tsv"
        result = runner.invoke(app, ["-r", str(tsv_path)])

        assert "chr1" in result.stdout
        assert "chr2" in result.stdout

    def test_run_shows_strands(self) -> None:
        """CLI output includes strand information."""
        tsv_path = Path(__file__).parent / "input" / "risearch_output.tsv"
        result = runner.invoke(app, ["-r", str(tsv_path)])

        assert "+" in result.stdout
        assert "-" in result.stdout
