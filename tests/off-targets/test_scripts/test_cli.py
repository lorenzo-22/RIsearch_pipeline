"""Tests for the CLI module."""

from pathlib import Path

from typer.testing import CliRunner

from RIsearch_pipeline.cli import app

runner = CliRunner()


class TestCLIMain:
    """Tests for the main command."""

    def test_run_valid_file(self) -> None:
        """CLI loads valid TSV and prints summary."""
        tsv_path = Path(__file__).parent.parent / "input" / "risearch_siRNAID.out"
        result = runner.invoke(app, ["off-targets", "-r", str(tsv_path)])

        assert result.exit_code == 0
        assert "60 predictions" in result.stdout
        assert "Energy range" in result.stdout

    def test_run_missing_file(self) -> None:
        """CLI exits with error for missing file."""
        result = runner.invoke(app, ["off-targets", "-r", "/nonexistent/file.tsv"])

        assert result.exit_code != 0

    def test_run_shows_chromosomes(self) -> None:
        """CLI output includes chromosome information."""
        # Use verbose to see the dataframe head where chromosomes/targets are listed
        tsv_path = Path(__file__).parent.parent / "input" / "risearch_siRNAID.out"
        result = runner.invoke(app, ["off-targets", "-r", str(tsv_path), "-v"])

        assert "transcript_3" in result.stdout

    def test_run_shows_strands(self) -> None:
        """CLI output includes strand information."""
        tsv_path = Path(__file__).parent.parent / "input" / "risearch_siRNAID.out"
        result = runner.invoke(app, ["off-targets", "-r", str(tsv_path), "-v"])

        assert "+" in result.stdout
