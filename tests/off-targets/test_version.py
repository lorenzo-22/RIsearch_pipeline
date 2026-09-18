"""`riot --version` and the single source of truth behind it.

Without this flag the only way to find out which RIOT is installed was
`python -c "import riot; print(riot.__version__)"`, which is a poor answer to
the first question on any bug report.

The version must also be *true*: a hardcoded `__version__` silently disagrees
with the installed distribution the moment one is bumped without the other, and
a version report that lies is worse than none.
"""

import subprocess
import sys

import pytest
from typer.testing import CliRunner

import riot
from riot.cli import app

runner = CliRunner()


class TestVersionFlag:
    def test_prints_the_version_and_exits_cleanly(self, plain):
        result = runner.invoke(app, ["--version"])

        assert result.exit_code == 0
        assert riot.__version__ in plain(result.stdout)

    def test_output_names_the_tool_not_just_a_bare_number(self, plain):
        """A bare '0.1.0' in a bug report is ambiguous about what produced it."""
        result = runner.invoke(app, ["--version"])

        assert "riot" in plain(result.stdout).lower()

    def test_works_without_a_subcommand(self, plain):
        """It must not require a subcommand, and must not print the help text."""
        result = runner.invoke(app, ["--version"])

        assert "Usage:" not in plain(result.stdout)

    def test_answers_even_when_another_option_would_fail_validation(self, plain):
        """`--config` has `exists=True`, so a missing path would fail parsing.

        Asking for the version must not be defeated by an unrelated bad argument
        — that is exactly the situation where someone is trying to report a bug.

        This asserts the behaviour, not its mechanism. `is_eager=True` is set on
        the option as the Typer idiom and as insurance against parameter-order
        changes, but it was measured to make no observable difference in this
        CLI's current shape, so no test here can meaningfully pin it.
        """
        result = runner.invoke(app, ["--version", "--config", "/does/not/exist.yaml"])

        assert result.exit_code == 0, plain(result.stdout)
        assert riot.__version__ in plain(result.stdout)

    def test_it_is_advertised_in_the_help(self, plain):
        result = runner.invoke(app, ["--help"])

        assert "--version" in plain(result.stdout)

    def test_reaches_the_real_console_script(self, plain):
        """Guards the installed entry point, not just the in-process app object."""
        result = subprocess.run(
            [sys.executable, "-m", "riot.cli", "--version"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, result.stderr
        assert riot.__version__ in plain(result.stdout)


class TestVersionIsTrue:
    def test_matches_the_installed_distribution_metadata(self, plain):
        """`__version__` must not drift from what was actually installed."""
        from importlib.metadata import PackageNotFoundError, version

        try:
            installed = version("riot-rna")
        except PackageNotFoundError:
            pytest.skip("riot-rna is not installed in this environment")

        assert riot.__version__ == installed

    def test_matches_the_version_declared_in_pyproject(self, plain):
        import re
        from pathlib import Path

        pyproject = (Path(__file__).resolve().parents[2] / "pyproject.toml").read_text()
        declared = re.search(r'(?m)^version\s*=\s*"([^"]+)"', pyproject)

        assert declared is not None, "pyproject.toml has no version"
        assert riot.__version__ == declared.group(1)
