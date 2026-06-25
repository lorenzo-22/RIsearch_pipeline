"""Tests for the public Python API (``import riot``).

These guard the Annotated-Typer refactor that makes the CLI command functions
directly callable from plain Python. The key regression is that parameter
defaults must be real values (e.g. ``None``, ``"1.0"``) rather than Typer
``OptionInfo`` / ``ArgumentInfo`` placeholder objects, so that calling the
functions with ordinary keyword arguments works without OptionInfo leakage.
"""

import inspect
from pathlib import Path

import polars as pl
import pytest

import riot

DATA_DIR = Path(__file__).parent / "data"
RISEARCH_FILE = DATA_DIR / "risearch_siRNAID.out"
GENOME_FASTA = DATA_DIR / "genome.fa"


def test_public_api_exposes_callables() -> None:
    """``import riot`` exposes off_targets, accessibility, index, search as callables."""
    for name in ("off_targets", "accessibility", "index", "search"):
        assert hasattr(riot, name), f"riot is missing public attribute {name!r}"
        assert callable(getattr(riot, name)), f"riot.{name} is not callable"


def test_off_targets_returns_dataframe_with_probabilities() -> None:
    """riot.off_targets(...) on a real single-file fixture returns a non-empty
    Polars DataFrame containing the probability column.

    This is the load-bearing regression: if any parameter default were still a
    Typer OptionInfo object, the call would not reach the energy-only probability
    path (it would raise typer.Exit) and there would be no P_off_target column.

    Energy-only path (no transcriptome): the fixture's chrom values
    (transcript_1..5) intersect neither the GTF transcript ids nor chr1, so a
    transcriptome join would zero out the frame. Passing only the predictions
    file exercises the single-siRNA probability calculation directly.
    """
    assert RISEARCH_FILE.exists(), f"Fixture not found: {RISEARCH_FILE}"

    df = riot.off_targets(risearch_file=RISEARCH_FILE)

    assert isinstance(df, pl.DataFrame), f"Expected pl.DataFrame, got {type(df)!r}"
    assert df.height > 0, "Expected a non-empty result DataFrame"
    assert "P_off_target" in df.columns, (
        f"Expected a P_off_target column; got columns: {df.columns}"
    )


def test_off_targets_signature_has_no_typer_placeholder_defaults() -> None:
    """No parameter default of riot.off_targets is a Typer OptionInfo/ArgumentInfo.

    With the modern Annotated idiom, the typer.Option(...) metadata lives in the
    annotation, and the parameter default is the real value. This guards against
    regressing back to ``param = typer.Option(default, ...)``.
    """
    from typer.models import ArgumentInfo, OptionInfo

    sig = inspect.signature(riot.off_targets)
    offenders = []
    for name, param in sig.parameters.items():
        if param.default is inspect.Parameter.empty:
            continue
        if isinstance(param.default, (OptionInfo, ArgumentInfo)):
            offenders.append(name)

    assert not offenders, (
        "These parameters still default to a Typer placeholder object "
        f"(OptionInfo/ArgumentInfo): {offenders}"
    )


def test_accessibility_returns_chrom_to_dataframe(tmp_path: Path) -> None:
    """riot.accessibility(...) returns a dict mapping chromosome -> in-memory
    DataFrame (schema [position, strand, u1..u{u}]) and writes NO files.

    Uses small RNAplfold parameters (W=10, L=5, u=3) suited to the short
    (~71 nt) genome fixture, matching the values used in test_accessibility.py.
    """
    if not GENOME_FASTA.exists():
        pytest.skip(f"Genome FASTA fixture not found: {GENOME_FASTA}")

    files_before = set(tmp_path.iterdir())
    result = riot.accessibility(
        genome=GENOME_FASTA,
        window_size=10,
        max_span=5,
        unpaired_prob=3,
    )

    assert isinstance(result, dict), f"Expected dict, got {type(result)!r}"
    assert result, "Expected at least one chromosome in the result mapping"
    for chrom, df in result.items():
        assert isinstance(chrom, str), f"Expected str chromosome key, got {chrom!r}"
        assert isinstance(df, pl.DataFrame), f"Expected pl.DataFrame, got {type(df)!r}"
        assert df.columns[:2] == ["position", "strand"], (
            f"Unexpected leading columns: {df.columns}"
        )
        assert {"u1", "u2", "u3"}.issubset(df.columns), (
            f"Expected u1..u3 columns; got {df.columns}"
        )
        assert set(df["strand"].unique()) <= {"+", "-"}
        assert df.height > 0

    # The in-memory API must not write anything to disk.
    assert set(tmp_path.iterdir()) == files_before, "accessibility() wrote files"


def test_off_targets_accepts_str_path() -> None:
    """Path parameters accept plain ``str`` (the README/scripting style), not only
    ``Path``. A direct Python call bypasses Typer's str->Path coercion, so the
    command body coerces str paths itself.
    """
    df = riot.off_targets(risearch_file=str(RISEARCH_FILE))
    assert isinstance(df, pl.DataFrame), f"Expected pl.DataFrame, got {type(df)!r}"
    assert df.height > 0, "Expected a non-empty result DataFrame"


def test_off_targets_directory_mode_yields_dataframes(tmp_path: Path) -> None:
    """A *directory* of per-siRNA prediction files returns a generator yielding
    one in-memory DataFrame per siRNA — and writes NO files (the in-memory API).

    This is the contract change from the earlier re-export, which returned a
    summary dict and streamed to disk. File writing is now the CLI's job.
    """
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / RISEARCH_FILE.name).write_bytes(RISEARCH_FILE.read_bytes())

    files_before = set(tmp_path.rglob("*"))
    gen = riot.off_targets(risearch_file=str(in_dir))

    import collections.abc

    assert isinstance(gen, collections.abc.Iterator), (
        f"Expected a generator/iterator, got {type(gen)!r}"
    )

    frames = list(gen)
    assert frames, "Expected at least one per-siRNA DataFrame"
    for df in frames:
        assert isinstance(df, pl.DataFrame), f"Expected pl.DataFrame, got {type(df)!r}"
        assert df.height > 0
        assert "P_off_target" in df.columns

    # The in-memory API must not write anything to disk.
    assert set(tmp_path.rglob("*")) == files_before, "directory mode wrote files"


def test_off_targets_no_input_raises_value_error() -> None:
    """Bad input raises a plain ValueError, not a Typer/Click Exit."""
    import typer

    with pytest.raises(ValueError):
        riot.off_targets()

    # Defensive: the raised exception must not be a Typer/Click Exit.
    try:
        riot.off_targets()
    except typer.Exit:  # pragma: no cover - should never hit
        pytest.fail("API leaked a typer.Exit instead of a plain exception")
    except ValueError:
        pass


def test_accessibility_missing_genome_raises(tmp_path: Path) -> None:
    """A missing genome raises FileNotFoundError (not typer.Exit)."""
    with pytest.raises(FileNotFoundError):
        riot.accessibility(genome=tmp_path / "does_not_exist.fa")
