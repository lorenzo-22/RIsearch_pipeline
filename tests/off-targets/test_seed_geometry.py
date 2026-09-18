"""Seed-geometry and energy-model options on `riot search`.

These expose the RIsearch2 knobs the CLI previously hid: the `-s n:m/l` seed
specification, `--noGUseed` and `-z`. They exist so a seed-geometry ablation can be
run from the CLI rather than by calling the PyO3 bindings directly.
"""

import pytest

from riot.core.risearch import parse_seed_spec
from riot.services.risearch_service import RIsearchError, RIsearchService


class TestParseSeedSpec:
    """`-s` accepts every form the RIsearch2 CLI accepts, so commands port over verbatim."""

    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("6", (None, None, 6)),  # length only — RIsearch2 default form
            ("7", (None, None, 7)),
            ("2:8/7", (2, 8, 7)),  # full form
            ("1:8/8", (1, 8, 8)),
            ("2:8", (2, 8, 7)),  # length defaults to the window width
            ("  2:8/7  ", (2, 8, 7)),  # tolerate surrounding whitespace
        ],
    )
    def test_accepted_forms(self, spec, expected):
        assert parse_seed_spec(spec) == expected

    @pytest.mark.parametrize("spec", ["", "   ", "abc", "2:", ":8", "2-8/7", "2:8/"])
    def test_rejects_malformed(self, spec):
        with pytest.raises(RIsearchError):
            parse_seed_spec(spec)


class TestSeedSpecValidation:
    """Validation happens before the search starts, with a message naming the cause."""

    def _run(self, **kw):
        # Validation precedes every filesystem and PyO3 touch, so the paths are
        # irrelevant here; any of these must raise before they are opened.
        return RIsearchService().run_search(
            query_path=__import__("pathlib").Path(__file__),
            index_path=__import__("pathlib").Path(__file__),
            target_fasta=__import__("pathlib").Path(__file__),
            **kw,
        )

    def test_rejects_seed_longer_than_window(self):
        """RIsearch2's own `-s 2:8/8`: positions 2-8 span 7 nt, so an 8-nt seed cannot fit.

        The C binary prints its complaint and then raises SIGABRT rather than exiting
        cleanly. We refuse up front instead.
        """
        with pytest.raises(RIsearchError, match="exceeds the 7-nt window"):
            self._run(seed_start=2, seed_end=8, seed_length=8)

    def test_accepts_seed_exactly_filling_window(self):
        """`-s 1:8/8` is the documented workaround and must stay valid."""
        with pytest.raises(RIsearchError) as exc:
            self._run(seed_start=1, seed_end=8, seed_length=8)
        assert "exceeds" not in str(exc.value)  # fails later, not on geometry

    def test_rejects_half_specified_window(self):
        with pytest.raises(RIsearchError, match="must be given together"):
            self._run(seed_start=2, seed_end=None, seed_length=7)

    def test_rejects_zero_based_start(self):
        with pytest.raises(RIsearchError, match="1-based"):
            self._run(seed_start=0, seed_end=8, seed_length=7)

    def test_rejects_inverted_window(self):
        with pytest.raises(RIsearchError, match="< seed_start"):
            self._run(seed_start=8, seed_end=2, seed_length=2)

    @pytest.mark.parametrize("bad", ["t05", "T04", "", "turner"])
    def test_rejects_unknown_matrix(self, bad):
        with pytest.raises(RIsearchError, match="must be 't04' or 't99'"):
            self._run(matrix=bad)
