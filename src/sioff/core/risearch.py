"""Pure core for RIsearch index/search.

Wraps :class:`RIsearchService` without any CLI side-effects (no stdout, no file
writes beyond the index binary itself). Returns in-memory objects.
"""

import re
from pathlib import Path
from typing import Optional, Tuple

import polars as pl

from sioff.services.risearch_service import RIsearchError, RIsearchService

_SEED_SPEC = re.compile(r"^\s*(?:(\d+):(\d+))?(?:/?(\d+))?\s*$")


def parse_seed_spec(spec: str) -> Tuple[Optional[int], Optional[int], int]:
    """Parse a RIsearch2 ``-s`` seed specification into ``(start, end, length)``.

    Accepts the two forms the RIsearch2 CLI accepts, so a command line can be
    copied between the two without translation:

    ``"7"``      length only, no positional constraint -> ``(None, None, 7)``
    ``"2:8/7"``  seed must lie in query positions 2..8 and be 7 nt -> ``(2, 8, 7)``
    ``"2:8"``    positions given, length defaults to the window width -> ``(2, 8, 7)``

    Positions are 1-based and inclusive, matching RIsearch2.
    """
    m = _SEED_SPEC.match(spec)
    if not m or spec.strip() == "":
        raise RIsearchError(
            f"Could not parse seed spec {spec!r}. Expected 'l', 'n:m' or 'n:m/l' "
            f"(e.g. '7', '2:8', '2:8/7')."
        )
    start, end, length = m.group(1), m.group(2), m.group(3)
    if start is None and length is None:
        raise RIsearchError(f"Could not parse seed spec {spec!r}: no length given.")
    if start is None:
        return None, None, int(length)
    s, e = int(start), int(end)
    return s, e, int(length) if length is not None else e - s + 1


def build_index(target: Path, output: Optional[Path] = None) -> Path:
    """Build (or reuse) a RIsearch index for a target FASTA.

    The index is a binary on-disk artifact consumed by :func:`run_search`, so this
    returns its :class:`Path` rather than in-memory data. Reuses an existing index
    if it is newer than the target.
    """
    target = Path(target)
    output = Path(output) if output is not None else None
    return RIsearchService().index_target(target, output)


def run_search(
    query: Path,
    index: Path,
    target: Optional[Path] = None,
    seed_length: int = 6,
    max_extension: int = 20,
    energy_threshold: float = -10.0,
    seed_start: Optional[int] = None,
    seed_end: Optional[int] = None,
    seed_wobble: bool = True,
    matrix: str = "t04",
) -> pl.DataFrame:
    """Run a RIsearch search and return hits as a Polars DataFrame.

    Columns: ``sirna_id, chrom, start, end, strand, energy``. ``target`` is required
    when the index was built outside this process (needed to resolve target names).

    ``seed_start``/``seed_end``/``seed_length`` are the ``-s n:m/l`` seed spec,
    ``seed_wobble=False`` is ``--noGUseed``, and ``matrix`` is ``-z``. See
    :meth:`RIsearchService.run_search` for the measured effect of each.
    """
    return RIsearchService().run_search(
        query_path=Path(query),
        index_path=Path(index),
        target_fasta=Path(target) if target is not None else None,
        seed_length=seed_length,
        max_extension=max_extension,
        energy_threshold=energy_threshold,
        seed_start=seed_start,
        seed_end=seed_end,
        seed_wobble=seed_wobble,
        matrix=matrix,
    )
