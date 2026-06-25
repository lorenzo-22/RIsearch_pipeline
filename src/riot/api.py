"""Public Python API for RIOT.

Importable functions that return results **in memory** and raise ordinary Python
exceptions — no files are written, nothing is printed, and no ``typer.Exit`` /
Click exceptions leak out. These wrap the pure :mod:`riot.core` layer; the
``riot`` CLI is a separate, file-writing wrapper over the same core.

Examples::

    import riot, polars as pl

    # Single predictions file -> one DataFrame
    df = riot.off_targets(risearch_file="predictions.tsv", gtf_file="ann.gtf")

    # A *directory* of per-siRNA files -> a generator of per-siRNA DataFrames
    for sirna_df in riot.off_targets(risearch_file="preds_dir/"):
        ...
    everything = pl.concat(list(riot.off_targets(risearch_file="preds_dir/")))

    # Accessibility profiles in memory, keyed by chromosome
    profiles = riot.accessibility(genome="genome.fa")   # dict[str, pl.DataFrame]

    # RIsearch index / search
    idx = riot.index("target.fa")                       # Path (binary artifact)
    hits = riot.search("query.fa", idx, target="target.fa")   # pl.DataFrame

Notes:
- ``index`` returns a :class:`~pathlib.Path`: a RIsearch index is a binary on-disk
  artifact, so the path (not in-memory data) is the natural result.
- ``riot.index`` / ``riot.search`` require the external ``risearch`` package, the
  same dependency the CLI's index/search commands need.
"""

from pathlib import Path
from typing import Iterator, Optional, Union

import polars as pl

from riot.core import accessibility as _accessibility
from riot.core import off_targets as _off_targets
from riot.core import risearch as _risearch

__all__ = ["off_targets", "accessibility", "index", "search"]


def _p(value: Optional[Union[str, Path]]) -> Optional[Path]:
    """Coerce a str path (scripting style) to Path; pass Path/None through."""
    return Path(value) if isinstance(value, str) else value


def off_targets(
    risearch_file: Optional[Union[str, Path]] = None,
    sirna_fasta: Optional[Union[str, Path]] = None,
    target_fasta: Optional[Union[str, Path]] = None,
    target_index: Optional[Union[str, Path]] = None,
    gtf_file: Optional[Union[str, Path]] = None,
    feature_type: str = "exon",
    expression_metric: str = "RPKM",
    transcriptome_format: str = "auto",
    accessibility_dir: Optional[Union[str, Path]] = None,
    genome_file: Optional[Union[str, Path]] = None,
    window_size: int = 80,
    max_span: int = 40,
    unpaired_prob: int = 30,
    temperature: float = 37.0,
    on_target_file: Optional[Union[str, Path]] = None,
    on_target_risearch_file: Optional[Union[str, Path]] = None,
    query_file: Optional[Union[str, Path]] = None,
    on_target_expression: float = 1000.0,
    on_target_accessibility: Optional[Union[str, Path]] = None,
    on_target_ids_file: Optional[Union[str, Path]] = None,
    alpha: str = "1.0",
    gamma: str = "1.0",
    theta: str = "",
    sense_only: bool = False,
    predictions_type: str = "gw",
    n_workers: int = 1,
) -> Union[pl.DataFrame, Iterator[pl.DataFrame]]:
    """Analyse siRNA off-target predictions, returning results in memory.

    - A single predictions file (or inline RIsearch via ``sirna_fasta`` +
      ``target_fasta``) returns one :class:`polars.DataFrame`.
    - A *directory* passed to ``risearch_file`` returns a **generator** yielding
      one :class:`polars.DataFrame` per siRNA (consume it fully, or wrap in
      ``contextlib.closing``, for prompt cleanup of the worker pool).

    Writes no files. Raises ``ValueError`` / ``FileNotFoundError`` on bad input.
    """
    rf = _p(risearch_file)

    if rf is not None and rf.is_dir():
        core_gen = _off_targets.compute_off_targets_directory(
            input_dir=rf,
            sirna_fasta=_p(sirna_fasta),
            gtf_file=_p(gtf_file),
            feature_type=feature_type,
            expression_metric=expression_metric,
            transcriptome_format=transcriptome_format,
            accessibility_dir=_p(accessibility_dir),
            temperature=temperature,
            on_target_ids_file=_p(on_target_ids_file),
            on_target_expression=on_target_expression,
            alpha=alpha,
            gamma=gamma,
            theta=theta,
            sense_only=sense_only,
            predictions_type=predictions_type,
            n_workers=n_workers,
        )
        # Yield bare per-siRNA DataFrames (the core yields (df, meta) tuples).
        # Closing this generator propagates GeneratorExit to core_gen, tearing
        # down its worker pool + temp-IPC dir.
        return (frame for frame, _meta in core_gen)

    df, _meta = _off_targets.compute_off_targets_single(
        risearch_file=rf,
        sirna_fasta=_p(sirna_fasta),
        target_fasta=_p(target_fasta),
        target_index=_p(target_index),
        gtf_file=_p(gtf_file),
        feature_type=feature_type,
        expression_metric=expression_metric,
        transcriptome_format=transcriptome_format,
        accessibility_dir=_p(accessibility_dir),
        genome_file=_p(genome_file),
        window_size=window_size,
        max_span=max_span,
        unpaired_prob=unpaired_prob,
        temperature=temperature,
        on_target_file=_p(on_target_file),
        on_target_risearch_file=_p(on_target_risearch_file),
        query_file=_p(query_file),
        on_target_expression=on_target_expression,
        on_target_accessibility=_p(on_target_accessibility),
        on_target_ids_file=_p(on_target_ids_file),
        alpha=alpha,
        gamma=gamma,
        theta=theta,
        sense_only=sense_only,
        predictions_type=predictions_type,
        n_workers=n_workers,
    )
    return df


def accessibility(
    genome: Union[str, Path],
    window_size: int = 80,
    max_span: int = 40,
    unpaired_prob: int = 30,
    temperature: float = 37.0,
) -> dict[str, pl.DataFrame]:
    """Compute per-chromosome accessibility profiles in memory.

    Returns a dict mapping chromosome name to a :class:`polars.DataFrame` with
    columns ``[position, strand, u1..u{unpaired_prob}]`` (both strands stacked).
    Writes no files. Raises ``FileNotFoundError`` if *genome* does not exist.
    """
    return _accessibility.compute_accessibility(
        _p(genome),
        window_size=window_size,
        max_span=max_span,
        unpaired_prob=unpaired_prob,
        temperature=temperature,
    )


def index(
    target: Union[str, Path],
    output: Optional[Union[str, Path]] = None,
) -> Path:
    """Build (or reuse) a RIsearch index and return its :class:`~pathlib.Path`."""
    return _risearch.build_index(_p(target), _p(output))


def search(
    query: Union[str, Path],
    index: Union[str, Path],
    target: Optional[Union[str, Path]] = None,
    seed_length: int = 6,
    max_extension: int = 20,
    energy_threshold: float = -10.0,
) -> pl.DataFrame:
    """Run a RIsearch search and return the hits as a :class:`polars.DataFrame`."""
    return _risearch.run_search(
        query=_p(query),
        index=_p(index),
        target=_p(target),
        seed_length=seed_length,
        max_extension=max_extension,
        energy_threshold=energy_threshold,
    )
