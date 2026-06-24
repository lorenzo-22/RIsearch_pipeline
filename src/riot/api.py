"""Public Python API for RIOT.

Mirrors the CLI commands as importable functions that return results in-memory
(Polars DataFrames / dicts) instead of only writing files.

Examples:
    import riot
    df = riot.off_targets(risearch_file="preds.tsv", gtf_file="ann.gtf")
    profiles = riot.accessibility(genome="genome.fa", output="acc/")

Notes:
- off-targets directory-mode input streams results to disk by design and returns a
  summary dict (paths + counts), not a DataFrame; single-file input returns a DataFrame.
- On error the underlying command raises typer.Exit (a Click exception).
- riot.index / riot.search require the external 'risearch' package to be installed
  (the same dependency the CLI's index/search commands need).
"""
from riot.commands.accessibility import run as accessibility
from riot.commands.off_targets import run as off_targets
from riot.commands.risearch import index, search

__all__ = ["off_targets", "accessibility", "index", "search"]
