"""siRNA off-target discovery pipeline."""

from importlib.metadata import PackageNotFoundError, version as _dist_version

from sioff.api import accessibility, index, off_targets, search

try:
    # Read the installed distribution so `sioff --version` reports what is
    # actually installed. A literal here silently disagrees with the wheel the
    # moment one is bumped without the other, and a version that lies is worse
    # than none.
    __version__ = _dist_version("sioff")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.1.0"

__all__ = ["off_targets", "accessibility", "index", "search", "__version__"]
