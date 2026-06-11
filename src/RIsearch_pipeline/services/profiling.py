"""Pipeline profiling utilities.

Lightweight wall-time + RSS-memory profiler that requires no third-party
dependencies (uses /proc/self/status on Linux, time.perf_counter everywhere).

Usage::

    from RIsearch_pipeline.services.profiling import PipelineProfiler

    profiler = PipelineProfiler(enabled=True)

    with profiler.stage("Intersection", rows_in=df.height) as stage:
        result = intersect(df, df_trans)
        stage.rows_out = result.height

    profiler.print_summary(console)
"""

import contextlib
import time
from dataclasses import dataclass
from typing import Generator

from rich.table import Table


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _rss_mb() -> float:
    """Return process RSS in MB via /proc/self/status (Linux).

    Returns 0.0 on platforms that don't have /proc (macOS, Windows).
    """
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except Exception:
        pass
    return 0.0


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class StageRecord:
    name: str
    elapsed: float = 0.0
    mem_before_mb: float = 0.0
    mem_after_mb: float = 0.0
    rows_in: int = 0
    rows_out: int = 0


class _StageContext:
    """Mutable handle yielded inside a ``with profiler.stage(...)`` block.

    The caller can set ``stage.rows_out`` after the work is done so the
    profiler can report throughput without requiring the caller to track timing.
    """

    def __init__(self, record: StageRecord) -> None:
        self._record = record

    @property
    def rows_out(self) -> int:
        return self._record.rows_out

    @rows_out.setter
    def rows_out(self, value: int) -> None:
        self._record.rows_out = value


# ---------------------------------------------------------------------------
# Main profiler
# ---------------------------------------------------------------------------

class PipelineProfiler:
    """Collects per-stage timing and memory deltas, then renders a Rich table.

    When ``enabled=False`` every call is a no-op with negligible overhead.
    """

    def __init__(self, enabled: bool = True) -> None:
        self.enabled = enabled
        self._stages: list[StageRecord] = []

    @contextlib.contextmanager
    def stage(
        self, name: str, rows_in: int = 0
    ) -> Generator[_StageContext, None, None]:
        """Context manager that records wall time and RSS delta for one stage.

        Example::

            with profiler.stage("Intersection", rows_in=df.height) as stage:
                result = intersect(df, df_trans)
                stage.rows_out = result.height
        """
        record = StageRecord(name=name, rows_in=rows_in)
        ctx = _StageContext(record)

        if not self.enabled:
            yield ctx
            return

        record.mem_before_mb = _rss_mb()
        t0 = time.perf_counter()
        try:
            yield ctx
        finally:
            record.elapsed = time.perf_counter() - t0
            record.mem_after_mb = _rss_mb()
            self._stages.append(record)

    def print_summary(self, console) -> None:
        """Print a Rich table summarising all recorded stages.

        Does nothing if the profiler is disabled or no stages were recorded.
        """
        if not self.enabled or not self._stages:
            return

        total_time = sum(s.elapsed for s in self._stages)

        table = Table(
            title="Pipeline Profile",
            show_header=True,
            header_style="bold cyan",
            show_lines=False,
        )
        table.add_column("Stage", style="bold", no_wrap=True)
        table.add_column("Time (s)", justify="right")
        table.add_column("% Total", justify="right")
        table.add_column("Rows in → out", justify="right")
        table.add_column("RSS (MB)", justify="right")
        table.add_column("Throughput (K rows/s)", justify="right")

        peak_rss = max((s.mem_after_mb for s in self._stages), default=0.0)

        for s in self._stages:
            pct = (s.elapsed / total_time * 100) if total_time > 0 else 0.0
            delta_mem = s.mem_after_mb - s.mem_before_mb
            # Throughput is based on whichever row count is non-zero
            ref_rows = s.rows_out if s.rows_out > 0 else s.rows_in
            throughput = (ref_rows / s.elapsed / 1_000) if s.elapsed > 0 and ref_rows > 0 else 0.0

            if s.rows_out > 0 and s.rows_in > 0:
                rows_str = f"{s.rows_in:,} → {s.rows_out:,}"
            elif s.rows_in > 0:
                rows_str = f"{s.rows_in:,}"
            elif s.rows_out > 0:
                rows_str = f"→ {s.rows_out:,}"
            else:
                rows_str = "–"

            # Show absolute RSS after stage + delta in parentheses.
            # RSS is always shown; delta shown only when ≥ 0.5 MB to avoid noise.
            rss_str = f"{s.mem_after_mb:.0f}"
            if abs(delta_mem) >= 0.5:
                sign = "+" if delta_mem > 0 else ""
                rss_str += f" ({sign}{delta_mem:.0f})"
            tput_str = f"{throughput:.1f}" if throughput > 0 else "–"

            table.add_row(
                s.name,
                f"{s.elapsed:.2f}",
                f"{pct:.1f}%",
                rows_str,
                rss_str,
                tput_str,
            )

        # Separator + total (peak RSS shown in the RSS column)
        table.add_section()
        table.add_row(
            "[bold]TOTAL[/bold]",
            f"[bold]{total_time:.2f}[/bold]",
            "[bold]100%[/bold]",
            "",
            f"[bold]peak {peak_rss:.0f}[/bold]",
            "",
        )

        console.print(table)

    def to_dict(self) -> dict:
        """Return profile data as a plain dict (useful for JSON export)."""
        return {
            "stages": [
                {
                    "name": s.name,
                    "elapsed_s": round(s.elapsed, 4),
                    "rss_before_mb": round(s.mem_before_mb, 1),
                    "rss_after_mb": round(s.mem_after_mb, 1),
                    "rss_delta_mb": round(s.mem_after_mb - s.mem_before_mb, 1),
                    "rows_in": s.rows_in,
                    "rows_out": s.rows_out,
                }
                for s in self._stages
            ],
            "total_elapsed_s": round(sum(s.elapsed for s in self._stages), 4),
        }
