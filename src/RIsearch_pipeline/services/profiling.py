import time
import contextlib
from collections.abc import Generator
from dataclasses import dataclass


def _rss_mb() -> float:
    """RSS in MB via /proc/self/status; returns 0.0 on macOS/Windows."""
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except Exception:
        pass
    return 0.0


@dataclass
class StageRecord:
    name: str
    elapsed: float = 0.0
    mem_before_mb: float = 0.0
    mem_after_mb: float = 0.0
    rows_in: int = 0
    rows_out: int = 0


class _StageContext:
    """Mutable handle yielded inside ``with profiler.stage(...)``.

    Caller sets ``stage.rows_out`` after work completes; profiler picks it up.
    """

    def __init__(self, record: StageRecord) -> None:
        self._record = record

    @property
    def rows_out(self) -> int:
        return self._record.rows_out

    @rows_out.setter
    def rows_out(self, value: int) -> None:
        self._record.rows_out = value


class PipelineProfiler:
    """Per-stage wall-time + RSS profiler. ``enabled=False`` → no-op."""

    def __init__(self, enabled: bool = True) -> None:
        self.enabled = enabled
        self._stages: list[StageRecord] = []

    @contextlib.contextmanager
    def stage(self, name: str, rows_in: int = 0) -> Generator[_StageContext, None, None]:
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
        if not self.enabled or not self._stages:
            return

        from rich.table import Table

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

            # delta suppressed below 0.5 MB to avoid noise from normal allocator fluctuation
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
