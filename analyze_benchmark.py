#!/usr/bin/env python3
"""Benchmark analysis: new vs old risearch pipeline scaling and speedup."""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.2)
PALETTE = {"new": "#2196F3", "old": "#E53935"}

# ── 1. Parse TSV ──────────────────────────────────────────────────────────────
df = pd.read_csv("benchmark_results.tsv", sep="\t")
main = df[df["pipeline"].isin(["new", "old"])].copy()

agg = (
    main.groupby(["subset_size", "pipeline"])
    .agg(mean_s=("wall_clock_s", "mean"), std_s=("wall_clock_s", "std"),
         mean_rss_kb=("peak_rss_kb", "mean"))
    .reset_index()
)

pivot_s   = agg.pivot(index="subset_size", columns="pipeline", values="mean_s")
pivot_std = agg.pivot(index="subset_size", columns="pipeline", values="std_s")
pivot_rss = agg.pivot(index="subset_size", columns="pipeline", values="mean_rss_kb")

speedup = pivot_s["old"] / pivot_s["new"]
rss_reduction_pct = (pivot_rss["old"] - pivot_rss["new"]) / pivot_rss["old"] * 100

# ── 2. Parse perf log (per-siRNA stage breakdowns for new pipeline) ────────────
# Sections: "=== new n=NNN runR ===" followed by "  stage_name: X.XXXs  (avg Y.Yms/siRNA)"
section_re = re.compile(r"=== new n=(\d+) run(\d+) ===")
stage_re   = re.compile(r"  (\w+):\s+([\d.]+)s\s+\(avg ([\d.]+)ms/siRNA\)")

perf_rows = []
cur_n, cur_run = None, None
with open("benchmark_perf.log") as fh:
    for line in fh:
        m = section_re.match(line)
        if m:
            cur_n, cur_run = int(m.group(1)), int(m.group(2))
            continue
        m = stage_re.match(line)
        if m and cur_n is not None:
            perf_rows.append(dict(n=cur_n, run=cur_run, stage=m.group(1),
                                  total_s=float(m.group(2)),
                                  avg_ms=float(m.group(3))))

perf = pd.DataFrame(perf_rows)
# n=2000 run1 has cold Parquet cache (load_s 10× higher); exclude from stage avg
perf_warm = perf[~((perf["n"] == 2000) & (perf["run"] == 1))]

# Top-level stages only (intersect_s is aggregate; drop its sub-stages + noise)
TOP_STAGES = ["load_s", "intersect_s", "prob_s", "serialize_s",
              "drain_deserialize_s", "drain_summary_s"]
stage_avg = (
    perf_warm[perf_warm["stage"].isin(TOP_STAGES)]
    .groupby(["n", "stage"])["avg_ms"]
    .mean()
    .reset_index()
)
STAGE_LABELS = {
    "load_s":             "Load\n(Parquet)",
    "intersect_s":        "Intersect\n(NCLS)",
    "prob_s":             "Probability\n(partition fn)",
    "serialize_s":        "Serialize\n(output)",
    "drain_deserialize_s":"Drain\ndeserialize",
    "drain_summary_s":    "Drain\nsummary",
}

# ── 3. Plot 1: Scaling curves ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
for pipeline, grp in agg.groupby("pipeline"):
    grp = grp.sort_values("subset_size")
    ax.plot(grp["subset_size"], grp["mean_s"], marker="o", linewidth=2.2,
            label=f"{pipeline.capitalize()} pipeline", color=PALETTE[pipeline])
    ax.fill_between(grp["subset_size"],
                    grp["mean_s"] - grp["std_s"],
                    grp["mean_s"] + grp["std_s"],
                    alpha=0.15, color=PALETTE[pipeline])
ax.set_xlabel("siRNA Dataset Size", fontsize=13)
ax.set_ylabel("Wall-Clock Time (s)", fontsize=13)
ax.set_title("Pipeline Scaling: Wall-Clock Time vs Dataset Size",
             fontsize=14, fontweight="bold")
ax.set_xticks(sorted(agg["subset_size"].unique()))
ax.legend(fontsize=12)
fig.tight_layout()
fig.savefig("benchmark_scaling.png", dpi=150)
plt.close(fig)

# ── 4. Plot 2: Speedup + memory comparison ─────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sizes = sorted(pivot_s.index)
x = np.arange(len(sizes))

bars = ax1.bar(x, [speedup[s] for s in sizes],
               color=["#43A047", "#66BB6A", "#81C784"], edgecolor="white")
for bar, s in zip(bars, sizes):
    ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
             f"{speedup[s]:.1f}×", ha="center", va="bottom",
             fontsize=12, fontweight="bold")
ax1.axhline(1, color="gray", linestyle="--", linewidth=1)
ax1.set_xticks(x)
ax1.set_xticklabels([str(s) for s in sizes])
ax1.set_xlabel("siRNA Dataset Size", fontsize=13)
ax1.set_ylabel("Speedup  (old ÷ new)", fontsize=13)
ax1.set_title("Wall-Clock Speedup", fontsize=14, fontweight="bold")
ax1.set_ylim(0, max(speedup) * 1.2)

w = 0.35
ax2.bar(x - w/2, [pivot_rss.loc[s, "old"] / 1e6 for s in sizes], w,
        label="Old", color="#E53935", alpha=0.85)
ax2.bar(x + w/2, [pivot_rss.loc[s, "new"] / 1e6 for s in sizes], w,
        label="New", color="#2196F3", alpha=0.85)
ax2.set_xticks(x)
ax2.set_xticklabels([str(s) for s in sizes])
ax2.set_xlabel("siRNA Dataset Size", fontsize=13)
ax2.set_ylabel("Peak RSS (GB)", fontsize=13)
ax2.set_title("Peak Memory Usage", fontsize=14, fontweight="bold")
ax2.legend(fontsize=12)

fig.suptitle("New vs Old Pipeline: Speedup & Memory", fontsize=15,
             fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("benchmark_comparison.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ── 5. Plot 3: Per-siRNA stage breakdown (new pipeline) ────────────────────────
stage_pivot = (
    stage_avg.pivot(index="stage", columns="n", values="avg_ms")
    .reindex([s for s in TOP_STAGES if s in stage_avg["stage"].values])
)
sorted_ns = sorted(stage_pivot.columns)
stage_labels = [STAGE_LABELS.get(s, s) for s in stage_pivot.index]

fig, ax = plt.subplots(figsize=(10, 5))
bar_w = 0.22
bar_colors = ["#90CAF9", "#42A5F5", "#1565C0"]
for i, (n, color) in enumerate(zip(sorted_ns, bar_colors)):
    positions = np.arange(len(stage_pivot)) + i * bar_w
    ax.bar(positions, stage_pivot[n], bar_w, label=f"N={n}", color=color)
ax.set_xticks(np.arange(len(stage_pivot)) + bar_w)
ax.set_xticklabels(stage_labels, fontsize=10)
ax.set_ylabel("Avg Time per siRNA (ms)", fontsize=13)
ax.set_title("New Pipeline: Per-siRNA Stage Breakdown (warm cache)",
             fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
fig.tight_layout()
fig.savefig("benchmark_stage_breakdown.png", dpi=150)
plt.close(fig)

# ── 6. Summary ────────────────────────────────────────────────────────────────
print("## Benchmark Results Summary\n")
print(f"{'Size':>6}  {'Old mean (s)':>12}  {'New mean (s)':>12}  "
      f"{'Speedup':>8}  {'Old RSS (GB)':>12}  {'New RSS (GB)':>12}  {'RSS −%':>7}")
print("─" * 80)
for s in sizes:
    print(f"{s:>6}  {pivot_s.loc[s,'old']:>12.1f}  {pivot_s.loc[s,'new']:>12.1f}  "
          f"{speedup[s]:>7.1f}×  {pivot_rss.loc[s,'old']/1e6:>12.2f}  "
          f"{pivot_rss.loc[s,'new']/1e6:>12.2f}  {rss_reduction_pct[s]:>6.1f}%")

max_n = max(sizes)
print(f"\nAt N={max_n}: {speedup[max_n]:.2f}× speedup, "
      f"{rss_reduction_pct[max_n]:.1f}% lower peak RSS.\n")

# Scaling exponents (log-log slope → O(N^k))
for pipeline in ["old", "new"]:
    log_n = np.log(sizes)
    log_t = np.log([pivot_s.loc[s, pipeline] for s in sizes])
    k = np.polyfit(log_n, log_t, 1)[0]
    print(f"{pipeline.capitalize()} pipeline scaling exponent: O(N^{k:.2f})")

print("\nSaved: benchmark_scaling.png  benchmark_comparison.png  "
      "benchmark_stage_breakdown.png")
