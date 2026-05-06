#!/usr/bin/env python3
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_theme(style="whitegrid", font_scale=1.2)
PALETTE = {"new": "#2196F3", "old": "#E53935"}

df = pd.read_csv("benchmark_results.csv")
df["peak_memory_gb"] = df["peak_memory_kb"] / 1_048_576

# t(0.975, df) critical values; fallback to 1.96 for large n
T_CRITS = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
           6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}

# Per-group stats for annotations and CI upper bounds
def ci95_upper(x):
    n = len(x)
    t = T_CRITS.get(n - 1, 1.96)
    return x.mean() + t * x.std(ddof=1) / np.sqrt(n)

agg = df.groupby(["pipeline_version", "siRNA_count"]).agg(
    mean_dur=("duration_seconds", "mean"),
    ci_upper_dur=("duration_seconds", ci95_upper),
    mean_mem=("peak_memory_gb", "mean"),
).reset_index()

means = agg.pivot(index="siRNA_count", columns="pipeline_version",
                  values=["mean_dur", "ci_upper_dur", "mean_mem"])
speedup = (means["mean_dur"]["old"] / means["mean_dur"]["new"]).to_dict()
ci_upper_new = means["ci_upper_dur"]["new"].to_dict()
mean_new_dur = means["mean_dur"]["new"].to_dict()
mean_new_mem = means["mean_mem"]["new"].to_dict()
mean_old_mem = means["mean_mem"]["old"].to_dict()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# ── Runtime plot ──────────────────────────────────────────────────────────────
sns.lineplot(
    data=df,
    x="siRNA_count",
    y="duration_seconds",
    hue="pipeline_version",
    marker="o",
    palette=PALETTE,
    errorbar=("ci", 95),
    err_style="bars",
    err_kws={"capsize": 5, "capthick": 1.5, "elinewidth": 1.5},
    ax=ax1,
)
ax1.set_xlabel("siRNA Count", fontsize=13)
ax1.set_ylabel("Duration (s)", fontsize=13)
ax1.set_title("Runtime vs. siRNA Count", fontsize=14, fontweight="bold")
ax1.set_xticks(sorted(df["siRNA_count"].unique()))
ax1.legend(title="Pipeline", fontsize=11)

# Speedup annotations above 'new' CI upper bound
y_pad = (ax1.get_ylim()[1] - ax1.get_ylim()[0]) * 0.02
for n in sorted(speedup):
    ax1.annotate(
        f"{speedup[n]:.1f}×",
        xy=(n, mean_new_dur[n]),
        xytext=(n, ci_upper_new[n] + y_pad),
        ha="center", va="bottom",
        fontsize=10, fontweight="bold", color="#1a6e2e",
        arrowprops=dict(arrowstyle="-", color="#1a6e2e", lw=0.8),
    )

# Expand ylim so annotations don't clip
ax1.set_ylim(bottom=0)
_, cur_top = ax1.get_ylim()
ax1.set_ylim(top=cur_top * 1.15)

# ── Memory plot ───────────────────────────────────────────────────────────────
sns.lineplot(
    data=df,
    x="siRNA_count",
    y="peak_memory_gb",
    hue="pipeline_version",
    marker="o",
    palette=PALETTE,
    errorbar=("ci", 95),
    err_style="bars",
    err_kws={"capsize": 5, "capthick": 1.5, "elinewidth": 1.5},
    ax=ax2,
)
ax2.set_xlabel("siRNA Count", fontsize=13)
ax2.set_ylabel("Peak Memory (GB)", fontsize=13)
ax2.set_title("Peak Memory vs. siRNA Count", fontsize=14, fontweight="bold")
ax2.set_xticks(sorted(df["siRNA_count"].unique()))
ax2.legend(title="Pipeline", fontsize=11)
ax2.set_ylim(bottom=0)

fig.suptitle("New vs Old Pipeline: Performance Comparison", fontsize=15,
             fontweight="bold")
fig.tight_layout()
fig.savefig("pipeline_comparison.png", dpi=150, bbox_inches="tight")
print("Saved: pipeline_comparison.png")
