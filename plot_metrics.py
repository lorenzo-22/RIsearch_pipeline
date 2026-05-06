#!/usr/bin/env python3
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.2)
PALETTE = {"new": "#2196F3", "old": "#E53935"}

df = pd.read_csv("benchmark_results.csv")
df["peak_memory_gb"] = df["peak_memory_kb"] / 1_048_576

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

sns.lineplot(
    data=df,
    x="siRNA_count",
    y="duration_seconds",
    hue="pipeline_version",
    marker="o",
    palette=PALETTE,
    errorbar="ci",
    ax=ax1,
)
ax1.set_xlabel("siRNA Count", fontsize=13)
ax1.set_ylabel("Duration (s)", fontsize=13)
ax1.set_title("Runtime vs. siRNA Count", fontsize=14, fontweight="bold")
ax1.set_xticks(sorted(df["siRNA_count"].unique()))
ax1.legend(title="Pipeline", fontsize=11)

sns.lineplot(
    data=df,
    x="siRNA_count",
    y="peak_memory_gb",
    hue="pipeline_version",
    marker="o",
    palette=PALETTE,
    errorbar="ci",
    ax=ax2,
)
ax2.set_xlabel("siRNA Count", fontsize=13)
ax2.set_ylabel("Peak Memory (GB)", fontsize=13)
ax2.set_title("Peak Memory vs. siRNA Count", fontsize=14, fontweight="bold")
ax2.set_xticks(sorted(df["siRNA_count"].unique()))
ax2.legend(title="Pipeline", fontsize=11)

fig.suptitle("New vs Old Pipeline: Performance Comparison", fontsize=15,
             fontweight="bold")
fig.tight_layout()
fig.savefig("pipeline_comparison.png", dpi=150, bbox_inches="tight")
print("Saved: pipeline_comparison.png")
