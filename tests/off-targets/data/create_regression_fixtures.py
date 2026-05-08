"""
One-time script to generate regression test fixtures from validated pipeline output.

Run from the repo root:
    uv run python tests/off-targets/data/create_regression_fixtures.py

Inputs (must exist on the machine that generated the validated outputs):
    /tmp/compare4/new_out3/output.tsv      – pipeline output for 4 siRNAs
    /tmp/compare4/new_out3/*.summary       – summary files (Z values)

Outputs (committed to the repo):
    tests/off-targets/data/regression_input.parquet    – input for ProbabilityService
    tests/off-targets/data/regression_reference.json   – expected Zoff_noacc values
"""

import json
import re
import pathlib

import polars as pl

INPUT_TSV = pathlib.Path("/tmp/compare4/new_out3/output.tsv")
SUMMARY_DIR = pathlib.Path("/tmp/compare4/new_out3")
OUT_DIR = pathlib.Path(__file__).parent

# Only the columns consumed by ProbabilityService.calculate_probabilities_per_sirna()
KEEP_COLS = [
    "sirna_id",
    "energy",
    "opening_energy",
    "exp_value",
    "raw_e_min",
    "is_on_target",
    "gene_id",
    "transcript_id",
]

print(f"Reading {INPUT_TSV} ...")
df = pl.read_csv(INPUT_TSV, separator="\t").select(KEEP_COLS)
print(f"  {df.shape[0]:,} rows, {df.shape[1]} columns")

out_parquet = OUT_DIR / "regression_input.parquet"
df.write_parquet(out_parquet)
print(f"Wrote {out_parquet} ({out_parquet.stat().st_size // 1024} KB)")

# Parse Zoff_noacc from .summary files
reference: dict[str, dict[str, float]] = {}
for summary_path in sorted(SUMMARY_DIR.glob("*.summary")):
    sid = summary_path.stem
    reference[sid] = {}
    with open(summary_path) as f:
        for line in f:
            m = re.search(
                r"alpha=(\S+) and gamma=(\S+);.*?Zoff_noacc: ([\d.e+\-]+)", line
            )
            if m:
                key = f"alpha={m.group(1)},gamma={m.group(2)}"
                reference[sid][key] = float(m.group(3))
            m2 = re.search(r"For theta=(\S+);.*?Zoff_noacc: ([\d.e+\-]+)", line)
            if m2:
                reference[sid][f"theta={m2.group(1)}"] = float(m2.group(2))
    print(f"  {sid}: {len(reference[sid])} parameter combos")

out_json = OUT_DIR / "regression_reference.json"
with open(out_json, "w") as f:
    json.dump(reference, f, indent=2)
print(f"Wrote {out_json}")
