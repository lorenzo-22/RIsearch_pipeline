# Pipeline Refactoring Plan

## Phase 1 — Local Validation [DONE]

- [x] Audit directory structure, CLI, and config schemas
- [x] Fix CPU compatibility: replace `polars>=1.0` with `polars[rtcompat]>=1.0`
- [x] Fix `risearch` git dep: pin to buildable commit `5242668c` (commit `69aa6d76` has Rust
      compilation bug — `from_fastas` removed from core but not from Python bindings)
- [x] Fix `config/off-targets.example.yaml`: update paths to `tests/off-targets/data/`
- [x] Verify CLI mode (`risearch-pipeline off-targets ...`) runs on test data
- [x] Verify config mode (`risearch-pipeline -c config/off-targets.example.yaml`) runs
- [x] 95/95 tests pass
- [x] Commit and push Phase 1 changes

## Phase 2 — Slurm Orchestrator [DONE]

Pipeline stages in dependency order:
1. (optional) `risearch-pipeline index` — build RIsearch index
2. (optional) `convert_risearch_to_parquet.py` — convert `.out.gz` → per-siRNA `.parquet`
3. `risearch-pipeline accessibility` — compute per-chromosome accessibility profiles
4. `risearch-pipeline off-targets` — main analysis (depends on step 3)

Tasks:

- [x] Create `run_pipeline.py` orchestrator script with `--slurm` flag
  - Local mode: run steps sequentially, capture stdout/stderr to structured log
  - Slurm mode: wrap each step in `sbatch --parsable`, chain deps via `--dependency=afterok:<job_id>`
- [x] Expose per-step Slurm resources at script top or via CLI flags:
  - `--partition`, `--time`, `--mem`, `--cpus-per-task`, `--account`
- [x] Support `--dry-run`: print sbatch commands without submitting
- [x] Support `--config` passthrough to underlying pipeline commands
- [x] Generate structured logs: `logs/<timestamp>/{step}.log` for local; job IDs + status for Slurm
- [x] Test `--dry-run` output locally (no Slurm scheduler available on dev machine)
- [ ] Commit and push Phase 2 changes

## Phase 3 — Documentation [TODO]

- [ ] Update `README.md`:
  - Local run instructions (CLI and config modes)
  - Cluster run instructions (`run_pipeline.py --slurm`)
  - Environment setup (`uv venv && uv sync`)
  - Note on risearch git dep and pinned commit
- [ ] Commit and push Phase 3 changes

---

## Notes

- `uv sync` now works cleanly. Fresh `uv venv && uv sync` will rebuild risearch from
  commit `5242668c` which compiles correctly.
- Smoke tests show 0 intersections in test data (chrom name mismatch between
  `risearch_siRNAID.out` and `expression_data.gtf` in gw mode) — functional tests rely
  on the regression fixture at `tests/off-targets/data/regression_input.parquet`.
- No Slurm scheduler on dev machine — Slurm mode tested via `--dry-run` only.
