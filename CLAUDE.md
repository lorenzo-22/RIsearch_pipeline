# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Setup
```bash
uv venv && source .venv/bin/activate
uv sync
```
Requires Python ≥3.14 and ViennaRNA 2.7.2 (installed via `uv sync`). The `risearch` Python package is a local file dependency pinned to `file:///dev/shm/src/risearch/risearch-python` (PyO3 bindings built from the `risearch/` git submodule).

### Run the pipeline
```bash
risearch-pipeline --help
risearch-pipeline off-targets --risearch-file predictions.tsv --transcriptome annotations.gtf ...
risearch-pipeline -c config/off-targets.example.yaml   # YAML config mode
```
The `gw` entrypoint is an alias for `risearch-pipeline`. `--risearch-file` can be a single TSV/`.out.gz` file **or a directory of `.parquet` files** (one per siRNA, produced by `convert_risearch_to_parquet.py`).

### Lint & type check
```bash
ruff check src/
ruff format src/
# Type checking (two options, both in dev dependencies)
ty check src/
pyrefly check src/
```

### Tests
```bash
pytest tests/
pytest tests/off-targets/test_scripts/test_probability.py   # single test file
```
Test data lives in `tests/off-targets/` (unit tests).

The golden regression tests (`test_regression_golden.py`) require pre-generated fixtures. If `tests/off-targets/fixtures/regression_input.parquet` is missing, regenerate them:
```bash
python tests/off-targets/fixtures/create_regression_fixtures.py
```

### Profile
```bash
scalene risearch-pipeline off-targets ...
```

### Convert RIsearch output to Parquet
```bash
python convert_risearch_to_parquet.py <results_dir> --ids-file sirna_ids.txt \
    --out-dir output_parquet/ --workers 32
```
Converts `.out.gz` files (one per siRNA) to per-siRNA `.parquet` files. The new pipeline reads these via memory-mapped Arrow IPC, avoiding re-parsing TSV on every run.

## Architecture

### Layer structure

```
cli.py                          ← Typer app, --config YAML dispatch
commands/{off_targets,accessibility}.py  ← one module per subcommand
services/                       ← stateless processing units
  risearch_service.py           ← wraps risearch PyO3 bindings (index + search, in-process)
  risearch_parser.py            ← parses RIsearch TSV/parquet output into Polars DataFrames
  transcriptome_parser.py       ← GTF/GFF/BED → Polars DataFrame with exp_value
  intersection_service.py       ← joins predictions to annotations (gw or tw mode)
  accessibility.py              ← RNAplfold wrapper + memmap binary profile lookup
  probability.py                ← partition-function probability calculation
  helpers.py / profiling.py
models.py                       ← shared dataclasses / typed structs
config.py                       ← OmegaConf schema + YAML loading + path resolution
```

### The `risearch` submodule

`risearch/` is a git submodule containing the Rust RIsearch binary **and** its PyO3 Python bindings (`risearch-python/`). `RIsearchService` calls the bindings **in-process** — there is no subprocess or intermediate TSV. The binary at `~/.cargo/bin/RIsearch` is only needed for the legacy CLI path.

### Data flow for `off-targets`

1. `RIsearchParser.load()` → Polars DataFrame of predictions
2. `TranscriptomeParser` → DataFrame with `transcript_id` + `exp_value`
3. `IntersectionService.intersect()` → join on chromosome (gw) or transcript ID (tw)
4. `AccessibilityService._annotate_opening_energy()` → add `opening_energy` via memmap lookup
5. `ProbabilityService.calculate_probabilities_per_sirna()` → add `P_off_target` columns
6. Write TSV (and optionally `.results` legacy format)

### Key design decisions

- **Polars throughout**: all DataFrames use Polars (not pandas) for vectorized, Rust-backed operations. Use lazy evaluation (`scan_csv`, `lazy()`) where possible.
- **Memory-mapped accessibility profiles**: binary `<chrom>.open.acc.bin` files store `uint8` values (energy × 10). Accessed via `np.memmap` — never load them fully into RAM.
- **Per-siRNA partition functions**: `Z_s` is computed independently for each siRNA. Do not mix predictions from different siRNAs into one partition function.
- **Alpha/gamma/theta parameter sweeps** are computed in a single Polars `group_by` pass — all parameter combinations produce separate `P_off_target:alpha=X,gamma=Y` columns in one aggregation rather than N separate passes.
- **Config system**: YAML configs are loaded via OmegaConf structured schemas (`PipelineConfig` → `OffTargetsConfig`/etc.). Relative paths in configs are resolved relative to the config file's directory, not the working directory.
- **CLI key mapping**: `config.py:config_to_kwargs()` translates YAML field names (e.g., `transcriptome`, `window`) to the Typer function parameter names (`gtf_file`, `window_size`). When adding new CLI params, update both `OffTargetsConfig` and `config_to_kwargs()`.
- **On-target detection**: pass `--on-target-id-map/-oi <tsv>` (columns: `sirna_id`, `gene_id`) to identify which rows in the predictions DataFrame are on-target hits. The on-target weight `W_on` is subtracted from `Z_s` to yield `Zoff`.

### Predictions modes
- **gw (genome-wide)**: RIsearch ran against a genome; predictions have genomic coordinates; intersection uses chromosome join.
- **tw (transcriptome-wide)**: RIsearch ran against transcript sequences; `chrom` column matches `transcript_id` directly.
