# RIsearch Pipeline

A bioinformatics pipeline for **siRNA off-target discovery and probability quantification**. Integrates RNA-RNA interaction predictions with transcriptome annotations, RNA accessibility profiling, and thermodynamic modeling to rank off-target binding sites.

---

## Pipeline Workflow

```mermaid
flowchart TD
    A([siRNA FASTA]) --> B[RIsearchService\nin-process PyO3 bindings]
    C([Target FASTA / Index]) --> B
    B --> D[(Predictions\nTSV / Parquet)]
    D --> E[RIsearchParser]

    E --> F[Predictions DataFrame]
    G([Annotation\nGTF / BED]) --> H[AnnotationParser]
    H --> I[Annotations + Expression]

    F --> J[IntersectionService\ngw: chrom join · tw: transcript_id match]
    I --> J

    K([Accessibility Profiles\n.parquet per chrom]) --> L[GenomeAccessibilityService\nRNA.pfl_fold_up lookup]
    J --> L

    L --> M[ProbabilityService\nper-siRNA partition function]
    N([On-target ID map\noptional]) --> M

    M --> O([TSV output\nP_off_target columns])
    M --> P([.results legacy format\noptional])

    style A fill:#d4edda,stroke:#28a745
    style C fill:#d4edda,stroke:#28a745
    style G fill:#d4edda,stroke:#28a745
    style K fill:#d4edda,stroke:#28a745
    style N fill:#d4edda,stroke:#28a745
    style O fill:#cce5ff,stroke:#004085
    style P fill:#cce5ff,stroke:#004085
```

---

## CLI Reference

### `off-targets`

| Flag | Description |
|------|-------------|
| `-r / --risearch-file` | Pre-computed predictions (TSV, `.out.gz`, or directory of Parquet files — directory triggers parallel per-siRNA mode) |
| `-s / --sirna-fasta` | siRNA FASTA — runs RIsearch in-process via PyO3 bindings |
| `--target-fasta / --genome` | Target FASTA for in-process RIsearch |
| `-t / --transcriptome` | GTF/GFF3 or BED annotation file |
| `-a / --accessibility-dir` | Directory of per-chromosome accessibility Parquet files (from `accessibility` command) |
| `--expression-metric` | GTF attribute for expression weighting (default: `RPKM`) |
| `--type` | `gw` (genome-wide) or `tw` (transcriptome-wide, default: `gw`) |
| `--alpha / --gamma / --theta` | Parameter sweep values (semicolon-separated, e.g. `0.8;1.0`) |
| `--on-target-ids / -oi` | TSV mapping `sirna_id → transcript_id` for on-target normalization |
| `-j / --workers` | Parallel worker processes (default: CPU count) |
| `-o / --output` | Output file path (TSV by default) |

### `accessibility`

| Flag | Description |
|------|-------------|
| `-f / --fasta` | Genome or transcriptome FASTA |
| `-o / --output` | Output directory (one `{chrom}.accessibility.parquet` per chromosome) |
| `-W / --window` | RNAplfold window size W (default: 80) |
| `-L / --span` | Max base-pair span L (default: 40) |
| `-u / --unpaired` | Unpaired probability length u (default: 30) |
| `-T / --temperature` | Folding temperature °C (default: 37.0) |
| `-j / --workers` | Parallel workers, one per chromosome (default: 1) |

---

## Installation

Requires **Python ≥ 3.14** and **ViennaRNA 2.7.2**.

```bash
git clone git@github.com:lorenzo-22/RIsearch_pipeline.git
cd RIsearch_pipeline

# Create virtual environment and install dependencies
uv venv && source .venv/bin/activate
uv sync
```

### The `risearch` dependency

The `risearch` PyO3 bindings are only required for the **in-process `index` and
`search`** commands (computing RNA-RNA interaction predictions in-process). The
core off-target analysis — `off-targets` and `accessibility` running on
**pre-computed** RIsearch output (TSV / `.out.gz` / Parquet) — works **without**
`risearch` installed; it is imported lazily.

`risearch` is currently fetched from a **private** repository over SSH and is
**not yet on PyPI**, so `uv sync` requires SSH access to that repo:

```
git+ssh://git@github.com/saiden89/risearch.git@5242668c…#subdirectory=risearch-python
```

The commit is pinned because a later commit (`69aa6d7`) removed `from_fastas`
from the Rust core without updating the Python bindings, breaking compilation.
Update the pin only when upstream fixes the mismatch.

> **PyPI:** Because this `git+ssh` dependency is a direct URL, the package
> cannot be published to PyPI as-is. Once `risearch` is released to PyPI, swap
> the dependency in `pyproject.toml` for a normal version pin and the package
> becomes publishable / `pip`-installable.

```bash
# Verify
risearch-pipeline --help
```

---

## Running the Pipeline

### Single-command mode (CLI)

```bash
# Pre-computed predictions file
risearch-pipeline off-targets \
  -r predictions.tsv \
  -t annotation.gtf \
  -a accessibility_profiles/ \
  --on-target-ids on_target_map.tsv \
  -o results.tsv

# Via YAML config (paths relative to config file)
risearch-pipeline -c config/off-targets.example.yaml
```

### Orchestrated multi-step mode (local)

> `run_pipeline.py` and `convert_risearch_to_parquet.py` live at the repo root —
> they are available when you **clone** the repo, but are **not** installed by
> `pip`/`uv` as console scripts. Run them with `python <script>.py` from a clone.

`run_pipeline.py` runs all pipeline stages in dependency order:

| Step | Command | Notes |
|------|---------|-------|
| `index` | `risearch-pipeline index` | Optional — build RIsearch index |
| `convert` | `convert_risearch_to_parquet.py` | Optional — `.out.gz` → Parquet |
| `accessibility` | `risearch-pipeline accessibility` | Compute RNA accessibility profiles |
| `off-targets` | `risearch-pipeline off-targets` | Main analysis |

`index`, `convert`, and `accessibility` are independent and run in parallel on Slurm. `off-targets` waits for both `accessibility` and `convert`.

```bash
# Dry-run (print commands without executing)
python run_pipeline.py --config config/run-pipeline.example.yaml --dry-run

# Run locally (sequential, logs to logs/<timestamp>/)
python run_pipeline.py --config config/run-pipeline.example.yaml

# Run a subset of steps
python run_pipeline.py --config config/run-pipeline.example.yaml --steps accessibility,off-targets
```

### Cluster mode (Slurm)

```bash
# Dry-run — shows full sbatch commands with dependency chain
python run_pipeline.py --config config/run-pipeline.example.yaml --slurm --dry-run

# Submit to Slurm
python run_pipeline.py --config config/run-pipeline.example.yaml --slurm

# Override resources for all steps
python run_pipeline.py --config config/run-pipeline.example.yaml --slurm \
  --partition gpu --mem 128G --cpus-per-task 32 --account mylab
```

Slurm resources are controlled from the YAML config under the `slurm:` key:

```yaml
slurm:
  partition: batch
  account: mylab
  accessibility:          # per-step overrides
    time: "08:00:00"
    mem: 32G
    cpus_per_task: 8
  off_targets:
    time: "04:00:00"
    mem: 128G
    cpus_per_task: 16
```

Each step gets its own `--output`/`--error` log under `logs/<timestamp>/`. CLI flags (`--partition`, `--time`, `--mem`, `--cpus-per-task`, `--account`) override YAML for all steps.

### Orchestrator config format

See `config/run-pipeline.example.yaml` for the full reference. All paths resolve relative to the config file's directory.

```yaml
steps: [accessibility, off-targets]   # default; add index/convert to enable

slurm:
  partition: batch
  account: mylab

accessibility:
  fasta: ../data/genome.fa
  output: ../data/accessibility/
  workers: 8

off_targets:
  risearch_file: ../data/parquet/
  transcriptome: ../data/annotations.gtf
  accessibility_dir: ../data/accessibility/
  output: ../results/off_targets.tsv
  type: gw
  alpha: "0.8;1.0"
  gamma: "1.0"
```

---

## Thermodynamic Model

Off-target probability is computed from a Boltzmann partition function over all predicted binding sites for each siRNA:

```
W_i  = Expression_i × exp(−ΔG_total_i / RT)

       ΔG_total = ΔG_hybridization + ΔG_opening

Z_s  = Σ W_i  (all off-targets of siRNA s)  +  W_on-target

P(off-target_i | siRNA_s) = W_i / Z_s
```

- **ΔG_hybridization**: RNA-RNA interaction energy from RIsearch (Turner 2004 parameters).
- **ΔG_opening**: Accessibility penalty — cost to unfold the target region, retrieved from pre-computed `RNA.pfl_fold_up` profiles.
- **Expression weighting**: annotation-derived RPKM/TPM values scale each site's contribution.
- **Per-siRNA normalization**: Partition functions are computed independently per siRNA; mixing them is biologically incorrect.

### Parameter sweeps

`--alpha`/`--gamma` clamp extremely favorable energies; `--theta` scales energy differences around a −10 kcal/mol reference. All combinations are computed in a single Polars `group_by` pass, producing separate `P_off_target:alpha=X,gamma=Y` columns.

---

## Output

### TSV (default)

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome or transcript ID |
| `start`, `end` | Binding site coordinates |
| `strand` | Strand orientation |
| `energy` | Hybridization energy (kcal/mol) |
| `opening_energy` | Accessibility penalty (kcal/mol) |
| `dG_total` | Combined free energy |
| `exp_value` | Expression level |
| `P_off_target` | Off-target probability (baseline α=γ=θ=1) |
| `P_off_target:alpha=X,gamma=Y` | Per-parameter-set probabilities |

### Legacy `.results` (optional)

```
# On-target info for siRNA #
# For alpha=1.0 and gamma=1.0; Pon: 0.847; Poff: 0.153; ...
## End of on-target info ##
```

---

## Performance

| Dataset | Time | Memory |
|---------|------|--------|
| 1 k predictions | ~0.2 s | ~80 MB |
| 10 k predictions | ~0.8 s | ~145 MB |
| 100 k predictions | ~6.5 s | ~200 MB |

Key optimizations:
- **Polars** throughout — Rust-backed columnar operations, lazy query planning.
- **Parquet accessibility profiles** — memory-mapped columnar lookup, no full-file reads.
- **ProcessPoolExecutor with Arrow IPC** — per-siRNA files processed in parallel; workers share transcriptome pages via OS page cache.
- **Single `group_by` pass** — all parameter-sweep columns computed in one aggregation.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `polars` | High-performance DataFrames |
| `pyarrow` | Parquet I/O and Arrow IPC |
| `viennaRNA 2.7.2` | RNA folding (`RNA.pfl_fold_up`) |
| `numpy` | Memory-mapped array operations |
| `biopython` | FASTA parsing |
| `typer` | CLI framework |
| `rich` | Progress bars and terminal output |
| `loguru` | Structured logging |
| `omegaconf` | YAML config loading |
| `ncls` | Interval tree for genomic intersection |

---

## Related: Rust RIsearch Core

`risearch` is a separate Rust project providing the RIsearch core and its PyO3
Python bindings, installed as the `risearch` dependency (see
[Installation](#installation)). The pipeline calls the bindings **in-process** —
no subprocess, no intermediate TSV. Features:

- Suffix-array based seed-and-extend search
- Turner 2004 thermodynamic parameters
- Multi-threaded parallel search via Rayon
- SIMD-optimized alignment kernels

---

## License

MIT
