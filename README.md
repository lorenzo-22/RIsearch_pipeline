# RIOT

RIOT — siRNA off-target discovery pipeline.

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

## Python API

RIOT can be used as a library — `import riot`, then call the commands as plain functions; they return their results in-memory.

```python
import riot

# Off-target analysis on a pre-computed predictions file → polars.DataFrame
df = riot.off_targets(risearch_file="predictions.tsv", gtf_file="annotations.gtf")

# Pre-compute per-chromosome accessibility profiles → dict[chrom -> Path]
profiles = riot.accessibility(genome="genome.fa", output="acc_dir/")

# Build a RIsearch index → Path (needs the external 'risearch' package)
idx = riot.index("target.fa")

# Run a RIsearch search → polars.DataFrame (needs the external 'risearch' package)
hits = riot.search("query.fa", "target.fa.idx", target="target.fa")
```

**Notes**

- `riot.off_targets` returns a `polars.DataFrame` when given a single predictions file. With a **directory** of per-siRNA Parquet files it streams results to disk and returns a summary `dict` (output paths + row/summary counts) instead.
- On error the underlying command raises `typer.Exit` (a Click exception).
- `riot.index` / `riot.search` require the external `risearch` package — the same dependency the CLI's `index`/`search` commands need.

---

## Installation

Requires **Python ≥ 3.14** and **ViennaRNA 2.7.2**.

```bash
git clone git@github.com:lorenzo-22/RIOT.git
cd RIOT

# Create virtual environment and install dependencies
uv venv && source .venv/bin/activate
uv sync
```

### Publishing / PyPI

The PyPI distribution name is **`riot-sirna`** (plain `riot` is already taken on
PyPI). Installing from PyPI with `pip install riot-sirna` is **pending
publication of the upstream `risearch` dependency** (currently a `git+ssh`
direct URL — see below) — until then, install from source / git as shown above.

> **Note:** the import name `riot` could collide with Datadog's PyPI `riot`
> package if both are installed in the same environment.

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
riot --help
```

---

## Running the Pipeline

### Single-command mode (CLI)

```bash
# Pre-computed predictions file
riot off-targets \
  -r predictions.tsv \
  -t annotation.gtf \
  -a accessibility_profiles/ \
  --on-target-ids on_target_map.tsv \
  -o results.tsv

# Via YAML config (paths relative to config file)
riot -c example_yaml/off-targets.example.yaml
```

### Orchestrated multi-step mode (local)

> `scripts/run_pipeline.py` and `scripts/convert_risearch_to_parquet.py` live in
> the `scripts/` directory — they are available when you **clone** the repo, but
> are **not** installed by `pip`/`uv` as console scripts. Run them with
> `python scripts/<script>.py` from a clone.

`scripts/run_pipeline.py` runs all pipeline stages in dependency order:

| Step | Command | Notes |
|------|---------|-------|
| `index` | `riot index` | Optional — build RIsearch index |
| `convert` | `scripts/convert_risearch_to_parquet.py` | Optional — `.out.gz` → Parquet |
| `accessibility` | `riot accessibility` | Compute RNA accessibility profiles |
| `off-targets` | `riot off-targets` | Main analysis |

`index`, `convert`, and `accessibility` are independent and run in parallel on Slurm. `off-targets` waits for both `accessibility` and `convert`.

```bash
# Dry-run (print commands without executing)
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml --dry-run

# Run locally (sequential, logs to logs/<timestamp>/)
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml

# Run a subset of steps
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml --steps accessibility,off-targets
```

### Cluster mode (Slurm)

```bash
# Dry-run — shows full sbatch commands with dependency chain
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml --slurm --dry-run

# Submit to Slurm
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml --slurm

# Override resources for all steps
python scripts/run_pipeline.py --config example_yaml/run-pipeline.example.yaml --slurm \
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

### Multiple transcriptomes (fan-out)

Add a top-level `transcriptomes:` list to analyze several genomes/transcriptomes in one launch — Slurm-native, **one transcriptome per node**, run in parallel. Each entry is an independent run (its own predictions, annotation, and output); groups never mix, so the off-target probability math (`Z_s`) is unchanged. The top-level `off_targets:`/`accessibility:`/`convert:` blocks act as shared defaults; each group overrides its per-group fields. `index` is never fanned out.

```yaml
steps: [off-targets]

off_targets:            # shared defaults for every group
  alpha: "0.8;1.0"
  type: gw

transcriptomes:
  - name: human         # required — job names, logs, default output
    risearch_file: ../data/human.out
    transcriptome: ../data/human.gtf
    accessibility_dir: ../data/human_acc/   # precomputed
    output: results/human.tsv               # must be unique per group
  - name: mouse
    risearch_file: ../data/mouse.out
    transcriptome: ../data/mouse.gtf
    accessibility_dir: ../data/mouse_acc/
    output: results/mouse.tsv
```

Each group submits its own Slurm job(s) (`rip_off_targets_human`, `rip_off_targets_mouse`, …); a group's `off-targets` waits only on its own upstream jobs. To compute accessibility per group, add `accessibility` to `steps` and give each group **both** a `fasta:` and an `accessibility_dir:` (the profiles are written there and read back by that group's off-targets; both are required and validated). Omitting a group `output` defaults it to `results/<name>.tsv`.

> `convert` also fans out per group when in `steps`, but its output is **not** auto-wired to that group's `risearch_file` — set each group's `input_dir`/`out_dir` and `risearch_file` explicitly.

### Orchestrator config format

See `example_yaml/run-pipeline.example.yaml` for the full reference. All paths resolve relative to the config file's directory.

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
