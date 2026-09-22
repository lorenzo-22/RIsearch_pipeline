# siOFF

siOFF — siRNA off-target discovery pipeline.

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

### Global options

Available on `sioff` itself, before any subcommand.

| Flag | Description |
|------|-------------|
| `-c / --config` | Path to a YAML config file; runs the command named in it. Top-level only — `sioff -c cfg.yaml`, not `sioff off-targets -c cfg.yaml` |
| `-v / --verbose` | Enable DEBUG-level logging |
| `--version` | Print the installed version and exit (long form only — `-v` is `--verbose`) |

### `off-targets`

| Flag | Description |
|------|-------------|
| `-r / --risearch-file` | Pre-computed predictions (TSV, `.out.gz`, or directory of Parquet files — directory triggers parallel per-siRNA mode) |
| `-s / --sirna-fasta` | siRNA FASTA — runs RIsearch in-process via PyO3 bindings |
| `--target-fasta / --genome` | Target FASTA for in-process RIsearch |
| `-t / --transcriptome` | GTF, GFF3, BED6 or BED7 annotation file |
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

siOFF can be used as a library — `import sioff`, then call the commands as plain functions; they return their results in-memory.

```python
import sioff

# Off-target analysis on a pre-computed predictions file → polars.DataFrame
df = sioff.off_targets(risearch_file="predictions.tsv", gtf_file="annotations.gtf")

# Pre-compute per-chromosome accessibility profiles → dict[chrom -> polars.DataFrame]
profiles = sioff.accessibility(genome="genome.fa")

# Build a RIsearch index → Path (needs the external 'risearch' package)
idx = sioff.index("target.fa")

# Run a RIsearch search → polars.DataFrame (needs the external 'risearch' package)
hits = sioff.search("query.fa", "target.fa.idx", target="target.fa")
```

**Notes**

- `sioff.off_targets` returns a `polars.DataFrame` when given a single predictions file. With a **directory** of per-siRNA Parquet files it returns a **generator** yielding one `polars.DataFrame` per siRNA; iterate it to consume the results. Neither form writes files — the API layer returns results in memory, and writing output is the CLI's job.
- `sioff.accessibility` likewise writes nothing: it returns `dict[chrom -> polars.DataFrame]`. Use the `accessibility` CLI command (or `sioff -c <config>`) if you want `{chrom}.accessibility.parquet` files on disk.
- On bad input the API functions raise ordinary Python exceptions — `FileNotFoundError` or `ValueError` — not `typer.Exit`. `typer.Exit` is raised only by the CLI layer for its own argument validation.
- `sioff.index` / `sioff.search` require the external `risearch` package — the same dependency the CLI's `index`/`search` commands need.

---

## Installation

Requires **Python ≥ 3.11** (tested on 3.11–3.14) and **ViennaRNA 2.7.2**.

```bash
git clone git@github.com:lorenzo-22/siOFF.git
cd siOFF

# Create virtual environment and install dependencies
uv venv && source .venv/bin/activate
uv sync    # full pipeline, including the in-process index/search commands
```

Plain `uv sync` installs the complete pipeline: the `risearch` dependency group
is part of `tool.uv.default-groups`, so the in-process `index` / `search`
commands work out of the box. Installing it requires SSH access to the private
`risearch` repository (see [The `risearch` dependency](#the-risearch-dependency));
without that access, `uv sync --no-group risearch` installs the core
`off-targets` / `accessibility` pipeline, which runs on pre-computed predictions.

### Publishing / PyPI

The PyPI distribution name and the import name are both **`sioff`**
(`pip install sioff`, `import sioff`).

`risearch` is declared as a [PEP 735](https://peps.python.org/pep-0735/)
dependency group rather than a normal dependency only because it is not yet on
PyPI: a `git+ssh` direct reference inside `[project.dependencies]` is copied
verbatim into the published `Requires-Dist` metadata, and PyPI rejects
distributions carrying a direct URL. A dependency group is resolver-only and
never reaches that metadata, so the package stays publishable while `risearch`
remains a private git repository. Once `risearch` is published on PyPI, it moves
into `[project.dependencies]` as a normal version pin so that
`pip install sioff` brings in the full pipeline, in-process `index` / `search`
included.

Until then, users installing `sioff` from PyPI get the `off-targets` /
`accessibility` pipeline; the in-process `index` / `search` commands additionally
need `risearch` installed from git (see below).

### The `risearch` dependency

The `risearch` PyO3 bindings power the **in-process `index` and `search`**
commands (computing RNA-RNA interaction predictions in-process) and are part of
the default install. The core off-target analysis — `off-targets` and
`accessibility` running on **pre-computed** RIsearch output (TSV / `.out.gz` /
Parquet) — works **without** `risearch` installed; it is imported lazily.

`risearch` is currently fetched from a **private** repository over SSH and is
**not on PyPI**, so the default `uv sync` requires SSH access to that repo
(use `uv sync --no-group risearch` without it):

```
git+ssh://git@github.com/saiden89/risearch.git@1a03a47…#subdirectory=bindings/python
```

The commit is pinned rather than tracking a branch so that installs are
reproducible: `risearch` is a fast-moving repository whose public API and search
results have both changed across releases, so the pin is bumped deliberately and
verified, not automatically.

> **PyPI:** this direct URL no longer blocks publication — `risearch` is declared
> as a PEP 735 dependency group, which never reaches published metadata, so
> `sioff` is publishable with its core dependencies. `risearch` itself still
> has to be installed from git; if it is ever released to PyPI, swap this for a
> normal version pin.

```bash
# Verify
sioff --help
```

---

## Running the Pipeline

### Single-command mode (CLI)

```bash
# Pre-computed predictions file
sioff off-targets \
  -r predictions.tsv \
  -t annotation.gtf \
  -a accessibility_profiles/ \
  --on-target-ids on_target_map.tsv \
  -o results.tsv

# Via YAML config (paths relative to config file)
sioff -c example_yaml/off-targets.example.yaml
```

### Orchestrated multi-step mode (local)

> `scripts/run_pipeline.py` lives in the `scripts/` directory — it is available
> when you **clone** the repo, but is **not** installed by `pip`/`uv` as a console
> script. Run it with `python scripts/run_pipeline.py` from a clone.

`scripts/run_pipeline.py` runs all pipeline stages in dependency order:

| Step | Command | Notes |
|------|---------|-------|
| `index` | `sioff index` | Optional — build RIsearch index |
| `accessibility` | `sioff accessibility` | Compute RNA accessibility profiles |
| `off-targets` | `sioff off-targets` | Main analysis |

`index` and `accessibility` are independent and run in parallel on Slurm. `off-targets` waits for `accessibility`.

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

Add a top-level `transcriptomes:` list to analyze several genomes/transcriptomes in one launch — Slurm-native, **one transcriptome per node**, run in parallel. Each entry is an independent run (its own predictions, annotation, and output); groups never mix, so the off-target probability math (`Z_s`) is unchanged. The top-level `off_targets:`/`accessibility:` blocks act as shared defaults; each group overrides its per-group fields. `index` is never fanned out.

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

### Orchestrator config format

See `example_yaml/run-pipeline.example.yaml` for the full reference. All paths resolve relative to the config file's directory.

```yaml
steps: [accessibility, off-targets]   # default; add index to enable

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

- **ΔG_hybridization**: duplex interaction energy from RIsearch. The nearest-neighbour parameter set is selectable with `sioff search -z/--matrix` (default Turner 2004).
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
- Selectable nearest-neighbour parameter sets: Turner 2004 (default) and 1999 for RNA-RNA, SantaLucia-Hicks 2004 for DNA-DNA, and Sugimoto 1995 for RNA/DNA hybrids
- Multi-threaded parallel search via Rayon
- SIMD-optimized alignment kernels

---

## Citation

If you use siOFF in published work, please cite:

> Roncelli S, Favaro L, Anthon C, Gorodkin J. *RIsearch and siOFF: An integrated,
> high-performance framework for RNA-RNA interaction and siRNA off-target
> prediction.* Bioinformatics.

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

Note that citation is not merely requested but a **term of the BUSL-1.1
licence** for any production use of siOFF or `risearch` — see [License](#license).

---

## License

siOFF (`sioff`) is released under the **Business Source License 1.1
(BUSL-1.1)** — see [LICENSE](LICENSE) — the same licence as its `risearch`
dependency. BUSL-1.1 is *source-available, not open source* and is not
OSI-approved. In brief:

- **Research and production use are free**, including commercial use, with one
  exception below.
- **No hosted services.** You may not use it to provide SaaS, PaaS or any other
  hosted or cloud-based service to third parties — commercial *or*
  non-commercial — except where such services are provided exclusively to
  accredited academic institutions, or individuals affiliated with them, for
  non-commercial research or educational purposes.
- **Citation is a licence term.** Any production use must cite the work; see
  [Citation](#citation).
- **Converts to Apache-2.0 on 14 September 2030.**
- Licensor: RTH, University of Copenhagen. Commercial licensing enquiries:
  <software@rth.dk>.

The summary above is provided for orientation only; [LICENSE](LICENSE) is the
authoritative text.
