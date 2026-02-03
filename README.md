# RIsearch Pipeline

A specialized bioinformatics pipeline for **siRNA off-target discovery and analysis**. This tool integrates RNA-RNA interaction predictions with transcriptome annotations, accessibility profiling, and thermodynamic modeling to quantify off-target binding probabilities.

---

## 🎯 What This Pipeline Does

siRNA (small interfering RNA) therapeutics are powerful tools for gene silencing, but they can inadvertently bind to unintended mRNA targets, causing **off-target effects**. This pipeline helps researchers:

1. **Identify** potential off-target binding sites across the transcriptome
2. **Quantify** the probability of off-target hybridization using thermodynamic models
3. **Prioritize** which off-targets are biologically relevant based on expression levels and RNA accessibility
4. **Compare** across orthologs in different species for conservation analysis

---

## 🧬 Commands

### `off-targets` — Core Analysis Command

Analyzes RIsearch2 predictions to calculate off-target probabilities.

```bash
risearch-pipeline off-targets \
  -r predictions.tsv \           # RIsearch2 output
  -t transcriptome.gtf \         # Annotated transcriptome
  -a accessibility_profiles/ \   # Pre-computed accessibility
  --on-target on_target.fa \     # On-target sequence
  -q sirna.fa \                  # siRNA query
  -o results.tsv
```

**Key Features:**

- Partition function-based probability calculation
- Integration of hybridization energy + RNA accessibility (opening energy)
- Expression-weighted Boltzmann statistics
- Legacy format output compatible with existing pipelines

### `accessibility` — Pre-compute Accessibility Profiles

Generates RNA accessibility profiles using ViennaRNA's `RNAplfold`.

```bash
risearch-pipeline accessibility \
  -i genome.fa \
  -o profiles/ \
  --window 80 \
  --span 40 \
  --unpaired 30
```

### `orthologs` — Cross-species Ortholog Analysis

Downloads orthologs from OrthoDB and transcriptomes from NCBI for comparative analysis.

```bash
risearch-pipeline orthologs \
  --target-gene PSMC2 \
  --species-list species.txt \
  --email user@example.com \
  --plot-lengths \
  --run-msa
```

---

## 📦 Dependencies

| Package                | Purpose                                                           |
| ---------------------- | ----------------------------------------------------------------- |
| **polars**             | High-performance DataFrames for processing large prediction files |
| **biopython**          | Sequence I/O (FASTA/FASTQ parsing)                                |
| **viennaRNA**          | RNA secondary structure and accessibility (`RNAplfold`)           |
| **numpy**              | Numerical operations and memory-mapped array storage              |
| **typer**              | CLI framework with automatic help generation                      |
| **rich**               | Beautiful terminal output (progress bars, tables)                 |
| **loguru**             | Structured logging                                                |
| **httpx**              | Async HTTP client for OrthoDB/NCBI APIs                           |
| **matplotlib/seaborn** | Visualization (ortholog length plots)                             |
| **muscle**             | Multiple sequence alignment                                       |
| **omegaconf**          | YAML configuration file support                                   |

---

## 🔬 Why This Is Interesting

### 1. **Thermodynamic-Based Probability Model**

Unlike simple sequence matching, this pipeline uses **free energy (ΔG)** calculations to predict binding affinity:

```
P(off-target) = W_off / Z_total

where:
  W = Expression × exp(-ΔG_total / RT)
  ΔG_total = ΔG_hybridization + ΔG_opening
  Z_total = Σ(all targets) + W_on-target
```

This partition function approach properly accounts for competition between all potential binding sites.

### 2. **RNA Accessibility Integration**

Not all binding sites are equally accessible—some are buried in secondary structures. The pipeline:

- Pre-computes accessibility profiles using `RNAplfold`
- Stores profiles as memory-mapped binary arrays for fast random access
- Looks up **opening energy** (cost to unfold the target region) for each prediction

### 3. **Expression-Weighted Analysis**

Highly expressed transcripts contribute more to off-target effects. The pipeline:

- Parses expression values from GTF attributes (RPKM, TPM, etc.)
- Weights Boltzmann probabilities by expression level
- Properly normalizes across the entire transcriptome

### 4. **Genome-Wide & Transcriptome-Wide Modes**

Supports both:

- **GW (Genome-Wide)**: Predictions mapped to genomic coordinates, intersected with annotations
- **TW (Transcriptome-Wide)**: Predictions directly on transcript sequences

---

## ⚡ Performance Optimizations

### Memory-Mapped Accessibility Profiles

```python
# Binary format: uint8 values (energy × 10) for compact storage
# Memory-mapped for constant-memory access to multi-GB profiles
profile = np.memmap(path, dtype=np.uint8, mode='r')
```

Genome-wide accessibility profiles can be gigabytes—memory mapping avoids loading them entirely.

### Polars Over Pandas

```python
# Lazy evaluation and parallel execution
df = pl.scan_csv(file).filter(...).collect()
```

Polars provides:

- Zero-copy operations where possible
- Automatic parallelization
- Lazy query optimization

### Efficient Coordinate Intersection

```python
# Group by chromosome, then binary search or interval trees
intersector.intersect(predictions_df, transcriptome_df, mode="gw")
```

Genomic interval operations are optimized for large-scale analysis.

### Caching & Deduplication

- Unique target coordinates are extracted before accessibility lookups
- Results are joined back, avoiding redundant queries
- Temporary files are cleaned up automatically

---

## 📊 Output Format

### Standard TSV Output

| Column           | Description                     |
| ---------------- | ------------------------------- |
| `chrom`          | Chromosome/Transcript ID        |
| `start`, `end`   | Binding site coordinates        |
| `strand`         | Strand orientation              |
| `energy`         | Hybridization energy (kcal/mol) |
| `opening_energy` | Accessibility penalty           |
| `dG_total`       | Combined free energy            |
| `exp_value`      | Expression level                |
| `P_off_target`   | Off-target probability          |

### Legacy Format (`.results`)

Compatible with older pipeline outputs:

```
# On-target info for siRNA #
# For alpha=1.0 and gamma=1.0; Pon: 0.847; Poff: 0.153; ...
## End of on-target info ##
```

---

## 🛠️ Installation

### Quick Start (Local)

```bash
# 1. Install Rust (for RIsearch binary)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# 2. Clone and build RIsearch
git clone https://github.com/your-org/RIsearch_pipeline.git
cd RIsearch_pipeline/risearch
cargo install --path .
cd ..

# 3. Install Python pipeline
uv venv && source .venv/bin/activate
uv sync

# 4. Verify
risearch-pipeline --help
```

Requires **Python ≥3.14**, **Rust 1.70+**, and **ViennaRNA** (for accessibility).

### Production Deployment

For SSH servers, HPC clusters, or containerized environments, see:

- **[DEPLOYMENT.md](DEPLOYMENT.md)** - Complete deployment guide with troubleshooting
- **Installation script**: Run `bash scripts/install.sh` for automated setup

---

## 📈 Performance Benchmarks

| Dataset Size        | Legacy Pipeline | New Pipeline | Speedup |
| ------------------- | --------------- | ------------ | ------- |
| 100 predictions     | 2.1s            | 0.15s        | **14x** |
| 1,000 predictions   | 8.5s            | 0.22s        | **38x** |
| 10,000 predictions  | 45s             | 0.8s         | **56x** |
| 100,000 predictions | ~7min           | 6.5s         | **65x** |

_Benchmarked on M1 Mac, 16 GB RAM, with accessibility enabled_

**Memory usage**: ~145 MB vs ~800 MB (legacy) for 10k predictions

---

## 📖 Documentation

- **[MIGRATION.md](MIGRATION.md)** - Migrating from legacy `old_pipeline.py`
- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Technical deep dive for developers
- **[DEPLOYMENT.md](DEPLOYMENT.md)** - Production deployment guide
- **[README.md](README.md)** - This file (overview and quick start)

---

## 📚 Related: Rust RIsearch Core

This Python pipeline is designed to work with predictions from [RIsearch](./risearch/), a high-performance RNA-RNA interaction predictor written in Rust. The Rust core provides:

- Suffix-array based seed-and-extend search
- Turner 2004 thermodynamic parameters
- Multi-threaded parallel search via Rayon
- SIMD-optimized alignment kernels (see `docs/SIMD_ARCHITECTURE.md`)

---

## License

MIT
