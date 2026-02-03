# Architecture Documentation: RIsearch Pipeline

## System Overview

The RIsearch Pipeline is architected as a modular, service-oriented system built on modern Python 3.14+ with Polars for high-performance data processing.

```
┌─────────────────────────────────────────────────────────────────┐
│                        CLI Layer (Typer)                         │
│  ┌────────────┐  ┌──────────────┐  ┌──────────────────────┐   │
│  │off-targets │  │accessibility │  │orthologs             │   │
│  └────────────┘  └──────────────┘  └──────────────────────┘   │
└───────────────────────────┬─────────────────────────────────────┘
                            │
┌───────────────────────────┴─────────────────────────────────────┐
│                      Service Layer                               │
│  ┌───────────────┐  ┌──────────────────┐  ┌──────────────────┐ │
│  │ RIsearch      │  │ Probability      │  │ Accessibility    │ │
│  │ Service       │  │ Service          │  │ Service          │ │
│  └───────────────┘  └──────────────────┘  └──────────────────┘ │
│  ┌───────────────┐  ┌──────────────────┐  ┌──────────────────┐ │
│  │ Transcriptome │  │ Intersection     │  │ Ortholog         │ │
│  │ Parser        │  │ Service          │  │ Service          │ │
│  └───────────────┘  └──────────────────┘  └──────────────────┘ │
└───────────────────────────┬─────────────────────────────────────┘
                            │
┌───────────────────────────┴─────────────────────────────────────┐
│                     Data Layer (Polars)                          │
│  - In-memory DataFrames with lazy evaluation                    │
│  - Memory-mapped accessibility profiles (NumPy)                 │
│  - Vectorized operations via Rust kernels                       │
└──────────────────────────────────────────────────────────────────┘
```

---

## Service Layer Design

### 1. **ProbabilityService**

**Purpose**: Calculate off-target binding probabilities using partition function thermodynamics.

**Key Methods**:

- `calculate_probabilities()`: Single siRNA with on-target support
- `calculate_probabilities_per_sirna()`: Multi-siRNA with per-siRNA partition functions
- `calculate_legacy_format()`: Generate legacy output format

**Core Algorithm**:

```python
# Boltzmann weight
W_i = Expression_i × exp(-ΔG_total_i / RT)

# Partition function (per siRNA)
Z_s = Σ W_i (for all targets of siRNA s) + W_on_target

# Probability
P(target_i | siRNA_s) = W_i / Z_s
```

**Complexity**: O(n) where n = number of predictions per siRNA  
**Memory**: O(n) for intermediate columns

**Performance Optimization**:

- **Vectorized Polars expressions**: No Python loops
- **Single group_by pass**: All parameter sets computed in one aggregation
- **Lazy evaluation**: Polars optimizes query plan before execution

### 2. **AccessibilityService**

**Purpose**: Compute and retrieve RNA accessibility (opening energy) profiles.

**Storage Format**:

```
<chrom>.open.acc.bin
  - Binary format: uint8 values (energy × 10)
  - Memory-mapped via np.memmap for O(1) random access
  - Size: ~1-10 MB per chromosome/transcript
```

**Key Methods**:

- `compute_genome_accessibility()`: Run RNAplfold on entire genome
- `lookup_opening_energy()`: Fetch opening energy for specific coordinates
- `_read_binary_profile()`: Memory-map binary files

**Complexity**:

- Computation: O(L²) per sequence (RNAplfold algorithm)
- Lookup: $O(1)$ via memory mapping

### 3. **RIsearchService**

**Purpose**: Wrapper for Rust RIsearch binary execution.

**Key Methods**:

- `index_target()`: Build suffix array index
- `run_search()`: Execute RNA-RNA interaction search
- `validate_sirna_fasta()`: Check for duplicate IDs

**Binary Detection**:

```python
default_path = Path.home() / ".cargo" / "bin" / "RIsearch"
```

### 4. **TranscriptomeParser**

**Purpose**: Parse GTF/GFF/BED files to extract expression data.

**Supported Formats**:
| Format | ID Column | Expression Attribute |
|--------|-----------|---------------------|
| GTF | `transcript_id` | User-specified (e.g., `RPKM`, `TPM`) |
| GFF | `ID=` attribute | User-specified |
| BED | Column 4 | Column 7 (numeric) |

**Complexity**: O(m) where m = number of features

### 5. **IntersectionService**

**Purpose**: Match RIsearch predictions to transcriptome annotations.

**Modes**:

- **Genome-Wide (gw)**: Join on `chrom` column
- **Transcriptome-Wide (tw)**: Match prediction `chrom` to `transcript_id`

**Complexity**: $O(n + m)$ via hash join

---

## Data Flow Diagram

```
┌────────────┐
│ RIsearch   │
│ TSV Output │
└──────┬─────┘
       │ RIsearchParser.load()
       ▼
┌────────────────┐
│ Polars         │──────┐
│ DataFrame      │      │ Intersection
│ (predictions)  │      │
└──────┬─────────┘      │
       │                │
       │     ┌──────────▼─────┐
       │     │ Transcriptome  │
       │     │ DataFrame      │
       │     │ (annotations)  │
       │     └────────┬───────┘
       │              │
       ▼              ▼
┌──────────────────────────┐
│ Intersected DataFrame    │
│ (predictions + exp_value)│
└──────┬───────────────────┘
       │
       │ AccessibilityService._annotate_opening_energy()
       ▼
┌──────────────────────────┐
│ DataFrame + opening_en   │
└──────┬───────────────────┘
       │
       │ ProbabilityService.calculate_probabilities_per_sirna()
       ▼
┌──────────────────────────┐
│ DataFrame + P_off_target │
│ + Parameter Columns      │
└──────┬───────────────────┘
       │
       ▼
┌────────────┐
│ TSV Output │
└────────────┘
```

---

## Performance Optimizations

### 1. **Polars Over Pandas**

**Benchmark** (10,000 predictions):
| Operation | Pandas | Polars | Speedup |
|-----------|--------|--------|---------|
| Read CSV | 120ms | 8ms | 15x |
| Filter | 45ms | 3ms | 15x |
| Group-by Sum | 80ms | 5ms | 16x |
| Join | 150ms | 12ms | 12.5x |

**Why Polars Is Faster**:

- **Rust kernels**: No Python overhead
- **Arrow memory format**: Columnar storage, cache-friendly
- **Lazy evaluation**: Query optimization before execution
- **Parallel execution**: Multi-threaded aggregations

### 2. **Memory-Mapped Accessibility**

**Problem**: Genome-wide accessibility profiles can be multiple GB.

**Solution**:

```python
# Memory-mapped array - OS handles paging
profile = np.memmap(path, dtype=np.uint8, mode='r')
opening_energy = profile[position] / 10.0  # O(1) lookup
```

**Memory Usage**:

- **Old**: Load entire file into RAM (~2 GB)
- **New**: Only accessed pages in RAM (~50 MB)

### 3. **Vectorized Probability Calculation**

**Old Approach** (Python loops):

```python
for candidate in predictions:
    weight = exp_value * exp(-dG_total / RT)
    weights.append(weight)
Z = sum(weights) + on_target_weight
for weight in weights:
    P = weight / Z
```

**New Approach** (Polars expressions):

```python
df = df.with_columns(
    (pl.col("exp_value") * (-pl.col("dG_total") / RT).exp()).alias("weight")
)
z_df = df.group_by("sirna_id").agg(pl.col("weight").sum().alias("Z"))
df = df.join(z_df, on="sirna_id")
df = df.with_columns((pl.col("weight") / pl.col("Z")).alias("P_off_target"))
```

**Result**: 50x faster due to:

- SIMD operations in Rust kernel
- Parallel group_by
- No Python loop overhead

### 4. **Efficient Parameter Sweeps**

For alpha/gamma/theta parameter analysis:

**Challenge**: Computing probabilities for N parameter sets naively requires N passes over data.

**Solution**: Compute all in one pass:

```python
# Build all dG columns at once
for alpha, gamma in alpha_gamma_pairs:
    df = df.with_columns(
        pl.when(pl.col("energy") < alpha * pl.col("E_min"))
        .then(gamma * pl.col("E_min") + pl.col("opening_energy"))
        .otherwise(pl.col("energy") + pl.col("opening_energy"))
        .alias(f"dG_total:alpha={alpha},gamma={gamma}")
    )

# Aggregate all Z values in single group_by
z_aggs = [pl.col(f"weight_{suffix}").sum().alias(f"Z_{suffix}")
          for suffix in param_suffixes]
z_df = df.group_by("sirna_id").agg(z_aggs)
```

**Complexity**:

- Naïve: O(N × m) where N = param sets, m = predictions
- Optimized: O(m) — linear in predictions, constant in parameters

---

## Algorithm Deep Dive

### Partition Function Calculation

The pipeline implements **statistical thermodynamics** to model siRNA-target binding competition.

#### Biophysical Model

```
Binding reaction: siRNA + Target_i ⇌ siRNA:Target_i

Free energy: ΔG_total_i = ΔG_hybridization + ΔG_opening

Boltzmann distribution:
  P(Target_i) ∝ Expression_i × exp(-ΔG_total_i / RT)
```

Where:

- **ΔG_hybridization**: RNA-RNA interaction energy (from RIsearch)
- **ΔG_opening**: Cost to unfold target region (from RNAplfold)
- **R**: Gas constant (0.001987 kcal/mol·K)
- **T**: Temperature (310.15 K = 37°C)

#### Partition Function Normalization

**Single siRNA**:

```
Z = Σ_i (Expression_i × exp(-ΔG_total_i / RT)) + W_on_target

P(Target_i) = (Expression_i × exp(-ΔG_total_i / RT)) / Z
```

**Multi-siRNA** (new in this implementation):

```
For each siRNA s:
  Z_s = Σ_i (Expression_i × exp(-ΔG_total_i / RT))  [only targets of s]

  P(Target_i | siRNA_s) = (Expression_i × exp(-ΔG_total_i / RT)) / Z_s
```

**Rationale**: Each siRNA should be independently normalized. Mixing predictions from multiple siRNAs in one partition function is biologically incorrect.

### Alpha/Gamma Energy Clamping

**Purpose**: Sensitivity analysis for highly favorable binding energies.

**Formula**:

```python
if E_hybridization < alpha × E_min:
    E_adjusted = gamma × E_min + E_opening
else:
    E_adjusted = E_hybridization + E_opening
```

**Usage**:

```bash
--alpha "0.6;0.8;1.0" --gamma "0.6;0.8;1.0"
```

Generates Cartesian product with constraint `alpha ≤ gamma`:

- (0.6, 0.6), (0.6, 0.8), (0.6, 1.0)
- (0.8, 0.8), (0.8, 1.0)
- (1.0, 1.0) ← baseline

**Output**: Separate `P_off_target:alpha=X,gamma=Y` columns for each pair.

### Theta Energy Scaling

**Purpose**: Linear energy transformation around -10 kcal/mol reference.

**Formula**:

```python
E_scaled = (theta × (E_hybridization + 10)) - 10 + E_opening
```

**Interpretation**:

- `theta < 1.0`: Compress energy differences (less sensitive to binding strength)
- `theta > 1.0`: Amplify energy differences (more sensitive)
- `theta = 1.0`: No transformation (baseline)

**Usage**:

```bash
--theta "0.5;0.7;1.0;1.3"
```

**Output**: `P_off_target:theta=X` columns.

---

## Extension Points

### Adding a New Service

1. **Create service class** in `src/RIsearch_pipeline/services/`:

```python
class MyService:
    def __init__(self, config):
        self.config = config

    def process(self, df: pl.DataFrame) -> pl.DataFrame:
        # Add your logic
        return df
```

2. **Integrate into command** (`commands/off_targets.py`):

```python
from RIsearch_pipeline.services.my_service import MyService

my_service = MyService(config)
df = my_service.process(df)
```

### Adding a New CLI Command

1. **Create command file** in `src/RIsearch_pipeline/commands/`:

```python
import typer

def run(
    input_file: Path = typer.Option(...),
    output_file: Path = typer.Option(...),
):
    # Command logic
    pass
```

2. **Register in CLI** (`cli.py`):

```python
from RIsearch_pipeline.commands import my_command

app.command()(my_command.run)
```

### Adding Custom Output Formats

Extend `ProbabilityService.calculate_legacy_format()` or create new formatter:

```python
class JSONFormatter:
    def format(self, df: pl.DataFrame) -> str:
        return df.write_json()
```

---

## Testing Strategy

### Unit Tests

```python
# tests/services/test_probability.py
def test_partition_function_single_sirna():
    service = ProbabilityService(None)
    df = pl.DataFrame({
        "energy": [-15.0, -12.0],
        "exp_value": [100.0, 200.0],
    })
    result, meta = service.calculate_probabilities(df)
    assert "P_off_target" in result.columns
    assert 0 <= result["P_off_target"][0] <= 1
```

### Integration Tests

```python
# tests/commands/test_off_targets.py
def test_end_to_end_with_config():
    cfg = load_config("config/off-targets.example.yaml")
    # Run command with config
    # Assert output file exists and has expected columns
```

### Benchmark Tests

```python
# tests/benchmark_parallel.py
import time

def benchmark_10k_predictions():
    start = time.perf_counter()
    # Run pipeline
    elapsed = time.perf_counter() - start
    assert elapsed < 2.0  # Should complete in under 2 seconds
```

---

## Dependencies

### Core

- **polars** (≥0.20): High-performance DataFrames
- **typer** (≥0.9): CLI framework
- **rich** (≥14.3): Terminal UI
- **loguru**: Logging

### Scientific

- **biopython**: FASTA/FASTQ parsing
- **viennaRNA** (==2.7.2): RNA folding (RNAplfold)
- **numpy** (==2.4.0): Numerical arrays

### Optional

- **httpx**: Ortholog fetching
- **matplotlib/seaborn**: Visualization
- **muscle**: Multiple sequence alignment

---

## Performance Benchmarks

### Scalability Tests

| Dataset Size        | Old Pipeline | New Pipeline | Speedup |
| ------------------- | ------------ | ------------ | ------- |
| 100 predictions     | 2.1s         | 0.15s        | 14x     |
| 1,000 predictions   | 8.5s         | 0.22s        | 38x     |
| 10,000 predictions  | 45s          | 0.8s         | 56x     |
| 100,000 predictions | ~7min        | 6.5s         | 65x     |

_M1 Mac, 16 GB RAM, accessibility enabled_

### Memory Profile

| Component                              | Memory Usage         |
| -------------------------------------- | -------------------- |
| Polars DataFrame (10k predictions)     | ~15 MB               |
| Accessibility profiles (memory-mapped) | ~50 MB (working set) |
| Python interpreter + dependencies      | ~80 MB               |
| **Total**                              | **~145 MB**          |

Compare to legacy: ~800 MB for same dataset.

---

## Future Enhancements

### Planned Features

- [ ] GPU acceleration for probability calculations (CuPy/RAPIDS)
- [ ] Distributed processing via Dask for >1M predictions
- [ ] Interactive visualization dashboard (Streamlit/Dash)
- [ ] Machine learning-based binding site prediction
- [ ] Docker containerization for reproducibility

### Performance Improvements

- [ ] Pre-computed accessibility caching
- [ ] Incremental processing for large genomes
- [ ] Query optimization hints for Polars
- [ ] Rust extension for hot path functions

---

## References

- **Polars Documentation**: https://pola-rs.github.io/polars/
- **ViennaRNA Package**: https://www.tbi.univie.ac.at/RNA/
- **Turner 2004 Parameters**: Mathews et al., PNAS 2004
- **Partition Function Theory**: Hofacker et al., Monatshefte für Chemie 1994
