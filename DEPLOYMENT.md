# Deployment Guide: RIsearch Pipeline

## Overview

This guide covers deploying the RIsearch Pipeline in production environments, with specific focus on SSH server deployment (the most common use case for bioinformatics pipelines).

---

## Prerequisites

### System Requirements

- **OS**: Linux (Ubuntu 20.04+, CentOS 8+) or macOS (11+)
- **CPU**: 4+ cores recommended for parallel processing
- **RAM**: 8 GB minimum, 16 GB recommended
- **Storage**: 50 GB+ for genome-wide accessibility profiles
- **Python**: 3.14+
- **Rust**: 1.70+ (for RIsearch binary)

### Network Requirements

- **Outbound**: Required for `uv` package downloads, OrthoDB/NCBI API access
- **Inbound**: Not required (CLI tool, not a service)

---

## Quick Start (Local Development)

```bash
# 1. Clone repository
git clone https://github.com/your-org/RIsearch_pipeline.git
cd RIsearch_pipeline

# 2. Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# 3. Build RIsearch binary
cd risearch
cargo install --path .  # Installs to ~/.cargo/bin/RIsearch
cd ..

# 4. Create Python environment
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# 5. Install pipeline
uv sync

# 6. Verify installation
risearch-pipeline --help
RIsearch --version
```

---

## SSH Server Deployment

### Step 1: Prepare Environment

```bash
# Connect to server
ssh user@server.example.com

# Install system dependencies
# Ubuntu/Debian
sudo apt update
sudo apt install -y build-essential git curl pkg-config

# CentOS/RHEL
sudo yum groupinstall "Development Tools"
sudo yum install -y git curl pkg-config
```

### Step 2: Install Rust Toolchain

```bash
# Install rustup (Rust installer)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Select option 1 (default installation)
# Add Rust to PATH for current session
source $HOME/.cargo/env

# Verify installation
rustc --version  # Should show 1.70+
cargo --version
```

### Step 3: Clone Repository

```bash
# Clone to desired location
cd ~
git clone https://github.com/your-org/RIsearch_pipeline.git
cd RIsearch_pipeline
```

### Step 4: Build RIsearch Binary

```bash
# Navigate to Rust subproject
cd risearch

# Build and install
cargo install --path .

# Verify binary is accessible
which RIsearch  # Should print: /home/user/.cargo/bin/RIsearch
RIsearch --help
```

### Step 5: Install Python Pipeline

#### Option A: Using `uv` (Recommended)

```bash
# Install uv if not already available
curl -LsSf https://astral.sh/uv/install.sh | sh

# Return to pipeline root
cd ~/RIsearch_pipeline

# Create virtual environment
uv venv .venv

# Activate environment
source .venv/bin/activate

# Install pipeline
uv sync
```

#### Option B: Using `pip`

```bash
# Install Python 3.14+ (if needed)
# Ubuntu
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.14 python3.14-venv

# Create venv
python3.14 -m venv .venv
source .venv/bin/activate

# Install pipeline
pip install -e .
```

### Step 6: Verify Installation

```bash
# Test CLI
risearch-pipeline --help

# Run example
risearch-pipeline off-targets \
  --risearch-file tests/off-targets/single/input/risearch_siRNAID.out \
  --transcriptome tests/off-targets/single/input/expression_data.bed \
  --accessibility-dir tests/off-targets/single/input/accessibility \
  --output test_output.tsv \
  --type tw \
  --sense-only
```

---

## Configuration Management

### Using YAML Configs

Create reusable configuration files:

```yaml
# config/production.yaml
command: off-targets

off_targets:
  risearch_file: /data/predictions/run_001.tsv
  transcriptome: /data/annotations/hg38.gtf
  accessibility_dir: /data/accessibility/hg38/
  output: /results/run_001_results.tsv
  type: gw
  sense_only: true
  alpha: "0.8;1.0"
  gamma: "0.8;1.0"
  on_target_expression: 1000.0
  verbose: true
```

Run with:

```bash
risearch-pipeline -c config/production.yaml
```

### Environment Variables

Set defaults via environment:

```bash
export RISEARCH_BINARY=/custom/path/to/RIsearch
export RISEARCH_ACCESSIBILITY_DIR=/data/accessibility/
```

---

## Performance Tuning

### Multi-Threading

Polars automatically uses all CPU cores. To limit:

```bash
# Limit to 8 threads
export POLARS_MAX_THREADS=8

risearch-pipeline off-targets ...
```

### Memory Optimization

For large datasets (>100k predictions):

1. **Streaming mode** (future feature):

```bash
risearch-pipeline off-targets --streaming ...
```

2. **Batch processing**:

```bash
# Split large prediction files
split -l 10000 large_predictions.tsv batch_

# Process in parallel
for batch in batch_*; do
  risearch-pipeline off-targets -r $batch ... &
done
wait
```

3. **Memory-mapped accessibility**:

```bash
# Already default behavior - no action needed
# Profiles are memory-mapped, not fully loaded
```

---

## Monitoring & Logging

### Log Levels

```bash
# Enable verbose logging
risearch-pipeline --verbose off-targets ...

# Capture logs to file
risearch-pipeline off-targets ... 2>pipeline.log
```

### Log Output Format (Loguru)

```
2026-02-03 21:00:00.123 | INFO     | probability:calculate_probabilities_per_sirna:215 - Processing 5 siRNAs
2026-02-03 21:00:00.456 | WARNING  | accessibility:lookup_opening_energy:89 - Profile not found for chr22
```

### Performance Profiling

```bash
# Install Scalene profiler
uv add --dev scalene

# Profile execution
scalene risearch-pipeline off-targets ...

# Or use built-in benchmark
python tests/off-targets/profile_benchmark.py
```

---

## Troubleshooting

### Common Issues

#### 1. "RIsearch binary not found"

**Symptom**:

```
WARNING: RIsearch binary not found at /home/user/.cargo/bin/RIsearch
```

**Solution**:

```bash
# Verify Rust installation
which cargo

# Rebuild RIsearch
cd risearch
cargo install --path .

# Check binary
ls -lh ~/.cargo/bin/RIsearch
```

#### 2. "ViennaRNA not found"

**Symptom**:

```
ModuleNotFoundError: No module named 'RNA'
```

**Solution A** (pip install failed to build):

```bash
# Install ViennaRNA system package first
# Ubuntu
sudo apt install libvienna-rna-dev

# macOS
brew install viennarna

# Then reinstall Python package
pip install --force-reinstall viennaRNA==2.7.2
```

**Solution B** (use conda):

```bash
conda install -c bioconda viennarna
```

#### 3. "Out of memory" for large genomes

**Symptom**:

```
MemoryError: Unable to allocate array
```

**Solution**:

```bash
# 1. Use memory-mapped profiles (already default)
# 2. Reduce worker count
export POLARS_MAX_THREADS=2

# 3. Process chromosomes separately
risearch-pipeline off-targets -r chr1_predictions.tsv ...
risearch-pipeline off-targets -r chr2_predictions.tsv ...
```

#### 4. "Permission denied" for accessibility directory

**Symptom**:

```
PermissionError: [Errno 13] Permission denied: '/data/accessibility/chr1.open.acc.bin'
```

**Solution**:

```bash
# Fix permissions
chmod -R u+rw /data/accessibility/

# Or copy to user-owned location
cp -r /data/accessibility ~/my_accessibility
risearch-pipeline off-targets -a ~/my_accessibility ...
```

---

## Production Best Practices

### 1. **Version Control Configuration**

```bash
# Store configs in Git
git add config/production.yaml
git commit -m "Add production config for hg38"

# Track pipeline version
risearch-pipeline --version > VERSION.txt
```

### 2. **Automated Pipelines**

Example Slurm job script:

```bash
#!/bin/bash
#SBATCH --job-name=risearch_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j.out

source ~/.cargo/env
source ~/RIsearch_pipeline/.venv/bin/activate

risearch-pipeline -c $CONFIG_FILE
```

Submit:

```bash
sbatch --export=CONFIG_FILE=config/run_001.yaml pipeline.slurm
```

### 3. **Data Management**

Directory structure:

```
/data/
├── genomes/
│   ├── hg38.fa
│   └── mm10.fa
├── annotations/
│   ├── hg38.gtf
│   └── mm10.gtf
├── accessibility/
│   ├── hg38/
│   │   ├── chr1.open.acc.bin
│   │   └── ...
│   └── mm10/
│       └── ...
├── risearch_outputs/
│   └── run_001.tsv
└── results/
    └── run_001_results.tsv
```

### 4. **Backup Critical Data**

```bash
# Accessibility profiles take hours to compute - back them up!
rsync -avz /data/accessibility/ backup:/data/accessibility/

# Version control configs
git add config/ && git commit -m "Update config"
```

---

## Scaling to High-Performance Computing

### Batch Job Arrays

Process multiple siRNAs in parallel:

```bash
#!/bin/bash
#SBATCH --array=1-100  # 100 siRNAs

SIRNA_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sirna_list.txt)

risearch-pipeline off-targets \
  -s sirnas/${SIRNA_ID}.fa \
  --target-fasta /data/genomes/hg38.fa \
  -t /data/annotations/hg38.gtf \
  -a /data/accessibility/hg38/ \
  -o results/${SIRNA_ID}_results.tsv
```

### Database Integration

Store results in PostgreSQL:

```python
import polars as pl
import psycopg2

# Read results
df = pl.read_csv("results.tsv", separator="\t")

# Write to database
df.write_database("off_targets", "postgresql://user:pass@host/db")
```

---

## Security Considerations

### 1. **Sensitive Data**

If working with proprietary sequences:

```bash
# Encrypt sensitive files
gpg --encrypt --recipient user@example.com sirna_sequences.fa

# Decrypt for processing
gpg --decrypt sirna_sequences.fa.gpg | risearch-pipeline ...
```

### 2. **Restricted Access**

```bash
# Limit file permissions
chmod 600 config/production.yaml  # Only owner can read/write
chmod 700 ~/RIsearch_pipeline/    # Only owner can access
```

### 3. **Audit Logging**

```bash
# Log all commands
risearch-pipeline off-targets ... | tee -a audit.log
```

---

## Containerization (Future)

### Docker (Planned)

```dockerfile
FROM python:3.14-slim

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install pipeline
WORKDIR /app
COPY . .
RUN cargo install --path risearch/
RUN pip install -e .

ENTRYPOINT ["risearch-pipeline"]
```

Build and run:

```bash
docker build -t risearch-pipeline .
docker run -v /data:/data risearch-pipeline off-targets -c /data/config.yaml
```

---

## Support & Maintenance

### Getting Help

- **GitHub Issues**: Report bugs or request features
- **Email**: [your-email@example.com]
- **Documentation**: [README.md](README.md), [MIGRATION.md](MIGRATION.md), [ARCHITECTURE.md](ARCHITECTURE.md)

### Updating the Pipeline

```bash
cd ~/RIsearch_pipeline

# Pull latest changes
git pull origin main

# Rebuild Rust binary if updated
cd risearch && cargo install --path . && cd ..

# Update Python dependencies
source .venv/bin/activate
uv sync
```

### Monitoring for Updates

```bash
# Check for new releases
git fetch --tags
git tag -l
```

---

## Appendix: Full Installation Script

Copy-paste ready deployment script:

```bash
#!/bin/bash
set -e  # Exit on error

echo "=== RIsearch Pipeline Installation ==="

# 1. Install Rust
if ! command -v cargo &> /dev/null; then
    echo "Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source $HOME/.cargo/env
fi

# 2. Clone repository
if [ ! -d "RIsearch_pipeline" ]; then
    echo "Cloning repository..."
    git clone https://github.com/your-org/RIsearch_pipeline.git
fi
cd RIsearch_pipeline

# 3. Build RIsearch binary
echo "Building RIsearch binary..."
cd risearch
cargo install --path .
cd ..

# 4. Install Python pipeline
echo "Installing Python pipeline..."
if ! command -v uv &> /dev/null; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv venv .venv
source .venv/bin/activate
uv sync

# 5. Verify
echo "Verifying installation..."
risearch-pipeline --version
RIsearch --version

echo "✅ Installation complete!"
echo "Activate environment with: source .venv/bin/activate"
```

Save as `install.sh`, then:

```bash
chmod +x install.sh
./install.sh
```
