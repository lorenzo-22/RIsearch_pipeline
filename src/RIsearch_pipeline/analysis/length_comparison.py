"""Ortholog length comparison analysis."""

from pathlib import Path
from typing import Dict, List
from Bio import SeqIO
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns
import json


def extract_species_from_header(description: str) -> str:
    """
    Extract species name from OrthoDB FASTA header.

    Header format: >id {"organism_name":"Species name",...}
    """
    try:
        # Find JSON portion in header
        json_start = description.find("{")
        if json_start != -1:
            json_str = description[json_start:]
            data = json.loads(json_str)
            return data.get("organism_name", "Unknown")
    except (json.JSONDecodeError, KeyError):
        pass
    return "Unknown"


def collect_ortholog_lengths(orthologs_dir: Path) -> Dict[str, List[int]]:
    """
    Collect sequence lengths grouped by species.

    Args:
        orthologs_dir: Directory containing ortholog FASTA files.

    Returns:
        Dict mapping species name to list of sequence lengths.
    """
    species_lengths: Dict[str, List[int]] = {}

    # Read all .fa files (excluding cluster files)
    for fasta_file in orthologs_dir.glob("*.fa"):
        if "cluster_" in fasta_file.name:
            continue

        for record in SeqIO.parse(fasta_file, "fasta"):
            species = extract_species_from_header(record.description)
            seq_len = len(record.seq)

            if species not in species_lengths:
                species_lengths[species] = []
            species_lengths[species].append(seq_len)

    logger.info(f"Collected lengths for {len(species_lengths)} species")
    return species_lengths


def plot_ortholog_lengths(orthologs_dir: Path, output_path: Path):
    """
    Generate box plot showing ortholog length distribution per species.

    Args:
        orthologs_dir: Directory containing ortholog FASTA files.
        output_path: Path to save the plot (PNG).
    """
    species_lengths = collect_ortholog_lengths(orthologs_dir)

    if not species_lengths:
        logger.warning("No ortholog data found for plotting")
        return

    # Prepare data for seaborn
    plot_data = []
    for species, lengths in species_lengths.items():
        for length in lengths:
            plot_data.append({"Species": species, "Length (bp)": length})

    import pandas as pd

    df = pd.DataFrame(plot_data)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.boxplot(
        data=df,
        x="Species",
        y="Length (bp)",
        hue="Species",
        legend=False,
        ax=ax,
        palette="viridis",
    )

    ax.set_title("Ortholog CDS Length Distribution by Species")
    ax.set_xlabel("Species")
    ax.set_ylabel("Sequence Length (bp)")

    # Force integer ticks on Y-axis
    from matplotlib.ticker import MaxNLocator

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Rotate x labels for readability
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Save
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Length comparison plot saved to {output_path}")
