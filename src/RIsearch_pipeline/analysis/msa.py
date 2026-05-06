"""Multiple Sequence Alignment (MSA) analysis."""

import subprocess
import shutil
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
import json


def extract_species_from_header(description: str) -> str:
    """
    Extract species name from OrthoDB FASTA header.
    Reused logic to ensure consistency.
    """
    try:
        json_start = description.find("{")
        if json_start != -1:
            json_str = description[json_start:]
            data = json.loads(json_str)
            return data.get("organism_name", "Unknown")
    except (json.JSONDecodeError, KeyError):
        pass
    return "Unknown"


def run_msa_per_species(orthologs_dir: Path, output_dir: Path):
    """
    Group orthologs by species and run MSA using MUSCLE.

    Args:
        orthologs_dir: Directory containing individual ortholog FASTA files.
        output_dir: Directory to save aligned FASTA files.
    """
    # Check if MUSCLE is installed
    muscle_exe = shutil.which("muscle")
    if not muscle_exe:
        logger.error("MUSCLE not found in PATH. Cannot run MSA.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Group sequences by species
    species_seqs: Dict[str, List[SeqRecord]] = {}

    for fasta_file in orthologs_dir.glob("*.fa"):
        if "cluster_" in fasta_file.name:
            continue

        for record in SeqIO.parse(fasta_file, "fasta"):
            species = extract_species_from_header(record.description)
            if species == "Unknown":
                continue

            if species not in species_seqs:
                species_seqs[species] = []

            # Ensure ID is unique if processing combined files,
            # though here files are individual.
            species_seqs[species].append(record)

    # Run MSA for each species
    for species, records in species_seqs.items():
        if len(records) < 2:
            logger.debug(
                f"Skipping MSA for {species}: only {len(records)} sequence(s)."
            )
            continue

        safe_species = species.replace(" ", "_")
        input_fa = output_dir / f"{safe_species}_unaligned.fa"
        output_aln = output_dir / f"{safe_species}_aligned.fa"

        # Write unaligned sequences
        SeqIO.write(records, input_fa, "fasta")

        logger.info(
            f"Running MUSCLE alignment for {species} ({len(records)} sequences)..."
        )
        try:
            # Run MUSCLE
            # muscle -in input.fa -out output.fa
            cmd = [muscle_exe, "-in", str(input_fa), "-out", str(output_aln)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Clean up input file
            input_fa.unlink()
            logger.info(f"Saved alignment to {output_aln}")

        except subprocess.CalledProcessError as e:
            logger.error(f"MSA failed for {species}: {e.stderr}")
        except Exception as e:
            logger.error(f"Error running MSA for {species}: {e}")
