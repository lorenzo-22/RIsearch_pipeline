import typer
from pathlib import Path
from typing import List, Optional
from Bio import SeqIO
from loguru import logger
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TimeElapsedColumn,
)

from RIsearch_pipeline.services.orthodb_client import OrthoDBClient
from RIsearch_pipeline.services.ncbi_client import NCBIClient

console = Console()


def run(
    target_gene_file: Path = typer.Option(
        ...,
        "--target-gene",
        "-g",
        help="Path to file containing target gene symbol (e.g. PSMC2)",
    ),
    species_list_file: Path = typer.Option(
        ...,
        "--species-list",
        "-s",
        help="Path to file containing list of species (scientific names)",
    ),
    output_dir: Path = typer.Option(
        Path("./"), "--output-dir", "-o", help="Directory to save output files"
    ),
    email: str = typer.Option(..., "--email", help="Email address for NCBI Entrez API"),
    insecta_level: int = typer.Option(
        50557, "--insecta-level", help="OrthoDB Taxonomy ID for Insecta level"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Enable verbose logging"
    ),
    plot_lengths: bool = typer.Option(
        False, "--plot-lengths", "-p", help="Generate ortholog length comparison plot"
    ),
    run_msa: bool = typer.Option(
        False, "--run-msa", "-m", help="Run multiple sequence alignment per species"
    ),
):
    """
    Download orthologs (OrthoDB) and transcriptomes (NCBI) for a list of species.
    """
    if verbose:
        logger.level("DEBUG")

    # Imports for analysis
    from RIsearch_pipeline.analysis.length_comparison import plot_ortholog_lengths
    from RIsearch_pipeline.analysis.msa import run_msa_per_species

    # 1. Setup Directories
    orthologs_dir = output_dir / "orthologs"
    transcriptomes_dir = output_dir / "transcriptomes"
    alignments_dir = output_dir / "alignments"

    orthologs_dir.mkdir(parents=True, exist_ok=True)
    transcriptomes_dir.mkdir(parents=True, exist_ok=True)
    if run_msa:
        alignments_dir.mkdir(parents=True, exist_ok=True)

    # 2. Read Input Files
    try:
        target_gene = target_gene_file.read_text().strip()
        species_list = [
            line.strip()
            for line in species_list_file.read_text().splitlines()
            if line.strip()
        ]

        console.print(f"[bold green]Target Gene:[/bold green] {target_gene}")
        console.print(
            f"[bold green]Species List:[/bold green] {len(species_list)} species loaded"
        )
    except Exception as e:
        console.print(f"[bold red]Error reading input files:[/bold red] {e}")
        raise typer.Exit(code=1)

    # 3. OrthoDB Operations
    orthodb = OrthoDBClient()

    with console.status(
        f"[bold cyan]Searching OrthoDB for {target_gene} (Level: {insecta_level})..."
    ) as status:
        ortholog_groups = orthodb.search_orthologs(target_gene, level=insecta_level)

    if not ortholog_groups:
        console.print(
            f"[yellow]No ortholog groups found for {target_gene}. Skipping ortholog download.[/yellow]"
        )
    else:
        # We assume the first group is the relevant one for now
        cluster_id = ortholog_groups[0].ortholog_id
        console.print(f"[green]Found Cluster ID:[/green] {cluster_id}")

        # Download Raw Cluster FASTA
        raw_fasta_path = orthologs_dir / f"cluster_{cluster_id}.fasta"
        if not raw_fasta_path.exists():
            with console.status(f"[cyan]Downloading cluster FASTA..."):
                orthodb.download_fasta(cluster_id, raw_fasta_path)
        else:
            console.print(f"[dim]Using cached cluster FASTA: {raw_fasta_path}[/dim]")

        # Filter and Split
        # OrthoDB Header format: >{pub_gene_id} {tax_id} {organism_name} ...
        # But actually format might vary. Let's inspect or handle robustly.
        # Based on docs: >pub_gene_id tax_id organism name

        with console.status("[cyan]Filtering and splitting orthologs...") as status:
            count_saved = 0
            for record in SeqIO.parse(raw_fasta_path, "fasta"):
                # Check if organism name matches any in our species list
                # Header description usually contains organism name
                # Simple containment check

                matched_species = None
                for species in species_list:
                    # Case insensitive check? Or exact?
                    # Normalized check
                    if species.lower() in record.description.lower():
                        matched_species = species
                        break

                if matched_species:
                    # Save to individual file
                    # File name: {species_sanitized}_{gene_id}.fa
                    safe_species = matched_species.replace(" ", "_")
                    safe_id = record.id.replace(":", "_")
                    out_file = orthologs_dir / f"{safe_species}_{safe_id}.fa"

                    SeqIO.write(record, out_file, "fasta")
                    count_saved += 1

            console.print(
                f"[green]Saved {count_saved} ortholog sequences for specified species[/green]"
            )

    # 4. NCBI Transcriptome Download
    ncbi = NCBIClient(email=email)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            "Downloading Transcriptomes...", total=len(species_list)
        )

        for species in species_list:
            progress.update(task, description=f"Processing {species}...")

            # Check if already exists
            safe_species = species.replace(" ", "_")
            out_file = transcriptomes_dir / f"{safe_species}_TSA.fa"

            if out_file.exists():
                console.print(f"  [dim]Skipping {species} (already exists)[/dim]")
                progress.advance(task)
                continue

            # Download complete transcriptome assembly
            try:
                ncbi.download_transcriptome_assembly(species, out_file)
                console.print(
                    f"  [green]Downloaded Transcriptome for {species}[/green]"
                )
            except Exception as e:
                console.print(
                    f"  [red]Failed to download transcriptome for {species}: {e}[/red]"
                )

            progress.advance(task)

    # 5. Analysis
    if plot_lengths:
        with console.status("[cyan]Generating length comparison plot..."):
            plot_path = output_dir / "ortholog_lengths.png"
            plot_ortholog_lengths(orthologs_dir, plot_path)
            console.print(f"[green]Length comparison plot saved to {plot_path}[/green]")

    if run_msa:
        with console.status("[cyan]Running Multiple Sequence Alignment..."):
            run_msa_per_species(orthologs_dir, alignments_dir)
            console.print(
                f"[green]MSA completed. Alignments saved to {alignments_dir}[/green]"
            )

    console.print(f"\n[bold green]Orthologs pipeline completed![/bold green]")
