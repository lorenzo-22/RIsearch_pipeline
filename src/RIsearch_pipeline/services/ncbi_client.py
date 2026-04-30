from pathlib import Path
from typing import List
from Bio import Entrez
from loguru import logger
import urllib.error


class NCBIClient:
    """Client for NCBI Entrez API."""

    def __init__(self, email: str):
        """
        Initialize NCBI client.

        Args:
            email: User email (required by NCBI).
        """
        Entrez.email = email
        # Start with default tool name
        Entrez.tool = "RIsearchPipeline"

    def search_transcriptomes(self, species_name: str) -> List[str]:
        """
        Search for Transcriptome Shotgun Assembly (TSA) records for a species.

        Args:
            species_name: Scientific name (e.g. "Drosophila melanogaster").

        Returns:
            List of GenBank Accession IDs (e.g. "JI445566").
        """
        # Query for mRNA sequences (RefSeq/GenBank) which represent the transcriptome
        # TSA keyword often returns Master records that cannot be downloaded directly as FASTA
        query = f'"{species_name}"[Organism] AND biomol_mrna[PROP]'

        # We might want to prioritize RefSeq or specific types, but TSA is requested.
        # Also could restrict to "biomol mrna"[Properties] if we want mRNA.

        logger.debug(f"Searching NCBI for transcriptomes of {species_name}...")

        try:
            handle = Entrez.esearch(
                db="nucleotide", term=query, retmax=50
            )  # Limit to 50 for now
            record = Entrez.read(handle)
            handle.close()

            id_list = record["IdList"]
            logger.info(f"Found {len(id_list)} TSA records for {species_name}")
            return id_list

        except Exception as e:
            logger.error(f"NCBI Search failed for {species_name}: {e}")
            return []

    def download_transcriptome(self, accession_id: str, output_path: Path):
        """
        Download sequence data for an accession.

        Args:
            accession_id: NCBI Accession ID.
            output_path: Path to save FASTA file.
        """
        logger.info(f"Downloading {accession_id} to {output_path}...")

        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession_id, rettype="fasta", retmode="text"
            )

            with open(output_path, "w") as f:
                f.write(handle.read())

            handle.close()
            logger.debug(f"Download complete: {output_path}")

        except urllib.error.HTTPError as e:
            if e.code == 429:
                logger.warning("NCBI Rate limit exceeded. Waiting...")
                # Biopython handles some retries, but we might need explicit sleep
                # For now, let's just raise
            logger.error(f"Failed to download {accession_id}: {e}")
            raise
        except Exception as e:
            logger.error(f"Error downloading {accession_id}: {e}")
            raise

    def download_transcriptome_assembly(self, species_name: str, output_path: Path):
        """
        Download complete transcriptome assembly for a species.

        Args:
            species_name: Scientific name (e.g. "Drosophila melanogaster").
            output_path: Path to save transcriptome FASTA file.
        """
        logger.info(f"Searching for transcriptome assembly for {species_name}...")

        try:
            # Query Assembly database - prioritize RefSeq (better annotated)
            query = f'"{species_name}"[Organism] AND "latest refseq"[filter]'
            handle = Entrez.esearch(db="assembly", term=query, retmax=5)
            record = Entrez.read(handle, validate=False)
            handle.close()

            id_list = record["IdList"]

            if not id_list:
                logger.warning(f"No assembly found for {species_name}")
                raise ValueError(f"No assembly found for {species_name}")

            # Get assembly details
            assembly_id = id_list[0]
            logger.debug(f"Found assembly ID: {assembly_id}")

            handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
            summary = Entrez.read(handle, validate=False)
            handle.close()

            # Extract FTP path - prefer RefSeq (better annotated)
            doc_sum = summary["DocumentSummarySet"]["DocumentSummary"][0]

            ftp_path = doc_sum.get("FtpPath_RefSeq")

            if not ftp_path:
                logger.warning(f"No FTP path available for {species_name}")
                raise ValueError(f"No FTP downloadable assembly for {species_name}")

            # Construct CDS FASTA filename
            # FTP structure: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../GCF_XXX_YYY/
            # Files: GCF_XXX_YYY_*.cds_from_genomic.fna.gz
            import gzip

            assembly_accession = doc_sum.get("AssemblyAccession", "")
            assembly_name = doc_sum.get("AssemblyName", "")
            logger.info(f"Downloading {assembly_accession} ({assembly_name})...")

            # List FTP directory to find CDS file
            from ftplib import FTP

            ftp_url_parts = ftp_path.replace("ftp://", "").split("/", 1)
            ftp_host = ftp_url_parts[0]
            ftp_dir = "/" + ftp_url_parts[1] if len(ftp_url_parts) > 1 else "/"

            ftp = FTP(ftp_host)
            ftp.login()
            ftp.cwd(ftp_dir)

            files = ftp.nlst()
            cds_file = None

            # Priority order: CDS > RNA > transcript
            search_patterns = [
                "cds_from_genomic.fna.gz",
                "rna_from_genomic.fna.gz",
                "rna.fna.gz",
                "_transcripts.fna.gz",
            ]

            for pattern in search_patterns:
                for f in files:
                    if pattern in f:
                        cds_file = f
                        break
                if cds_file:
                    break

            if not cds_file:
                logger.warning(f"No CDS file found in {ftp_dir}")
                ftp.quit()
                raise ValueError(f"No CDS FASTA available for {species_name}")

            # Download
            logger.info(f"Downloading {cds_file}...")
            temp_gz = output_path.with_suffix(".fna.gz")

            with open(temp_gz, "wb") as local_file:
                ftp.retrbinary(f"RETR {cds_file}", local_file.write)

            ftp.quit()

            # Decompress
            with gzip.open(temp_gz, "rb") as gz_file:
                with open(output_path, "wb") as out_file:
                    out_file.write(gz_file.read())

            temp_gz.unlink()  # Remove .gz
            logger.info(f"Transcriptome saved to {output_path}")

        except Exception as e:
            logger.error(f"Failed to download transcriptome for {species_name}: {e}")
            raise
