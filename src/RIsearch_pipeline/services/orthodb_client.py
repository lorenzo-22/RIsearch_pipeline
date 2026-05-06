import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import httpx
from loguru import logger


@dataclass
class OrthologInfo:
    """OrthoDB ortholog information."""

    ortholog_id: str
    organism_name: str
    ncbi_taxid: str
    gene_id: str


class OrthoDBClient:
    """Client for OrthoDB API."""

    BASE_URL = "https://data.orthodb.org/current"
    INSECTA_LEVEL = 50557

    def __init__(self, rate_limit_delay: float = 1.0):
        """
        Initialize OrthoDB client.

        Args:
            rate_limit_delay: Seconds to wait between requests (default: 1.0).
        """
        self.rate_limit_delay = rate_limit_delay
        self.last_request_time = 0.0
        self.client = httpx.Client(base_url=self.BASE_URL, timeout=30.0)

    def _wait_for_rate_limit(self):
        """Enforce rate limiting."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit_delay:
            time.sleep(self.rate_limit_delay - elapsed)
        self.last_request_time = time.time()

    def search_orthologs(
        self, gene_id: str, level: int = INSECTA_LEVEL
    ) -> List[OrthologInfo]:
        """
        Search for orthologs of a gene at a specific taxonomic level.

        Args:
            gene_id: Gene symbol or ID (e.g. "PSMC2", "FBgn0003023").
            level: Taxonomy ID for the level (default: 50557 for Insecta).

        Returns:
            List of OrthologInfo objects.
        """
        self._wait_for_rate_limit()

        params = {
            "query": gene_id,
            "level": level,
        }

        logger.debug(f"Searching OrthoDB for {gene_id} at level {level}...")
        try:
            response = self.client.get("/search", params=params)
            response.raise_for_status()
            data = response.json()

            results = []

            # If data is a list of strings, these are cluster IDs
            clusters = data.get("data", [])
            if not clusters:
                logger.warning(f"No orthologs found for {gene_id}")
                return []

            cluster_id = clusters[0]

            return [
                OrthologInfo(
                    ortholog_id=cluster_id,
                    organism_name="Unknown",
                    ncbi_taxid="0",
                    gene_id=gene_id,
                )
            ]

        except httpx.HTTPError as e:
            logger.error(f"OrthoDB API error: {e}")
            raise

    def download_fasta(
        self,
        cluster_id: str,
        output_path: Path,
        species_filter: Optional[List[str]] = None,
    ):
        """
        Download FASTA sequences for an ortholog cluster.

        Args:
            cluster_id: OrthoDB Cluster ID.
            output_path: Path to save FASTA.
            species_filter: Optional list of species names to filter.
        """
        self._wait_for_rate_limit()

        logger.info(f"Downloading CDS FASTA for cluster {cluster_id}...")
        params = {
            "id": cluster_id,
            "seqtype": "cds",  # Request coding DNA sequences
        }

        try:
            with self.client.stream("GET", "/fasta", params=params) as response:
                response.raise_for_status()

                with open(output_path, "wb") as f:
                    for chunk in response.iter_bytes():
                        f.write(chunk)

            logger.info(f"Saved orthologs to {output_path}")

        except httpx.HTTPError as e:
            logger.error(f"Failed to download FASTA: {e}")
            raise
