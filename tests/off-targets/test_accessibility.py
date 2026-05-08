import shutil
import unittest
from pathlib import Path
import numpy as np
from RIsearch_pipeline.services.accessibility import (
    GenomeAccessibilityService,
    AccessibilityError,
)
from RIsearch_pipeline.services.helpers import read_fasta

# Dummy FASTA content
DUMMY_FASTA = """>chr1
ACGTACGTACGT
>chr2
GGGGCCCC
"""


class TestAccessibility(unittest.TestCase):
    def setUp(self):
        self.test_dir = Path("tests/output/accessibility_test")
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
        self.test_dir.mkdir(parents=True)

        self.fasta_path = self.test_dir / "genome.fa"
        with open(self.fasta_path, "w") as f:
            f.write(DUMMY_FASTA)

    def tearDown(self):
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)

    def test_read_fasta(self):
        """Test the helper FASTA reader."""
        entries = list(read_fasta(self.fasta_path))
        self.assertEqual(len(entries), 2)
        self.assertEqual(entries[0][0], "chr1")
        self.assertEqual(entries[0][1], "ACGTACGTACGT")
        self.assertEqual(entries[1][0], "chr2")
        self.assertEqual(entries[1][1], "GGGGCCCC")

    def test_compute_accessibility_real(self):
        """Test compute logic with REAL ViennaRNA."""
        import polars as pl

        service = GenomeAccessibilityService(self.test_dir)
        # Use small window/span for short sequence
        # Sequence: ACGT... length 12
        # Use W=10, L=5, u=3
        results = service.compute_genome_accessibility(
            self.fasta_path, window_size=10, max_span=5, unpaired_prob=3
        )

        # Results keyed by chromosome name
        self.assertIn("chr1", results)

        # Verify Parquet file created
        chr1_path = results["chr1"]
        self.assertTrue(chr1_path.exists())
        self.assertEqual(chr1_path.suffix, ".parquet")

        # Read and validate content
        df = pl.read_parquet(chr1_path)
        # Should have 12 positions × 2 strands = 24 rows
        self.assertEqual(df.height, 24)
        # Must contain position, strand, u1, u2, u3
        self.assertIn("position", df.columns)
        self.assertIn("strand", df.columns)
        self.assertIn("u1", df.columns)
        self.assertIn("u2", df.columns)
        self.assertIn("u3", df.columns)
        # Check some non-default values exist (accessibility computed)
        u1_vals = df.filter(pl.col("strand") == "+")["u1"].to_list()
        self.assertTrue(
            any(v != 25.5 for v in u1_vals),
            "Should have non-default accessibility values",
        )

    def test_query_service(self):
        """Test querying the created profiles."""
        service = GenomeAccessibilityService(self.test_dir)

        # Manually create a profile for testing query
        chrom = "chr1"
        data = np.array([0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float32)
        path = service.get_profile_path(chrom)
        np.save(path, data)

        # Query
        res = service.query(chrom, 1, 4)  # indices 1, 2, 3 -> 0.2, 0.3, 0.4
        self.assertEqual(len(res), 3)
        np.testing.assert_array_almost_equal(
            res, np.array([0.2, 0.3, 0.4], dtype=np.float32)
        )

    def test_lru_eviction(self):
        """Test LRU eviction with max_cached=2."""
        service = GenomeAccessibilityService(self.test_dir, max_cached=2)

        # Create 3 chromosome profiles on disk
        for chrom in ["chrA", "chrB", "chrC"]:
            data = np.array([1.0, 2.0, 3.0], dtype=np.float32)
            path = service.get_profile_path(chrom, "+")
            np.save(path, data)

        # Load chrA, then chrB -> cache: [chrA_+, chrB_+]
        service.query("chrA", 0, 1)
        service.query("chrB", 0, 1)
        self.assertEqual(len(service._profiles), 2)
        self.assertIn("chrA_+", service._profiles)
        self.assertIn("chrB_+", service._profiles)

        # Load chrC -> evicts chrA (LRU), cache: [chrB_+, chrC_+]
        service.query("chrC", 0, 1)
        self.assertEqual(len(service._profiles), 2)
        self.assertNotIn("chrA_+", service._profiles)
        self.assertIn("chrB_+", service._profiles)
        self.assertIn("chrC_+", service._profiles)

        # Re-query chrA -> reloads from disk, evicts chrB
        res = service.query("chrA", 0, 3)
        self.assertEqual(len(service._profiles), 2)
        self.assertIn("chrA_+", service._profiles)
        self.assertNotIn("chrB_+", service._profiles)
        np.testing.assert_array_almost_equal(
            res, np.array([1.0, 2.0, 3.0], dtype=np.float32)
        )

    def test_query_single(self):
        """Test query_single returns a single quantized float."""
        service = GenomeAccessibilityService(self.test_dir)

        # 1D profile
        data = np.array([0.15, 0.25, 0.35, 0.45, 0.55], dtype=np.float32)
        path = service.get_profile_path("chr1", "+")
        np.save(path, data)

        # query_single uses 1-based coords: start=2, end=4 -> 0-based [1,4)
        # strand "+": picks end0-1 = index 3 -> 0.45 -> quantized round(4.5)/10 = 0.4
        val = service.query_single("chr1", 2, 4, "+")
        self.assertAlmostEqual(val, 0.4, places=1)


if __name__ == "__main__":
    unittest.main()
