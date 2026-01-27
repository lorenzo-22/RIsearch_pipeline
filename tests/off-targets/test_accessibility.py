import shutil
import unittest
from unittest.mock import patch
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

        service = GenomeAccessibilityService(self.test_dir)
        # Use small window/span for short sequence
        # Sequence: ACGT... length 12
        # Use W=10, L=5, u=3
        results = service.compute_genome_accessibility(
            self.fasta_path, window_size=10, max_span=5, unpaired_prob=3
        )

        self.assertIn("chr1", results)

        # Verify files created
        chr1_path = results["chr1"]
        self.assertTrue(chr1_path.exists())

        # Verify content
        arr = np.load(chr1_path)
        self.assertEqual(len(arr), 12)
        self.assertTrue(arr.dtype == np.float32)
        # Check that we have some non-zero values (accessibility > 0)
        self.assertTrue(np.any(arr > 0), "Should have some non-zero accessibility")

    @patch("RIsearch_pipeline.services.accessibility.HAS_VIENNA_BINDINGS", False)
    def test_missing_bindings_raises(self):
        """Test that missing bindings raise error."""
        service = GenomeAccessibilityService(self.test_dir)
        with self.assertRaises(AccessibilityError):
            service.compute_genome_accessibility(self.fasta_path)

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


if __name__ == "__main__":
    unittest.main()
