import shutil
import unittest
from pathlib import Path
import numpy as np
import polars as pl
from RIsearch_pipeline.services.helpers import merge_intervals
from RIsearch_pipeline.services.accessibility import (
    GenomeAccessibilityService,
    AccessibilityError,
)
from RIsearch_pipeline.services.probability import ProbabilityService


class TestMergeIntervals(unittest.TestCase):
    """Test the merge_intervals utility."""

    def test_empty(self):
        self.assertEqual(merge_intervals([]), [])

    def test_single(self):
        self.assertEqual(merge_intervals([(10, 20)]), [(10, 20)])

    def test_no_overlap(self):
        result = merge_intervals([(10, 20), (50, 60), (100, 110)])
        self.assertEqual(result, [(10, 20), (50, 60), (100, 110)])

    def test_overlapping(self):
        result = merge_intervals([(10, 30), (20, 50), (45, 60)])
        self.assertEqual(result, [(10, 60)])

    def test_adjacent_no_padding(self):
        result = merge_intervals([(10, 20), (20, 30)], padding=0)
        self.assertEqual(result, [(10, 30)])

    def test_with_padding(self):
        """Intervals within padding distance should merge."""
        result = merge_intervals([(10, 20), (50, 60)], padding=30)
        self.assertEqual(result, [(10, 60)])

    def test_with_padding_no_merge(self):
        """Intervals beyond padding distance should not merge."""
        result = merge_intervals([(10, 20), (100, 110)], padding=30)
        self.assertEqual(result, [(10, 20), (100, 110)])

    def test_unsorted_input(self):
        """Input order should not matter."""
        result = merge_intervals([(50, 60), (10, 20), (30, 40)], padding=15)
        # gap 20→30 = 10 ≤ 15, gap 40→50 = 10 ≤ 15 → all merge
        self.assertEqual(result, [(10, 60)])

    def test_nested_intervals(self):
        """Nested intervals should be handled correctly."""
        result = merge_intervals([(10, 100), (20, 30), (50, 60)])
        self.assertEqual(result, [(10, 100)])


class TestComputeBindingSiteAccessibility(unittest.TestCase):
    """Test binding-site-only accessibility computation."""

    def setUp(self):
        self.test_dir = Path("tests/output/binding_site_acc_test")
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
        self.test_dir.mkdir(parents=True)

        # Create a short genome FASTA
        self.fasta_path = self.test_dir / "genome.fa"
        # Sequence long enough for pfl_fold_up (100 bp)
        seq = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU" * 2
        # Use DNA alphabet — the code should convert T→U
        seq_dna = seq.replace("U", "T")
        with open(self.fasta_path, "w") as f:
            f.write(f">chr1\n{seq_dna}\n")

        # Create a mock RIsearch output file (per-siRNA)
        self.risearch_dir = self.test_dir / "risearch"
        self.risearch_dir.mkdir()
        risearch_file = self.risearch_dir / "siRNA1.tsv"
        # RIsearch2 format: siRNA_id, ?, ?, chrom, start, end, strand, energy
        with open(risearch_file, "w") as f:
            f.write("siRNA1\tseq1\t0\tchr1\t10\t20\t+\t-5.0\n")
            f.write("siRNA1\tseq1\t0\tchr1\t30\t40\t+\t-3.0\n")
            f.write("siRNA1\tseq1\t0\tchr1\t50\t60\t-\t-4.0\n")

    def tearDown(self):
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)

    def test_compute_binding_sites(self):
        """Test that binding-site accessibility produces valid Parquet output."""
        service = GenomeAccessibilityService(self.test_dir)
        output_path = self.test_dir / "output.parquet"

        result_path = service.compute_binding_site_accessibility(
            genome_path=self.fasta_path,
            risearch_dir=self.risearch_dir,
            output_path=output_path,
            window_size=10,
            max_span=5,
            unpaired_prob=3,
        )

        self.assertTrue(result_path.exists())

        # Load and verify structure
        df = pl.read_parquet(result_path)
        self.assertIn("chrom", df.columns)
        self.assertIn("start", df.columns)
        self.assertIn("end", df.columns)
        self.assertIn("strand", df.columns)
        self.assertIn("opening_energy", df.columns)

        # Should have 3 rows (3 unique binding sites)
        self.assertEqual(df.height, 3)

        # Opening energies should be non-negative
        self.assertTrue((df["opening_energy"] >= 0).all())


class TestParquetJoinAnnotation(unittest.TestCase):
    """Test that precomputed Parquet accessibility joins correctly."""

    def test_join_path(self):
        """Verify Parquet join produces correct opening_energy column."""
        # Create a mock precomputed accessibility DataFrame
        acc_df = pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "start": [10, 20, 30],
                "end": [20, 30, 40],
                "strand": ["+", "+", "-"],
                "opening_energy": [1.5, 2.0, 3.0],
            }
        )

        # Create ProbabilityService with precomputed accessibility
        service = ProbabilityService(
            accessibility_service=None,
            precomputed_accessibility=acc_df,
        )

        # Create prediction DataFrame
        predictions = pl.DataFrame(
            {
                "energy": [-10.0, -5.0, -3.0, -2.0],
                "chrom": ["chr1", "chr1", "chr2", "chr3"],
                "start": [10, 20, 30, 50],
                "end": [20, 30, 40, 60],
                "strand": ["+", "+", "-", "+"],
            }
        )

        result = service._annotate_opening_energy(predictions)

        self.assertIn("opening_energy", result.columns)
        self.assertEqual(result.height, 4)

        # Check matched values
        oe = result["opening_energy"].to_list()
        self.assertAlmostEqual(oe[0], 1.5)
        self.assertAlmostEqual(oe[1], 2.0)
        self.assertAlmostEqual(oe[2], 3.0)
        # chr3 not in precomputed → default 10.0
        self.assertAlmostEqual(oe[3], 10.0)


if __name__ == "__main__":
    unittest.main()
