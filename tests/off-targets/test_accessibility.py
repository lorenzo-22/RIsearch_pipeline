import shutil
import sys
import unittest
from pathlib import Path
import numpy as np
import polars as pl
from riot.services.accessibility import (
    GenomeAccessibilityService,
    AccessibilityError,
)
from riot.services.helpers import read_fasta

# Dummy FASTA content
DUMMY_FASTA = """>chr1
ACGTACGTACGT
>chr2
GGGGCCCC
"""


def _write_test_parquet(
    path: Path,
    positions: list[int],
    strands: list[str],
    u_values: list[list[float]],
) -> None:
    """Write a minimal accessibility Parquet for testing."""
    n_u = len(u_values[0])
    data = {"position": positions, "strand": strands}
    for i in range(1, n_u + 1):
        data[f"u{i}"] = [row[i - 1] for row in u_values]
    schema = {
        "position": pl.Int32,
        "strand": pl.Utf8,
        **{f"u{i}": pl.Float32 for i in range(1, n_u + 1)},
    }
    pl.DataFrame(data).cast(schema).write_parquet(path)


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
        # Use small window/span for short sequence (length 12, W=10, L=5, u=3)
        results = service.compute_genome_accessibility(
            self.fasta_path, window_size=10, max_span=5, unpaired_prob=3
        )

        self.assertIn("chr1", results)

        chr1_path = results["chr1"]
        self.assertTrue(chr1_path.exists())
        self.assertEqual(chr1_path.suffix, ".parquet")

        df = pl.read_parquet(chr1_path)
        # Should have 12 positions × 2 strands = 24 rows
        self.assertEqual(df.height, 24)
        self.assertIn("position", df.columns)
        self.assertIn("strand", df.columns)
        self.assertIn("u1", df.columns)
        self.assertIn("u2", df.columns)
        self.assertIn("u3", df.columns)
        u1_vals = df.filter(pl.col("strand") == "+")["u1"].to_list()
        self.assertTrue(
            any(v != 25.5 for v in u1_vals),
            "Should have non-default accessibility values",
        )

    def test_missing_bindings_raises(self):
        """Test that missing ViennaRNA raises AccessibilityError."""
        # Temporarily hide RNA from sys.modules to simulate missing bindings
        saved = sys.modules.get("RNA", None)
        sys.modules["RNA"] = None  # type: ignore[assignment]
        try:
            service = GenomeAccessibilityService(self.test_dir)
            with self.assertRaises((AccessibilityError, ImportError)):
                service.compute_genome_accessibility(self.fasta_path)
        finally:
            if saved is None:
                sys.modules.pop("RNA", None)
            else:
                sys.modules["RNA"] = saved


    def test_query_service(self):
        """Test querying a pre-computed profile."""
        service = GenomeAccessibilityService(self.test_dir)

        # Create parquet: 5 positions on + strand, one u-column
        parquet_path = self.test_dir / "chr1.accessibility.parquet"
        _write_test_parquet(
            parquet_path,
            positions=[1, 2, 3, 4, 5],
            strands=["+", "+", "+", "+", "+"],
            u_values=[[0.1], [0.2], [0.3], [0.4], [0.5]],
        )

        # query(chrom, start, end, strand) is 0-based half-open
        # Positions 1-3 (0-based) → values 0.2, 0.3, 0.4
        res = service.query("chr1", 1, 4, "+")
        self.assertEqual(len(res), 3)
        np.testing.assert_array_almost_equal(
            res[:, 0], np.array([0.2, 0.3, 0.4], dtype=np.float32)
        )

    def test_lru_eviction(self):
        """Test LRU eviction with max_cached=2."""
        service = GenomeAccessibilityService(self.test_dir, max_cached=2)

        # Create 3 chromosome profiles on disk
        for chrom in ["chrA", "chrB", "chrC"]:
            _write_test_parquet(
                self.test_dir / f"{chrom}.accessibility.parquet",
                positions=[1, 2, 3],
                strands=["+", "+", "+"],
                u_values=[[1.0], [2.0], [3.0]],
            )

        # Load chrA, then chrB -> cache: [chrA_+, chrB_+]
        service.query("chrA", 0, 1, "+")
        service.query("chrB", 0, 1, "+")
        self.assertEqual(len(service._profiles), 2)
        self.assertIn("chrA_+", service._profiles)
        self.assertIn("chrB_+", service._profiles)

        # Load chrC -> evicts chrA (LRU), cache: [chrB_+, chrC_+]
        service.query("chrC", 0, 1, "+")
        self.assertEqual(len(service._profiles), 2)
        self.assertNotIn("chrA_+", service._profiles)
        self.assertIn("chrB_+", service._profiles)
        self.assertIn("chrC_+", service._profiles)

        # Re-query chrA -> reloads from disk, evicts chrB
        res = service.query("chrA", 0, 3, "+")
        self.assertEqual(len(service._profiles), 2)
        self.assertIn("chrA_+", service._profiles)
        self.assertNotIn("chrB_+", service._profiles)
        np.testing.assert_array_almost_equal(
            res[:, 0], np.array([1.0, 2.0, 3.0], dtype=np.float32)
        )

    def test_query_single(self):
        """Test query_single returns a single quantized float."""
        service = GenomeAccessibilityService(self.test_dir)

        # 5 positions, one u-column
        _write_test_parquet(
            self.test_dir / "chr1.accessibility.parquet",
            positions=[1, 2, 3, 4, 5],
            strands=["+", "+", "+", "+", "+"],
            u_values=[[0.15], [0.25], [0.35], [0.45], [0.55]],
        )

        # query_single uses 1-based coords: start=2, end=4
        # start0=1, end0=4, interaction_len=3, col_idx=min(3,1)-1=0
        # strand "+": row_idx = end0 - 1 = 3 → value 0.45
        # quantized: round(0.45 * 10) / 10 = 0.4
        val = service.query_single("chr1", 2, 4, "+")
        self.assertAlmostEqual(val, 0.4, places=1)

    def test_compute_genome_accessibility_temp37_matches_default(self):
        """Passing temperature=37.0 must give identical output to the default (no temp arg)."""
        service_default = GenomeAccessibilityService(self.test_dir / "default")
        service_t37    = GenomeAccessibilityService(self.test_dir / "t37")

        results_default = service_default.compute_genome_accessibility(
            self.fasta_path, window_size=10, max_span=5, unpaired_prob=3
        )
        results_t37 = service_t37.compute_genome_accessibility(
            self.fasta_path, window_size=10, max_span=5, unpaired_prob=3,
            temperature=37.0,
        )

        for chrom in results_default:
            df_default = pl.read_parquet(results_default[chrom])
            df_t37     = pl.read_parquet(results_t37[chrom])
            self.assertEqual(df_default.shape, df_t37.shape)
            for col in df_default.columns:
                np.testing.assert_array_equal(
                    df_default[col].to_numpy(),
                    df_t37[col].to_numpy(),
                    err_msg=f"Column {col!r} of {chrom} differs at temperature=37",
                )

    def test_compute_sequence_accessibility_temp37_matches_default(self):
        """Passing temperature=37.0 to compute_sequence_accessibility must match default."""
        service = GenomeAccessibilityService(self.test_dir)
        seq = "ACGUACGUACGU"

        profile_default = service.compute_sequence_accessibility(
            seq, window_size=10, max_span=5, unpaired_prob=3
        )
        profile_t37 = service.compute_sequence_accessibility(
            seq, window_size=10, max_span=5, unpaired_prob=3,
            temperature=37.0,
        )

        np.testing.assert_array_equal(
            profile_default, profile_t37,
            err_msg="Profiles differ at temperature=37 vs default",
        )


if __name__ == "__main__":
    unittest.main()
