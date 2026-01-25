#!/usr/bin/env python3
"""
Compare accessibility profiles between old and new pipeline.

Usage:
    python compare_accessibility.py <old_dir> <new_dir> <chrom>
"""

import sys
import numpy as np
from pathlib import Path


def load_old_profile(
    old_dir: Path, chrom: str, strand: str, u: int = 30, u_minus: int = None
):
    """Load old pipeline accessibility profile.

    Args:
        u_minus: If specified, use this u-index for minus strand (to test if reversal flips u values)
    """
    # The old pipeline uses .open.acc.bin and .rev.open.acc.bin
    # These files are seq_len * 30 bytes (uint8)
    if strand == "+":
        names = [f"{chrom}.open.acc.bin", f"{chrom}_plus.access.npy"]
        actual_u = u
    else:
        names = [f"{chrom}.rev.open.acc.bin", f"{chrom}_minus.access.npy"]
        # If u_minus is specified, use it for minus strand (for testing reversal theory)
        actual_u = u_minus if u_minus is not None else u

    for name in names:
        path = old_dir / name
        if path.exists():
            if path.suffix == ".npy":
                return np.load(path)
            elif path.suffix == ".bin":
                # 30 bytes per position (one for each u=1..30)
                data = np.fromfile(path, dtype=np.uint8)
                num_pos = len(data) // 30
                if len(data) % 30 != 0:
                    print(
                        f"⚠ Warning: {name} size is not a multiple of 30. Truncating."
                    )
                    data = data[: num_pos * 30]

                matrix = data.reshape((num_pos, 30))
                # Use actual_u for extraction (allows testing minus strand u-flip)
                return matrix[:, actual_u - 1].astype(np.float32) / 10.0
            else:
                try:
                    return np.loadtxt(path)
                except Exception:
                    pass

    raise FileNotFoundError(
        f"Could not find old profile for {chrom} {strand} in {old_dir}"
    )


def load_new_profile(new_dir: Path, chrom: str, strand: str):
    """Load new pipeline accessibility profile."""
    suffix = "plus" if strand == "+" else "minus"
    path = new_dir / f"{chrom}_{suffix}.access.npy"
    return np.load(path)


def compare_profiles(old_profile, new_profile):
    """Compare two accessibility profiles."""
    if len(old_profile) != len(new_profile):
        print(f"⚠ Length mismatch: old={len(old_profile)}, new={len(new_profile)}")
        min_len = min(len(old_profile), len(new_profile))
        old_profile = old_profile[:min_len]
        new_profile = new_profile[:min_len]

    # Calculate statistics
    diff = new_profile - old_profile
    abs_diff = np.abs(diff)

    print("\nComparison Results:")
    print(f"  Length: {len(old_profile)}")
    print(f"  Mean absolute difference: {np.mean(abs_diff):.6f}")
    print(f"  Max absolute difference: {np.max(abs_diff):.6f}")
    print(f"  Correlation: {np.corrcoef(old_profile, new_profile)[0, 1]:.6f}")

    # Count exact matches and close matches (tolerance based on 0.1 quantization)
    exact_matches = np.sum(old_profile == new_profile)
    close_matches = np.sum(np.isclose(old_profile, new_profile, atol=0.05, rtol=0.01))
    print(f"  Exact matches: {exact_matches} / {len(old_profile)}")
    print(
        f"  Close matches (atol=0.05): {close_matches} / {len(old_profile)} ({100 * close_matches / len(old_profile):.1f}%)"
    )

    # Show distribution of differences
    percentiles = [0, 25, 50, 75, 95, 99, 100]
    print("\n  Difference percentiles:")
    for p in percentiles:
        val = np.percentile(abs_diff, p)
        print(f"    {p:3d}%: {val:.6f}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    old_dir = Path(sys.argv[1])
    new_dir = Path(sys.argv[2])
    chrom = sys.argv[3]

    for strand in ["+", "-"]:
        print(f"\n{'=' * 60}")
        print(f"Comparing {chrom} strand {strand}")
        print("=" * 60)

        try:
            # For minus strand, test if using u=1 improves correlation
            # (the old pipeline reverses the flat array which may flip u-indices)
            u_minus_test = 1 if strand == "-" else None
            old = load_old_profile(old_dir, chrom, strand, u=30, u_minus=u_minus_test)
            new = load_new_profile(new_dir, chrom, strand)
            compare_profiles(old, new)
        except FileNotFoundError as e:
            print(f"⚠ {e}")
        except Exception as e:
            print(f"✗ Error: {e}")
