from pathlib import Path
from typing import Generator, List, Tuple


def merge_intervals(
    intervals: List[Tuple[int, int]], padding: int = 0
) -> List[Tuple[int, int]]:
    """
    Merge overlapping or nearby intervals.

    Intervals within `padding` bp of each other are merged into one.
    Input intervals are (start, end) with start < end (0-based, half-open).

    Args:
        intervals: List of (start, end) tuples.
        padding: Maximum gap between intervals to still merge them.

    Returns:
        Sorted list of non-overlapping merged (start, end) tuples.

    Complexity: O(N log N) time (dominated by sort), O(N) space.
    """
    if not intervals:
        return []

    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged: List[Tuple[int, int]] = [sorted_ivs[0]]

    for start, end in sorted_ivs[1:]:
        prev_start, prev_end = merged[-1]
        # Merge if overlapping or within padding distance
        if start <= prev_end + padding:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged


def read_fasta(path: Path) -> Generator[Tuple[str, str], None, None]:
    """
    Simple FASTA reader.
    Yields (header, sequence).
    """
    header = None
    sequence = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]  # Take first word as ID
                sequence = []
            else:
                sequence.append(line)

        if header:
            yield header, "".join(sequence)


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    complement_map = str.maketrans("ATCGatcg", "TAGCtagc")
    return sequence.translate(complement_map)[::-1]
