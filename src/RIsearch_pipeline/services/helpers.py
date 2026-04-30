from pathlib import Path
from typing import Generator


def merge_intervals(intervals: list[tuple[int, int]], padding: int = 0) -> list[tuple[int, int]]:
    """Merge overlapping or nearby intervals (within `padding` bp). O(N log N)."""
    if not intervals:
        return []

    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_ivs[0]]

    for start, end in sorted_ivs[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + padding:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged


def read_fasta(path: Path) -> Generator[tuple[str, str], None, None]:
    """Yield (id, sequence) pairs from a FASTA file."""
    header = None
    sequence: list[str] = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if header:
        yield header, "".join(sequence)


def reverse_complement(sequence: str) -> str:
    complement_map = str.maketrans("ATCGatcg", "TAGCtagc")
    return sequence.translate(complement_map)[::-1]
