from pathlib import Path
from typing import Generator, Tuple


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
