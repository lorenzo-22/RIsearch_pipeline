#!/usr/bin/env python3
"""
Generate 10,000 siRNA FASTA and RIsearch output files for testing scalability.
"""

from pathlib import Path
import shutil


def generate_10k_dataset():
    # Input paths (using single siRNA data as template)
    base_dir = Path("tests/off-targets/single/input")
    sirna_fa = base_dir / "sirna.fa"
    risearch_out = base_dir / "risearch_siRNAID.out"

    # Output paths
    out_dir = Path("tests/off-targets/10k_dataset")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_fa = out_dir / "sirnas_10k.fa"
    out_risearch = out_dir / "risearch_10k.out"

    print(f"Generating 10k dataset in {out_dir}...")

    # 1. Read templates
    with open(sirna_fa, "r") as f:
        # Assume valid FASTA: >header\nSEQUENCE
        lines = f.readlines()
        seq = lines[1].strip()

    with open(risearch_out, "r") as f:
        risearch_lines = f.readlines()
        # Filter comments and keep data lines
        data_lines = [l for l in risearch_lines if not l.startswith("#") and l.strip()]

    # 2. Generate FASTA
    print(f"Writing {out_fa}...")
    with open(out_fa, "w") as f:
        for i in range(1, 10001):
            f.write(f">siRNAID{i}\n{seq}\n")

    # 3. Generate RIsearch output
    print(f"Writing {out_risearch}...")
    with open(out_risearch, "w") as f:
        # Write header roughly matching RIsearch format if needed,
        # but our parser handles no header. We'll skip complex headers.

        for i in range(1, 10001):
            id_str = f"siRNAID{i}"
            for line in data_lines:
                parts = line.split()
                if not parts:
                    continue

                # Replace the first column (QueryID)
                # RIsearch output: QueryID TargetID ...
                # But wait, let's check the format in 'risearch_siRNAID.out'
                # Standard RIsearch: QueryID QStart QEnd TargetID TStart TEnd Strand Energy ...

                # We need to construct the line carefully to preserve tab separation
                # The file 'risearch_siRNAID.out' typically has:
                # siRNAID	transcript_22	433	444	+	-11.81
                # Wait, that's not standard 8-col format. Let's check the file content first.
                pass

    # Actually, simpler approach: read line, replace "siRNAID" (literal) with "siRNAID{i}"
    # This assumes the template file literally uses "siRNAID" as the identifier.

    with open(out_risearch, "w") as f:
        for i in range(1, 10001):
            new_id = f"siRNAID{i}"
            for line in data_lines:
                # Replace the literal string "siRNAID" with the new ID
                # Be careful not to replace it if it appears elsewhere, but it's likely just the first col.
                # Just replace the first occurrence.
                new_line = line.replace("siRNAID", new_id, 1)
                f.write(new_line)

    print("Done.")
    print(f"Generated {out_fa} ({out_fa.stat().st_size / 1024 / 1024:.2f} MB)")
    print(
        f"Generated {out_risearch} ({out_risearch.stat().st_size / 1024 / 1024:.2f} MB)"
    )


if __name__ == "__main__":
    generate_10k_dataset()
