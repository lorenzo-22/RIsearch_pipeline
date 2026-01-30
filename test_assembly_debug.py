from Bio import Entrez
from ftplib import FTP

Entrez.email = "test@example.com"

species = "Drosophila melanogaster"

# Query Assembly
query = (
    f'"{species}"[Organism] AND ("latest refseq"[filter] OR "latest genbank"[filter])'
)
handle = Entrez.esearch(db="assembly", term=query, retmax=5)
record = Entrez.read(handle, validate=False)
handle.close()

print(f"Found {len(record['IdList'])} assemblies")

for assembly_id in record["IdList"][:3]:
    handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
    summary = Entrez.read(handle, validate=False)
    handle.close()

    doc_sum = summary["DocumentSummarySet"]["DocumentSummary"][0]

    print(f"\nAssembly: {doc_sum.get('AssemblyAccession')}")
    print(f"  Name: {doc_sum.get('AssemblyName')}")
    print(f"  Status: {doc_sum.get('AssemblyStatus')}")

    ftp_path = doc_sum.get("FtpPath_RefSeq") or doc_sum.get("FtpPath_GenBank")
    if ftp_path:
        print(f"  FTP: {ftp_path}")

        # List files
        ftp_url_parts = ftp_path.replace("ftp://", "").split("/", 1)
        ftp_host = ftp_url_parts[0]
        ftp_dir = "/" + ftp_url_parts[1]

        try:
            ftp = FTP(ftp_host)
            ftp.login()
            ftp.cwd(ftp_dir)
            files = ftp.nlst()
            ftp.quit()

            rna_files = [
                f for f in files if "rna" in f or "cds" in f or "transcript" in f
            ]
            print(f"  RNA/CDS files: {rna_files}")
        except Exception as e:
            print(f"  FTP error: {e}")
