"""
download_baselines.py - Download baseline sequences for entropy comparison.

    Conserved:  RplA (50S ribosomal protein L1) - reviewed SwissProt bacteria
    Divergent:  MFS transporters (500 reviewed) - highly variable family

Output:
    data/external/conservation_baseline/rplA.fasta
    data/external/conservation_baseline/mfs_transporter.fasta
"""

import sys
import requests
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

OUT_DIR = Path(PATHS["data_external"]) / "conservation_baseline"
OUT_DIR.mkdir(parents=True, exist_ok=True)

URLS = {
    "rplA.fasta": (
        "https://rest.uniprot.org/uniprotkb/stream"
        "?query=gene:rplA+AND+reviewed:true+AND+taxonomy_id:2"
        "&format=fasta"
    ),
    "mfs_transporter.fasta": (
        "https://rest.uniprot.org/uniprotkb/stream"
        "?query=family:%22major+facilitator+superfamily%22+AND+reviewed:true"
        "&format=fasta&size=500"
    ),
}

for filename, url in URLS.items():
    out_path = OUT_DIR / filename
    if out_path.exists():
        print(f"[SKIP] {filename} already exists")
        continue
    print(f"[DOWNLOAD] {filename} ...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    out_path.write_text(r.text)
    n = r.text.count(">")
    print(f"[DONE] {filename} -> {n} sequences -> {out_path}")

print("Done. Next: align with MAFFT, then run entropy.py --baseline-conserved / --baseline-divergent")