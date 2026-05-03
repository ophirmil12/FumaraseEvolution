"""
generate_tables.py - Produce Table 1 (Class I) and Table 2 (Class II).

Inputs:
  - proteome_metadata.tsv  (UniProt reference proteomes with Taxonomic lineage)
  - class1_sequences.fasta
  - class2_sequences.fasta

Logic:
  - Total organisms per group  → from metadata (Taxonomic lineage column)
  - Harboring a class          → proteome ID present in that FASTA
  - % harboring / % lacking    → computed from the above

Taxonomic groups:
  All, Bacteria, Archaea, Fungi, Animals, Plants, Protozoal, Algae,
  Alphaproteobacteria, SAR

NOTE - Group overlaps:
  Alphaproteobacteria is a subclass of Bacteria - proteomes in this group
  are counted in BOTH the "Bacteria" row and the "Alphaproteobacteria" row.
  Similarly, SAR is a subset of Eukaryota and overlaps with "Algae".

Usage:
    python scripts/01_search/generate_tables.py
    python scripts/01_search/generate_tables.py \\
        --metadata data/processed/proteome_metadata.tsv \\
        --class1   data/processed/class1_sequences.fasta \\
        --class2   data/processed/class2_sequences.fasta \\
        --outdir   results/stats/

Output:
    results/stats/table1_class1.csv
    results/stats/table2_class2.csv
"""

import sys
import argparse
import logging
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

LIN = "taxonomic_lineage"


# ---------------------------------------------------------------------------
# Taxonomic group definitions
# ---------------------------------------------------------------------------

def define_groups(df: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    """
    Return (label, boolean_mask) pairs for each taxonomic group.
    """
    lin = df[LIN]
    
    # Pre-calculate masks for logic reuse
    mask_bacteria = lin.str.contains(r",\s*Bacteria,", regex=True, na=False)
    mask_alpha    = lin.str.contains("Alphaproteobacteria", na=False)

    return [
        ("All",
            pd.Series(True, index=df.index)),

        ("Other Bacteria",
            mask_bacteria & ~mask_alpha),

        ("Alphaproteobacteria",
            mask_alpha),

        ("Archaea",
            lin.str.contains(r",\s*Archaea,", regex=True, na=False)),

        ("Eukaryota",
            lin.str.contains(r",\s*Eukaryota,", regex=True, na=False)),

        ("Fungi",
            lin.str.contains("Fungi", na=False)),

        ("Animals",
            lin.str.contains("Metazoa", na=False)),

        ("Plants",
            lin.str.contains("Viridiplantae", na=False)),

        ("Protozoal",
            lin.str.contains(
                "Discoba|Metamonada|Amoebozoa|Parabasalia|Fornicata",
                na=False)),

        ("Algae",
            (lin.str.contains(r"\bSar\b", na=False) & ~lin.str.contains("Metazoa", na=False)) |
            lin.str.contains("Rhodophyta|Haptophyta|Cryptophyceae", na=False) |
            (lin.str.contains("Viridiplantae", na=False) & ~lin.str.contains("Embryophyta", na=False))),

        ("SAR",
            lin.str.contains(r"\bSar\b", na=False)),
    ]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_metadata(path: Path) -> pd.DataFrame:
    log.info(f"Loading metadata from {path}")
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = df.columns.str.strip()
    required = {"proteome_id", LIN}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Metadata TSV missing columns: {missing}\n"
            f"Found: {list(df.columns)}"
        )
    df = df.set_index("proteome_id")
    log.info(f"  {len(df):,} proteomes loaded")
    return df


def extract_proteome_ids(fasta_path: Path) -> set[str]:
    """Return the set of proteome IDs present in a FASTA (from header prefix)."""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split("|")[0].strip())
    log.info(f"  {fasta_path.name}: {len(ids):,} unique proteome IDs")
    return ids


def build_table(metadata: pd.DataFrame,
                harboring_ids: set[str],
                fumarase_class: int) -> pd.DataFrame:
    """
    Build the distribution table:
        Taxonomic group | # Proteomes | # Harboring | % Harboring | % Lacking
    """
    rows = []
    for label, mask in define_groups(metadata):
        group     = metadata[mask]
        total     = len(group)
        harboring = int(group.index.isin(harboring_ids).sum())
        lacking   = total - harboring
        pct_h     = round(harboring / total * 100, 2) if total > 0 else 0.0
        pct_l     = round(lacking   / total * 100, 2) if total > 0 else 0.0

        rows.append({
            "Taxonomic group": label,
            "# Proteomes":     total,
            "# Harboring":     harboring,
            "% Harboring":     pct_h,
            "% Lacking":       pct_l,       # In the paper, we only report % harboring, as asked in review
        })

    df_out = pd.DataFrame(rows)
    log.info(
        f"Class {fumarase_class}: "
        f"{df_out['# Harboring'].iloc[0]:,} / {df_out['# Proteomes'].iloc[0]:,} "
        f"proteomes harbor this class"
    )
    return df_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate Table 1 (Class I) and Table 2 (Class II)."
    )
    parser.add_argument(
        "--metadata", type=Path,
        default=Path(PATHS["data_processed"]) / "proteome_metadata.tsv",
    )
    parser.add_argument(
        "--class1", type=Path,
        default=Path(PATHS["data_processed"]) / "class1_sequences.fasta",
    )
    parser.add_argument(
        "--class2", type=Path,
        default=Path(PATHS["data_processed"]) / "class2_sequences.fasta",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["stats"]),
    )
    args = parser.parse_args()

    for p in [args.metadata, args.class1, args.class2]:
        if not p.exists():
            log.error(f"File not found: {p}")
            sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata   = load_metadata(args.metadata)
    class1_ids = extract_proteome_ids(args.class1)
    class2_ids = extract_proteome_ids(args.class2)

    table1 = build_table(metadata, class1_ids, fumarase_class=1)
    table2 = build_table(metadata, class2_ids, fumarase_class=2)

    out1 = args.outdir / "table1_class1.csv"
    out2 = args.outdir / "table2_class2.csv"
    table1.to_csv(out1, index=False)
    table2.to_csv(out2, index=False)

    log.info(f"Table 1 → {out1}")
    log.info(f"Table 2 → {out2}")

    for label, tbl in [("Table 1 - Class I", table1), ("Table 2 - Class II", table2)]:
        print(f"\n{label}")
        print(tbl.to_string(index=False))


if __name__ == "__main__":
    main()
