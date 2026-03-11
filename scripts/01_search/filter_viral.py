"""
filter_viral.py - Remove viral proteomes from the UniProt reference proteome list.

Downloads (or reads locally) the full UniProt reference proteome metadata and
removes any entry whose taxonomic lineage includes "Viruses". Writes a clean
TSV to data/processed/ ready for mmseq_run.sh.

Usage:
    # Download fresh from UniProt and filter:
    python scripts/01_search/filter_viral.py

    # Use a locally cached TSV instead of downloading:
    python scripts/01_search/filter_viral.py --input data/raw/uniprot_proteomes_raw.tsv

Output:
    data/processed/proteome_metadata.tsv   (viral entries removed, all columns kept)
"""

import sys
import argparse
import logging
from pathlib import Path

import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_PROTEOMES_URL = (
    "https://rest.uniprot.org/proteomes/stream"
    "?query=proteome_type%3A1&format=tsv"
    "&fields=upid,organism,organism_id,protein_count,busco,cpd,lineage"
)

# Columns that must be present in the downloaded/cached TSV.
# "Taxonomic lineage" is used for viral filtering AND kept in the output
# so that downstream scripts (filter_hits.py, generate_tables.py) can
# assign taxonomy without a separate lookup.
REQUIRED_COLS = {
    "Proteome Id",
    "Organism",
    "Organism Id",
    "Protein count",
    "BUSCO",
    "CPD",
    "Taxonomic lineage",
}

# Snake-case rename map applied before writing.
# All columns are renamed and written; none are dropped.
RENAME = {
    "Proteome Id":       "proteome_id",
    "Organism":          "organism",
    "Organism Id":       "organism_id",
    "Protein count":     "protein_count",
    "BUSCO":             "busco",
    "CPD":               "cpd",
    "Taxonomic lineage": "taxonomic_lineage",
    "Taxon mnemonic":    "taxon_mnemonic",   # present in real downloads; kept if available
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def download_proteome_list(raw_path: Path) -> pd.DataFrame:
    """Fetch full proteome list from UniProt and cache locally."""
    log.info("Downloading UniProt reference proteome list ...")
    raw_path.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(UNIPROT_PROTEOMES_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        raw_path.write_bytes(r.content)

    log.info(f"Cached raw proteome list -> {raw_path}")
    return pd.read_csv(raw_path, sep="\t", dtype=str)


def load_local(path: Path) -> pd.DataFrame:
    log.info(f"Loading local proteome list from {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def validate_columns(df: pd.DataFrame) -> None:
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(
            f"Input TSV is missing required columns: {missing}\n"
            f"Found: {list(df.columns)}"
        )


def filter_viral(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
    Drop rows where 'Taxonomic lineage' contains 'Viruses' or 'Organism' is empty/NaN.
    Returns (filtered_df, n_removed).
    """
    # Create masks for criteria we want to EXCLUDE
    is_viral = df["Taxonomic lineage"].str.contains("Viruses", case=False, na=False)
    is_empty_org = df["Organism"].isna() | (df["Organism"] == "")

    # Keep only the rows that satisfy both conditions
    df_clean = df[~(is_viral | is_empty_org)].copy()

    n_removed = len(df) - len(df_clean)
    return df_clean, n_removed

def rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename UniProt column headers to snake_case.
    Unknown extra columns are kept as-is.
    """
    return df.rename(columns=RENAME)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Remove viral proteomes from the UniProt reference proteome list."
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help="Path to a locally cached raw proteome TSV. "
             "If omitted, downloads fresh from UniProt."
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(PATHS["data_processed"]) / "proteome_metadata.tsv",
        help="Output TSV path (default: data/processed/proteome_metadata.tsv)"
    )
    args = parser.parse_args()

    raw_cache = Path(PATHS["data_raw"]) / "uniprot_proteomes_raw.tsv"

    # Load
    if args.input:
        df = load_local(args.input)
    else:
        df = download_proteome_list(raw_cache)

    log.info(f"Loaded {len(df):,} proteomes")

    # Validate
    validate_columns(df)

    # Filter viral
    df_clean, n_removed = filter_viral(df)
    log.info(f"Removed {n_removed:,} viral proteomes -> {len(df_clean):,} remaining")

    # Rename and write — all columns retained
    df_out = rename_columns(df_clean)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.out, sep="\t", index=False)
    log.info(f"Written -> {args.out}")
    log.info(f"Columns: {list(df_out.columns)}")


if __name__ == "__main__":
    main()