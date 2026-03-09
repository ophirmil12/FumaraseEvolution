"""
filter_viral.py — Remove viral proteomes from the UniProt reference proteome list.

Downloads (or reads locally) the full UniProt reference proteome metadata and
removes any entry whose taxonomic lineage includes "Viruses". Writes a clean
TSV to data/processed/ ready for mmseq_run.sh.

Usage:
    # Download fresh from UniProt and filter:
    python scripts/01_search/filter_viral.py

    # Use a locally cached TSV instead of downloading:
    python scripts/01_search/filter_viral.py --input data/raw/uniprot_proteomes_raw.tsv

Output:
    data/processed/proteomes_filtered.tsv
"""

import sys
import argparse
import logging
from pathlib import Path

import pandas as pd
import requests

# --- Project root on path ---
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_PROTEOMES_URL = (
    "https://rest.uniprot.org/proteomes/stream"
    "?query=*&format=tsv"
    "&fields=upid,organism,organism_id,protein_count,busco,cpd,lineage"
)

EXPECTED_COLS = {
    "Proteome Id", "Organism", "Organism Id",
    "Protein count", "BUSCO", "CPD", "Lineage",
}

OUTPUT_COLS = [
    "proteome_id", "organism", "organism_id",
    "protein_count", "busco", "cpd",
]

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
    log.info("Downloading UniProt reference proteome list...")
    raw_path.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(UNIPROT_PROTEOMES_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        raw_path.write_bytes(r.content)

    log.info(f"Cached raw proteome list → {raw_path}")
    return pd.read_csv(raw_path, sep="\t", dtype=str)


def load_local(path: Path) -> pd.DataFrame:
    log.info(f"Loading local proteome list from {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def validate_columns(df: pd.DataFrame) -> None:
    missing = EXPECTED_COLS - set(df.columns)
    if missing:
        raise ValueError(
            f"Input TSV is missing expected columns: {missing}\n"
            f"Found: {list(df.columns)}"
        )


def filter_viral(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
    Drop rows whose Lineage contains 'Viruses'.
    Returns (filtered_df, n_removed).
    """
    before = len(df)
    is_viral = df["Lineage"].str.contains("Viruses", case=False, na=False)
    df_clean = df[~is_viral].copy()
    return df_clean, before - len(df_clean)


def rename_and_select(df: pd.DataFrame) -> pd.DataFrame:
    """Rename UniProt column headers to snake_case and drop Lineage."""
    col_map = {
        "Proteome Id":   "proteome_id",
        "Organism":      "organism",
        "Organism Id":   "organism_id",
        "Protein count": "protein_count",
        "BUSCO":         "busco",
        "CPD":           "cpd",
    }
    return df.rename(columns=col_map)[OUTPUT_COLS]


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
        default=Path(PATHS["data_processed"]) / "proteomes_filtered.tsv",
        help="Output TSV path (default: data/processed/proteomes_filtered.tsv)"
    )
    args = parser.parse_args()

    raw_cache = Path(PATHS["data_raw"]) / "uniprot_proteomes_raw.tsv"

    # Load
    if args.input:
        df = load_local(args.input)
    else:
        df = download_proteome_list(raw_cache)

    log.info(f"Loaded {len(df):,} proteomes total")

    # Validate
    validate_columns(df)

    # Filter
    df_clean, n_removed = filter_viral(df)
    log.info(f"Removed {n_removed:,} viral proteomes → {len(df_clean):,} remaining")

    # Rename and write
    df_out = rename_and_select(df_clean)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.out, sep="\t", index=False)
    log.info(f"Written → {args.out}")

    # Also write a metadata-only TSV for filter_hits.py
    metadata_path = Path(PATHS["data_processed"]) / "proteome_metadata.tsv"
    df_out.to_csv(metadata_path, sep="\t", index=False)
    log.info(f"Metadata copy written → {metadata_path}")


if __name__ == "__main__":
    main()