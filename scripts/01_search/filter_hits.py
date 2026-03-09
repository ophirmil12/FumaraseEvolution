"""
filter_hits.py - Filter MMseqs2 .m8 hits and merge into a single JSON file.

Applies thresholds from config.py (e-value, identity, coverage) across all
.m8 files in data/raw/m8s/, cross-references organism metadata, and writes
a single merged JSON to data/processed/filtered_hits.json.

Usage:
    python scripts/01_search/filter_hits.py --class1
    python scripts/01_search/filter_hits.py --class2
    python scripts/01_search/filter_hits.py --class1 --class2   # both

Output JSON structure (one entry per passing hit):
    [
        {
            "proteome_id":   "UP000000625",
            "query_id":      "P0AC33",
            "fumarase_class": 1,
            "hit_id":        "sp|P0AC33|FUMA_ECOLI",
            "identity":      0.98,
            "alignment_len": 540,
            "mismatches":    10,
            "gap_opens":     2,
            "q_start":       1,   "q_end":   548,
            "s_start":       1,   "s_end":   548,
            "evalue":        1e-200,
            "bitscore":      980.0,
            "coverage":      0.985,
            "organism":      "Escherichia coli K-12",
            "organism_id":   "83333",
            "kingdom":       "Bacteria",
            "phylum":        "Proteobacteria"
        },
        ...
    ]
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd

# --- Project root on path ---
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS, MMSEQS2, QUERIES

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Standard MMseqs2 .m8 column names (BLAST tabular format)
M8_COLS = [
    "query_id", "hit_id",
    "identity", "alignment_len", "mismatches", "gap_opens",
    "q_start", "q_end", "s_start", "s_end",
    "evalue", "bitscore",
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

def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """
    Load organism metadata TSV produced by the search step.
    Expected columns: proteome_id, organism, organism_id, kingdom, phylum.
    """
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Metadata file not found: {metadata_path}\n"
            "Run the MMseqs search step first to generate organism metadata."
        )
    df = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"proteome_id", "organism", "organism_id", "kingdom", "phylum"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Metadata TSV is missing columns: {missing}")
    return df.set_index("proteome_id")


def parse_m8(m8_path: Path) -> pd.DataFrame:
    """Read a single .m8 file into a DataFrame. Returns empty DF if file is empty."""
    if m8_path.stat().st_size == 0:
        return pd.DataFrame(columns=M8_COLS)
    return pd.read_csv(m8_path, sep="\t", header=None, names=M8_COLS)


def compute_coverage(df: pd.DataFrame, query_length: int) -> pd.Series:
    """Coverage = alignment length / query length."""
    return df["alignment_len"] / query_length


def apply_filters(df: pd.DataFrame, query_length: int) -> pd.DataFrame:
    """Apply e-value, identity, and coverage filters from config."""
    df = df.copy()
    df["identity"] = df["identity"] / 100.0          # MMseqs2 outputs 0–100
    df["coverage"] = compute_coverage(df, query_length)

    before = len(df)
    df = df[
        (df["evalue"]   <  MMSEQS2["evalue_cutoff"])   &
        (df["identity"] >= MMSEQS2["identity_cutoff"]) &
        (df["coverage"] >= MMSEQS2["coverage_cutoff"])
    ]
    after = len(df)
    log.info(f"  Filter: {before} → {after} hits "
             f"(removed {before - after})")
    return df


def flag_manual_review(hits_by_proteome: dict) -> set:
    """Return proteome IDs with <= min_hits_for_auto hits (need manual review)."""
    threshold = MMSEQS2["min_hits_for_auto"]
    return {
        pid for pid, hits in hits_by_proteome.items()
        if len(hits) <= threshold
    }


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def process_class(fumarase_class: int, metadata: pd.DataFrame) -> list[dict]:
    """
    Process all .m8 files for a given fumarase class.
    Returns a list of hit dicts ready for JSON serialisation.
    """
    query_key  = f"class{fumarase_class}"
    query_id   = QUERIES[query_key]["uniprot_id"]
    query_len  = QUERIES[query_key]["length"]
    m8_dir     = Path(PATHS["data_raw"]) / "m8s"

    m8_files = sorted(m8_dir.glob("*_output.m8"))
    if not m8_files:
        log.warning(f"No .m8 files found in {m8_dir}")
        return []

    log.info(f"Class {fumarase_class}: processing {len(m8_files)} .m8 files "
             f"(query={query_id}, length={query_len})")

    hits_by_proteome: dict[str, list] = defaultdict(list)
    total_raw = 0

    for m8_path in m8_files:
        proteome_id = m8_path.name.replace("_output.m8", "")
        df = parse_m8(m8_path)
        total_raw += len(df)

        if df.empty:
            continue

        df = apply_filters(df, query_len)

        if df.empty:
            continue

        # Attach metadata
        meta = metadata.loc[proteome_id] if proteome_id in metadata.index else {}

        for _, row in df.iterrows():
            hits_by_proteome[proteome_id].append({
                "proteome_id":    proteome_id,
                "query_id":       query_id,
                "fumarase_class": fumarase_class,
                "hit_id":         row["hit_id"],
                "identity":       round(float(row["identity"]), 4),
                "alignment_len":  int(row["alignment_len"]),
                "mismatches":     int(row["mismatches"]),
                "gap_opens":      int(row["gap_opens"]),
                "q_start":        int(row["q_start"]),
                "q_end":          int(row["q_end"]),
                "s_start":        int(row["s_start"]),
                "s_end":          int(row["s_end"]),
                "evalue":         float(row["evalue"]),
                "bitscore":       float(row["bitscore"]),
                "coverage":       round(float(row["coverage"]), 4),
                "organism":       str(meta.get("organism",    "NA")),
                "organism_id":    str(meta.get("organism_id", "NA")),
                "kingdom":        str(meta.get("kingdom",     "NA")),
                "phylum":         str(meta.get("phylum",      "NA")),
            })

    # Flatten to list
    all_hits = [hit for hits in hits_by_proteome.values() for hit in hits]

    # Flag low-hit proteomes
    manual = flag_manual_review(hits_by_proteome)
    if manual:
        log.warning(
            f"  {len(manual)} proteomes have <= {MMSEQS2['min_hits_for_auto']} hits "
            f"and require manual verification: {sorted(manual)}"
        )

    log.info(
        f"Class {fumarase_class} summary: "
        f"{total_raw} raw hits → {len(all_hits)} passing hits "
        f"across {len(hits_by_proteome)} proteomes"
    )
    return all_hits


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter MMseqs2 .m8 hits and write merged JSON."
    )
    parser.add_argument("--class1", action="store_true", help="Process Class I fumarase")
    parser.add_argument("--class2", action="store_true", help="Process Class II fumarase")
    parser.add_argument(
        "--metadata", type=Path,
        default=Path(PATHS["data_processed"]) / "proteome_metadata.tsv",
        help="Path to organism metadata TSV (default: data/processed/proteome_metadata.tsv)"
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path(PATHS["data_processed"]) / "filtered_hits.json",
        help="Output JSON path (default: data/processed/filtered_hits.json)"
    )
    args = parser.parse_args()

    if not args.class1 and not args.class2:
        parser.error("Specify at least one of --class1, --class2")

    log.info(f"Loading metadata from {args.metadata}")
    metadata = load_metadata(args.metadata)

    all_hits = []
    if args.class1:
        all_hits.extend(process_class(1, metadata))
    if args.class2:
        all_hits.extend(process_class(2, metadata))

    # Write output
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(all_hits, f, indent=2)

    log.info(f"Written {len(all_hits)} total hits → {args.out}")


if __name__ == "__main__":
    main()