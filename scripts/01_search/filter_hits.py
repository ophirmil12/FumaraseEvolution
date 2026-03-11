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

# Extended M8 column names — col 13 is query_length, present in our pipeline output
M8_COLS_EXT = M8_COLS + ["query_length"]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_lineage(lineage: str) -> tuple[str, str]:
    """
    Derive kingdom and phylum from a UniProt taxonomic lineage string.

    Lineage format (comma-separated, most general to most specific):
        "cellular organisms, Bacteria, Pseudomonadota, Gammaproteobacteria, ..."
        "cellular organisms, Eukaryota, Opisthokonta, Fungi, ..."
        "cellular organisms, Archaea, TACK group, Thermoproteota, ..."

    Returns (kingdom, phylum); "Unknown" if not determinable.
    """
    if not isinstance(lineage, str) or not lineage.strip():
        return "Unknown", "Unknown"

    parts = [p.strip() for p in lineage.split(",")]

    # Intermediate grouping nodes that are not informative phyla
    _SKIP = {
        "Opisthokonta", "Archaeplastida", "Discoba", "Metamonada",
        "Sar", "SAR", "Alveolata", "Rhizaria", "Stramenopiles",
        "Amoebozoa", "Apusozoa", "Cryptophyceae", "Haptophyta",
        "TACK group", "DPANN group", "Asgard group", "FCB group",
        "PVC group", "Terrabacteria group", "Acidobacteriota group",
    }

    for kingdom_name in ("Eukaryota", "Bacteria", "Archaea", "Viruses"):
        if kingdom_name in parts:
            idx = parts.index(kingdom_name)
            phylum = "Unknown"
            for p in parts[idx + 1:]:
                if p and p not in _SKIP and "group" not in p.lower():
                    phylum = p
                    break
            return kingdom_name, phylum

    return "Unknown", "Unknown"


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """
    Load organism metadata TSV (produced by filter_viral.py) and derive
    kingdom and phylum from the taxonomic_lineage column.

    Required columns: proteome_id, organism, organism_id, taxonomic_lineage.
    Adds derived columns: kingdom, phylum.
    """
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Metadata file not found: {metadata_path}\n"
            "Run filter_viral.py first to generate proteome_metadata.tsv."
        )
    df = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"proteome_id", "organism", "organism_id", "taxonomic_lineage"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Metadata TSV is missing columns: {missing}\n"
            "Re-run filter_viral.py — the output must include taxonomic_lineage."
        )

    parsed        = df["taxonomic_lineage"].apply(parse_lineage)
    df["kingdom"] = parsed.apply(lambda x: x[0])
    df["phylum"]  = parsed.apply(lambda x: x[1])

    log.info(
        f"Metadata: {len(df):,} proteomes | "
        f"kingdom counts: {df['kingdom'].value_counts().to_dict()}"
    )
    return df.set_index("proteome_id")


def parse_m8(m8_path: Path) -> pd.DataFrame:
    """
    Read a single .m8 file into a DataFrame.

    Column layout produced by mmseq_run.sh:
        query_id   : proteome protein accession, prefixed "class1_" or "class2_"
        hit_id     : fumarase reference accession (P0AC33 or P05042)
        cols 3-12  : standard BLAST tabular fields
        col 13     : query_length (aa length of the proteome protein)

    Identity values are already in 0-1 scale.
    An optional text header row (starting with "Query_ID") is skipped if present.
    Returns empty DataFrame if file is empty or unreadable.
    """
    if m8_path.stat().st_size == 0:
        return pd.DataFrame(columns=M8_COLS_EXT)

    with open(m8_path) as fh:
        first = fh.readline()
    skip = 1 if first.lower().startswith("query") else 0

    df = pd.read_csv(m8_path, sep="\t", header=None, skiprows=skip)

    if df.shape[1] >= 13:
        df = df.iloc[:, :13]
        df.columns = M8_COLS_EXT
    else:
        df = df.iloc[:, :12]
        df.columns = M8_COLS
        df["query_length"] = None   # coverage will be skipped if missing

    # Coerce numeric columns
    for col in ("identity", "evalue", "bitscore", "alignment_len", "query_length"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def apply_filters(df: pd.DataFrame, class_prefix: str) -> pd.DataFrame:
    """
    Filter rows for a given class and apply quality thresholds.

    - Keeps only rows whose query_id starts with class_prefix ("class1_" / "class2_")
    - Identity is already 0-1 (no division needed)
    - Coverage = alignment_len / query_length  (using per-row query_length from col 13)
    - Applies evalue, identity, coverage cutoffs from config
    - Returns one row per unique hit_id (best bitscore), so each proteome protein
      appears at most once in the output
    """
    df = df[df["query_id"].str.startswith(class_prefix)].copy()
    if df.empty:
        return df

    # Coverage per row using the proteome protein's own length
    if df["query_length"].notna().all():
        df["coverage"] = df["alignment_len"] / df["query_length"]
    else:
        df["coverage"] = 0.0   # can't compute — will be filtered out

    before = len(df)
    df = df[
        (df["evalue"]   <  MMSEQS2["evalue_cutoff"])   &
        (df["identity"] >= MMSEQS2["identity_cutoff"]) &
        (df["coverage"] >= MMSEQS2["coverage_cutoff"])
    ]

    # Keep best hit per proteome protein (highest bitscore)
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("query_id")

    after = len(df)
    log.debug(f"  {class_prefix}: {before} -> {after} passing hits")
    return df


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def process_class(fumarase_class: int, metadata: pd.DataFrame) -> list[dict]:
    """
    Process all .m8 files for a given fumarase class.

    Each .m8 file corresponds to one proteome (named UP000XXXXX_output.m8)
    and contains one row per proteome protein that matched the fumarase query.
    We keep the single best-scoring protein per proteome as the representative hit.

    Returns a list of hit dicts ready for JSON serialisation.
    """
    class_prefix = f"class{fumarase_class}_"
    query_id     = QUERIES[f"class{fumarase_class}"]["uniprot_id"]
    m8_dir       = Path(PATHS["data_raw"]) / "m8s"

    m8_files = sorted(m8_dir.glob("*_output.m8"))
    if not m8_files:
        log.warning(f"No .m8 files found in {m8_dir}")
        return []

    log.info(f"Class {fumarase_class}: scanning {len(m8_files)} .m8 files "
             f"for prefix '{class_prefix}' (reference={query_id})")

    all_hits   = []
    n_files_hit = 0
    n_raw_total = 0

    for m8_path in m8_files:
        proteome_id = m8_path.name.replace("_output.m8", "")
        df = parse_m8(m8_path)

        if df.empty:
            continue

        n_raw_total += df["query_id"].str.startswith(class_prefix).sum()
        df = apply_filters(df, class_prefix)

        if df.empty:
            continue

        n_files_hit += 1
        meta = metadata.loc[proteome_id] if proteome_id in metadata.index else {}

        # One representative hit per proteome = top bitscore row
        row = df.sort_values("bitscore", ascending=False).iloc[0]

        all_hits.append({
            "proteome_id":    proteome_id,
            "query_id":       query_id,
            "fumarase_class": fumarase_class,
            "hit_id":         str(row["hit_id"]),
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

    log.info(
        f"Class {fumarase_class}: {n_raw_total:,} raw hits across all files "
        f"-> {len(all_hits):,} proteomes with a passing representative hit"
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

    log.info(f"Written {len(all_hits)} total hits -> {args.out}")


if __name__ == "__main__":
    main()