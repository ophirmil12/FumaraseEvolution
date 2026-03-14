"""
stats_identity.py - Pairwise sequence identity statistics for Class II fumarase.

Loads the all-vs-all MMseqs2 output (class2_allvsall.tsv), annotates each
sequence pair with taxonomic group metadata, filters to biologically meaningful
cross-domain comparisons, and computes descriptive + inferential statistics.

Identity metrics:
    local_identity  = fident (fraction identical over aligned region)
    global_identity = alnlen * fident / max(qlen, tlen)
                    = identical matches / max sequence length
                    (more conservative, penalises partial alignments)

Comparisons (one sequence from group A paired with one from group B):
    Euk–Alphaproteobacteria   (mitochondrial endosymbiont origin hypothesis)
    Euk–OtherBacteria         (other bacterial donors)
    Euk–Archaea               (archaeal homologs)
    Alpha–OtherBacteria       (within-prokaryote reference)

Statistics per comparison group:
    Descriptive:
        n                     number of pairs
        median_identity_pct   median local identity %
        mean_identity_pct     mean local identity %
        std_identity_pct      standard deviation
        q25, q75              interquartile range
        median_global_pct     median global identity %
        mean_global_pct       mean global identity %
        pct_above_30          % of pairs with local identity > 30%
        pct_above_50          % of pairs with local identity > 50%

    Inferential:
        Kruskal-Wallis        omnibus test across all groups (global)
        Mann-Whitney U        pairwise, BH-corrected

Usage:
    python scripts/04_identity/stats_identity.py

    # Explicit paths:
    python scripts/04_identity/stats_identity.py \\
        --allvsall  data/processed/class2_allvsall.tsv \\
        --metadata  data/processed/proteome_metadata.tsv \\
        --outdir    results/stats/

Output:
    results/stats/identity_descriptive.csv    per-group descriptive stats
    results/stats/identity_pairwise.csv       pairwise Mann-Whitney U + BH correction
    results/stats/identity_kruskal.csv        Kruskal-Wallis omnibus result
    results/stats/identity_pairs_annotated.tsv  full annotated pairs (for plotting)
"""

import sys
import argparse
import logging
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Taxonomic group assignment
# ---------------------------------------------------------------------------

def assign_group(lineage: str) -> str:
    """
    Assign a broad taxonomic group from a UniProt lineage string.
    Returns the most specific relevant group for our comparisons.
    """
    if not isinstance(lineage, str) or not lineage.strip():
        return "Unknown"

    if "Alphaproteobacteria" in lineage:
        return "Alphaproteobacteria"
    if re.search(r",\s*Bacteria,", lineage):
        return "OtherBacteria"
    if re.search(r",\s*Archaea,", lineage):
        return "Archaea"
    if "Eukaryota" in lineage:
        return "Eukaryota"

    return "Unknown"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_allvsall(path: Path) -> pd.DataFrame:
    """
    Load MMseqs2 all-vs-all TSV.
    Columns: query, target, fident, alnlen, qlen, tlen
    Computes global_identity = (fident * alnlen) / max(qlen, tlen)
    """
    log.info(f"Loading all-vs-all from {path} ...")
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["query", "target", "fident", "alnlen", "qlen", "tlen"],
        dtype={"query": str, "target": str, "fident": float,
               "alnlen": int, "qlen": int, "tlen": int},
    )
    log.info(f"  Loaded {len(df):,} pairs")

    # Compute global identity
    df["global_identity"] = (df["fident"] * df["alnlen"]) / df[["qlen", "tlen"]].max(axis=1)

    # Convert to %
    df["local_pct"]  = df["fident"]         * 100
    df["global_pct"] = df["global_identity"] * 100

    return df


def load_metadata(path: Path) -> pd.DataFrame:
    """Load proteome metadata, assign taxonomic groups."""
    import re as _re
    global re
    re = _re

    log.info(f"Loading metadata from {path} ...")
    meta = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    meta.columns = meta.columns.str.strip()
    meta = meta.set_index("proteome_id")
    meta["group"] = meta["taxonomic_lineage"].apply(assign_group)

    counts = meta["group"].value_counts().to_dict()
    log.info(f"  Group counts: {counts}")
    return meta


def extract_proteome_id(seq_id: str) -> str:
    """
    Extract proteome ID from sequence header.
    Expected format: UP000XXXXXX|...  (pipe-delimited, proteome ID first)
    """
    return seq_id.split("|")[0].strip()


# ---------------------------------------------------------------------------
# Pair annotation
# ---------------------------------------------------------------------------

def annotate_pairs(df: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate each pair with the taxonomic group of query and target.
    Drops pairs where either sequence cannot be mapped to metadata.
    """
    log.info("Annotating pairs with taxonomic groups ...")

    df["query_proteome"]  = df["query"].apply(extract_proteome_id)
    df["target_proteome"] = df["target"].apply(extract_proteome_id)

    df["query_group"]  = df["query_proteome"].map(meta["group"])
    df["target_group"] = df["target_proteome"].map(meta["group"])

    before = len(df)
    df = df.dropna(subset=["query_group", "target_group"])
    df = df[df["query_group"] != "Unknown"]
    df = df[df["target_group"] != "Unknown"]
    after = len(df)
    log.info(f"  Dropped {before - after:,} pairs with unmappable/unknown groups")
    log.info(f"  Remaining: {after:,} pairs")

    return df


# ---------------------------------------------------------------------------
# Comparison group filtering
# ---------------------------------------------------------------------------

COMPARISONS = [
    ("Euk-Alpha",         "Eukaryota",         "Alphaproteobacteria"),
    ("Euk-OtherBacteria", "Eukaryota",         "OtherBacteria"),
    ("Euk-Archaea",       "Eukaryota",         "Archaea"),
    ("Alpha-OtherBac",    "Alphaproteobacteria","OtherBacteria"),
]


def extract_comparison_pairs(df: pd.DataFrame) -> dict[str, pd.Series]:
    """
    For each comparison, extract local_pct values for all cross-group pairs.
    A pair (A, B) is included if {query_group, target_group} == {groupA, groupB}.
    Returns dict: label -> Series of local_pct values.
    """
    result = {}
    for label, g1, g2 in COMPARISONS:
        mask = (
            ((df["query_group"] == g1) & (df["target_group"] == g2)) |
            ((df["query_group"] == g2) & (df["target_group"] == g1))
        )
        subset = df[mask]["local_pct"]
        result[label] = subset.reset_index(drop=True)
        log.info(f"  {label}: {len(subset):,} pairs")

    return result


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def descriptive_stats(groups: dict[str, pd.Series]) -> pd.DataFrame:
    """Compute descriptive statistics per comparison group."""
    rows = []
    for label, values in groups.items():
        if values.empty:
            log.warning(f"  {label}: no pairs — skipping")
            continue
        rows.append({
            "group":              label,
            "n":                  len(values),
            "median_identity_pct": round(values.median(), 4),
            "mean_identity_pct":   round(values.mean(),   4),
            "std_identity_pct":    round(values.std(),    4),
            "q25":                 round(values.quantile(0.25), 4),
            "q75":                 round(values.quantile(0.75), 4),
            "iqr":                 round(values.quantile(0.75) - values.quantile(0.25), 4),
            "min_pct":             round(values.min(), 4),
            "max_pct":             round(values.max(), 4),
            "pct_above_30":        round((values > 30).mean() * 100, 2),
            "pct_above_50":        round((values > 50).mean() * 100, 2),
            "pct_above_70":        round((values > 70).mean() * 100, 2),
        })

    return pd.DataFrame(rows)


def kruskal_wallis(groups: dict[str, pd.Series]) -> pd.DataFrame:
    """Run Kruskal-Wallis omnibus test across all comparison groups."""
    valid = {k: v for k, v in groups.items() if not v.empty}
    if len(valid) < 2:
        log.warning("  Fewer than 2 groups for Kruskal-Wallis — skipping")
        return pd.DataFrame()

    stat, p = stats.kruskal(*valid.values())
    log.info(f"  Kruskal-Wallis: H={stat:.4f}, p={p:.4e}")

    return pd.DataFrame([{
        "test":    "Kruskal-Wallis",
        "groups":  ", ".join(valid.keys()),
        "n_groups": len(valid),
        "H_stat":  round(stat, 4),
        "p_value": p,
        "significant_0.05": p < 0.05,
    }])


def pairwise_mannwhitney(groups: dict[str, pd.Series]) -> pd.DataFrame:
    """
    Run pairwise Mann-Whitney U tests with Benjamini-Hochberg FDR correction.
    """
    valid = {k: v for k, v in groups.items() if not v.empty}
    pairs = list(combinations(valid.keys(), 2))

    rows = []
    for g1, g2 in pairs:
        u_stat, p_raw = stats.mannwhitneyu(
            valid[g1], valid[g2], alternative="two-sided"
        )
        # Effect size: rank-biserial correlation
        n1, n2 = len(valid[g1]), len(valid[g2])
        r = 1 - (2 * u_stat) / (n1 * n2)

        rows.append({
            "group_1":    g1,
            "group_2":    g2,
            "n_1":        n1,
            "n_2":        n2,
            "U_stat":     round(u_stat, 2),
            "p_raw":      p_raw,
            "effect_r":   round(r, 4),   # rank-biserial: -1 to 1
        })

    df = pd.DataFrame(rows)

    # BH correction
    _, p_adj, _, _ = multipletests(df["p_raw"], method="fdr_bh")
    df["p_adj_BH"]          = p_adj
    df["significant_0.05"]  = p_adj < 0.05
    df["significant_0.01"]  = p_adj < 0.01

    # Sort by adjusted p-value
    df = df.sort_values("p_adj_BH").reset_index(drop=True)

    log.info(f"  {len(df)} pairwise tests, "
             f"{df['significant_0.05'].sum()} significant at FDR<0.05")

    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Sequence identity statistics for Class II fumarase."
    )
    parser.add_argument(
        "--allvsall", type=Path,
        default=Path(PATHS["data_processed"]) / "class2_allvsall.tsv",
    )
    parser.add_argument(
        "--metadata", type=Path,
        default=Path(PATHS["data_processed"]) / "proteome_metadata.tsv",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["stats"]),
    )
    args = parser.parse_args()

    for p in [args.allvsall, args.metadata]:
        if not p.exists():
            log.error(f"File not found: {p}")
            sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load
    df   = load_allvsall(args.allvsall)
    meta = load_metadata(args.metadata)

    # Annotate
    df = annotate_pairs(df, meta)

    # Save annotated pairs for downstream plotting
    annotated_path = args.outdir / "identity_pairs_annotated.tsv"
    df.to_csv(annotated_path, sep="\t", index=False)
    log.info(f"Annotated pairs -> {annotated_path}")

    # Extract comparison groups
    groups = extract_comparison_pairs(df)

    # Descriptive stats
    df_desc = descriptive_stats(groups)
    desc_path = args.outdir / "identity_descriptive.csv"
    df_desc.to_csv(desc_path, index=False)
    log.info(f"Descriptive stats -> {desc_path}")
    print("\n--- Descriptive Statistics ---")
    print(df_desc.to_string(index=False))

    # Kruskal-Wallis
    df_kw = kruskal_wallis(groups)
    kw_path = args.outdir / "identity_kruskal.csv"
    df_kw.to_csv(kw_path, index=False)
    log.info(f"Kruskal-Wallis -> {kw_path}")
    print("\n--- Kruskal-Wallis ---")
    print(df_kw.to_string(index=False))

    # Pairwise Mann-Whitney U + BH
    df_pw = pairwise_mannwhitney(groups)
    pw_path = args.outdir / "identity_pairwise.csv"
    df_pw.to_csv(pw_path, index=False)
    log.info(f"Pairwise tests -> {pw_path}")
    print("\n--- Pairwise Mann-Whitney U (BH corrected) ---")
    print(df_pw.to_string(index=False))


if __name__ == "__main__":
    import re
    main()