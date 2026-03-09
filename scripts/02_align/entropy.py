"""
entropy.py - Compute per-position Shannon entropy conservation scores from
             MAFFT-aligned Class I and II fumarase sequences.

Shannon entropy H at position i:
    H(i) = -sum( f_a * log2(f_a) ) for each amino acid a with frequency f_a > 0

Lower H = higher conservation (H=0 means perfectly conserved).
Scale: 0 (perfect conservation) to ~4.32 bits (log2(20), maximum variability).

Also computes baseline entropy for a ribosomal protein (highly conserved) and
a membrane transporter (highly divergent) for context, if provided.

Usage:
    # Both classes (default paths from config):
    python scripts/02_align/entropy.py

    # Explicit paths:
    python scripts/02_align/entropy.py \\
        --class1 results/alignments/class1_aligned.fasta \\
        --class2 results/alignments/class2_aligned.fasta

    # With baselines:
    python scripts/02_align/entropy.py \\
        --baseline-conserved  data/external/rplA_aligned.fasta \\
        --baseline-divergent  data/external/membrane_transporter_aligned.fasta

Output:
    results/stats/entropy_class1.csv     per-position scores
    results/stats/entropy_class2.csv     per-position scores
    results/stats/entropy_summary.csv    mean scores per alignment + baselines
"""

import sys
import argparse
import logging
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS, ENTROPY, QUERIES

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Gap characters to exclude from frequency calculation
GAP_CHARS = {"-", ".", "X"}


# ---------------------------------------------------------------------------
# Core entropy calculation
# ---------------------------------------------------------------------------

def parse_aligned_fasta(path: Path) -> tuple[list[str], list[str]]:
    """
    Parse an aligned FASTA into (headers, sequences).
    All sequences must be the same length (alignment columns).
    """
    headers, sequences = [], []
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                headers.append(line[1:])
            else:
                current_seq.append(line.upper())

    if current_seq:
        sequences.append("".join(current_seq))

    if not sequences:
        raise ValueError(f"No sequences found in {path}")

    lengths = set(len(s) for s in sequences)
    if len(lengths) > 1:
        raise ValueError(
            f"Sequences in {path} have inconsistent lengths: {lengths}\n"
            "Is this file actually aligned?"
        )

    log.info(f"  Parsed {len(sequences):,} sequences, alignment length {len(sequences[0]):,}")
    return headers, sequences


def column_entropy(column: str) -> float:
    """
    Compute Shannon entropy (bits) for a single alignment column.
    Gap characters are excluded from frequency calculation.
    Returns 0.0 if the column is all gaps.
    """
    residues = [aa for aa in column if aa not in GAP_CHARS]
    if not residues:
        return 0.0

    counts = Counter(residues)
    total = len(residues)
    entropy = 0.0
    for count in counts.values():
        freq = count / total
        entropy -= freq * np.log2(freq)

    return entropy


def compute_entropy(sequences: list[str]) -> np.ndarray:
    """
    Compute per-position Shannon entropy across all sequences.
    Returns array of shape (alignment_length,).
    """
    aln_len = len(sequences[0])
    scores = np.zeros(aln_len)

    for i in range(aln_len):
        column = "".join(seq[i] for seq in sequences)
        scores[i] = column_entropy(column)

    return scores


# ---------------------------------------------------------------------------
# Functional site mapping
# ---------------------------------------------------------------------------

# Known catalytic residues from PDB annotations (0-indexed in reference sequence)
# Class I  (P0AC33, FumA): [4Fe-4S] cluster coordinating cysteines
# Class II (P05042, FumC): active site residues
FUNCTIONAL_SITES = {
    "class1": {
        "name":     "Class I (P0AC33)",
        "residues": [104, 223, 317],          # Cys105, Cys224, Cys318
    },
    "class2": {
        "name":     "Class II (P05042)",
        "residues": [187, 317, 323, 330],     # His188, Ser318, Lys324, Glu331
        # "substrate_binding": [99, 138, 139, 140, 186, 325]    # the others are the catalytic residues
    },
}


def map_sites_to_alignment(sequences: list[str],
                           ref_header_fragment: str,
                           residue_positions: list[int]) -> list[int]:
    """
    Find the reference sequence in the alignment and map ungapped residue
    positions to alignment column indices.

    Returns list of alignment column indices (one per functional site).
    Returns empty list if reference sequence is not found.
    """
    ref_seq = None
    for seq in sequences:
        # Heuristic: reference sequence has fewest gaps
        if seq.count("-") < len(seq) * 0.05:
            ref_seq = seq
            break

    if ref_seq is None:
        log.warning("Could not identify reference sequence for site mapping.")
        return []

    # Build map: ungapped position → alignment column
    ungapped_to_col = {}
    ungapped_pos = 0
    for col_idx, aa in enumerate(ref_seq):
        if aa not in GAP_CHARS:
            ungapped_to_col[ungapped_pos] = col_idx
            ungapped_pos += 1

    mapped = []
    for res_pos in residue_positions:
        if res_pos in ungapped_to_col:
            mapped.append(ungapped_to_col[res_pos])
        else:
            log.warning(f"Residue position {res_pos} not found in reference sequence.")

    return mapped


# ---------------------------------------------------------------------------
# Per-alignment processing
# ---------------------------------------------------------------------------

def process_alignment(fasta_path: Path,
                      label: str,
                      class_key: str | None = None) -> dict:
    """
    Parse alignment, compute entropy, optionally map functional sites.
    Returns summary dict and per-position DataFrame.
    """
    log.info(f"Processing {label}: {fasta_path.name}")
    headers, sequences = parse_aligned_fasta(fasta_path)
    scores = compute_entropy(sequences)

    mean_entropy = float(np.mean(scores))
    log.info(f"  Mean entropy: {mean_entropy:.4f} bits")

    # Per-position DataFrame
    df = pd.DataFrame({
        "position":  np.arange(1, len(scores) + 1),   # 1-indexed
        "entropy":   scores,
        "conserved": scores < ENTROPY["scale_max"] * 0.25,  # bottom quartile = highly conserved
    })

    # Map functional sites if this is a fumarase class
    site_col_indices = []
    if class_key and class_key in FUNCTIONAL_SITES:
        sites = FUNCTIONAL_SITES[class_key]
        site_col_indices = map_sites_to_alignment(
            sequences, sites["name"], sites["residues"]
        )
        df["is_functional_site"] = df.index.isin(site_col_indices)

        if site_col_indices:
            site_entropy = scores[site_col_indices]
            mean_site_entropy = float(np.mean(site_entropy))
            log.info(
                f"  Functional site mean entropy: {mean_site_entropy:.4f} bits "
                f"({len(site_col_indices)} sites mapped)"
            )
        else:
            mean_site_entropy = None
    else:
        df["is_functional_site"] = False
        mean_site_entropy = None

    summary = {
        "label":               label,
        "n_sequences":         len(sequences),
        "alignment_length":    len(scores),
        "mean_entropy":        round(mean_entropy, 4),
        "mean_site_entropy":   round(mean_site_entropy, 4) if mean_site_entropy else None,
        "n_functional_sites":  len(site_col_indices),
        "pct_highly_conserved": round((scores < 0.5).mean() * 100, 2),
    }

    return summary, df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute per-position Shannon entropy for fumarase alignments."
    )
    parser.add_argument(
        "--class1", type=Path,
        default=Path(PATHS["alignments"]) / "class1_aligned.fasta",
    )
    parser.add_argument(
        "--class2", type=Path,
        default=Path(PATHS["alignments"]) / "class2_aligned.fasta",
    )
    parser.add_argument(
        "--baseline-conserved", type=Path, default=None,
        metavar="FASTA",
        help="Aligned FASTA of a highly conserved reference protein (e.g. RplA)"
    )
    parser.add_argument(
        "--baseline-divergent", type=Path, default=None,
        metavar="FASTA",
        help="Aligned FASTA of a divergent reference protein (e.g. membrane transporter)"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["stats"]),
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    alignments_to_run = [
        (args.class1, "Class I",  "class1"),
        (args.class2, "Class II", "class2"),
    ]
    if args.baseline_conserved:
        alignments_to_run.append((args.baseline_conserved, "Baseline (conserved)", None))
    if args.baseline_divergent:
        alignments_to_run.append((args.baseline_divergent, "Baseline (divergent)", None))

    # Validate inputs
    for path, label, _ in alignments_to_run:
        if not path.exists():
            log.error(f"{label} alignment not found: {path}")
            if "class" in label.lower():
                log.error("Run scripts/02_align/run_mafft.sh first.")
            sys.exit(1)

    summaries = []
    for path, label, class_key in alignments_to_run:
        summary, df_pos = process_alignment(path, label, class_key)
        summaries.append(summary)

        # Write per-position scores
        safe_label = label.lower().replace(" ", "_").replace("(", "").replace(")", "")
        out_path = args.outdir / f"entropy_{safe_label}.csv"
        df_pos.to_csv(out_path, index=False)
        log.info(f"  Per-position scores → {out_path}")

    # Write summary
    df_summary = pd.DataFrame(summaries)
    summary_path = args.outdir / "entropy_summary.csv"
    df_summary.to_csv(summary_path, index=False)
    log.info(f"\nSummary → {summary_path}")

    print("\n" + df_summary.to_string(index=False))


if __name__ == "__main__":
    main()