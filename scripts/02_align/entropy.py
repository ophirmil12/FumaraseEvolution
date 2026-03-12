"""
entropy.py - Compute per-position Shannon entropy conservation scores from
             MAFFT-aligned Class I and II fumarase sequences.

Shannon entropy H at position i:
    H(i) = -sum( f_a * log2(f_a) ) for each amino acid a with frequency f_a > 0

Lower H = higher conservation (H=0 means perfectly conserved).
Scale: 0 (perfect conservation) to ~4.32 bits (log2(20), maximum variability).

Scores are computed over ALL alignment columns, then filtered to only the
columns where the E. coli reference sequence (P0AC33 / P05042) has a residue.
This maps scores back to reference sequence positions (0-548 for Class I,
0-467 for Class II), matching the x-axis in Figure 5.

For baseline proteins (no reference defined), the first sequence in the
alignment is used as the reference for filtering.

Also computes baseline entropy for a ribosomal protein (highly conserved) and
a membrane transporter (highly divergent) for context, if provided.

# RplA (conserved) - reviewed Swiss-Prot bacteria:
https://rest.uniprot.org/uniprotkb/stream?query=gene:rplA+AND+reviewed:true&format=fasta

# MFS transporter (divergent):
https://rest.uniprot.org/uniprotkb/stream?query=family:"major+facilitator+superfamily"+AND+reviewed:true&format=fasta&size=500

Usage:
    # Both classes (default paths from config):
    python scripts/02_align/entropy.py

    # Explicit paths:
    python scripts/02_align/entropy.py \\
        --class1 results/alignments/class1_aligned.fasta \\
        --class2 results/alignments/class2_aligned.fasta

    # With baselines:
    python scripts/02_align/entropy.py \\
        --baseline-conserved  data/external/conservation_baseline/rplA_aligned.fasta \\
        --baseline-divergent  data/external/conservation_baseline/mfs_transporter_aligned.fasta

    # Subsample all groups to at most N sequences:
    python scripts/02_align/entropy.py --max-seqs 2000

    # Reproducible subsampling with custom seed:
    python scripts/02_align/entropy.py --max-seqs 2000 --seed 123

Output:
    results/stats/entropy_class_i.csv          per-position scores (reference positions only)
    results/stats/entropy_class_ii.csv         per-position scores (reference positions only)
    results/stats/entropy_baseline_conserved.csv
    results/stats/entropy_baseline_divergent.csv
    results/stats/entropy_summary.csv          mean scores per alignment + baselines
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
from config import PATHS, ENTROPY, QUERIES, FUNCTIONAL_SITES

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
# Subsampling
# ---------------------------------------------------------------------------

def subsample_sequences(headers: list[str],
                        sequences: list[str],
                        max_seqs: int,
                        seed: int = 42) -> tuple[list[str], list[str]]:
    """
    Randomly subsample sequences to max_seqs. Deterministic given the same seed.
    Returns unchanged lists if len(sequences) <= max_seqs.
    """
    if len(sequences) <= max_seqs:
        return headers, sequences

    rng = np.random.default_rng(seed)
    indices = sorted(rng.choice(len(sequences), size=max_seqs, replace=False))
    log.info(f"  Subsampling {len(sequences):,} -> {max_seqs:,} sequences (seed={seed})")
    return [headers[i] for i in indices], [sequences[i] for i in indices]


# ---------------------------------------------------------------------------
# Reference column filtering
# ---------------------------------------------------------------------------

def get_reference_sequence(headers: list[str],
                            sequences: list[str],
                            uniprot_id: str) -> str | None:
    """
    Find and return the reference sequence in the alignment by UniProt ID.
    Returns None if not found.
    """
    for header, seq in zip(headers, sequences):
        if uniprot_id in header:
            log.info(f"  Reference sequence found: {header[:80]}")
            return seq
    log.warning(f"  Reference sequence {uniprot_id} not found in alignment.")
    return None


def get_reference_columns(ref_seq: str) -> list[int]:
    """
    Return alignment column indices where the reference has a residue (not a gap).
    These are the columns that map to reference sequence positions.
    """
    return [i for i, aa in enumerate(ref_seq) if aa not in GAP_CHARS]


def filter_to_reference(scores: np.ndarray, ref_seq: str) -> np.ndarray:
    """
    Filter full-MSA entropy scores to only the columns where the reference
    sequence has a residue. Returns array of length = ungapped reference length.
    """
    ref_cols = get_reference_columns(ref_seq)
    return scores[ref_cols]


# ---------------------------------------------------------------------------
# Functional site mapping
# ---------------------------------------------------------------------------


def map_sites_to_reference(residue_positions: list[int],
                            ref_seq: str,
                            uniprot_id: str) -> list[int]:
    """
    Map ungapped residue positions (0-indexed) to reference-sequence indices.

    Since scores are already filtered to reference columns, the functional site
    positions ARE the indices into the filtered scores array — we just need to
    verify they fall within the ungapped reference length.

    Returns list of 0-indexed positions into the reference-filtered scores array.
    """
    ref_len = sum(1 for aa in ref_seq if aa not in GAP_CHARS)
    mapped = []
    for pos in residue_positions:
        if pos < ref_len:
            mapped.append(pos)
        else:
            log.warning(
                f"  Residue position {pos} exceeds reference length {ref_len} "
                f"for {uniprot_id} — skipping."
            )
    return mapped


# ---------------------------------------------------------------------------
# Per-alignment processing
# ---------------------------------------------------------------------------

def process_alignment(fasta_path: Path,
                      label: str,
                      class_key: str | None = None,
                      max_seqs: int | None = None,
                      seed: int = 42) -> tuple[dict, pd.DataFrame]:
    """
    Parse alignment, optionally subsample, compute entropy over all columns,
    then filter to reference sequence positions.

    For fumarase classes: reference is the E. coli UniProt sequence (P0AC33/P05042).
      The reference sequence is guaranteed to survive subsampling.
    For baselines: reference is the first sequence in the alignment (before subsampling).

    Returns (summary_dict, per_position_DataFrame).
    """
    log.info(f"Processing {label}: {fasta_path.name}")
    headers, sequences = parse_aligned_fasta(fasta_path)

    # --- Subsample if requested ---
    if max_seqs is not None:
        if class_key and class_key in FUNCTIONAL_SITES:
            # Save reference before subsampling so we can re-insert if dropped
            uniprot_id   = QUERIES[class_key]["uniprot_id"]
            ref_idx      = next((i for i, h in enumerate(headers) if uniprot_id in h), None)
            ref_header   = headers[ref_idx]   if ref_idx is not None else None
            ref_sequence = sequences[ref_idx] if ref_idx is not None else None

            headers, sequences = subsample_sequences(headers, sequences, max_seqs, seed)

            # Re-insert reference if it was dropped by subsampling
            if ref_header is not None and not any(uniprot_id in h for h in headers):
                headers.append(ref_header)
                sequences.append(ref_sequence)
                log.info(f"  Reference {uniprot_id} re-inserted after subsampling")
        else:
            # Baseline: record first sequence before subsampling as reference
            baseline_ref_header   = headers[0]
            baseline_ref_sequence = sequences[0]

            headers, sequences = subsample_sequences(headers, sequences, max_seqs, seed)

            # Re-insert baseline reference if dropped
            if not any(baseline_ref_header in h for h in headers):
                headers.insert(0, baseline_ref_header)
                sequences.insert(0, baseline_ref_sequence)
                log.info("  Baseline reference re-inserted after subsampling")

    # --- Compute entropy over all MSA columns ---
    scores_full = compute_entropy(sequences)
    log.info(f"  Full alignment entropy computed: {len(scores_full):,} columns")

    # --- Find reference sequence and filter to its columns ---
    if class_key and class_key in FUNCTIONAL_SITES:
        uniprot_id = QUERIES[class_key]["uniprot_id"]
        ref_seq = get_reference_sequence(headers, sequences, uniprot_id)

        if ref_seq is not None:
            scores = filter_to_reference(scores_full, ref_seq)
            log.info(
                f"  Filtered to reference columns: {len(scores)} positions "
                f"(= ungapped {uniprot_id} length)"
            )
        else:
            # Should not happen after re-insertion guard above
            scores  = scores_full
            ref_seq = None
            log.warning("  Falling back to full alignment scores.")
    else:
        # Baseline: use first sequence as reference
        ref_seq = sequences[0]
        scores  = filter_to_reference(scores_full, ref_seq)
        log.info(
            f"  Baseline: using first sequence as reference ({headers[0][:60]}), "
            f"{len(scores)} positions"
        )

    mean_entropy = float(np.mean(scores))
    log.info(f"  Mean entropy: {mean_entropy:.4f} bits")

    # --- Per-position DataFrame ---
    df = pd.DataFrame({
        "position":  np.arange(1, len(scores) + 1),   # 1-indexed
        "entropy":   scores,
        "conserved": scores < ENTROPY["scale_max"] * 0.25,  # bottom quartile
    })

    # --- Map functional sites ---
    site_indices = []
    if class_key and class_key in FUNCTIONAL_SITES and ref_seq is not None:
        sites        = FUNCTIONAL_SITES[class_key]
        uniprot_id   = QUERIES[class_key]["uniprot_id"]
        site_indices = map_sites_to_reference(sites["residues"], ref_seq, uniprot_id)
        df["is_functional_site"] = df.index.isin(site_indices)

        if site_indices:
            site_entropy      = scores[site_indices]
            mean_site_entropy = float(np.mean(site_entropy))
            log.info(
                f"  Functional site mean entropy: {mean_site_entropy:.4f} bits "
                f"({len(site_indices)} sites mapped)"
            )
        else:
            mean_site_entropy = None
    else:
        df["is_functional_site"] = False
        mean_site_entropy        = None

    summary = {
        "label":                label,
        "n_sequences":          len(sequences),
        "subsampled":           max_seqs is not None and max_seqs < len(sequences),
        "alignment_length":     len(scores_full),
        "reference_length":     len(scores),
        "mean_entropy":         round(mean_entropy, 4),
        "mean_site_entropy":    round(mean_site_entropy, 4) if mean_site_entropy is not None else None,
        "n_functional_sites":   len(site_indices),
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
        help="Aligned FASTA of a divergent reference protein (e.g. MFS transporter)"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["stats"]),
    )
    parser.add_argument(
        "--max-seqs", type=int, default=None,
        metavar="N",
        help="Subsample each alignment to at most N sequences before computing entropy"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for subsampling (default: 42)"
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

    if args.max_seqs is not None:
        log.info(
            f"Subsampling enabled: max {args.max_seqs:,} sequences per alignment "
            f"(seed={args.seed})"
        )

    summaries = []
    for path, label, class_key in alignments_to_run:
        summary, df_pos = process_alignment(
            path, label, class_key,
            max_seqs=args.max_seqs,
            seed=args.seed,
        )
        summaries.append(summary)

        # Write per-position scores
        safe_label = label.lower().replace(" ", "_").replace("(", "").replace(")", "")
        out_path = args.outdir / f"entropy_{safe_label}.csv"
        df_pos.to_csv(out_path, index=False)
        log.info(f"  Per-position scores -> {out_path}")

    # Write summary
    df_summary   = pd.DataFrame(summaries)
    summary_path = args.outdir / "entropy_summary.csv"
    df_summary.to_csv(summary_path, index=False)
    log.info(f"\nSummary -> {summary_path}")

    print("\n" + df_summary.to_string(index=False))


if __name__ == "__main__":
    main()