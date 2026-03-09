"""
clean_sequences.py - Remove fragments and contaminants from fumarase FASTA files.

Operates IN-PLACE on:
    data/processed/class1_sequences.fasta
    data/processed/class2_sequences.fasta

Filters applied:
    Both classes:
        - Remove fragments  (header contains "(Fragment)")
        - Remove sequences shorter than min_length (default: 300 aa)

    Class 1 only:
        - Remove TRZ/ATZ family proteins  (MMseqs2 false positives)

    Class 2 only:
        - Remove aspartate ammonia-lyase / AspA  (paralogues, same superfamily)
        - Remove UPI-only headers  (no UniProt accession, broken format)

A backup is written before any modification:
    data/processed/class1_sequences.fasta.bak
    data/processed/class2_sequences.fasta.bak

Usage:
    python scripts/01_search/clean_sequences.py           # both classes
    python scripts/01_search/clean_sequences.py --class1  # class 1 only
    python scripts/01_search/clean_sequences.py --class2  # class 2 only
    python scripts/01_search/clean_sequences.py --dry-run # report only, no writes

Output:
    Prints a per-class summary of sequences removed and reason.
    Overwrites the input FASTA in-place (after backup).
"""

import sys
import argparse
import logging
import shutil
from pathlib import Path
from collections import Counter

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
# Contaminant / false-positive patterns
# ---------------------------------------------------------------------------

# Applied to both classes
FRAGMENT_TAG = "(Fragment)"

BOTH_EXCLUDE = [
    # EC number placeholder with no gene name annotation
    "E4.2.1.2A",
]

# Class 1 specific exclusions (case-insensitive substring match on header)
CLASS1_EXCLUDE = [
    "TRZ/ATZ family protein",      # triazine hydrolases - wrong enzyme
]

# Class 2 specific exclusions
CLASS2_EXCLUDE = [
    "Aspartate ammonia-lyase",     # AspA paralogues - same superfamily, wrong function
    "aspartate ammonia-lyase",
    " aspA ",                      # gene name in description (space-padded to avoid partial matches)
    "AspA2",
]

MIN_LENGTH = 300   # below this almost certainly a fragment or wrong hit

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header_line, sequence) tuples. Header includes '>'."""
    records = []
    current_header = None
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                    current_seq = []
                current_header = line
            else:
                current_seq.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_seq)))

    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for header, seq in records:
            f.write(header + "\n")
            # wrap at 60 chars
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def is_upi_only(header: str) -> bool:
    """
    Detect headers with no proper UniProt accession.
    Standard format: >UPID|db|accession|...
    Broken format:   >UPID|UPI... status=active
    """
    parts = header.lstrip(">").split("|")
    if len(parts) < 2:
        return True
    db_field = parts[1]
    return db_field.startswith("UPI") or "status=active" in db_field


def filter_records(
    records: list[tuple[str, str]],
    class_key: str,
    min_length: int = MIN_LENGTH,
) -> tuple[list[tuple[str, str]], Counter]:
    """
    Apply all filters for the given class.
    Returns (kept_records, removal_reason_counter).
    """
    kept = []
    removed = Counter()

    exclude_patterns = list(BOTH_EXCLUDE)
    if class_key == "class1":
        exclude_patterns += CLASS1_EXCLUDE
    elif class_key == "class2":
        exclude_patterns += CLASS2_EXCLUDE

    for header, seq in records:
        # Fragment tag
        if FRAGMENT_TAG in header:
            removed["fragment_tag"] += 1
            continue

        # Too short
        seq_clean = seq.replace("-", "").replace(".", "")
        if len(seq_clean) < min_length:
            removed["too_short"] += 1
            continue

        # UPI-only headers (class2 only)
        if class_key == "class2" and is_upi_only(header):
            removed["upi_only_header"] += 1
            continue

        # Contaminant patterns
        matched = False
        for pattern in exclude_patterns:
            if pattern in header:
                removed[f"contaminant:{pattern.strip()[:40]}"] += 1
                matched = True
                break
        if matched:
            continue

        kept.append((header, seq))

    return kept, removed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_class(fasta_path: Path, class_key: str, dry_run: bool, min_length: int) -> None:
    if not fasta_path.exists():
        log.error(f"File not found: {fasta_path}")
        log.error("Rename your FASTA files to class1_sequences.fasta / class2_sequences.fasta")
        sys.exit(1)

    records = parse_fasta(fasta_path)
    n_before = len(records)

    kept, removed = filter_records(records, class_key, min_length)
    n_removed = n_before - len(kept)

    log.info(f"{class_key}: {n_before:,} -> {len(kept):,} sequences ({n_removed} removed)")
    for reason, count in sorted(removed.items(), key=lambda x: -x[1]):
        log.info(f"  {count:>5}  {reason}")

    if dry_run:
        log.info("  [DRY RUN] no files written")
        return

    # Backup
    backup_path = fasta_path.with_suffix(fasta_path.suffix + ".bak")
    shutil.copy2(fasta_path, backup_path)
    log.info(f"  Backup -> {backup_path}")

    # Write in-place
    write_fasta(fasta_path, kept)
    log.info(f"  Written -> {fasta_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Clean fumarase FASTA files in-place."
    )
    parser.add_argument("--class1",   action="store_true")
    parser.add_argument("--class2",   action="store_true")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Report removals without writing any files")
    parser.add_argument("--min-length", type=int, default=MIN_LENGTH,
                        help=f"Minimum sequence length (default: {MIN_LENGTH})")
    args = parser.parse_args()

    # Default: both classes
    if not args.class1 and not args.class2:
        args.class1 = True
        args.class2 = True

    processed = Path(PATHS["data_processed"])

    if args.class1:
        process_class(
            processed / "class1_sequences.fasta",
            "class1",
            args.dry_run,
            args.min_length,
        )

    if args.class2:
        process_class(
            processed / "class2_sequences.fasta",
            "class2",
            args.dry_run,
            args.min_length,
        )

    if not args.dry_run:
        log.info("Done. Run scripts/02_align/run_mafft.sh next.")


if __name__ == "__main__":
    main()
