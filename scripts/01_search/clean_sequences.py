"""
clean_sequences.py - Remove fragments and contaminants from fumarase FASTA files.

Operates IN-PLACE on:
    data/processed/class1_sequences.fasta
    data/processed/class2_sequences.fasta

Filters applied to BOTH classes:
    1. Fragment tag       - header contains "(Fragment)"
    2. Too short         - sequence < MIN_LENGTH aa (default 100)
    3. Domain-only hits  - N-terminal or C-terminal domain fragments annotated
                           as such ("N-terminal domain-containing", "C-terminal
                           domain-containing") but lacking a "(Fragment)" tag

Class I only:
    4. TRZ/ATZ family   - triazine hydrolases, same Fe-S superfamily, wrong enzyme
    5. Alpha subunit     - tartrate dehydratase alpha-type catalytic domain;
                           these are the small alpha subunit of TtdAB (not FumA)
    6. EC 4.2.1.2A      - placeholder annotation with no gene name

Class II only:
    7. AspA paralogues   - Aspartate ammonia-lyase / AspA; same superfamily,
                           different reaction (EC 4.3.1.1 vs EC 4.2.1.2)
    8. Completely unrelated hits (INO80, arginase, succinate dehydrogenase)
    9. UPI-only headers  - no UniProt accession, broken format from extraction

A .bak backup is written before any modification.

Usage:
    python scripts/01_search/clean_sequences.py           # both classes
    python scripts/01_search/clean_sequences.py --class1  # Class I only
    python scripts/01_search/clean_sequences.py --class2  # Class II only
    python scripts/01_search/clean_sequences.py --dry-run # report only, no writes
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
# Filter definitions
# ---------------------------------------------------------------------------

MIN_LENGTH = 100   # aa; below this almost certainly a fragment or domain-only hit

# Applied to both classes
SHARED_EXCLUDE = [
    # Domain-only annotations that lack a (Fragment) tag — these are partial
    # sequences that would distort alignment and entropy calculations.
    # Full fumarases are never annotated as "N-terminal domain-containing" or
    # "C-terminal domain-containing" in UniProt.
    "N-terminal domain-containing protein",
    "C-terminal domain-containing protein",
    # EC number placeholder with no gene-name annotation
    "E4.2.1.2A",
]

# Class I — Fe-S type (FumA/FumB in bacteria)
CLASS1_EXCLUDE = [
    # Triazine hydrolases (TrzA/AtzA family): same Fe-S hydro-lyase superfamily
    # as Class I fumarase but catalyse ring-opening of triazines (EC 3.5.99.3).
    # MMseqs2 picks them up due to structural similarity in the Fe-S domain.
    "TRZ/ATZ family protein",
    # Tartrate dehydratase alpha subunit (TtdA): this is the small, iron-free
    # alpha subunit of the heterodimeric TtdAB complex. The beta subunit (TtdB)
    # is the Class I fumarase homologue; the alpha subunit is not.
    "tartrate dehydratase alpha-type",
    "tartrate dehydratase alpha type",
]

# Class II — FumC-type (eukaryotes + most bacteria)
CLASS2_EXCLUDE = [
    # Aspartate ammonia-lyase (AspA / EC 4.3.1.1): the closest paralogue to
    # Class II fumarase. Same (beta/alpha)8 barrel, ~25% identity to FumC.
    # Distinguished by a diagnostic loop that positions the amine-abstracting
    # Lys instead of the fumarase Ser318. MMseqs2 cannot distinguish them at
    # 40% identity cutoff.
    "Aspartate ammonia-lyase",
    "aspartate ammonia-lyase",
    "Putative aspartate ammonia-lyase",
    "L-Aspartase-like protein",
    "L-aspartase-like protein",
    " aspA ",           # gene name embedded in description (space-padded)
    "AspA2",
    # Completely unrelated proteins — false positives from low-complexity
    # regions or mis-annotated proteomes
    "Chromatin-remodeling ATPase INO80",
    "arginase",         # manganese metalloenzyme, unrelated fold
    "Succinate dehydrogenase",  # FAD-dependent, unrelated fold
]

# ---------------------------------------------------------------------------
# Parsing / writing
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header_line, sequence) tuples. Header retains '>'."""
    records = []
    current_header = None
    current_seq: list[str] = []

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
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def is_upi_only(header: str) -> bool:
    """
    Detect headers with no proper UniProt accession.
    Standard: >UPID|sp|ACCESSION|...  or  >UPID|tr|ACCESSION|...
    Broken:   >UPID|UPI...  or  >UPID|... status=active
    """
    parts = header.lstrip(">").split("|")
    if len(parts) < 2:
        return True
    return parts[1].startswith("UPI") or "status=active" in parts[1]


def filter_records(
    records: list[tuple[str, str]],
    class_key: str,
    min_length: int,
) -> tuple[list[tuple[str, str]], Counter]:
    """
    Apply all filters for the given class_key ("class1" or "class2").
    Returns (kept_records, removal_counter).
    """
    kept: list[tuple[str, str]] = []
    removed: Counter = Counter()

    # Build exclusion list for this class
    exclude = list(SHARED_EXCLUDE)
    if class_key == "class1":
        exclude += CLASS1_EXCLUDE
    elif class_key == "class2":
        exclude += CLASS2_EXCLUDE

    for header, seq in records:
        seq_clean = seq.replace("-", "").replace(".", "")

        # 1. Fragment tag in header
        if "(Fragment)" in header:
            removed["(Fragment) tag"] += 1
            continue

        # 2. Sequence too short
        if len(seq_clean) < min_length:
            removed[f"too short (<{min_length} aa)"] += 1
            continue

        # 3. UPI-only header (class2 only — class1 headers already validated
        #    by extract_sequences.py; class2 may still carry stragglers)
        if class_key == "class2" and is_upi_only(header):
            removed["UPI-only header"] += 1
            continue

        # 4–9. Exclusion patterns
        matched_pattern = None
        for pattern in exclude:
            if pattern in header:
                matched_pattern = pattern.strip()[:50]
                break
        if matched_pattern:
            removed[f"contaminant: {matched_pattern}"] += 1
            continue

        kept.append((header, seq))

    return kept, removed


# ---------------------------------------------------------------------------
# Per-class driver
# ---------------------------------------------------------------------------

def process_class(
    fasta_path: Path,
    class_key: str,
    dry_run: bool,
    min_length: int,
) -> None:
    if not fasta_path.exists():
        log.error(f"File not found: {fasta_path}")
        sys.exit(1)

    records  = parse_fasta(fasta_path)
    n_before = len(records)

    kept, removed = filter_records(records, class_key, min_length)
    n_removed = n_before - len(kept)

    label = "Class I" if class_key == "class1" else "Class II"
    log.info(f"[{label}] {n_before:,} -> {len(kept):,} ({n_removed} removed)")
    for reason, count in sorted(removed.items(), key=lambda x: -x[1]):
        log.info(f"    {count:>5}  {reason}")

    if dry_run:
        log.info(f"[{label}] dry-run: no files written")
        return

    backup = fasta_path.with_suffix(fasta_path.suffix + ".bak")
    shutil.copy2(fasta_path, backup)
    log.info(f"[{label}] backup  -> {backup.name}")

    write_fasta(fasta_path, kept)
    log.info(f"[{label}] written -> {fasta_path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Clean fumarase FASTA files in-place."
    )
    parser.add_argument("--class1",     action="store_true")
    parser.add_argument("--class2",     action="store_true")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Report removals without writing any files")
    parser.add_argument("--min-length", type=int, default=MIN_LENGTH,
                        help=f"Minimum sequence length in aa (default: {MIN_LENGTH})")
    args = parser.parse_args()

    if not args.class1 and not args.class2:
        args.class1 = True
        args.class2 = True

    processed = Path(PATHS["data_processed"])

    if args.class1:
        process_class(processed / "class1_sequences.fasta", "class1",
                      args.dry_run, args.min_length)
    if args.class2:
        process_class(processed / "class2_sequences.fasta", "class2",
                      args.dry_run, args.min_length)

    if not args.dry_run:
        log.info("Done. Run scripts/02_align/run_mafft.sh next.")


if __name__ == "__main__":
    main()