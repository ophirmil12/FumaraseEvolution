"""
extract_sequences.py - Validate and register the final representative FASTA files.

At this point in the pipeline the representative sequences have already been
produced by MMseqs2 and filtered. This script:

  1. Validates both FASTA files (header format, sequence content, duplicates)
  2. Parses and logs a summary per taxonomic kingdom (from the header)
  3. Copies them into data/processed/ under canonical names:
       class1_sequences.fasta
       class2_sequences.fasta

These canonical files are the single input source for all downstream steps
(02_align, 03_phylogeny, 04_embeddings).

Usage:
    python scripts/01_search/extract_sequences.py \\
        --class1 path/to/class1.fasta \\
        --class2 path/to/class2.fasta

    # Validate only, no copy:
    python scripts/01_search/extract_sequences.py \\
        --class1 path/to/class1.fasta --class2 path/to/class2.fasta --dry-run
"""

import sys
import re
import shutil
import argparse
import logging
from pathlib import Path
from collections import Counter

# --- Project root on path ---
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Expected header format: >PROTEOME_ID|db|UNIPROT_ID|GENE_NAME description
HEADER_RE = re.compile(r"^>(UP\d+)\|(\w+)\|(\w+)\|(\S+)")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> list[dict]:
    """
    Parse a FASTA file into a list of dicts:
        {header, proteome_id, uniprot_id, sequence}
    Raises on malformed headers or empty sequences.
    """
    records = []
    current_header = None
    current_seq_lines = []

    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous record
                if current_header is not None:
                    seq = "".join(current_seq_lines)
                    if not seq:
                        raise ValueError(
                            f"Empty sequence for header: {current_header}"
                        )
                    records.append({**current_header, "sequence": seq})
                    current_seq_lines = []

                m = HEADER_RE.match(line)
                if not m:
                    raise ValueError(
                        f"Unexpected header format at line {lineno}:\n  {line}\n"
                        f"Expected: >PROTEOME_ID|db|UNIPROT_ID|GENE_NAME ..."
                    )
                current_header = {
                    "header":      line[1:],
                    "proteome_id": m.group(1),
                    "uniprot_id":  m.group(3),
                }
            else:
                if current_header is None:
                    raise ValueError(f"Sequence data before first header at line {lineno}")
                current_seq_lines.append(line)

    # Save last record
    if current_header is not None:
        seq = "".join(current_seq_lines)
        if not seq:
            raise ValueError(f"Empty sequence for header: {current_header}")
        records.append({**current_header, "sequence": seq})

    return records


def check_duplicates(records: list[dict], label: str) -> None:
    """Warn if any UniProt IDs appear more than once."""
    counts = Counter(r["uniprot_id"] for r in records)
    dupes = {uid: n for uid, n in counts.items() if n > 1}
    if dupes:
        log.warning(
            f"[{label}] {len(dupes)} duplicate UniProt IDs found "
            f"(showing first 5): {list(dupes.items())[:5]}"
        )
    else:
        log.info(f"[{label}] No duplicate UniProt IDs.")


def summarise(records: list[dict], label: str) -> None:
    """Log sequence count and length distribution."""
    lengths = [len(r["sequence"]) for r in records]
    proteomes = len(set(r["proteome_id"] for r in records))
    log.info(
        f"[{label}] {len(records):,} sequences | "
        f"{proteomes:,} unique proteomes | "
        f"length: min={min(lengths)}, max={max(lengths)}, "
        f"mean={sum(lengths)//len(lengths)}"
    )


def validate(path: Path, label: str) -> list[dict]:
    log.info(f"[{label}] Validating {path.name} ...")
    records = parse_fasta(path)
    if not records:
        raise ValueError(f"[{label}] No sequences found in {path}")
    check_duplicates(records, label)
    summarise(records, label)
    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validate and register representative FASTA files."
    )
    parser.add_argument("--class1", type=Path, required=True,
                        help="Class I representative FASTA")
    parser.add_argument("--class2", type=Path, required=True,
                        help="Class II representative FASTA")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate only, do not copy files")
    args = parser.parse_args()

    # Validate inputs exist
    for p, label in [(args.class1, "Class I"), (args.class2, "Class II")]:
        if not p.exists():
            log.error(f"{label} FASTA not found: {p}")
            sys.exit(1)

    # Validate content
    validate(args.class1, "Class I")
    validate(args.class2, "Class II")

    if args.dry_run:
        log.info("Dry-run complete. No files copied.")
        return

    # Copy to canonical paths in data/processed/
    out_dir = Path(PATHS["data_processed"])
    out_dir.mkdir(parents=True, exist_ok=True)

    destinations = {
        args.class1: out_dir / "class1_sequences.fasta",
        args.class2: out_dir / "class2_sequences.fasta",
    }

    for src, dst in destinations.items():
        shutil.copy2(src, dst)
        log.info(f"Copied {src.name} → {dst}")

    log.info("01_search complete. Canonical FASTAs ready for 02_align:")
    log.info(f"  {destinations[args.class1]}")
    log.info(f"  {destinations[args.class2]}")


if __name__ == "__main__":
    main()