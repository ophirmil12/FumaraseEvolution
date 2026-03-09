"""
extract_sequences.py - Fetch representative FASTA files from UniProt.

Reads data/processed/filtered_hits.json (produced by filter_hits.py), collects
the unique UniProt accessions per fumarase class, fetches their sequences from
the UniProt REST API in batches, validates the results, and writes the canonical
FASTA files used by all downstream steps.

Output:
    data/processed/class1_sequences.fasta
    data/processed/class2_sequences.fasta

Usage:
    python scripts/01_search/extract_sequences.py             # both classes
    python scripts/01_search/extract_sequences.py --class1    # Class I only
    python scripts/01_search/extract_sequences.py --class2    # Class II only
    python scripts/01_search/extract_sequences.py --dry-run   # count only, no fetch
    python scripts/01_search/extract_sequences.py --force     # re-fetch if file exists

Note:
    Requires internet access to rest.uniprot.org.
    Large datasets (10,000+ sequences) may take several minutes.
"""

import sys
import json
import time
import argparse
import logging
from pathlib import Path
from collections import defaultdict, Counter

import requests

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
# Constants
# ---------------------------------------------------------------------------

UNIPROT_BATCH_URL = "https://rest.uniprot.org/uniprotkb/accessions"
BATCH_SIZE    = 200   # UniProt recommends <= 200 accessions per request
RETRY_MAX     = 3
RETRY_DELAY   = 5     # seconds between retries
MIN_SEQ_LEN   = 100   # shorter sequences are flagged (fragments)


# ---------------------------------------------------------------------------
# Load hits
# ---------------------------------------------------------------------------

def load_hits(json_path: Path) -> dict[int, dict[str, str]]:
    """
    Read filtered_hits.json and return per-class mapping:
        {fumarase_class: {hit_id: proteome_id}}

    One hit_id may appear in multiple proteomes; the first occurrence is kept
    (highest-scoring, since filter_hits.py processes files in sorted order).
    """
    if not json_path.exists():
        log.error(f"filtered_hits.json not found: {json_path}")
        log.error("Run filter_hits.py first.")
        sys.exit(1)

    with open(json_path) as f:
        hits = json.load(f)

    by_class: dict[int, dict[str, str]] = defaultdict(dict)
    for hit in hits:
        cls      = hit["fumarase_class"]
        acc      = hit["hit_id"]
        proteome = hit["proteome_id"]
        if acc not in by_class[cls]:
            by_class[cls][acc] = proteome

    for cls, acc_map in by_class.items():
        log.info(f"Class {cls}: {len(acc_map):,} unique accessions to fetch")

    return dict(by_class)


# ---------------------------------------------------------------------------
# Fetch from UniProt
# ---------------------------------------------------------------------------

def fetch_batch(accessions: list[str]) -> str:
    """Fetch a batch of UniProt accessions as FASTA text, with retries."""
    params = {"accessions": ",".join(accessions), "format": "fasta"}
    for attempt in range(1, RETRY_MAX + 1):
        try:
            r = requests.get(UNIPROT_BATCH_URL, params=params, timeout=60)
            r.raise_for_status()
            return r.text
        except requests.RequestException as e:
            log.warning(f"  Attempt {attempt}/{RETRY_MAX} failed: {e}")
            if attempt < RETRY_MAX:
                time.sleep(RETRY_DELAY)
    log.error(f"Failed to fetch batch after {RETRY_MAX} attempts.")
    return ""


def parse_uniprot_fasta(fasta_text: str,
                        acc_to_proteome: dict[str, str]) -> list[tuple[str, str]]:
    """
    Parse a raw UniProt FASTA response and prepend the proteome ID to each header.

    UniProt headers look like:
        >sp|P0AC33|FUMA_ECOLI Fumarate hydratase class I OS=...
        >tr|A0A066VRU3|A0A066VRU3_9PROT Fumarate hydratase OS=...

    Output header format used by all downstream scripts:
        >UP000000625|sp|P0AC33|FUMA_ECOLI Fumarate hydratase class I OS=...
    """
    records: list[tuple[str, str]] = []
    current_header: str | None = None
    current_seq: list[str] = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header and current_seq:
                records.append((current_header, "".join(current_seq)))
                current_seq = []
            parts     = line[1:].split("|")
            accession = parts[1] if len(parts) >= 2 else ""
            proteome  = acc_to_proteome.get(accession, "UNKNOWN")
            current_header = f"{proteome}|{line[1:]}"
        else:
            current_seq.append(line)

    if current_header and current_seq:
        records.append((current_header, "".join(current_seq)))

    return records


def fetch_all(acc_to_proteome: dict[str, str]) -> list[tuple[str, str]]:
    """
    Fetch all accessions in batches and return (header, sequence) pairs.
    Logs any accessions that UniProt did not return.
    """
    accessions = list(acc_to_proteome.keys())
    total      = len(accessions)
    records: list[tuple[str, str]] = []
    missing: list[str] = []

    for i in range(0, total, BATCH_SIZE):
        batch = accessions[i : i + BATCH_SIZE]
        end   = min(i + len(batch), total)
        log.info(f"  Fetching {i + 1}-{end} / {total} ...")
        fasta_text = fetch_batch(batch)

        if not fasta_text.strip():
            log.warning(f"  Empty response for batch starting at index {i}")
            missing.extend(batch)
            continue

        batch_records = parse_uniprot_fasta(fasta_text, acc_to_proteome)
        records.extend(batch_records)

        returned     = {r[0].split("|")[2] for r in batch_records}
        not_returned = set(batch) - returned
        if not_returned:
            sample = sorted(not_returned)[:5]
            suffix = "..." if len(not_returned) > 5 else ""
            log.warning(f"  {len(not_returned)} accessions not returned: {sample}{suffix}")
            missing.extend(not_returned)

        time.sleep(0.5)   # be polite to the API

    if missing:
        log.warning(f"{len(missing)} accessions total could not be fetched from UniProt.")

    log.info(f"Fetched {len(records):,} sequences ({len(missing)} missing)")
    return records


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

def validate(records: list[tuple[str, str]], label: str) -> None:
    """Log summary stats and warn about duplicates or short sequences."""
    if not records:
        log.error(f"[{label}] No sequences — nothing to write.")
        return

    lengths    = [len(seq) for _, seq in records]
    accessions = []
    for header, _ in records:
        parts = header.split("|")
        accessions.append(parts[2] if len(parts) >= 3 else header)

    proteomes = len({h.split("|")[0] for h, _ in records})

    log.info(
        f"[{label}] {len(records):,} sequences | "
        f"{proteomes:,} proteomes | "
        f"length min={min(lengths)} max={max(lengths)} mean={sum(lengths)//len(lengths)}"
    )

    dupes = {acc: n for acc, n in Counter(accessions).items() if n > 1}
    if dupes:
        log.warning(f"[{label}] {len(dupes)} duplicate accessions: {list(dupes.items())[:5]}")
    else:
        log.info(f"[{label}] No duplicate accessions.")

    short = sum(1 for l in lengths if l < MIN_SEQ_LEN)
    if short:
        log.warning(f"[{label}] {short} sequences shorter than {MIN_SEQ_LEN} aa "
                    f"(will be removed by clean_sequences.py)")


# ---------------------------------------------------------------------------
# Write FASTA
# ---------------------------------------------------------------------------

def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write (header, sequence) pairs to a FASTA file, wrapping sequence at 60 chars."""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
    log.info(f"Written {len(records):,} sequences -> {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fetch fumarase sequences from UniProt and write canonical FASTAs."
    )
    parser.add_argument("--class1",  action="store_true", help="Fetch Class I only")
    parser.add_argument("--class2",  action="store_true", help="Fetch Class II only")
    parser.add_argument(
        "--hits", type=Path,
        default=Path(PATHS["data_processed"]) / "filtered_hits.json",
        help="Path to filtered_hits.json (default: data/processed/filtered_hits.json)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report accession counts only, do not fetch or write"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-fetch even if output file already exists"
    )
    args = parser.parse_args()

    # Default: both classes
    if not args.class1 and not args.class2:
        args.class1 = True
        args.class2 = True

    out_dir = Path(PATHS["data_processed"])
    out_dir.mkdir(parents=True, exist_ok=True)

    by_class = load_hits(args.hits)

    class_cfg = {
        1: (out_dir / "class1_sequences.fasta", args.class1),
        2: (out_dir / "class2_sequences.fasta", args.class2),
    }

    for cls, (out_path, should_run) in class_cfg.items():
        if not should_run or cls not in by_class:
            continue

        if args.dry_run:
            log.info(f"Class {cls}: would fetch {len(by_class[cls]):,} sequences [dry-run]")
            continue

        if out_path.exists() and not args.force:
            log.info(f"Class {cls}: {out_path.name} already exists — skipping "
                     "(use --force to re-fetch)")
            continue

        log.info(f"Class {cls}: fetching from UniProt ...")
        records = fetch_all(by_class[cls])

        if not records:
            log.error(f"Class {cls}: no sequences retrieved. Check network and accessions.")
            continue

        validate(records, f"Class {cls}")
        write_fasta(out_path, records)

    if not args.dry_run:
        log.info("Done. Next step: python scripts/01_search/clean_sequences.py")


if __name__ == "__main__":
    main()