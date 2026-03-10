#!/bin/bash
# =============================================================================
#  mmseq_run.sh - Download proteome FASTAs and run MMseqs2 easy-search
#
#  Usage:
#    bash mmseq_run.sh <tsv_file> <query_fasta> [--jobs N]
#
#  Arguments:
#    tsv_file      TSV of proteomes (columns: proteome_id, organism,
#                  organism_id, protein_count, busco, cpd)
#    query_fasta   Query FASTA (e.g. E. coli FumA or FumC)
#    --jobs N      Number of parallel jobs (default: 4)
#
#  Paths are resolved from project config (config.py via resolve_paths.py):
#    FASTAs  → data/raw/fastas/
#    .m8s    → data/raw/m8s/
#    tmp     → data/raw/tmp/
#
#  Example:
#    bash scripts/01_search/mmseq_run.sh \
#         data/processed/proteomes_class1.tsv \
#         data/external/ecoli_fumA.fasta \
#         --jobs 8
# =============================================================================

set -euo pipefail

# --- Resolve project root (two levels up from this script) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# --- Parse arguments ---
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <tsv_file> <query_fasta> [--jobs N]"
  exit 1
fi

TSV_FILE="$1"
QUERY_FASTA="$2"
N_JOBS=4

shift 2
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) N_JOBS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

# --- Resolve output paths from config.py ---
RAW_DIR=$(python3 -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
from config import PATHS
print(PATHS['data_raw'])
")

FASTA_DIR="$RAW_DIR/fastas"
M8_DIR="$RAW_DIR/m8s"
TMP_DIR="$RAW_DIR/tmp"

mkdir -p "$FASTA_DIR" "$M8_DIR" "$TMP_DIR"

# --- Validate inputs ---
if [[ ! -f "$TSV_FILE" ]]; then
  echo "Error: TSV file '$TSV_FILE' not found."
  exit 1
fi
if [[ ! -f "$QUERY_FASTA" ]]; then
  echo "Error: Query FASTA '$QUERY_FASTA' not found."
  exit 1
fi

# --- Derive class prefix from query filename ---
QUERY_BASENAME=$(basename "$QUERY_FASTA" .fasta)
if [[ "$QUERY_BASENAME" == *"fumA"* || "$QUERY_BASENAME" == *"fuma"* ]]; then
  CLASS_PREFIX="class1"
elif [[ "$QUERY_BASENAME" == *"fumC"* || "$QUERY_BASENAME" == *"fumc"* ]]; then
  CLASS_PREFIX="class2"
else
  echo "Error: cannot determine class from query filename '$QUERY_BASENAME'."
  echo "Filename must contain 'fumA' or 'fumC'."
  exit 1
fi

# --- Rewrite query header with class prefix into a temp file ---
# Original: >sp|P0AC33|FUMA_ECOLI ...
# Becomes:  >class1_P0AC33
PREFIXED_QUERY="$TMP_DIR/query_${CLASS_PREFIX}.fasta"
awk -v prefix="$CLASS_PREFIX" '
  /^>/ { split($0, a, "|"); print ">" prefix "_" a[2]; next }
  { print }
' "$QUERY_FASTA" > "$PREFIXED_QUERY"

echo "Query prefix: $CLASS_PREFIX ($(head -1 "$PREFIXED_QUERY"))"

# --- Per-proteome worker function (called by GNU parallel) ---
process_proteome() {
  local proteome_id="$1"
  local query_fasta="$2"
  local fasta_dir="$3"
  local m8_dir="$4"
  local tmp_dir="$5"

  local fasta_file="$fasta_dir/${proteome_id}.fasta"
  local m8_file="$m8_dir/${proteome_id}_output.m8"

  # Skip if .m8 already exists (allows safe re-runs)
  if [[ -f "$m8_file" ]]; then
    echo "[SKIP] $proteome_id - .m8 already exists."
    return 0
  fi

  # Download FASTA
  local url="https://rest.uniprot.org/uniprotkb/stream?query=proteome:${proteome_id}&format=fasta"
  echo "[DOWNLOAD] $proteome_id"
  if ! curl -s --retry 3 --retry-delay 5 -o "$fasta_file" "$url"; then
    echo "[ERROR] Failed to download FASTA for $proteome_id."
    return 1
  fi

  # Sanity-check: non-empty FASTA
  if [[ ! -s "$fasta_file" ]]; then
    echo "[WARN] Empty FASTA for $proteome_id - skipping MMseqs."
    return 1
  fi

  # Run MMseqs2
  echo "[MMSEQS] $proteome_id"
  if mmseqs easy-search \
      "$(realpath "$query_fasta")" \
      "$fasta_file" \
      "$m8_file" \
      "$tmp_dir/${proteome_id}_tmp" \
      --threads 1 \
      -s 7.5 \
      --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen" \
      2>>"$m8_dir/${proteome_id}.log"; then
    echo "[DONE] $proteome_id → $m8_file"
  else
    echo "[ERROR] MMseqs failed for $proteome_id. See $m8_dir/${proteome_id}.log"
    return 1
  fi
}

export -f process_proteome

# --- Extract proteome IDs from TSV (skip header) and run in parallel ---
echo "========================================"
echo " MMseqs2 search"
echo " TSV:    $TSV_FILE"
echo " Query:  $QUERY_FASTA"
echo " Class:  $CLASS_PREFIX"
echo " FASTAs: $FASTA_DIR"
echo " .m8s:   $M8_DIR"
echo " Jobs:   $N_JOBS"
echo "========================================"

tail -n +2 "$TSV_FILE" \
  | cut -f1 \
  | parallel \
      --jobs "$N_JOBS" \
      --bar \
      --halt soon,fail=1 \
      process_proteome {} \
        "$PREFIXED_QUERY" \
        "$FASTA_DIR" \
        "$M8_DIR" \
        "$TMP_DIR"

echo ""
echo "All proteomes processed."
echo "  FASTAs saved to: $FASTA_DIR"
echo "  .m8 files saved to: $M8_DIR"