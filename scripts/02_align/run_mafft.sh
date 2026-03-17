#!/bin/bash
# =============================================================================
#  run_mafft.sh - Align a FASTA file using MAFFT auto mode
#
#  Usage:
#    bash scripts/02_align/run_mafft.sh <input.fasta> <output.fasta> [--threads N] [--aamatrix MATRIX]
#
#  Example:
#    bash scripts/02_align/run_mafft.sh \
#         data/processed/class1_sequences.fasta \
#         results/alignments/class1_aligned.fasta \
#         --threads 8
# =============================================================================

set -euo pipefail

# --- Parse arguments ---
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.fasta> <output.fasta> [--threads N] [--aamatrix file.mat]"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"
THREADS=4             
AAMATRIX=""           # Default is empty (allows MAFFT to use internal JTT defaults)

shift 2
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --aamatrix) AAMATRIX="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

# --- Validate ---
if ! command -v mafft &>/dev/null; then
  echo "ERROR: mafft not found. Install via: conda install -c bioconda mafft"
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: Input file not found: $INPUT"
  exit 1
fi

if [[ -f "$OUTPUT" ]]; then
  echo "[SKIP] Output already exists: $OUTPUT"
  exit 0
fi

# If a custom matrix was provided, ensure the file actually exists
if [[ -n "$AAMATRIX" ]] && [[ ! -f "$AAMATRIX" ]]; then
  echo "ERROR: Custom amino acid matrix file not found: $AAMATRIX"
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"
LOG="${OUTPUT%.fasta}.log"

N_SEQS=$(grep -c "^>" "$INPUT")
echo "========================================"
echo " MAFFT alignment"
echo " Input:   $INPUT ($N_SEQS sequences)"
echo " Output:  $OUTPUT"
echo " Threads: $THREADS"
if [[ -n "$AAMATRIX" ]]; then
    echo " Matrix:  $AAMATRIX"
else
    echo " Matrix:  Default (Auto JTT)"
fi
echo "========================================"

# --- Build the MAFFT command dynamically ---
# Start with the base command
MAFFT_CMD=(mafft --auto --thread "$THREADS" --reorder)

# Add the aamatrix flag only if it was provided
if [[ -n "$AAMATRIX" ]]; then
    MAFFT_CMD+=(--aamatrix "$AAMATRIX")
fi

# Add the input file
MAFFT_CMD+=("$INPUT")

# Execute the array command
"${MAFFT_CMD[@]}" > "$OUTPUT" 2> "$LOG"

N_ALIGNED=$(grep -c "^>" "$OUTPUT")
echo "[DONE] $N_ALIGNED sequences aligned -> $OUTPUT"
echo "       Log: $LOG"