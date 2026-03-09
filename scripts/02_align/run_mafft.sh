#!/bin/bash
# =============================================================================
#  run_mafft.sh - Run MAFFT multiple sequence alignments for Class I and II
#
#  Aligns class1_sequences.fasta and class2_sequences.fasta separately using
#  MAFFT auto mode. Outputs aligned FASTAs to results/alignments/.
#
#  Usage:
#    bash scripts/02_align/run_mafft.sh            # both classes
#    bash scripts/02_align/run_mafft.sh --class1   # Class I only
#    bash scripts/02_align/run_mafft.sh --class2   # Class II only
#    bash scripts/02_align/run_mafft.sh --threads 8
#
#  Output:
#    results/alignments/class1_aligned.fasta
#    results/alignments/class2_aligned.fasta
#    results/alignments/class1_mafft.log
#    results/alignments/class2_mafft.log
#
# =============================================================================
#
#  INSTALLATION - MAFFT
#  ---------------------
#
#  Option 1: Conda (recommended - manages all bioinformatics dependencies)
#  ------------------------------------------------------------------------
#  If you don't have conda, install Miniconda first:
#    https://docs.conda.io/en/latest/miniconda.html
#
#  Then create the project environment (first time only):
#    conda env create -f envs/environment.yml
#    conda activate fumaraseevo
#
#  To add MAFFT to an existing environment:
#    conda activate fumaraseevo
#    conda install -c bioconda mafft
#
#  Verify installation:
#    mafft --version
#
#  Option 2: Homebrew (macOS)
#  --------------------------
#    brew install mafft
#
#  Option 3: apt (Ubuntu/Debian)
#  -----------------------------
#    sudo apt update && sudo apt install mafft
#
#  Option 4: Manual download
#  -------------------------
#    https://mafft.cbrc.jp/alignment/software/
#
#  Option 5: HPC cluster (HUJI-CSE)
#  ---------------------------------
#    module load mafft        # check available: module avail mafft
#    # or use the conda env (preferred):
#    conda activate fumaraseevo
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# --- Resolve paths from config.py ---
DATA_PROCESSED=$(python3 -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
from config import PATHS; print(PATHS['data_processed'])
")
ALIGNMENTS_DIR=$(python3 -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
from config import PATHS; print(PATHS['alignments'])
")

# --- Defaults ---
RUN_CLASS1=false
RUN_CLASS2=false
THREADS=4

# --- Parse arguments ---
if [[ $# -eq 0 ]]; then
  RUN_CLASS1=true
  RUN_CLASS2=true
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --class1)  RUN_CLASS1=true; shift ;;
    --class2)  RUN_CLASS2=true; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

# --- Check MAFFT is available ---
if ! command -v mafft &>/dev/null; then
  echo ""
  echo "ERROR: mafft not found in PATH."
  echo ""
  echo "Install via conda (recommended):"
  echo "  conda env create -f envs/environment.yml"
  echo "  conda activate fumaraseevo"
  echo ""
  echo "Or: conda install -c bioconda mafft"
  echo "See full install notes at the top of this script."
  exit 1
fi

mkdir -p "$ALIGNMENTS_DIR"

# --- Alignment function ---
run_alignment() {
  local class_label="$1"
  local input_fasta="$DATA_PROCESSED/${class_label}_sequences.fasta"
  local output_fasta="$ALIGNMENTS_DIR/${class_label}_aligned.fasta"
  local log_file="$ALIGNMENTS_DIR/${class_label}_mafft.log"

  if [[ ! -f "$input_fasta" ]]; then
    echo "ERROR: Input FASTA not found: $input_fasta"
    echo "Run scripts/01_search/extract_sequences.py first."
    exit 1
  fi

  # Skip if alignment already exists
  if [[ -f "$output_fasta" ]]; then
    echo "[SKIP] $class_label alignment already exists: $output_fasta"
    return 0
  fi

  local n_seqs
  n_seqs=$(grep -c "^>" "$input_fasta")
  echo ""
  echo "========================================"
  echo " Aligning $class_label ($n_seqs sequences)"
  echo " Input:   $input_fasta"
  echo " Output:  $output_fasta"
  echo " Threads: $THREADS"
  echo "========================================"

  # MAFFT auto mode:
  #   --auto     selects the best strategy based on dataset size
  #              (L-INS-i for <200 seqs, FFT-NS-2 for larger datasets)
  #   --thread   number of CPU threads
  #   --reorder  reorder sequences so similar ones are adjacent in output
  #   2>         redirect progress to log file
  mafft \
    --auto \
    --thread "$THREADS" \
    --reorder \
    "$input_fasta" \
    > "$output_fasta" \
    2> "$log_file"

  local n_aligned
  n_aligned=$(grep -c "^>" "$output_fasta")
  echo "[DONE] $class_label: $n_aligned sequences aligned -> $output_fasta"
  echo "       Log: $log_file"
}

# --- Run ---
$RUN_CLASS1 && run_alignment "class1"
$RUN_CLASS2 && run_alignment "class2"

echo ""
echo "Alignments complete. Next step:"
echo "  python scripts/02_align/entropy.py"