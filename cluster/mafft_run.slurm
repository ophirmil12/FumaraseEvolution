#!/bin/bash
#SBATCH --job-name=mafft_align
#SBATCH --out=slurm_out/mafft_%j.out
#SBATCH --time=09:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --killable

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

MAFFT="$PROJECT_ROOT/scripts/02_align/run_mafft.sh"
DATA="$PROJECT_ROOT/data"
RESULTS="$PROJECT_ROOT/results/alignments"
THREADS="$SLURM_CPUS_PER_TASK"

# --- Class I ---
bash "$MAFFT" \
    "$DATA/processed/class1_sequences.fasta" \
    "$RESULTS/class1_aligned.fasta" \
    --threads "$THREADS"

# --- Class II ---
bash "$MAFFT" \
    "$DATA/processed/class2_sequences.fasta" \
    "$RESULTS/class2_aligned.fasta" \
    --threads "$THREADS"

# --- Baseline: conserved (RplA) ---
bash "$MAFFT" \
    "$DATA/external/conservation_baseline/rplA.fasta" \
    "$DATA/external/conservation_baseline/rplA_aligned.fasta" \
    --threads "$THREADS"

# --- Baseline: divergent (MFS transporter) ---
bash "$MAFFT" \
    "$DATA/external/conservation_baseline/mfs_transporter.fasta" \
    "$DATA/external/conservation_baseline/mfs_transporter_aligned.fasta" \
    --threads "$THREADS"

# Note: need to run predict_3di.py first to generate the 3Di FASTAs [predict_3di.slurm], then align them with MAFFT
# --- 3Di sequences (class I) ---
bash "$MAFFT" \
    "$DATA/processed/class1_3di.fasta" \
    "$DATA/processed/class1_3di_aligned.fasta" \
    --threads "$THREADS"

# --- 3Di sequences (class II) ---
bash "$MAFFT" \
    "$DATA/processed/class2_3di.fasta" \
    "$DATA/processed/class2_3di_aligned.fasta" \
    --threads "$THREADS"

echo "All alignments complete."