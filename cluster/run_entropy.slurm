#!/bin/bash
#SBATCH --job-name=entropy
#SBATCH --out=slurm_out/entropy_%j.out
#SBATCH --killable
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

python "$PROJECT_ROOT/scripts/02_align/entropy.py" \
    --class1 "$PROJECT_ROOT/results/alignments/class1_aligned.fasta" \
    --class2 "$PROJECT_ROOT/results/alignments/class2_aligned.fasta" \
    --baseline-conserved "$PROJECT_ROOT/data/external/conservation_baseline/rplA_aligned.fasta" \
    --baseline-divergent "$PROJECT_ROOT/data/external/conservation_baseline/mfs_transporter_aligned.fasta" \
    --seed 42