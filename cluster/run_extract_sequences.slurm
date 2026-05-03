#!/bin/bash
#SBATCH --job-name=extract_seqs
#SBATCH --out=slurm_out/extract_seqs_%j.out
#SBATCH --killable
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

# Download sequences of hits
#python "$PROJECT_ROOT/scripts/01_search/extract_sequences.py"

echo "-------------------------------------------------------------------"

# Clean sequences
#python "$PROJECT_ROOT/scripts/01_search/clean_sequences.py"

echo "-------------------------------------------------------------------"

# Create tables
python "$PROJECT_ROOT/scripts/01_search/generate_tables.py"
