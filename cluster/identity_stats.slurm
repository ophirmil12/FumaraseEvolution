#!/bin/bash
#SBATCH --job-name=identity_stats
#SBATCH --out=slurm_out/identity_stats_%j.out
#SBATCH --killable
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

echo "========================================"
echo " Sequence Identity Statistics — Class II"
echo " $(date)"
echo "========================================"

python "$PROJECT_ROOT/scripts/04_identity/stats_identity.py" \
    --allvsall  "$PROJECT_ROOT/data/processed/class2_allvsall.tsv" \
    --metadata  "$PROJECT_ROOT/data/processed/proteome_metadata.tsv" \
    --outdir    "$PROJECT_ROOT/results/stats"

echo ""
echo "========================================"
echo " Done. $(date)"
echo " Outputs: $PROJECT_ROOT/results/stats/"
echo "========================================"
