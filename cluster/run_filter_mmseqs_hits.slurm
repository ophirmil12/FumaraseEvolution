#!/bin/bash
#SBATCH --job-name=filter_hits
#SBATCH --out=slurm_out/filter_hits_%j.out
#SBATCH --killable
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

python "$PROJECT_ROOT/scripts/01_search/filter_hits.py" --class1 --class2
