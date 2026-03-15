#!/bin/bash
#SBATCH --job-name=plot_identity
#SBATCH --out=slurm_out/plot_identity_%j.out
#SBATCH --killable
#SBATCH --time=01:00:00
#SBATCH --mem=72G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

python "$PROJECT_ROOT/scripts/08_plots/identity_plot.py" \
    --statsdir "$PROJECT_ROOT/results/stats" \
    --outdir   "$PROJECT_ROOT/results/figures" \
    --sample   5000 \
    --scatter  500

echo "Done. Outputs in $PROJECT_ROOT/results/figures/"