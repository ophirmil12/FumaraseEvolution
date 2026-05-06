#!/bin/bash
#SBATCH --job-name=tree_plotting_class1
#SBATCH --out=slurm_out/tree_plotting_class1_%j.out
#SBATCH --cpus-per-task=8            
#SBATCH --mem=32G
#SBATCH --time=48:00:00    

echo "SLURM script Starting on $(hostname) at $(date)..."

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

RESULT="$PROJECT_ROOT/results"
SCRIPT="$PROJECT_ROOT/scripts"
DATA="$PROJECT_ROOT/data"

echo "Starting iqtree script..."


#unrooted trees

Rscript "$SCRIPT/08_plots/tree_plot_unrooted.R" \
		"$DATA/processed/proteome_metadata.tsv" \
		"$RESULT/trees/class1/AA_clustered_2k.treefile" \
		"$RESULT/figures/unrooted/unrooted_tree_class1_AA_2k_final.png"


#circular trees
Rscript "$SCRIPT/08_plots/tree_plot_circular.R" \
		"$DATA/processed/proteome_metadata.tsv" \
		"$RESULT/trees/class1/AA_clustered_2k.treefile" \
		"$RESULT/figures/circular/circular_tree_class1_AA_clustered_2k_final.png"
