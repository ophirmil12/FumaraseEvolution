#!/bin/bash
#SBATCH --job-name=tree_data_extraction
#SBATCH --out=slurm_out/extracting_tree_data%j.out
#SBATCH --time=2:00:00
#SBATCH --mem=12
#SBATCH --cpus-per-task=2


set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env


echo "Starting tree information extraction..."
echo "Date: $(date)"


echo
mkdir -p "$PROJECT_ROOT/results/trees/tree_data"

#Syntax: 
#treefile processing:

#python "$PROJECT_ROOT/scripts/07_trees/extract_iqtree_data_final.py" \
#    <input treefile path>.treefile \
#    <output summary file path>.csv

#additional iqtree information scavanging:

#python "$PROJECT_ROOT/scripts_backup/07_trees/tree_stats_extended.py" \
#    <input summary file path>.csv \
#    <iqtree logs file path>.iqtree

python "$PROJECT_ROOT/scripts/07_trees/extract_iqtree_data_final.py" \
  "$PROJECT_ROOT/results/trees/class1/subtrees/class1_subtree1_ver2.treefile" \
  "$PROJECT_ROOT/results/trees/class1/tree_data/class1_subtree1_ver2_summary.csv" 

python "$PROJECT_ROOT/scripts_backup/07_trees/tree_stats_extended.py" \
  "$PROJECT_ROOT/results/trees/class1/tree_data/class1_subtree1_ver2_summary.csv" \
  "$PROJECT_ROOT/results/trees/class1/subtrees/class1_subtree1_ver2.iqtree"



echo "Extraction Script Finished at $(date)."

