#!/bin/bash
#SBATCH --job-name=data_extraction
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/extracting_tree_data%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=64
#SBATCH --cpus-per-task=8

echo "SLURM script Starting on $(hostname) at $(date)..."

#set -euo pipefail

PROJECT_ROOT=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution
cd $PROJECT_ROOT

umask 003

mkdir -p slurm_out

source /sci/labs/itamarsi/dror.jacobi/miniforge3/etc/profile.d/conda.sh
conda activate /sci/labs/itamarsi/dror.jacobi/miniforge3/envs/fumaraseevo

echo "Starting script..."


mkdir -p results/trees/tree_data

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class1/subtrees/class1_subtree1_ver2.treefile results/trees/class1/tree_data/class1_subtree1_ver2_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class1/tree_data/class1_subtree1_ver2_summary.csv results/trees/class1/subtrees/class1_subtree1_ver2.iqtree

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class1/subtrees/class1_subtree2.treefile results/trees/class1/tree_data/class1_subtree2_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class1/tree_data/class1_subtree2_summary.csv results/trees/class1/subtrees/class1_subtree2.iqtree


#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree1_thinned_300.treefile results/trees/class2/tree_data/class2_subtree1_thinned_300_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree1_thinned_300_summary.csv results/trees/class2/subtrees/class2_subtree1_thinned_300.iqtree

python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree2A_ver3.treefile results/trees/class2/tree_data/class2_subtree2A_ver3_summary.csv
python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree2A_ver3_summary.csv results/trees/class2/subtrees/class2_subtree2A_ver3.iqtree

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree2B.treefile results/trees/class2/tree_data/class2_subtree2B_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree2B_summary.csv results/trees/class2/subtrees/class2_subtree2B.iqtree

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree3.treefile results/trees/class2/tree_data/class2_subtree3_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree3_summary.csv results/trees/class2/subtrees/class2_subtree3.iqtree

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree3_ver2.treefile results/trees/class2/tree_data/class2_subtree3_ver2_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree3_ver2_summary.csv results/trees/class2/subtrees/class2_subtree3_ver2.iqtree

#python scripts_backup/07_trees/extract_iqtree_data_final.py results/trees/class2/subtrees/class2_subtree4.treefile results/trees/class2/tree_data/class2_subtree4_summary.csv
#python scripts_backup/07_trees/tree_stats_extended.py results/trees/class2/tree_data/class2_subtree4_summary.csv results/trees/class2/subtrees/class2_subtree4.iqtree


echo "Extraction Script Finished at $(date)."

