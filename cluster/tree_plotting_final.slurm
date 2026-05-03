#!/bin/bash
#SBATCH --job-name=tree_plotting_class2_circular_only
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/tree_plotting_class2_circ_only_%j.out
#SBATCH --cpus-per-task=8            
#SBATCH --mem=64G
#SBATCH --time=48:00:00    

#set -euo pipefail

PROJECT_ROOT=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution
cd $PROJECT_ROOT

umask 003

mkdir -p slurm_out

source /sci/labs/itamarsi/dror.jacobi/miniforge3/etc/profile.d/conda.sh
conda activate /sci/labs/itamarsi/dror.jacobi/miniforge3/envs/fumaraseevo

	




#unrooted trees

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class1/AA_clustered_2k.treefile results/figures/unrooted/unrooted_tree_class1_AA_2k_final.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class1/3di_clustered_2k_bb.treefile results/figures/unrooted/unrooted_tree_class1_3di_clustered_2k_bb_final.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered_qpfam_r8_fast.treefile results/figures/unrooted/unrooted_tree_class2_AA__clustered6k_R8_fast.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered2k_qpfam_r8.treefile results/figures/unrooted/unrooted_tree_class2_AA_clustered2k_R8_final.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/3di_clustered_2k_bb.treefile results/figures/unrooted/unrooted_tree_class2_3di_clustered_2k_bb_final.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered_clustered_qpfam_r8_2k_copy_it347_963.treefile
results/figures/unrooted/unrooted_tree_class2_AA_clustered_qpfam_r8_2k_copy_it347_963.treefile.png

# Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered_qpfam_r8_assisted_copy_it190_956.treefile results/figures/unrooted/unrooted_tree_class2_AA_clustered_qpfam_r8_assisted_copy_it190_956.png


#circular trees

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class1/AA_clustered_2k.treefile results/figures/circular/circular_tree_class1_AA_clustered_2k_final.png

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class1/3di_clustered_2k_bb.treefile results/figures/circular/circular_tree_class1_3di_clustered_2k_bb_final.png

Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered2k_qpfam_r8.treefile results/figures/circular/circular_tree_class2_AA_clustered2k_R8_final.png

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered2k_qpfam_r8_fast.treefile results/figures/circular/circular_tree_class2_AA_clustered2k_R8_fast.png

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/3di_clustered_2k_bb.treefile results/figures/circular/circular_tree_class2_3di_clustered_2k_bb_final.png

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered_clustered_qpfam_r8_2k_copy_it347_963.treefile
results/figures/circular/circular_tree_class2_AA_clustered_qpfam_r8_2k_copy_it347_963.png

# Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/AA_clustered_qpfam_r8_assisted_copy_it190_956.treefile results/figures/circular/circular_tree_class2_AA_clustered_qpfam_r8_assisted_copy_it190_956.png