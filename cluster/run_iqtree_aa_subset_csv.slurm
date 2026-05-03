#!/bin/bash
#SBATCH --job-name=subtree23B
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/subtree23B_%j.out
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=END
#SBATCH --mail-user=dror.jacobi@mail.huji.ac.il


echo "SLURM script Starting on $(hostname) at $(date)..."

#set -euo pipefail

PROJECT_ROOT=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution
cd $PROJECT_ROOT

umask 003

mkdir -p slurm_out

source /sci/labs/itamarsi/dror.jacobi/miniforge3/etc/profile.d/conda.sh
conda activate /sci/labs/itamarsi/dror.jacobi/miniforge3/envs/fumaraseevo

echo "Starting iqtree script..."

#------------ create csvs from text if needed ---------------------

#python scripts_backup/08_plots/generate_metadata_csv.py data/processed/subsets/class1_sequences_subset.txt data/processed/proteome_metadata.tsv data/processed/subsets/class1_sequences_subset.csv

#python scripts_backup/08_plots/generate_metadata.py \
	 data/processed/subsets/HGT_subtree_organisms.csv \
	 data/processed/proteome_metadata.tsv \
	 data/processed/subsets/class2_HGT_subset.csv

#------------ create subset fasta ---------------------

#python scripts_backup/08_plots/subset_fasta.py \
	 data/processed/subsets/class2_HGT_subset.csv \
	 results/alignments/class2_aligned.bak.fasta \
	 results/alignments/subset_alignments/class2_aligned_subset_HGT.fasta

#python scripts_backup/08_plots/subset_fasta.py data/processed/subsets/HGT_subtree_organisms.csv results/alignments/class1_aligned.bak.fasta results/alignments/subset_alignments/class1_aligned_subset_HGT.fasta


#----------- does all of the above together -----------------------------------------------

python scripts_backup/08_plots/super_subset.py data/processed/subsets/class2_subtree3B.txt \
	data/processed/proteome_metadata.tsv \
	results/alignments/class2_aligned.bak.fasta \
	results/alignments/subset_alignments/class2_aligned_subtree3B.fasta \
	data/processed/subsets/csv/class2_subtree3B.csv
#------------ run iqtree -------------------------------------------------------------


iqtree -s results/alignments/subset_alignments/class2_aligned_subtree3B.fasta \
        --prefix results/trees/class2/subtrees/class2_subtree3B \
	-mem 64G \
        -ntmax $SLURM_CPUS_PER_TASK \
        -nt AUTO \
        -mset Q.pfam -bb 1000 -alrt 1000 



#sub trees -------------------------------

Rscript scripts_backup/08_plots/tree_plot_unrooted.R data/processed/proteome_metadata.tsv results/trees/class2/subtrees/class2_subtree3B.treefile results/figures/unrooted/subtrees/unrooted_tree_class2_subtree3B.png


Rscript scripts_backup/08_plots/tree_plot_circular.R data/processed/proteome_metadata.tsv results/trees/class2/subtrees/class2_subtree3B.treefile results/figures/circular/subtrees/circular_tree_class2_subtree3B.png



echo "IQTree creation Script Finished at $(date)."
