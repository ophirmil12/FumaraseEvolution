#!/bin/bash
#SBATCH --job-name=subtree23B
#SBATCH --out=slurm_out/subtree23B_%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=6




echo "SLURM script Starting on $(hostname) at $(date)..."

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

RESULT="$PROJECT_ROOT/results"
DATA="$PROJECT_ROOT/data""
SCRIPTS="$PROJECT_ROOT/scripts"

echo "Starting script..."



#----------- creates subtree fro -----------------------------------------------
echo "creating alignment subset file..."
#for syntax, see super_subset.py file in scripts.

python "$SCRIPTS/08_plots/super_subset.py" \
	"$DATA/processed/subsets/class2_subtree3B.txt" \
	"$DATA/processed/proteome_metadata.tsv" \
	"$RESULT/alignments/class2_aligned.bak.fasta" \
	"$RESULT/alignments/subset_alignments/class2_aligned_subtree3B.fasta" \
	"$DATA/processed/subsets/csv/class2_subtree3B.csv"
#------------ run iqtree -------------------------------------------------------------

echo "Starting iqtree script..."
#

iqtree -s "$RESULT/alignments/subset_alignments/class2_aligned_subtree3B.fasta" \
        --prefix "$RESULT/trees/class2/subtrees/class2_subtree3B" \
		-mem 64G \
        -ntmax $SLURM_CPUS_PER_TASK \
        -nt AUTO \
        -mset Q.pfam -bb 1000 -alrt 1000 



#sub trees -------------------------------

echo "Plotting subtrees..."

Rscript "$SCRIPTS/08_plots/tree_plot_unrooted.R" \
		"$DATA/processed/proteome_metadata.tsv" \
		"$RESULT/trees/class2/subtrees/class2_subtree3B.treefile" \
		"$RESULT/figures/unrooted/subtrees/unrooted_tree_class2_subtree3B.png"

Rscript "$SCRIPTS/08_plots/tree_plot_circular.R" \
		"$DATA/processed/proteome_metadata.tsv" \
		"$RESULT/trees/class2/subtrees/class2_subtree3B.treefile" \
		"$RESULT/figures/circular/subtrees/circular_tree_class2_subtree3B.png"



echo "Tree and Plot creation finished at $(date)."
