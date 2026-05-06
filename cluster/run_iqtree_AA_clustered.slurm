#!/bin/bash
#SBATCH --job-name=AA_iqtree_class2_qpfram_r8
#SBATCH --out=slurm_out/iqtree_class2_AA_clustered2k_qpfam_r8%j.out
#SBATCH --time=:700:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12



echo "SLURM script Starting on $(hostname) at $(date)..."

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

RESULT="$PROJECT_ROOT/results"

echo "Starting iqtree script..."


iqtree2 -s "$RESULT/alignments/class2_aa_clustered_aligned_80.fasta" \
--prefix "$RESULT/trees/class2/backups/AA_clustered2k_qpfam_r8" \
-mem 128G -ntmax 12 -nt AUTO -bb 1000 -m Q.pfam+R8 -redo

echo "IQTree creation Script Finished at $(date)."
