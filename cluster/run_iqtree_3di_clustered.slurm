#!/bin/bash
#SBATCH --job-name=3di_iqtree
#SBATCH --out=slurm_out/iqtree_3di_tree_%j.out
#SBATCH --time=7-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=10


echo "SLURM script Starting on $(hostname) at $(date)..."

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

RESULT="$PROJECT_ROOT/results"

echo "Starting iqtree script..."


iqtree -s "$RESULT/alignments/class2_3di_clustered_aligned.fasta" \
        --prefix "$RESULT/trees/class2/3di_clustered" \
        -st AA \
		-mem 120G \
        -m "$PROJECT_ROOT/data/external/Q.3Di.LLM.m+G4" \
	 	-nt AUTO \
        -ntmax $SLURM_CPUS_PER_TASK \
	 	-bb 1000 

echo "IQTree creation Script Finished at $(date)."
