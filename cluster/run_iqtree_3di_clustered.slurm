#!/bin/bash
#SBATCH --job-name=1_3Di_2k_bb
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/iqtree_1_3Di_2k_bb_%j.out
#SBATCH --time=7-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=10
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


iqtree -s results/alignments/class2_3di_clustered_aligned.fasta \
         --prefix results/trees/class2/3di_clustered \
         -st AA \
	 -mem 128G \
         -m data/external/Q.3Di.LLM.m+G4 \
	 -nt AUTO \
         -ntmax $SLURM_CPUS_PER_TASK \
	 -bb 1000 

echo "IQTree creation Script Finished at $(date)."
