#!/bin/bash
#SBATCH --job-name=fumarase_iqtree_dummy
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/iqtree_3di_%j.out
#SBATCH --time=0:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
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


iqtree -s data/processed/class1_3di_aligned_dummy.fasta \
         --prefix results/trees/class1/3Di_ \
         -st AA \
         -m scripts/07_trees/Q.3Di.LLM.m+G4 \
         -ntmax $SLURM_CPUS_PER_TASK \
         -bb 1000

#iqtree -s data/processed/class2_3di_aligned.fasta \
#        --prefix results/trees/class2/3Di_ \
#        -st AA \
#        -m scripts/07_trees/Q.3Di.LLM.m+G4 \
#        -ntmax $SLURM_CPUS_PER_TASK \
#        -bb 1000


echo "IQTree creation Script Finished at $(date)."
