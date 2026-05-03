#!/bin/bash
#SBATCH --job-name=2_AA_clustered2k_qpfam_r8_finishing_1100
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/iqtree_class2_AA_clustered2k_qpfam_r8_finishing_1100_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12
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

iqtree2 -s results/alignments/class2_aa_clustered_aligned_80.fasta \
--prefix results/trees/class2/backups/AA_clustered2k_qpfam_r8 \
-mem 128G -ntmax 12 -nt AUTO -bb 1000 -m Q.pfam+R8

echo "IQTree creation Script Finished at $(date)."