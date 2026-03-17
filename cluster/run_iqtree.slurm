#!/bin/bash
#SBATCH --job-name=fumarase_iqtree_dummy
#SBATCH --out=/sci/labs/itamarsi/dror.jacobi/FumaraseEvolution/slurm_out/iqtree_AA_%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
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


iqtree -s results/alignments/class1_aligned.fasta \
        --prefix results/trees/class1/AA_ \
        -nt AUTO \
        -bb 1000 \
        -m MFP


#iqtree -s results/alignment/class2_aligned.fasta \
#        --prefix results/trees/class2/AA_ \
#        -nt AUTO \
#        -bb 1000 \
#        -m MFP


echo "IQTree creation Script Finished at $(date)."
