#!/bin/bash
#SBATCH --job-name=predict_3di
#SBATCH --out=slurm_out/predict_3di_%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

#SBATCH --gres=gg:g4:1			# Phoenix
# #SBATCH --gres=gpu:1

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

export HF_HOME="$PROJECT_ROOT/.hf_cache"
export HF_HUB_DISABLE_XET=1

python "$PROJECT_ROOT/scripts/03_3di/predict_3di.py" \
    --class1 "$PROJECT_ROOT/data/processed/class1_sequences.fasta" \
    --class2 "$PROJECT_ROOT/data/processed/class2_sequences.fasta" \
    --outdir "$PROJECT_ROOT/data/processed" \
    --batch-size 16

echo "3Di prediction complete."
echo "Next: align the 3Di FASTAs with run_mafft.sh"
