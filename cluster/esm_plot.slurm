#!/bin/bash
#SBATCH --job-name=esm_embeddings
#SBATCH --out=slurm_out/esm_embeddings_%j.out
#SBATCH --killable
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gg:g4:1
#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

export HF_HOME="$PROJECT_ROOT/.hf_cache"
export HF_HUB_DISABLE_XET=1
export TORCH_HOME="$PROJECT_ROOT/.torch_cache"

SCRIPTS="$PROJECT_ROOT/scripts"
DATA="$PROJECT_ROOT/data/processed"
EMB="$PROJECT_ROOT/data/processed/embeddings"
mkdir -p "$EMB"

echo "========================================"
echo " ESM-2 Embeddings + Dimensionality Reduction"
echo " $(date)"
echo "========================================"

# -------------------------------------------------------------------
# Step 1: Compute mean embeddings for both classes and save as .pt
# -------------------------------------------------------------------
echo ""
echo "[1/3] Computing ESM-2 mean embeddings..."

python "$SCRIPTS/06_esm/esm_factory.py" \
    --fasta  "$DATA/class1_sequences.fasta" \
    --output "$EMB/class1_embeddings.pt" \
    --batch-size 32

python "$SCRIPTS/06_esm/esm_factory.py" \
    --fasta  "$DATA/class2_sequences.fasta" \
    --output "$EMB/class2_embeddings.pt" \
    --batch-size 32

echo "  Embeddings saved."

# -------------------------------------------------------------------
# Step 2: Dimensionality reduction + plotting (all methods, 2D)
# -------------------------------------------------------------------
echo ""
echo "[2/3] Plotting — Class I..."

for METHOD in tsne umap pca; do
    echo "  Class I | $METHOD | 2D"
    python "$SCRIPTS/08_plots/esm_plot.py" \
        --embeddings "$EMB/class1_embeddings.pt" \
        --fasta      "$DATA/class1_sequences.fasta" \
        --method     "$METHOD" \
        --dims       2 \
        --label      "class1" \
        --outdir     "$PROJECT_ROOT/results/embeddings"
done

echo ""
echo "[3/3] Plotting — Class II..."

for METHOD in tsne umap pca; do
    echo "  Class II | $METHOD | 2D"
    python "$SCRIPTS/08_plots/esm_plot.py" \
        --embeddings "$EMB/class2_embeddings.pt" \
        --fasta      "$DATA/class2_sequences.fasta" \
        --method     "$METHOD" \
        --dims       2 \
        --label      "class2" \
        --outdir     "$PROJECT_ROOT/results/embeddings"
done

echo ""
echo "========================================"
echo " Done. $(date)"
echo " Outputs: $PROJECT_ROOT/results/embeddings/"
echo "========================================"
