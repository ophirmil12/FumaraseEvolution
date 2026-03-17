#!/bin/bash
#SBATCH --job-name=mafft_align
#SBATCH --out=slurm_out/mafft_%j.out
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --killable

# Note that this is not the regular MAFFT run, but a script that runs it on the clustered representative sequences.
# Here we also align the 3Di sequences using the custom 3Di matrix.


set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

MAFFT="$PROJECT_ROOT/scripts/02_align/run_mafft.sh"
DATA="$PROJECT_ROOT/data"
RESULTS="$PROJECT_ROOT/results/alignments"
THREADS="$SLURM_CPUS_PER_TASK"

MAFFT_3DI_MATRIX="$PROJECT_ROOT/data/external/mat3di.out"       # from https://www.nature.com/articles/s41587-023-01773-0#MOESM1

# --- Class I 3Di clustered Alignment ---
echo "Running MAFFT alignment for Class I 3Di sequences..."

bash "$MAFFT" \
    "$DATA/processed/clustered/class1_3di_clustered_rep_seq.fasta" \
    "$RESULTS/class1_3di_clustered_aligned.fasta" \
    --threads "$THREADS" \
    --aamatrix "$MAFFT_3DI_MATRIX"

# --- Class II 3Di clustered Alignment ---
echo "Running MAFFT alignment for Class II 3Di sequences..."
bash "$MAFFT" \
    "$DATA/processed/clustered/class2_3di_clustered_rep_seq.fasta" \
    "$RESULTS/class2_3di_clustered_aligned.fasta" \
    --threads "$THREADS" \
    --aamatrix "$MAFFT_3DI_MATRIX"

# -- Class I AA clustered Alignment ---
echo "Running MAFFT alignment for Class I AA sequences..."
bash "$MAFFT" \
    "$DATA/processed/clustered/class1_aa_clustered_rep_seq.fasta" \
    "$RESULTS/class1_aa_clustered_aligned.fasta" \
    --threads "$THREADS"

# -- Class II AA clustered Alignment ---
echo "Running MAFFT alignment for Class II AA sequences..."
bash "$MAFFT" \
    "$DATA/processed/clustered/class2_aa_clustered_rep_seq.fasta" \
    "$RESULTS/class2_aa_clustered_aligned.fasta" \
    --threads "$THREADS"


echo "All alignments complete."