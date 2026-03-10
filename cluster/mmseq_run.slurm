#!/bin/bash
#SBATCH --job-name=fumarase_mmseqs
#SBATCH --output=slurm_out/mmseqs_%j.out
#SBATCH --error=slurm_out/mmseqs_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

set -euo pipefail

mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"

source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

# --- Run class I ---
bash "$PROJECT_ROOT/scripts/01_search/mmseq_run.sh" \
    "$PROJECT_ROOT/data/processed/proteome_metadata.tsv" \
    "$PROJECT_ROOT/data/external/ecoli_fumA.fasta" \
    --jobs "$SLURM_CPUS_PER_TASK"

# Move class I results to a separate directory
mkdir -p "$PROJECT_ROOT/data/raw/m8s_class1"
mv "$PROJECT_ROOT/data/raw/m8s/"*_output.m8 "$PROJECT_ROOT/data/raw/m8s_class1/"

# --- Run class II (m8s/ is now empty, nothing gets skipped) ---
bash "$PROJECT_ROOT/scripts/01_search/mmseq_run.sh" \
    "$PROJECT_ROOT/data/processed/proteome_metadata.tsv" \
    "$PROJECT_ROOT/data/external/ecoli_fumC.fasta" \
    --jobs "$SLURM_CPUS_PER_TASK"

mkdir -p "$PROJECT_ROOT/data/raw/m8s_class2"
mv "$PROJECT_ROOT/data/raw/m8s/"*_output.m8 "$PROJECT_ROOT/data/raw/m8s_class2/"

# Merge both into m8s/ — append class2 results into class1 files
for m8 in "$PROJECT_ROOT/data/raw/m8s_class2/"*_output.m8; do
    proteome_id=$(basename "$m8" _output.m8)
    cat "$PROJECT_ROOT/data/raw/m8s_class1/${proteome_id}_output.m8" "$m8" \
        > "$PROJECT_ROOT/data/raw/m8s/${proteome_id}_output.m8"
done

echo "MMseqs2 search complete."
