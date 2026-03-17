#!/bin/bash
#SBATCH --job-name=allvsall_class2
#SBATCH --out=slurm_out/allvsall_%j.out
#SBATCH --killable
#SBATCH --constraint=avx2           # Ensure we run on nodes with AVX2 support for MMseqs2
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

# --- MMseqs2 ---
MMSEQS="$PROJECT_ROOT/tools/mmseqs_avx2/bin/mmseqs"

if [ ! -f "$MMSEQS" ]; then
    mkdir -p "$PROJECT_ROOT/tools"
    wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
         -O /tmp/mmseqs-avx2.tar.gz
    tar xf /tmp/mmseqs-avx2.tar.gz -C "$PROJECT_ROOT/tools" \
        --transform 's/mmseqs/mmseqs_avx2/'
    rm /tmp/mmseqs-avx2.tar.gz
    echo "MMseqs2 AVX2 binary installed."
fi

echo "Using MMseqs2: $("$MMSEQS" version)"
echo "Compute node CPU: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"

# --- Paths ---
FASTA="$PROJECT_ROOT/data/processed/class2_sequences.fasta"
OUT="$PROJECT_ROOT/data/processed/class2_allvsall.tsv"
TMP="$PROJECT_ROOT/data/raw/tmp/allvsall_tmp"
N_SEQS=$(grep -c "^>" "$FASTA")

echo "========================================"
echo " MMseqs2 all-vs-all Class II"
echo " Sequences: $N_SEQS"
echo " Output:    $OUT"
echo " Threads:   $SLURM_CPUS_PER_TASK"
echo "========================================"

mkdir -p "$TMP"

MERGED="$PROJECT_ROOT/data/processed/class2_allvsall_full.tsv"

"$MMSEQS" easy-search \
    "$FASTA" \
    "$FASTA" \
    "$MERGED" \
    "$TMP" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --max-seqs "$N_SEQS" \
    -s 7.5 \
    --k-score 80 \
    -e 100 \
    --min-seq-id 0 \
    --format-output "query,target,fident,alnlen,qlen,tlen"

# Filter diagonal and duplicates, write final output
awk -F'\t' '$1 < $2' "$MERGED" > "$OUT"

# Remove unfiltered file
rm "$MERGED"

echo ""
echo "Done. Lines written: $(wc -l < "$OUT")"
echo "File size:           $(du -sh "$OUT" | cut -f1)"

rm -rf "$TMP"
echo "Tmp cleaned."