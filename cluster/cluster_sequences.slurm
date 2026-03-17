#!/bin/bash
#SBATCH --job-name=cluster_all
#SBATCH --out=slurm_out/clustering_%j.out
#SBATCH --killable
#SBATCH --constraint=avx2           # Ensure we run on nodes with AVX2 support for MMseqs2
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END
#SBATCH --mail-user=ophir.miller@mail.huji.ac.il

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="/sci/labs/itamarsi/ophirmil12/FumaraseEvolution"
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

# --- MMseqs2 Setup ---
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

# --- File Definitions ---
declare -A INPUTS=(
    ["class1_aa"]="$PROJECT_ROOT/data/processed/class1_sequences.fasta"
    ["class2_aa"]="$PROJECT_ROOT/data/processed/class2_sequences.fasta"
    ["class1_3di"]="$PROJECT_ROOT/data/processed/class1_3di.fasta"
    ["class2_3di"]="$PROJECT_ROOT/data/processed/class2_3di.fasta"
)

declare -A SEQ_ID_THRESHOLDS=(
    ["class1_aa"]=0.94
    ["class2_aa"]=0.94
    ["class1_3di"]=0.96
    ["class2_3di"]=0.96
)

# --- Run Clustering ---
for prefix in "${!INPUTS[@]}"; do
    FASTA="${INPUTS[$prefix]}"
    
    if [ ! -f "$FASTA" ]; then
        echo "Skipping $prefix: File $FASTA not found."
        continue
    fi

    OUT_PREFIX="$PROJECT_ROOT/data/processed/${prefix}_clustered"
    TMP_DIR="$PROJECT_ROOT/data/raw/tmp/cluster_tmp_${prefix}"
    mkdir -p "$TMP_DIR"

    echo "========================================"
    echo " Clustering: $prefix"
    echo " Input:      $FASTA"
    echo " Output:     ${OUT_PREFIX}_rep_seq.fasta"
    echo "========================================"

    # Run easy-cluster
    # --min-seq-id 0.95 ensures representative sequences are highly unique
    # -c 0.8 ensures the alignment covers at least 80% of both sequences
    "$MMSEQS" easy-cluster \
        "$FASTA" \
        "$OUT_PREFIX" \
        "$TMP_DIR" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --min-seq-id "${SEQ_ID_THRESHOLDS[$prefix]}" \
        -c 0.8

    # Clean up temp files for this run
    rm -rf "$TMP_DIR"
    echo "$prefix clustering finished."
    echo ""
done


# --- Summary Table ---
cd ./data/processed/

printf "%-20s %-15s %-15s %-15s\n" "Dataset" "Before" "After" "% Retained"
echo "-----------------------------------------------------------------"

for rep in *_clustered_rep_seq.fasta; do
    # Extract the prefix (e.g., class1_aa or class2_3di)
    prefix=${rep%_clustered_rep_seq.fasta}
    
    # Map the prefix back to your exact original filenames
    if [[ "$prefix" == *_aa ]]; then
        # Converts "class1_aa" to "class1_sequences.fasta"
        base=${prefix%_aa}
        orig="${base}_sequences.fasta"
    else
        # Converts "class1_3di" to "class1_3di.fasta"
        orig="${prefix}.fasta"
    fi

    # Skip if for some reason the file still isn't found
    if [ ! -f "$orig" ]; then
        continue
    fi

    # Count the sequences
    before=$(grep -c "^>" "$orig")
    after=$(grep -c "^>" "$rep")
    
    # Calculate percentage
    pct=$(awk "BEGIN {printf \"%.1f%%\", ($after/$before)*100}")
    
    # Print row
    printf "%-20s %-15s %-15s %-15s\n" "$prefix" "$before" "$after" "$pct"
done


echo "All clustering runs completed successfully."