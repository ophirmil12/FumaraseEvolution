#!/bin/bash
#SBATCH --job-name=cluster_all
#SBATCH --out=slurm_out/clustering_%j.out
#SBATCH --killable
#SBATCH --constraint=avx2           # Ensure we run on nodes with AVX2 support for MMseqs2
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

set -euo pipefail
mkdir -p slurm_out

PROJECT_ROOT="./FumaraseEvolution"
source ./miniforge3/etc/profile.d/conda.sh
conda activate project_env

mkdir -p "$PROJECT_ROOT/data/processed/clustered"

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
    ["class1_aa"]=0.90
    ["class2_aa"]=0.90
    ["class1_3di"]=0.96
    ["class2_3di"]=0.96
)

# Required organisms
REQUIRED_ORGANISMS=(
    "Escherichia coli"
    "Bacillus subtilis"
    "Mycobacterium tuberculosis"
    "Giardia intestinalis"
    "Trichomonas vaginalis"
    "Entamoeba histolytica"
    "Saccharomyces cerevisiae"
    "Arabidopsis thaliana"
    "Caenorhabditis elegans"
    "Homo sapiens"
    "Rattus norvegicus"
)


# --- Run Clustering ---
for prefix in "${!INPUTS[@]}"; do
    FASTA="${INPUTS[$prefix]}"
    
    if [ ! -f "$FASTA" ]; then
        echo "Skipping $prefix: File $FASTA not found."
        continue
    fi

    OUT_PREFIX="$PROJECT_ROOT/data/processed/clustered/${prefix}_clustered"
    REP_FASTA="${OUT_PREFIX}_rep_seq.fasta"
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

    # --- Ensure Required Organisms are included ---
    echo "Ensuring required organisms are in $prefix representative set..."
    for species in "${REQUIRED_ORGANISMS[@]}"; do
        # Check if species scientific name exists in the representative file
        if ! grep -qi "OS=$species" "$REP_FASTA"; then
            # If missing, extract it from the original file and append
            # This awk script extracts the full FASTA record for the first match
            echo "  [+] $species missing from reps. Adding from original..."
            awk -v sp="OS=$species" 'BEGIN{RS=">"; IGNORECASE=1} $0 ~ sp {print ">"$0; exit}' "$FASTA" >> "$REP_FASTA"
        else
            echo "  [✓] $species already present."
        fi
    done

    # Clean up temp files for this run
    rm -rf "$TMP_DIR"
    echo "$prefix clustering finished."
    echo ""
done


# --- Summary Table ---
cd "$PROJECT_ROOT/data/processed/"

mkdir -p "$PROJECT_ROOT/results/stats"
OUTPUT_FILE="$PROJECT_ROOT/results/stats/clustering_summary.txt"

{
    printf "# Clustering the sequences (3Di and AA) into representative sequences for each class.\n"
    printf "# This summary shows the number of sequences before and after clustering, and the percentage retained.\n\n"
    printf "%-20s %-15s %-15s %-15s\n" "Dataset" "Before" "After" "% Retained"
    echo "-----------------------------------------------------------------"

    # Note: Added 'clustered/' to the search path here
    for rep in clustered/*_clustered_rep_seq.fasta; do
        
        # Strip the directory path and the suffix to get the raw prefix
        filename=$(basename "$rep")
        prefix=${filename%_clustered_rep_seq.fasta}
        
        # Map the prefix back to your exact original filenames located in data/processed/
        if [[ "$prefix" == *_aa ]]; then
            base=${prefix%_aa}
            orig="${base}_sequences.fasta"
        else
            orig="${prefix}.fasta"
        fi

        # Skip if for some reason the file still isn't found
        if [ ! -f "$orig" ]; then
            echo "Warning: Could not find original file $orig for $prefix"
            continue
        fi

        # Count the sequences
        before=$(grep -c "^>" "$orig")
        after=$(grep -c "^>" "$rep")
        
        # Calculate percentage
        pct=$(awk "BEGIN {printf \"%.2f%%\", ($after/$before)*100}")
        
        # Print row
        printf "%-20s %-15s %-15s %-15s\n" "$prefix" "$before" "$after" "$pct"
    done
} | tee "$OUTPUT_FILE"

echo "Summary saved to: $OUTPUT_FILE"


echo "All clustering runs completed successfully."
echo "Next steps: Run MAFFT alignments on the clustered representative sequences (don't forget to use the --aamatrix for 3Di datasets!)."
echo "Then create the trees with the Q matrix for 3Di alignments."
