"""
Universal Fumarase Subsetter-
----------------------------------------
Accepts UIDs (numbers), UPIDs (UP...), or a CSV.
Outputs a clean FASTA subset.

Usage: 
    python super_subset.py <mixed_ids.txt> <master_metadata.tsv> <input.fasta> <output.fasta>
"""
import sys
import os
import pandas as pd
from Bio import SeqIO

def resolve_ids(input_path, master_df):
    # 1. Load the raw input
    if input_path.endswith('.csv'):
        df = pd.read_csv(input_path)
        raw_list = df['organism_id'].astype(str).tolist()
    else:
        with open(input_path, 'r') as f:
            raw_list = [line.strip() for line in f if line.strip()]

    # 2. Separate UP IDs and Taxonomy IDs
    up_targets = {i for i in raw_list if i.startswith('UP')}
    tax_targets = {i for i in raw_list if not i.startswith('UP')}

    # 3. Map Taxonomy IDs to UP IDs using the master TSV
    if tax_targets:
        print(f"Mapping {len(tax_targets)} Taxonomy IDs to Proteomes...")
        master_df['organism_id'] = master_df['organism_id'].astype(str)
        matches = master_df[master_df['organism_id'].isin(tax_targets)]
        mapped_ups = set(matches['proteome_id'].dropna().astype(str))
        up_targets.update(mapped_ups)

    print(f"Total unique Proteome IDs to extract: {len(up_targets)}")
    return up_targets

def get_rank(lineage, index):
    """Parses the lineage string to extract Domain (1) or Phylum (2)."""
    if pd.isna(lineage): return "Unknown"
    parts = [p.strip() for p in lineage.split(',')]
    return parts[index] if len(parts) > index else "Unknown"

def save_final_metadata(found_ids, master_df, id_file, output_csv):
    """Filters master metadata for the IDs found in the FASTA and saves to CSV."""
    print(f"Creating final metadata for {len(found_ids)} sequences...")
    
    # Filter for the IDs that made the cut
    final_df = master_df[master_df['proteome_id'].isin(found_ids)].copy()
    
    # Parse lineage for Domain and Phylum
    final_df['domain'] = final_df['taxonomic_lineage'].apply(lambda x: get_rank(x, 1))
    final_df['phylum'] = final_df['taxonomic_lineage'].apply(lambda x: get_rank(x, 2))
    
    # If your original input was a CSV, try to bring back the 'role' column
    if id_file.endswith('.csv'):
        user_df = pd.read_csv(id_file)
        if 'role' in user_df.columns:
            user_df['organism_id'] = user_df['organism_id'].astype(str)
            final_df['organism_id'] = final_df['organism_id'].astype(str)
            final_df = pd.merge(final_df, user_df[['organism_id', 'role']], on='organism_id', how='left')
    
    # Standardize columns and save
    cols = ['organism_id', 'proteome_id', 'organism', 'domain', 'phylum', 'role']
    # Ensure columns exist before selecting
    existing_cols = [c for c in cols if c in final_df.columns]
    final_df[existing_cols].to_csv(output_csv, index=False)


def main(id_file, master_file, fasta_in, fasta_out, metadata_out):
    # Create directory if missing
    os.makedirs(os.path.dirname(fasta_out) or '.', exist_ok=True)

    # Load master metadata
    master = pd.read_csv(master_file, sep='\t')
    
    # Get the unified list of UP IDs
    final_up_ids = resolve_ids(id_file, master)

    # Filter the FASTA
    found_ids = set()
    count = 0
    with open(fasta_out, 'w') as out_h:
        for record in SeqIO.parse(fasta_in, "fasta"):
            matched_id = next((up for up in final_up_ids if up in record.description), None)
            if matched_id:
                SeqIO.write(record, out_h, "fasta")
                found_ids.add(matched_id)
                count += 1
    
    save_final_metadata(found_ids, master, id_file, metadata_out)
    print(f"Success! Created {fasta_out} with {count} sequences.")

if __name__ == "__main__":
    if len(sys.argv) == 6:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print("Usage: python super_subset.py <ids> <master_tsv> <in_fasta> <out_fasta> <out_metadata_csv>")