import re
import numpy as np
import sys
import os
import csv  # Added for CSV support

def analyze_tree(input_path, output_path=None):
    try:
        with open(input_path, 'r') as f:
            content = f.read()

        # Extract Newick string
        tree_match = re.search(r'(\(.*\);)', content, re.DOTALL)
        if not tree_match:
            print(f"Error: No Newick tree found in {input_path}")
            return
        
        tree_string = tree_match.group(1)

        # 1. Extract Taxa (Leaf) Count
        taxa_names = re.findall(r'([(),]|[A-Za-z0-9_|.-]+):', tree_string)
        taxa_count = len([x for x in taxa_names if x not in '(),'])

        # 2. Extract Support Values
        raw_supports = re.findall(r'\)([0-9./]+):', tree_string)
        bootstraps = []
        for b in raw_supports:
            val = b.split('/')[-1]
            if val:
                bootstraps.append(float(val))

        # 3. Extract Branch Lengths 
        branch_lengths = [float(l) for l in re.findall(r'(?<!\))[:]([0-9.eE-]+)', tree_string)]

        if not bootstraps:
            print(f"Error: No support values detected. Verify the tree format.")
            return

        # Statistical Calculations
        total = len(bootstraps)
        high = sum(1 for x in bootstraps if x >= 95)
        mod  = sum(1 for x in bootstraps if 70 <= x < 95)
        low  = sum(1 for x in bootstraps if x < 70)

        # Organized data as a list of dictionaries for clean CSV writing
        stats = [
            {"Metric": "Source File", "Value": os.path.basename(input_path)},
            {"Metric": "Total Taxa (Leaves)", "Value": taxa_count},
            {"Metric": "Total Internal Nodes", "Value": total},
            {"Metric": "Mean Support", "Value": round(np.mean(bootstraps), 2)},
            {"Metric": "Median Support", "Value": round(np.median(bootstraps), 2)},
            {"Metric": "Std Deviation", "Value": round(np.std(bootstraps), 2)},
            {"Metric": "Min Support", "Value": np.min(bootstraps)},
            {"Metric": "Max Support", "Value": np.max(bootstraps)},
            {"Metric": "Mean Branch Length", "Value": round(np.mean(branch_lengths), 4)},
            {"Metric": "Total Tree Depth", "Value": round(np.sum(branch_lengths), 4)},
            {"Metric": "High Support (>=95)", "Value": f"{high} ({(high/total)*100:.1f}%)"},
            {"Metric": "Moderate Support (70-94)", "Value": f"{mod} ({(mod/total)*100:.1f}%)"},
            {"Metric": "Low Support (<70)", "Value": f"{low} ({(low/total)*100:.1f}%)"},
        ]

        # Determine output location (changed default extension to .csv)
        if not output_path:
            output_path = os.path.splitext(input_path)[0] + "_stats.csv"

        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Write CSV File
        with open(output_path, 'w', newline='') as csvfile:
            fieldnames = ['Metric', 'Value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            writer.writerows(stats)
        
        print(f"\nAnalysis Complete.")
        print(f"Input:  {input_path}")
        print(f"Output: {output_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tree_exporter.py <input_file> [output_path]")
    elif len(sys.argv) == 2:
        analyze_tree(sys.argv[1])
    else:
        analyze_tree(sys.argv[1], sys.argv[2])