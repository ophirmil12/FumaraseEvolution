import re
import csv
import os
import sys


def extract_iqtree_metadata(iqtree_path):
    """
    Parses the .iqtree file using specific formatting for LogL and Information Criteria.
    """
    metadata = []
    if not os.path.exists(iqtree_path):
        print(f"Warning: Metadata file {iqtree_path} not found.")
        return metadata

    try:
        with open(iqtree_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
            # 1. Model Selection
            model = re.search(r"Best-fit model: (.*) chosen", content)
            if model:
                metadata.append({"Metric": "Best-fit Model", "Value": model.group(1)})

            # 2. Likelihood & Information Criteria
            # Matches the number but stops before the opening parenthesis (s.e.)
            logl = re.search(r"Log-likelihood of the tree: ([-0-9.]+) ", content)
            
            # Captures values for AIC, AICc, and BIC
            aic  = re.search(r"Akaike information criterion \(AIC\) score: ([-0-9.]+)", content)
            aicc = re.search(r"Corrected Akaike information criterion \(AICc\) score: ([-0-9.]+)", content)
            bic  = re.search(r"Bayesian information criterion \(BIC\) score: ([-0-9.]+)", content)
            
            # Number of free parameters
            params = re.search(r"Number of free parameters.*: (\d+)", content)

            if logl:   metadata.append({"Metric": "Log-Likelihood (LogL)", "Value": logl.group(1)})
            if params: metadata.append({"Metric": "Free Parameters", "Value": params.group(1)})
            if aic:    metadata.append({"Metric": "AIC Score", "Value": aic.group(1)})
            if aicc:   metadata.append({"Metric": "AICc Score", "Value": aicc.group(1)})
            if bic:    metadata.append({"Metric": "BIC Score", "Value": bic.group(1)})

            # 3. Check for the Consensus Tree Note
            if "Consensus tree has higher likelihood than ML tree" in content:
                metadata.append({"Metric": "Note", "Value": "Consensus tree used (higher likelihood than ML)"})

            # 4. Alignment Informativeness (Still useful for SI)
            const_sites = re.search(r"Number of constant sites: (\d+)", content)
            info_sites  = re.search(r"Number of parsimony informative sites: (\d+)", content)
            
            if const_sites: metadata.append({"Metric": "Constant Sites", "Value": const_sites.group(1)})
            if info_sites:  metadata.append({"Metric": "Parsimony Informative Sites", "Value": info_sites.group(1)})

    except Exception as e:
        print(f"Error reading .iqtree file: {e}")
        
    return metadata

def append_metadata_to_csv(csv_path, metadata):
    """
    Appends the metadata rows to the existing CSV file.
    """
    if not metadata:
        return

    try:
        with open(csv_path, 'a', newline='') as csvfile:
            fieldnames = ['Metric', 'Value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            # Add a separator row for clarity
            writer.writerow({"Metric": "---", "Value": "---"})
            writer.writerow({"Metric": "IQ-TREE RUN METADATA", "Value": ""})
            
            writer.writerows(metadata)
        print(f"Metadata successfully appended to {csv_path}")
    except Exception as e:
        print(f"Error writing to CSV: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python tree_stats_extended.py <stats_csv_file> <iqtree_report_file>")
    else:
        csv_input = sys.argv[1]
        iqtree_input = sys.argv[2]
        
        run_data = extract_iqtree_metadata(iqtree_input)
        append_metadata_to_csv(csv_input, run_data)