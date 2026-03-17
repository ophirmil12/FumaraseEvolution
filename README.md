# FumaraseEvolution

This repository contains the code, analysis pipeline, and partial results for the paper: **"Evolution of Fumarases (Class I and II): origin, function and subcellular targeting"** by O. Miller, L. Carmel, and O. Pines.

TODO: add Dror when Ophry adds him to the author list.


## 📂 Repository Structure

The repository is organized into scripts, configuration files, and generated results:

### Code & Pipeline (`scripts/`)
The analysis pipeline is broken down into sequential steps:
* **`scripts/01_search/`**: Searches for fumarase sequences in whole proteomes from the UniProt reference database using MMseqs2.
* **`scripts/02_align/`**: Aligns the sequences using MAFFT and handles conservation calculations (entropy).
* **`scripts/03_3di/`**: Predicts 3Di sequences via ProstT5's encoder-decoder translation (AA $\rightarrow$ 3Di).
* **`scripts/04_identity/`**: Calculates the all-vs-all sequence identity matrix for class II and generates descriptive statistics.
* **`scripts/05_structure/`**: General explanations and scripts on how to recreate the structural figures.
* **`scripts/06_esm/`**: Handles generating embeddings and mean-pooling of the sequences using the ESM-2 PLM.
* **`scripts/07_trees/`**: Parse FASTA files, runs IQ-TREE2 for regular and 3Di sequences, and handles the phylogenetic tree results.
* **`scripts/08_plots/`**: Scripts for plotting the results of the analysis, including trees, ESM dimensionality reductions (PCA, t-SNE, UMAP), and conservation scores.

### Results (`results/`)
Contains the generated outputs from the pipeline:
* **`embeddings/`**: HTML and PNG files of PCA, t-SNE, and UMAP dimensionality reductions for Class I and Class II sequences.
* **`figures/`**: Generated plots, including identity violin plots.
* **`stats/`**: Tabular data (`.csv`) including entropy baselines, identity statistics, Kruskal-Wallis test results, and summary tables.
* **`structure/`**: ChimeraX session files (`.cxs`) for protein structure visualization.

### High-Performance Computing (`cluster/`)
General scripts for running the computationally heavy steps of the analysis on a cluster are available in the `cluster/` directory. These `.slurm` scripts are designed to be submitted to a SLURM job scheduler and will execute the corresponding steps of the analysis pipeline. Specific required package versions are noted within the scripts.

---

## ⚙️ Installation & Environment setup

In the `envs/` directory, we provide a conda environment specification (`environment.yml`) that includes all the necessary packages and exact versions to reproduce the analysis. 

To create the environment, run the following command from the root of the repository:

```bash
conda env create -f envs/environment.yml
```

This will build a conda environment with all the required dependencies. Make sure to activate the environment before running any of the Python or Shell scripts:

```bash
conda activate fumaraseevo
```

---

## 📝 Citation
If you use this code or data in your research, please cite our paper:
```bibtex
TODO: Add the full citation for the paper here once available.
```

## 📄 License

TODO: Add license information here.