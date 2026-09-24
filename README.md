# FumaraseEvolution

Code, analysis pipeline, and selected results for the paper:

> Miller O, Jacobi D, Carmel L, Pines O. **Evolution and metabolite signaling of fumarases (class I and II): origin, function, and subcellular targeting.** *Molecular Biology and Evolution*. 2026. https://doi.org/10.1093/molbev/msag240

- **Paper:** https://doi.org/10.1093/molbev/msag240
- **Supplementary data (full trees and additional material):** https://data.mendeley.com/datasets/4kvht65rtb/1

## Overview

We searched ~21,000 UniProt reference proteomes for class I and class II fumarases using MMseqs2, then combined sequence-based phylogenetics, 3Di structure-based phylogenetics, conservation analysis, pairwise identity, and ESM-2 protein language model embeddings to study their origins across the tree of life.

This repository contains the scripts used at each step, the SLURM job scripts for the heavy computations, and the generated figures and statistics. Large intermediate files and the complete set of trees are hosted on Mendeley Data (link above).

## Repository structure

```
FumaraseEvolution/
├── scripts/      # Analysis pipeline (steps 01-08)
├── cluster/      # SLURM scripts for the heavy steps
├── envs/         # Conda environment specification
└── results/      # Generated outputs
```

### Pipeline (`scripts/`)

Run the steps in numerical order.

| Step | Directory | Description |
|------|-----------|-------------|
| 1 | `scripts/01_search/` | Search whole proteomes from the UniProt reference database for fumarase sequences with MMseqs2. |
| 2 | `scripts/02_align/` | Align sequences with MAFFT and compute per-position conservation (Shannon entropy). |
| 3 | `scripts/03_3di/` | Translate amino-acid sequences to 3Di structural sequences with the ProstT5 encoder-decoder model. |
| 4 | `scripts/04_identity/` | Compute the all-vs-all sequence identity matrix for class II and generate descriptive statistics. |
| 5 | `scripts/05_structure/` | Explanations and scripts for recreating the structural figures (ChimeraX). |
| 6 | `scripts/06_esm/` | Generate ESM-2 embeddings and mean-pool them per sequence. |
| 7 | `scripts/07_trees/` | Parse FASTA files, run IQ-TREE2 on amino-acid and 3Di sequences, and process the resulting trees. |
| 8 | `scripts/08_plots/` | Plot trees, ESM dimensionality reductions (PCA, t-SNE, UMAP), and conservation scores. |

### Cluster scripts (`cluster/`)

`.slurm` scripts for the computationally heavy steps. Submit them to a SLURM scheduler (e.g. `sbatch cluster/<script>.slurm`); each runs the corresponding pipeline step. Required package versions are noted inside the scripts.

### Results (`results/`)

| Directory | Contents |
|-----------|----------|
| `embeddings/` | HTML and PNG files of PCA, t-SNE, and UMAP projections for class I and class II sequences. |
| `figures/` | Generated plots, including the identity violin plots. |
| `stats/` | CSV tables: entropy baselines, identity statistics, Kruskal-Wallis test results, and summary tables. |
| `structure/` | ChimeraX session files (`.cxs`) for structure visualization. |

The complete supplementary material (including the full set of trees) is available on [Mendeley Data](https://data.mendeley.com/datasets/4kvht65rtb/1).

## Installation

A conda environment with all required packages and exact versions is provided in `envs/environment.yml`. From the repository root:

```bash
conda env create -f envs/environment.yml
conda activate fumaraseevo
```

Activate the environment before running any of the Python or shell scripts.

## Data

All proteomes were obtained from the [UniProt reference proteomes](https://www.uniprot.org/proteomes) database. Reference sequences for the searches are the *E. coli* K-12 fumarases: **P0AC33** (FumA, class I) and **P05042** (FumC, class II). No new experimental data were generated in this study.

## Citation

If you use this code or data, please cite the paper:

```bibtex
@article{miller2026fumarases,
  author  = {Miller, Ophir and Jacobi, Dror and Carmel, Liran and Pines, Ophry},
  title   = {Evolution and metabolite signaling of fumarases (class {I} and {II}): origin, function, and subcellular targeting},
  journal = {Molecular Biology and Evolution},
  year    = {2026},
  pages   = {msag240},
  doi     = {10.1093/molbev/msag240}
}
```

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
