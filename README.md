# FumaraseEvolution

This repository contains code and partial results for the paper "Evolution of Fumarases (Class I and II): origin, function and subcellular
targeting", by O. Miller, L. Carmel, and O. Pines.

## Code
The code is organized as follows:
- config.py: contains configuration parameters for the analysis.
- scripts/01_search - search for fumarase sequences in whole proteomes from the UniProt reference database, using MMseqs2.
- scripts/02_align - align the sequences using MAFFT, and handles conservation calculations.
- scripts/03_3di - predicts 3Di sequences via ProstT5's encoder-decoder translation (AA $\rightarrow$ 3Di).
- ...   TODO: finish

## Cluster scripts
General scripts for running the analysis on a cluster are available in the `cluster` directory. These scripts are designed to be submitted to a job scheduler (e.g., SLURM) and will execute the corresponding steps of the analysis pipeline. We noted some specific package versions that are required for the scripts to run correctly.

## envs
In enviroment.yml, we provide a conda environment specification that includes the necessary packages and their versions to run the analysis. To create the environment, use the following command:

```bash
conda env create -f enviroment.yml
```
This will create a conda environment with all the required dependencies installed. Make sure to activate the environment before running any of the scripts:

```bash
conda activate fumaraseevo
```

