from pathlib import Path

# ============================================================
#  Central Configuration
# ============================================================

# The root directory of the project
ROOT = Path(__file__).parent

# --- Paths ---
PATHS = {
    "data_raw":       ROOT / "data" / "raw",
    "data_processed": ROOT / "data" / "processed",
    "data_external":  ROOT / "data" / "external",
    "results":        ROOT / "results",
    "alignments":     ROOT / "results" / "alignments",
    "trees":          ROOT / "results" / "trees",
    "figures":        ROOT / "results" / "figures",
    "stats": ROOT / "results" / "stats",

}

# --- Query sequences (E. coli K-12) ---
QUERIES = {
    "class1": {
        "uniprot_id": "P0AC33",
        "gene":       "fumA",
        "length":     548,
    },
    "class2": {
        "uniprot_id": "P05042",
        "gene":       "fumC",
        "length":     467,
    },
}

# --- MMseqs2 filtering thresholds ---
MMSEQS2 = {
    "evalue_cutoff":     1e-3,   # 0.001
    "identity_cutoff":   0.40,   # 40%
    "coverage_cutoff":   0.40,   # 40% of query length
}

# --- Phylogeny (IQ-TREE) ---
IQTREE = {
    "model_selection":  "TEST",  # ModelFinder auto-selects best model
    "bootstraps":       1000,    # ultrafast bootstrap replicates
    "threads":          "AUTO",
    "gamma_categories": 4,       # +G4
}

# --- Conservation (Shannon entropy) ---
ENTROPY = {
    "scale_min": 0.0,    # perfect conservation
    "scale_max": 4.32,   # maximum variability (log2 of 20 AA)
}

# --- Embeddings ---
# TODO: Revise that (if wants to import esm, maybe use ESM2)
ESM = {
    "model": "esm1_t6_43M_UR50S",  # ESM-1b
    "layer": 6,                     # mean pooling layer
}

TSNE = {
    "perplexity":   30,
    "n_iter":       1000,
    "random_seed":  42,
}

UMAP = {
    "n_neighbors": 15,
    "min_dist":    0.1,
    "random_seed": 42,
}

# --- Figures ---
FIGURES = {
    "dpi":     600,     # Or higher for publication-quality
    "formats": ["png", "tiff", "pdf"],
}

# --- Taxonomic group color palette ---
COLORS = {
    "Eukaryote":            "#1f77b4",
    "Alpha-proteobacteria": "#ff7f0e",
    "Bacteria":             "#7f7f7f",
    "Archaea":              "#d62728",
}
