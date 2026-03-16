from pathlib import Path

# ============================================================
#  Central Configuration
# ============================================================

# TODO: search all files for "TODO", then remove environment-specific settings
#  e.g. paths (ophirmil12, dror...), mails, etc.

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

# Known catalytic residues from PDB annotations (0-indexed in reference sequence)
# Class I  (P0AC33, FumA): [4Fe-4S] cluster coordinating cysteines
# Class II (P05042, FumC): active site residues
FUNCTIONAL_SITES = {
    "class1": {
        "name":     "Class I (P0AC33)",
        "residues": [104, 223, 317],          # Cys105, Cys224, Cys318
    },
    "class2": {
        "name":     "Class II (P05042)",
        "residues": [187, 317, 323, 330],     # His188, Ser318, Lys324, Glu331
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
    "model_selection":  "TEST",  # ModelFinder auto-selects best model     TODO: identify the specific model chosen by ModelFinder (e.g., LG+G4, JTT+I+G) for both the sequence-based and 3Di-based trees and update the text accordingly.
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
ESM = {
    "model": "esm2_t33_650M_UR50D",   # ESM-2 model variant
    "layer": 33,                      # mean pooling layer
    "embedding_size": 1280,           # dimensionality of the chosen ESM model
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
    "Alphaproteobacteria": "#ff7f0e",
    "Bacteria":             "#A6CEE3",
    "Archaea":              "#d62728",
    # #49306B, purple
}
