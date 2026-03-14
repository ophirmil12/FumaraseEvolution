"""
plot_embeddings.py - Dimensionality reduction and visualization of ESM-2 embeddings.

Takes a dict of {seq_id: mean_embedding_tensor}, reduces to 2D or 3D using
t-SNE, UMAP, or PCA, and produces:
    - Static PNG (600 dpi) colored by taxonomic group
    - Interactive HTML (Plotly) with hover showing organism UP ID + name

Taxonomic coloring is defined in config.py COLORS.
Metadata (organism name) is loaded from proteome_metadata.tsv.

Usage:
    # Typically called from a higher-level script after ESMFactory.get_mean_embeddings()
    from scripts.08_plots.esm_plot import plot_embeddings

    embeddings = factory.get_mean_embeddings(sequence_data)
    plot_embeddings(
        embeddings  = embeddings,
        fasta_path  = Path("data/processed/class1_sequences.fasta"),
        method      = "tsne",       # "tsne" | "umap" | "pca"
        n_dims      = 2,            # 2 | 3
        label       = "class1",
        outdir      = Path("results/embeddings/"),
    )

    # Or run directly:
    python scripts/08_plots/esm_plot.py \\
        --embeddings data/processed/class1_embeddings.pt \\
        --fasta      data/processed/class1_sequences.fasta \\
        --method     tsne \\
        --dims       2 \\
        --label      class1

Output:
    results/embeddings/class1_tsne_2d.png
    results/embeddings/class1_tsne_2d.html
"""

import sys
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS, COLORS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Fallback color for groups not in COLORS
DEFAULT_COLOR = "#aaaaaa"


# ---------------------------------------------------------------------------
# Taxonomic group assignment
# ---------------------------------------------------------------------------

import re as _re

def assign_group(lineage: str) -> str:
    if not isinstance(lineage, str) or not lineage.strip():
        return "Other"
    if "Alphaproteobacteria" in lineage:
        return "Alpha-proteobacteria"
    if _re.search(r",\s*Bacteria,", lineage):
        return "Bacteria"
    if _re.search(r",\s*Archaea,", lineage):
        return "Archaea"
    if "Eukaryota" in lineage:
        return "Eukaryote"
    return "Other"


# ---------------------------------------------------------------------------
# FASTA parsing (to get seq_id -> proteome_id mapping)
# ---------------------------------------------------------------------------

def parse_fasta_headers(fasta_path: Path) -> dict[str, str]:
    """
    Return {header: proteome_id} from a FASTA file.
    Header format: UP000XXXXXX|sp|ACCESSION|... → proteome_id = UP000XXXXXX
    """
    mapping = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                proteome_id = header.split("|")[0].strip()
                mapping[header] = proteome_id
    log.info(f"  Parsed {len(mapping):,} headers from {fasta_path.name}")
    return mapping


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load proteome metadata, assign taxonomic groups and organism names."""
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
    meta.columns = meta.columns.str.strip()
    meta = meta.set_index("proteome_id")
    meta["group"]    = meta["taxonomic_lineage"].apply(assign_group)
    meta["org_name"] = meta.get("organism", meta.index.to_series())
    log.info(f"  Loaded metadata for {len(meta):,} proteomes")
    return meta


# ---------------------------------------------------------------------------
# Dimensionality reduction
# ---------------------------------------------------------------------------

def reduce_embeddings(matrix: np.ndarray,
                      method: str,
                      n_dims: int,
                      seed: int = 42) -> np.ndarray:
    """
    Reduce (N, D) embedding matrix to (N, n_dims) using method.
    method: "tsne" | "umap" | "pca"
    n_dims: 2 | 3
    """
    log.info(f"  Running {method.upper()} (n={len(matrix)}, D={matrix.shape[1]}) -> {n_dims}D ...")

    if method == "pca":
        from sklearn.decomposition import PCA
        reducer = PCA(n_components=n_dims, random_state=seed)
        return reducer.fit_transform(matrix)

    elif method == "tsne":
        from sklearn.manifold import TSNE
        # PCA pre-reduction for speed when D > 50
        if matrix.shape[1] > 50:
            from sklearn.decomposition import PCA
            matrix = PCA(n_components=50, random_state=seed).fit_transform(matrix)
        reducer = TSNE(
            n_components=n_dims,
            perplexity=min(30, len(matrix) - 1),
            max_iter=1000,
            random_state=seed,
            init="pca",
            learning_rate="auto",
        )
        return reducer.fit_transform(matrix)

    elif method == "umap":
        from umap import UMAP
        reducer = UMAP(
            n_components=n_dims,
            n_neighbors=15,
            min_dist=0.1,
            random_state=seed,
        )
        return reducer.fit_transform(matrix)

    else:
        raise ValueError(f"Unknown method: {method}. Choose from: tsne, umap, pca")


# ---------------------------------------------------------------------------
# Build annotation dataframe
# ---------------------------------------------------------------------------

def build_plot_df(embeddings: dict[str, torch.Tensor],
                  fasta_path: Path,
                  metadata: pd.DataFrame,
                  method: str,
                  n_dims: int) -> pd.DataFrame:
    """
    Align embeddings with metadata, run reduction, return annotated DataFrame.
    """
    # Build ordered lists
    seq_ids = list(embeddings.keys())
    matrix  = np.stack([embeddings[sid].numpy() for sid in seq_ids])

    # Header -> proteome_id mapping
    header_to_proteome = parse_fasta_headers(fasta_path)

    # Run reduction
    coords = reduce_embeddings(matrix, method, n_dims)

    rows = []
    for i, seq_id in enumerate(seq_ids):
        proteome_id = header_to_proteome.get(seq_id, seq_id.split("|")[0])
        meta_row    = metadata.loc[proteome_id] if proteome_id in metadata.index else {}

        group    = meta_row.get("group",    "Other") if hasattr(meta_row, "get") else "Other"
        org_name = meta_row.get("org_name", "")      if hasattr(meta_row, "get") else ""

        row = {
            "seq_id":      seq_id,
            "proteome_id": proteome_id,
            "org_name":    org_name,
            "group":       group,
            "color":       COLORS.get(group, DEFAULT_COLOR),
            "x":           coords[i, 0],
            "y":           coords[i, 1],
        }
        if n_dims == 3:
            row["z"] = coords[i, 2]
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Static plot (matplotlib)
# ---------------------------------------------------------------------------

def plot_static(df: pd.DataFrame,
                method: str,
                n_dims: int,
                label: str,
                outdir: Path) -> Path:
    """Produce a publication-quality static PNG."""

    fig = plt.figure(figsize=(10, 8), facecolor="white")

    if n_dims == 3:
        ax = fig.add_subplot(111, projection="3d")
    else:
        ax = fig.add_subplot(111)

    # Plot each group separately for legend
    groups_present = df["group"].unique()
    for group in sorted(groups_present):
        sub   = df[df["group"] == group]
        color = COLORS.get(group, DEFAULT_COLOR)
        kwargs = dict(
            c     = color,
            s     = 8,
            alpha = 0.7,
            linewidths = 0,
            label = group,
            zorder = 2,
        )
        if n_dims == 3:
            ax.scatter(sub["x"], sub["y"], sub["z"], **kwargs)
        else:
            ax.scatter(sub["x"], sub["y"], **kwargs)

    # Axes labels
    axis_label = method.upper()
    if n_dims == 3:
        ax.set_xlabel(f"{axis_label} 1", fontsize=9, labelpad=4)
        ax.set_ylabel(f"{axis_label} 2", fontsize=9, labelpad=4)
        ax.set_zlabel(f"{axis_label} 3", fontsize=9, labelpad=4)
    else:
        ax.set_xlabel(f"{axis_label} 1", fontsize=10)
        ax.set_ylabel(f"{axis_label} 2", fontsize=10)

    ax.tick_params(labelsize=7)
    ax.set_facecolor("#f9f9f9")

    # Legend
    handles = [
        mpatches.Patch(color=COLORS.get(g, DEFAULT_COLOR), label=g)
        for g in sorted(groups_present)
    ]
    ax.legend(
        handles=handles,
        title="Taxonomic group",
        title_fontsize=9,
        fontsize=8,
        loc="best",
        framealpha=0.9,
        edgecolor="#cccccc",
    )

    class_label = label.replace("_", " ").title()
    plt.title(
        f"ESM-2 Embeddings — {class_label} ({method.upper()} {n_dims}D)",
        fontsize=12, pad=12,
    )
    plt.tight_layout()

    out_path = outdir / f"{label}_{method}_{n_dims}d.png"
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()
    log.info(f"  Static PNG -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Interactive plot (Plotly)
# ---------------------------------------------------------------------------

def plot_interactive(df: pd.DataFrame,
                     method: str,
                     n_dims: int,
                     label: str,
                     outdir: Path) -> Path:
    """Produce an interactive HTML plot with hover labels."""

    traces = []
    for group in sorted(df["group"].unique()):
        sub   = df[df["group"] == group]
        color = COLORS.get(group, DEFAULT_COLOR)

        hover_text = (
            "<b>" + sub["org_name"] + "</b><br>" +
            "Proteome: " + sub["proteome_id"] + "<br>" +
            "Group: " + sub["group"]
        )

        common = dict(
            name        = group,
            mode        = "markers",
            marker      = dict(size=5, color=color, opacity=0.75,
                               line=dict(width=0)),
            text        = hover_text,
            hoverinfo   = "text",
        )

        if n_dims == 3:
            trace = go.Scatter3d(
                x=sub["x"], y=sub["y"], z=sub["z"],
                **common,
            )
        else:
            trace = go.Scatter(
                x=sub["x"], y=sub["y"],
                **common,
            )
        traces.append(trace)

    axis_label = method.upper()
    layout_common = dict(
        title=dict(
            text=f"ESM-2 Embeddings — {label.replace('_', ' ').title()} "
                 f"({method.upper()} {n_dims}D)",
            font=dict(size=15),
        ),
        paper_bgcolor="white",
        plot_bgcolor="#f9f9f9",
        legend=dict(title="Taxonomic group", font=dict(size=11)),
        hoverlabel=dict(bgcolor="white", font_size=12),
    )

    if n_dims == 3:
        layout = go.Layout(
            **layout_common,
            scene=dict(
                xaxis_title=f"{axis_label} 1",
                yaxis_title=f"{axis_label} 2",
                zaxis_title=f"{axis_label} 3",
            ),
        )
    else:
        layout = go.Layout(
            **layout_common,
            xaxis=dict(title=f"{axis_label} 1", showgrid=True, gridcolor="#e0e0e0"),
            yaxis=dict(title=f"{axis_label} 2", showgrid=True, gridcolor="#e0e0e0"),
        )

    fig = go.Figure(data=traces, layout=layout)

    out_path = outdir / f"{label}_{method}_{n_dims}d.html"
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  Interactive HTML -> {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_embeddings(embeddings: dict[str, torch.Tensor],
                    fasta_path: Path,
                    method: str = "tsne",
                    n_dims: int = 2,
                    label: str = "class",
                    outdir: Path | None = None,
                    seed: int = 42) -> tuple[Path, Path]:
    """
    Main entry point. Produces static PNG + interactive HTML.

    Args:
        embeddings : {seq_id: mean_embedding_tensor} from ESMFactory
        fasta_path : path to the input FASTA (for header -> proteome_id mapping)
        method     : "tsne" | "umap" | "pca"
        n_dims     : 2 | 3
        label      : output filename prefix (e.g. "class1", "class2")
        outdir     : output directory (default: results/embeddings/)
        seed       : random seed for reproducibility

    Returns:
        (png_path, html_path)
    """
    if outdir is None:
        outdir = Path(PATHS["results"]) / "embeddings"
    outdir.mkdir(parents=True, exist_ok=True)

    metadata_path = Path(PATHS["data_processed"]) / "proteome_metadata.tsv"
    metadata = load_metadata(metadata_path)

    log.info(f"Building plot dataframe ({label}, {method}, {n_dims}D) ...")
    df = build_plot_df(embeddings, fasta_path, metadata, method, n_dims)

    log.info(f"Group distribution:\n{df['group'].value_counts().to_string()}")

    png_path  = plot_static(df, method, n_dims, label, outdir)
    html_path = plot_interactive(df, method, n_dims, label, outdir)

    return png_path, html_path


# ---------------------------------------------------------------------------
# CLI (for standalone use with pre-saved embeddings)
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot ESM-2 embeddings with dimensionality reduction."
    )
    parser.add_argument(
        "--embeddings", type=Path, required=True,
        help="Path to .pt file containing {seq_id: tensor} dict saved with torch.save()"
    )
    parser.add_argument(
        "--fasta", type=Path, required=True,
        help="Input FASTA (used for header -> proteome_id mapping)"
    )
    parser.add_argument(
        "--method", choices=["tsne", "umap", "pca"], default="tsne",
    )
    parser.add_argument(
        "--dims", type=int, choices=[2, 3], default=2,
    )
    parser.add_argument(
        "--label", type=str, default="class",
        help="Output filename prefix"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["results"]) / "embeddings",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    for p in [args.embeddings, args.fasta]:
        if not p.exists():
            log.error(f"File not found: {p}")
            sys.exit(1)

    log.info(f"Loading embeddings from {args.embeddings} ...")
    embeddings = torch.load(args.embeddings, map_location="cpu")
    log.info(f"  Loaded {len(embeddings):,} embeddings")

    png_path, html_path = plot_embeddings(
        embeddings = embeddings,
        fasta_path = args.fasta,
        method     = args.method,
        n_dims     = args.dims,
        label      = args.label,
        outdir     = args.outdir,
        seed       = args.seed,
    )

    print(f"\nOutputs:")
    print(f"  Static:      {png_path}")
    print(f"  Interactive: {html_path}")


if __name__ == "__main__":
    main()