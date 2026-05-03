"""
plot_entropy.py - Entropy profile visualization for Class I and II fumarases.

Generates a publication-style 2-panel figure showing sequence entropy across
alignment positions for Class I (Panel A) and Class II (Panel B).

Features:
    - Raw entropy scatter/line (faint)
    - Rolling median (smoothed trend)
    - Interquartile range background highlighting
    - Functional sites marked with red stars/lines

Inputs:
    results/stats/entropy_class_i.csv       (Position, entropy, and functional site flags)
    results/stats/entropy_class_ii.csv      (Position, entropy, and functional site flags)
    results/stats/entropy_summary.csv       (High-level summary stats for annotations)

Output:
    results/figures/figure_entropy_profiles.png (and/or pdf, svg based on config)

Usage:
    # Run with default paths:
    python scripts/08_plots/plot_entropy.py

    # Run with custom inputs/outputs:
    python scripts/08_plots/plot_entropy.py \
        --class-i  data/custom_class_i.csv \
        --class-ii data/custom_class_ii.csv \
        --summary  data/custom_summary.csv \
        --outdir   results/custom_figures/
"""

import sys
import argparse
import logging
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS, FIGURES

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data Loading & Processing
# ---------------------------------------------------------------------------

def to_bool(series: pd.Series) -> pd.Series:
    """Safely convert a series of various truthy string values to boolean."""
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y", "t"})
    )


def load_profile(path: Path) -> pd.DataFrame:
    """Load and validate an entropy profile CSV."""
    log.info(f"Loading entropy profile from {path} ...")
    df = pd.read_csv(path)
    required = {"position", "entropy", "is_functional_site"}
    missing = required.difference(df.columns)
    
    if missing:
        log.error(f"Missing required columns in {path}: {sorted(missing)}")
        sys.exit(1)

    out = df.copy()
    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    out["entropy"] = pd.to_numeric(out["entropy"], errors="coerce")
    out["is_functional_site"] = to_bool(out["is_functional_site"])
    
    before = len(out)
    out = out.dropna(subset=["position", "entropy"]).sort_values("position")
    after = len(out)
    
    if before != after:
        log.info(f"  Dropped {before - after} rows with NaN values")
        
    log.info(f"  Loaded {after:,} positions")
    return out


def get_rolling_window(n_rows: int) -> int:
    """Calculate an appropriate odd-numbered rolling window size based on sequence length."""
    window = max(9, int(round(n_rows * 0.04)))
    if window % 2 == 0:
        window += 1
    return window


def load_summary(path: Path) -> pd.DataFrame:
    """Load the summary stats CSV if it exists."""
    if not path.exists():
        log.warning(f"Summary file not found at {path}. Proceeding without annotations.")
        return pd.DataFrame()
    log.info(f"Loading summary from {path} ...")
    return pd.read_csv(path)


def get_summary_note(summary: pd.DataFrame, label: str) -> str:
    """Format a summary note string for the plot panel."""
    if summary.empty or "label" not in summary.columns:
        return ""
        
    row = summary[summary["label"] == label]
    if row.empty:
        return ""
        
    r = row.iloc[0]
    n_sequences = int(r["n_sequences"]) if "n_sequences" in row.columns else 0
    reference_length = int(r["reference_length"]) if "reference_length" in row.columns else 0
    mean_entropy = float(r["mean_entropy"]) if "mean_entropy" in row.columns else 0.0
    
    return f"n={n_sequences:,}; Lref={reference_length}; mean H={mean_entropy:.3f}"


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_panel(ax: plt.Axes,
               panel_tag: str,
               label: str,
               profile: pd.DataFrame,
               summary_note: str,
               class_color: str) -> None:
    """Plot a single class entropy profile onto the provided Axes."""
    x = profile["position"]
    y = profile["entropy"]
    q1 = float(y.quantile(0.25))
    q3 = float(y.quantile(0.75))
    y_top = max(float(y.max()), q3) * 1.10

    window = get_rolling_window(len(profile))
    y_smooth = y.rolling(window=window, center=True, min_periods=1).median()

    # Raw and smoothed data
    ax.plot(x, y, color=class_color, linewidth=0.9, alpha=0.30, label="Raw entropy")
    ax.plot(
        x, y_smooth,
        color=class_color,
        linewidth=1.9, alpha=1.0,
        label=f"Rolling median (w={window})",
    )

    # Quartile spans and lines
    ax.axhspan(0, q1, color="#2f855a", alpha=0.08, zorder=0)
    ax.axhspan(q3, y_top, color="#c53030", alpha=0.06, zorder=0)
    ax.axhline(
        q1, color="#2f855a", linestyle="--", linewidth=1.2,
        label=f"25th percentile = {q1:.2f}",
    )
    ax.axhline(
        q3, color="#c53030", linestyle="--", linewidth=1.2,
        label=f"75th percentile = {q3:.2f}",
    )

    # Functional sites
    functional = profile[profile["is_functional_site"]]
    if not functional.empty:
        ax.vlines(
            x=functional["position"],
            ymin=0, ymax=functional["entropy"],
            color="#e11d48", linewidth=0.8, alpha=0.30, zorder=2,
        )
        ax.scatter(
            functional["position"],
            [0.0] * len(functional),
            marker="*", s=150,
            color="#e11d48", edgecolors="white", linewidths=0.6,
            clip_on=False, zorder=5, label="Functional sites",
        )

    # Aesthetics
    ax.set_ylim(0, y_top)
    ax.set_ylabel("Entropy (bits)", fontsize=12, labelpad=8)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.35)
    
    # Summary Box
    if summary_note:
        ax.text(
            0.995, 0.94, summary_note,
            transform=ax.transAxes, ha="right", va="top",
            fontsize=12, color="#334155",
            bbox={"facecolor": "white", "edgecolor": "#cbd5e1", 
                  "boxstyle": "round,pad=0.3", "alpha": 0.85},
            zorder=6,
        )

    ax.legend(loc="upper left", frameon=False, fontsize=12, ncols=2)


# ---------------------------------------------------------------------------
# Main Routine
# ---------------------------------------------------------------------------

def plot_entropy_profiles(class_i_csv: Path,
                          class_ii_csv: Path,
                          summary_csv: Path,
                          outdir: Path,
                          dpi: int) -> None:
    """Main plotting orchestration routine."""
    
    profile_i  = load_profile(class_i_csv)
    profile_ii = load_profile(class_ii_csv)
    summary    = load_summary(summary_csv)

    log.info("Generating figure ...")
    
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })

    fig, axes = plt.subplots(2, 1, figsize=(11.4, 7.4), constrained_layout=True)
    
    plot_panel(
        axes[0], panel_tag="A", label="Class I",
        profile=profile_i,
        summary_note=get_summary_note(summary, "Class I"),
        class_color="#b91c1c",
    )
    
    plot_panel(
        axes[1], panel_tag="B", label="Class II",
        profile=profile_ii,
        summary_note=get_summary_note(summary, "Class II"),
        class_color="#1d4ed8",
    )

    axes[1].set_xlabel("Position in reference sequence", fontsize=12, labelpad=8)
    
    # Save standard formats
    outdir.mkdir(parents=True, exist_ok=True)
    formats = FIGURES.get("formats", ["png", "pdf", "svg"])
    
    for fmt in formats:
        out_path = outdir / f"figure_entropy_profiles.{fmt}"
        save_dpi = dpi if fmt == "png" else FIGURES.get("dpi", 600)
        fig.savefig(out_path, dpi=save_dpi, bbox_inches="tight")
        log.info(f"  Saved -> {out_path} (dpi={save_dpi})")

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a publication-style entropy profile figure for Class I and Class II."
    )
    parser.add_argument(
        "--class-i", type=Path,
        default=Path(PATHS["results"]) / "stats" / "entropy_class_i.csv",
        help="CSV with Class I entropy profile.",
    )
    parser.add_argument(
        "--class-ii", type=Path,
        default=Path(PATHS["results"]) / "stats" / "entropy_class_ii.csv",
        help="CSV with Class II entropy profile.",
    )
    parser.add_argument(
        "--summary", type=Path,
        default=Path(PATHS["results"]) / "stats" / "entropy_summary.csv",
        help="Summary CSV used for panel annotation.",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["figures"]),
        help="Output directory for the figures.",
    )
    parser.add_argument(
        "--dpi", type=int, default=900,
        help="Raster DPI used for the PNG export.",
    )
    args = parser.parse_args()

    # Verify input existence
    for p in [args.class_i, args.class_ii]:
        if not p.exists():
            log.error(f"Input file not found: {p}")
            sys.exit(1)

    plot_entropy_profiles(
        class_i_csv=args.class_i,
        class_ii_csv=args.class_ii,
        summary_csv=args.summary,
        outdir=args.outdir,
        dpi=args.dpi,
    )
    
    log.info("\nDone.")


if __name__ == "__main__":
    main()