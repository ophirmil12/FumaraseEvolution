"""
plot_identity.py - Pairwise sequence identity boxplot with jitter.

Reproduces the style in the target figure:
    - Boxplot with overlaid jittered scatter per comparison group
    - Significance brackets with p-values from identity_pairwise.csv
    - Median labels inside each box
    - N= labels below x-axis
    - Colors from config.py COLORS

Inputs (from stats_identity.py outputs):
    results/stats/identity_descriptive.csv      medians, N per group
    results/stats/identity_pairwise.csv         BH-corrected p-values
    results/stats/identity_pairs_annotated.tsv  full data (sampled for scatter)

Output:
    results/figures/figure_identity.png   (600 dpi)
    results/figures/figure_identity.pdf

Usage:
    python scripts/08_plots/plot_identity.py

    # Custom sample size for jitter (default 5000 per group):
    python scripts/08_plots/plot_identity.py --sample 3000

    # Skip specific groups:
    python scripts/08_plots/plot_identity.py --groups Euk-Alpha Euk-OtherBacteria Euk-Archaea
"""

import sys
import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.ticker as ticker

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS, COLORS, FIGURES

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Group display config
# ---------------------------------------------------------------------------

# Order and display labels for x-axis (matching the target figure)
GROUP_ORDER = [
    "Euk-Alpha",
    "Euk-OtherBacteria",
    "Euk-Archaea",
]

GROUP_LABELS = {
    "Euk-Alpha":          "Alpha-proteobacteria\u2013Eukaryote",
    "Euk-OtherBacteria":  "Other Bacteria\u2013Eukaryote",
    "Euk-Archaea":        "Eukaryote\u2013Archaea",
}

# Map comparison group names -> COLORS key
GROUP_COLOR_KEY = {
    "Euk-Alpha":         "Alpha-proteobacteria",
    "Euk-OtherBacteria": "Bacteria",
    "Euk-Archaea":       "Archaea",
}

DEFAULT_COLOR = "#aaaaaa"


def get_color(group: str) -> str:
    return COLORS.get(GROUP_COLOR_KEY.get(group, ""), DEFAULT_COLOR)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_inputs(stats_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    desc_path = stats_dir / "identity_descriptive.csv"
    pw_path   = stats_dir / "identity_pairwise.csv"
    ann_path  = stats_dir / "identity_pairs_annotated.tsv"

    for p in [desc_path, pw_path, ann_path]:
        if not p.exists():
            log.error(f"Missing input: {p}")
            sys.exit(1)

    desc = pd.read_csv(desc_path)
    pw   = pd.read_csv(pw_path)
    log.info(f"Loading annotated pairs from {ann_path} ...")
    ann  = pd.read_csv(ann_path, sep="\t",
                       usecols=["query_group", "target_group", "local_pct"])
    log.info(f"  Loaded {len(ann):,} pairs")
    return desc, pw, ann


def build_group_series(ann: pd.DataFrame,
                       sample: int,
                       seed: int = 42) -> dict[str, pd.Series]:
    """
    Reconstruct per-comparison-group local_pct Series from annotated pairs.
    Samples up to `sample` rows per group for the jitter scatter.
    Returns full series for boxplot stats, sampled for scatter.
    """
    rng = np.random.default_rng(seed)
    result_full    = {}
    result_sampled = {}

    for label, g1, g2 in [
        ("Euk-Alpha",         "Eukaryota", "Alphaproteobacteria"),
        ("Euk-OtherBacteria", "Eukaryota", "OtherBacteria"),
        ("Euk-Archaea",       "Eukaryota", "Archaea"),
    ]:
        mask = (
            ((ann["query_group"] == g1) & (ann["target_group"] == g2)) |
            ((ann["query_group"] == g2) & (ann["target_group"] == g1))
        )
        full = ann[mask]["local_pct"].dropna().values
        result_full[label] = full

        n = min(sample, len(full))
        idx = rng.choice(len(full), size=n, replace=False)
        result_sampled[label] = full[idx]

        log.info(f"  {label}: {len(full):,} pairs, {n:,} sampled for scatter")

    return result_full, result_sampled


# ---------------------------------------------------------------------------
# Significance formatting
# ---------------------------------------------------------------------------

def format_pvalue(p: float) -> str:
    """Format p-value for display on significance bracket."""
    if p < 2.22e-16:
        return "p < 2.22e-16"
    elif p < 0.001:
        return f"{p:.2e}"
    elif p < 0.05:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.3f} (ns)"


def get_pvalue(pw: pd.DataFrame, g1: str, g2: str) -> float | None:
    """Look up BH-corrected p-value for a group pair."""
    mask = (
        ((pw["group_1"] == g1) & (pw["group_2"] == g2)) |
        ((pw["group_1"] == g2) & (pw["group_2"] == g1))
    )
    row = pw[mask]
    if row.empty:
        return None
    return float(row["p_adj_BH"].iloc[0])


# ---------------------------------------------------------------------------
# Significance bracket drawing
# ---------------------------------------------------------------------------

def draw_bracket(ax, x1: float, x2: float, y: float,
                 label: str, fontsize: float = 9) -> None:
    """Draw a significance bracket between positions x1 and x2 at height y."""
    tick_h = 1.0   # vertical tick height in data units

    ax.plot([x1, x1, x2, x2],
            [y, y + tick_h, y + tick_h, y],
            lw=1.0, color="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + tick_h + 0.3, label,
            ha="center", va="bottom", fontsize=fontsize, color="black")


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot_figure(desc: pd.DataFrame,
                pw: pd.DataFrame,
                data_full: dict[str, np.ndarray],
                data_sampled: dict[str, np.ndarray],
                outdir: Path,
                groups: list[str],
                scatter: int = 500) -> None:

    fig, ax = plt.subplots(figsize=(10, 7), facecolor="white")
    ax.set_facecolor("white")

    positions    = list(range(len(groups)))
    violin_width = 0.6
    rng          = np.random.default_rng(42)

    # Precompute all_vals and y_n before loop
    all_vals = np.concatenate([v for v in data_full.values() if len(v) > 0])
    y_n      = max(0, np.percentile(all_vals, 0.01) - 4)

    # --- Draw violin + scatter + labels per group ---
    for pos, group in zip(positions, groups):
        color     = get_color(group)
        full_vals = data_full.get(group, np.array([]))
        samp_vals = data_sampled.get(group, np.array([]))

        if len(full_vals) == 0:
            log.warning(f"  {group}: no data — skipping")
            continue

        # Violin
        parts = ax.violinplot(
            full_vals,
            positions=[pos],
            widths=violin_width,
            showmedians=True,
            showextrema=True,
            showmeans=False,
        )
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_edgecolor("#333333")
            pc.set_alpha(0.6)
            pc.set_linewidth(0.8)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(2.5)
        parts["cmins"].set_color("#555555")
        parts["cmaxes"].set_color("#555555")
        parts["cbars"].set_color("#555555")
        parts["cbars"].set_linewidth(0.8)

        # Jitter scatter
        n_scatter   = min(scatter, len(samp_vals))
        scatter_idx = rng.choice(len(samp_vals), size=n_scatter, replace=False)
        scatter_vals = samp_vals[scatter_idx]
        jitter = rng.uniform(-0.12, 0.12, size=n_scatter)
        ax.scatter(
            pos + jitter, scatter_vals,
            color=color, alpha=0.15, s=3, linewidths=0,
            zorder=2, rasterized=True,
        )

        # Median label
        median_val = np.median(full_vals)
        ax.text(pos, median_val + 0.8, f"Median: {median_val:.2f}%",
                ha="center", va="bottom", fontsize=9,
                fontweight="bold", zorder=5)

        # N= label
        ax.text(pos, y_n, f"N = {len(full_vals):,}",
                ha="center", va="top", fontsize=9, color="#333333")

    # --- Significance brackets ---
    y_max        = np.percentile(all_vals, 99.9)
    bracket_base = y_max + 3

    bracket_pairs   = [(0, 1), (0, 2), (1, 2)]
    bracket_heights = [bracket_base, bracket_base + 6, bracket_base + 12]

    for (i, j), y in zip(bracket_pairs, bracket_heights):
        g1, g2 = groups[i], groups[j]
        p = get_pvalue(pw, g1, g2)
        if p is not None:
            draw_bracket(ax, i, j, y, format_pvalue(p), fontsize=9)

    # --- Axes formatting ---
    ax.set_xticks(positions)
    ax.set_xticklabels([GROUP_LABELS.get(g, g) for g in groups], fontsize=10)
    ax.set_ylabel("Pairwise Sequence Identity (%)", fontsize=11, labelpad=8)
    ax.set_xlim(-0.6, len(groups) - 0.4)

    y_top = bracket_base + 12 + 6
    ax.set_ylim(
        max(0, np.percentile(all_vals, 0.01) - 5),
        y_top,
    )

    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.1, right=0.95, top=0.92, bottom=0.12)

    # --- Save ---
    outdir.mkdir(parents=True, exist_ok=True)
    for fmt in FIGURES.get("formats", ["png", "pdf"]):
        out_path = outdir / f"figure_identity_violin.{fmt}"
        dpi = 300 if fmt == "png" else FIGURES.get("dpi", 600)
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        log.info(f"  Saved -> {out_path}, dpi={dpi}")

    plt.close()

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot Figure: pairwise identity boxplot with jitter."
    )
    parser.add_argument(
        "--statsdir", type=Path,
        default=Path(PATHS["results"]) / "stats",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["figures"]),
    )
    parser.add_argument(
        "--sample", type=int, default=5000,
        help="Max jitter points per group (default: 5000)"
    )
    parser.add_argument(
        "--groups", nargs="+",
        default=GROUP_ORDER,
        help="Groups to plot (default: all three Euk comparisons)"
    )
    parser.add_argument(
        "--scatter", type=int, default=500,
        help="Points per group for jitter overlay (default: 500, independent of --sample)"
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    desc, pw, ann = load_inputs(args.statsdir)

    data_full, data_sampled = build_group_series(ann, args.sample, args.seed)

    plot_figure(desc, pw, data_full, data_sampled, args.outdir, args.groups, args.scatter)

    log.info(f"\nDone. Outputs in {args.outdir}")


if __name__ == "__main__":
    main()