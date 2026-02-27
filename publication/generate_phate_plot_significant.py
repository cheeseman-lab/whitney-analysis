#!/usr/bin/env python3
"""
Generate PHATE map plot for significant perturbations with high-confidence clusters colored.

Usage:
    python generate_phate_plot_significant.py [--output-dir figures/]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

# Base paths
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
WORKING_DIR = Path("/mnt/data/blainey/whitney-analysis")

# Input files
PHATE_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "phate_leiden_clustering.tsv"
CLUSTER_SUMMARIES = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "mozzarellm" / "claude-sonnet-4-5-20250929_results_summaries.tsv"

# Standard figure size
FIGURE_SIZE = (8, 7)

# Nontargeting control colors (grey shades)
CONTROL_COLORS = {
    "nontargeting_intergenic": "#888888",  # Medium grey
    "nontargeting_noncutting": "#aaaaaa",  # Light grey
    "nontargeting_or": "#666666",  # Dark grey
}

CONTROL_LABELS = {
    "nontargeting_intergenic": "Intergenic Controls",
    "nontargeting_noncutting": "Noncutting Controls",
    "nontargeting_or": "Olfactory Receptor Controls",
}


def setup_style():
    """Set up publication-quality matplotlib style with sans-serif font."""
    plt.style.use("seaborn-v0_8-white")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "Helvetica"],
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.bottom": False,
        "axes.spines.left": False,
        "axes.grid": False,
    })


def load_phate_data():
    """Load PHATE clustering data."""
    df = pd.read_csv(PHATE_FILE, sep="\t")
    return df


def load_high_confidence_clusters():
    """Load cluster summaries and return set of high-confidence cluster IDs."""
    df = pd.read_csv(CLUSTER_SUMMARIES, sep="\t")
    high_conf = df[df["pathway_confidence"] == "High"]["cluster_id"].tolist()
    return set(high_conf)


def get_cluster_colormap(high_conf_clusters):
    """Create a colormap for high-confidence clusters."""
    n_clusters = len(high_conf_clusters)
    # Use a qualitative colormap with good distinction
    cmap = cm.get_cmap("tab20", n_clusters)
    cluster_colors = {}
    for i, cluster_id in enumerate(sorted(high_conf_clusters)):
        cluster_colors[cluster_id] = cmap(i)
    return cluster_colors


def categorize_genes(df, high_conf_clusters):
    """
    Categorize genes into different groups for plotting.
    """
    categories = {}

    # Identify nontargeting controls
    for prefix in ["nontargeting_intergenic", "nontargeting_noncutting", "nontargeting_or"]:
        mask = df["gene_symbol_0"].str.startswith(prefix, na=False)
        categories[prefix] = df[mask].copy()

    # Identify genes in high-confidence clusters (excluding nontargeting)
    all_nontargeting = df["gene_symbol_0"].str.startswith("nontargeting", na=False)
    high_conf_mask = df["cluster"].isin(high_conf_clusters) & ~all_nontargeting
    categories["high_confidence"] = df[high_conf_mask].copy()

    # Background genes (low/medium confidence, not nontargeting)
    background_mask = ~high_conf_mask & ~all_nontargeting
    categories["background"] = df[background_mask].copy()

    return categories


def plot_phate_map(df, categories, high_conf_clusters, cluster_colors, output_path):
    """
    Create PHATE map plot with high-confidence clusters colored.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Plot background genes first (black)
    if len(categories["background"]) > 0:
        ax.scatter(
            categories["background"]["PHATE_0"],
            categories["background"]["PHATE_1"],
            c="#000000",
            s=6,
            alpha=0.3,
            edgecolors="none",
            label="Other Clusters",
            rasterized=True,
        )

    # Plot high-confidence cluster genes by cluster (only first one gets legend label)
    high_conf_df = categories["high_confidence"]
    first_high_conf = True
    cluster_centroids = {}
    for cluster_id in sorted(high_conf_clusters):
        cluster_data = high_conf_df[high_conf_df["cluster"] == cluster_id]
        if len(cluster_data) > 0:
            ax.scatter(
                cluster_data["PHATE_0"],
                cluster_data["PHATE_1"],
                c=[cluster_colors[cluster_id]],
                s=14,
                alpha=0.7,
                edgecolors="none",
                label="LLM High Confidence Clusters" if first_high_conf else None,
                rasterized=True,
                zorder=10,
            )
            # Store centroid for labeling
            cluster_centroids[cluster_id] = (
                cluster_data["PHATE_0"].mean(),
                cluster_data["PHATE_1"].mean()
            )
            first_high_conf = False

    # Add cluster number labels at centroids with jittering to avoid overlap
    import random
    random.seed(42)  # For reproducibility
    for cluster_id, (cx, cy) in cluster_centroids.items():
        # Add small random jitter to label position
        jitter_x = random.uniform(-0.008, 0.008)
        jitter_y = random.uniform(-0.008, 0.008)
        ax.annotate(
            str(cluster_id),
            (cx + jitter_x, cy + jitter_y),
            fontsize=10,
            fontweight='bold',
            ha='center',
            va='center',
            color='black',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.7),
            zorder=20,
        )

    # Plot nontargeting controls last (on top) - only intergenic for this dataset
    intergenic_df = categories["nontargeting_intergenic"]
    if len(intergenic_df) > 0:
        ax.scatter(
            intergenic_df["PHATE_0"],
            intergenic_df["PHATE_1"],
            c=CONTROL_COLORS["nontargeting_intergenic"],
            s=6,
            alpha=0.7,
            edgecolors="none",
            label="Intergenic Controls",
            rasterized=True,
        )

    # Formatting - remove ticks but keep labels
    ax.set_xlabel("PHATE 1")
    ax.set_ylabel("PHATE 2")
    ax.set_title("Clustering: Significant Perturbations", fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])

    # Simplified legend with rainbow gradient for high confidence clusters
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.legend_handler import HandlerBase
    import matplotlib.patches as mpatches

    # Custom handler for rainbow rectangle
    class RainbowHandler(HandlerBase):
        def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
            rainbow_colors = ['#e41a1c', '#ff7f00', '#ffff33', '#4daf4a', '#377eb8', '#984ea3']
            n = len(rainbow_colors)
            patches = []
            w = width / n
            for i, color in enumerate(rainbow_colors):
                patch = mpatches.Rectangle([xdescent + i * w, ydescent], w, height,
                                          facecolor=color, edgecolor='none', transform=trans)
                patches.append(patch)
            return patches

    # Create legend handles
    other_handle = Line2D([0], [0], marker='o', color='w', markerfacecolor='#000000',
                          markersize=5, alpha=0.5, label='Other Clusters')
    rainbow_handle = Rectangle((0, 0), 1, 1, label='LLM High Confidence')
    intergenic_handle = Line2D([0], [0], marker='o', color='w',
                               markerfacecolor=CONTROL_COLORS["nontargeting_intergenic"],
                               markersize=5, label='Intergenic Controls')

    ax.legend(
        handles=[other_handle, rainbow_handle, intergenic_handle],
        handler_map={Rectangle: RainbowHandler()},
        loc="lower right",
        framealpha=0.9,
        edgecolor="none",
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # Also save as PDF
    pdf_path = output_path.with_suffix(".pdf")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")


def print_summary(categories, high_conf_clusters):
    """Print summary of gene categories."""
    print("\nGene Category Summary:")
    print("-" * 40)
    print(f"Background (low/medium conf): {len(categories['background'])}")
    print(f"High-confidence clusters: {len(categories['high_confidence'])}")
    print(f"  Clusters: {sorted(high_conf_clusters)}")
    for prefix in ["nontargeting_intergenic", "nontargeting_noncutting", "nontargeting_or"]:
        print(f"{CONTROL_LABELS[prefix]}: {len(categories[prefix])}")
    print("-" * 40)
    total = sum(len(v) for v in categories.values())
    print(f"Total: {total}")


def main():
    parser = argparse.ArgumentParser(description="Generate PHATE map plot for significant perturbations")
    parser.add_argument("--output-dir", type=str, default="figures",
                       help="Output directory for figures")
    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    setup_style()

    print("=" * 60)
    print("Generating PHATE Map Plot (Significant Perturbations)")
    print("=" * 60)

    # Load data
    print("\nLoading PHATE data...")
    phate_df = load_phate_data()
    print(f"  Loaded {len(phate_df)} perturbations")

    print("Loading cluster summaries...")
    high_conf_clusters = load_high_confidence_clusters()
    print(f"  Found {len(high_conf_clusters)} high-confidence clusters: {sorted(high_conf_clusters)}")

    # Get cluster colors
    cluster_colors = get_cluster_colormap(high_conf_clusters)

    # Categorize genes
    print("\nCategorizing genes...")
    categories = categorize_genes(phate_df, high_conf_clusters)
    print_summary(categories, high_conf_clusters)

    # Generate plot
    print("\nGenerating PHATE map...")
    output_path = output_dir / "phate_map_significant.png"
    plot_phate_map(phate_df, categories, high_conf_clusters, cluster_colors, output_path)

    print("\n" + "=" * 60)
    print(f"Plot saved to: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
