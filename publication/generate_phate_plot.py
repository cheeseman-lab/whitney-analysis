#!/usr/bin/env python3
"""
Generate PHATE map plot highlighting nontargeting controls and essential genes.

Creates a publication-quality PHATE visualization where:
- All points are shown in grey by default
- Three nontargeting control categories are highlighted with distinct colors
- Essential genes are highlighted in a separate color

Usage:
    python generate_phate_plot.py [--output-dir figures/]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Base paths
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
WORKING_DIR = Path("/mnt/data/blainey/whitney-analysis")

# Input files
PHATE_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "20" / "phate_leiden_clustering.tsv"
ESSENTIAL_GENES_FILE = WORKING_DIR / "essential_genes.tsv"

# Standard figure size (matching QC panels style)
FIGURE_SIZE = (8, 7)

# Color scheme
COLORS = {
    "background": "#000000",  # Black for non-essential genes
    "nontargeting_intergenic": "#888888",  # Medium grey
    "nontargeting_noncutting": "#aaaaaa",  # Light grey
    "nontargeting_or": "#666666",  # Dark grey
    "essential": "#cc4444",  # Red
}

LABELS = {
    "background": "Non-Essential Genes",
    "nontargeting_intergenic": "Intergenic Controls",
    "nontargeting_noncutting": "Noncutting Controls",
    "nontargeting_or": "Olfactory Receptor Controls",
    "essential": "Essential Genes",
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
        "legend.fontsize": 10,
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


def load_essential_genes():
    """Load list of essential genes."""
    df = pd.read_csv(ESSENTIAL_GENES_FILE, sep="\t")
    return set(df["gene_symbol_0"].dropna().unique())


def categorize_genes(df, essential_genes):
    """
    Categorize genes into different groups for plotting.

    Returns a dictionary mapping category names to DataFrames.
    """
    categories = {}

    # Identify nontargeting controls
    for prefix in ["nontargeting_intergenic", "nontargeting_noncutting", "nontargeting_or"]:
        mask = df["gene_symbol_0"].str.startswith(prefix, na=False)
        categories[prefix] = df[mask].copy()

    # Identify essential genes (excluding nontargeting controls)
    all_nontargeting = df["gene_symbol_0"].str.startswith("nontargeting", na=False)
    essential_mask = df["gene_symbol_0"].isin(essential_genes) & ~all_nontargeting
    categories["essential"] = df[essential_mask].copy()

    # Background genes (everything else)
    highlighted = all_nontargeting | essential_mask
    categories["background"] = df[~highlighted].copy()

    return categories


def plot_phate_map(df, categories, output_path):
    """
    Create PHATE map plot with highlighted categories.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Plot background genes first (grey)
    if len(categories["background"]) > 0:
        ax.scatter(
            categories["background"]["PHATE_0"],
            categories["background"]["PHATE_1"],
            c=COLORS["background"],
            s=6,
            alpha=0.3,
            edgecolors="none",
            label=f"{LABELS['background']} (n={len(categories['background'])})",
            rasterized=True,
        )

    # Plot essential genes (with transparency)
    if len(categories["essential"]) > 0:
        ax.scatter(
            categories["essential"]["PHATE_0"],
            categories["essential"]["PHATE_1"],
            c=COLORS["essential"],
            s=6,
            alpha=0.5,
            edgecolors="none",
            label=f"{LABELS['essential']} (n={len(categories['essential'])})",
            rasterized=True,
        )

    # Plot nontargeting controls last (on top, small, muted colors)
    for prefix in ["nontargeting_intergenic", "nontargeting_noncutting", "nontargeting_or"]:
        cat_df = categories[prefix]
        if len(cat_df) > 0:
            ax.scatter(
                cat_df["PHATE_0"],
                cat_df["PHATE_1"],
                c=COLORS[prefix],
                s=6,
                alpha=0.7,
                edgecolors="none",
                label=f"{LABELS[prefix]} (n={len(cat_df)})",
                rasterized=True,
            )

    # Highlight specific genes with black circles
    highlight_genes = ["RPL36A-HNRNPH2", "nontargeting_intergenic_GAGGGTCGCTTC"]
    for gene in highlight_genes:
        gene_data = df[df["gene_symbol_0"] == gene]
        if len(gene_data) > 0:
            ax.scatter(
                gene_data["PHATE_0"],
                gene_data["PHATE_1"],
                s=80,
                facecolors="none",
                edgecolors="black",
                linewidths=2,
                zorder=100,
            )

    # Formatting - remove ticks but keep labels
    ax.set_xlabel("PHATE 1")
    ax.set_ylabel("PHATE 2")
    ax.set_title("Clustering: All Perturbations", fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])

    # Legend
    ax.legend(
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


def print_summary(categories):
    """Print summary of gene categories."""
    print("\nGene Category Summary:")
    print("-" * 40)
    print(f"Background genes: {len(categories['background'])}")
    print(f"Essential genes: {len(categories['essential'])}")
    for prefix in ["nontargeting_intergenic", "nontargeting_noncutting", "nontargeting_or"]:
        print(f"{LABELS[prefix]}: {len(categories[prefix])}")
    print("-" * 40)
    total = sum(len(v) for v in categories.values())
    print(f"Total: {total}")


def main():
    parser = argparse.ArgumentParser(description="Generate PHATE map plot")
    parser.add_argument("--output-dir", type=str, default="figures",
                       help="Output directory for figures")
    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    setup_style()

    print("=" * 60)
    print("Generating PHATE Map Plot")
    print("=" * 60)

    # Load data
    print("\nLoading PHATE data...")
    phate_df = load_phate_data()
    print(f"  Loaded {len(phate_df)} perturbations")

    print("Loading essential genes list...")
    essential_genes = load_essential_genes()
    print(f"  Loaded {len(essential_genes)} essential genes")

    # Categorize genes
    print("\nCategorizing genes...")
    categories = categorize_genes(phate_df, essential_genes)
    print_summary(categories)

    # Generate plot
    print("\nGenerating PHATE map...")
    output_path = output_dir / "phate_map_controls.png"
    plot_phate_map(phate_df, categories, output_path)

    print("\n" + "=" * 60)
    print(f"Plot saved to: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
