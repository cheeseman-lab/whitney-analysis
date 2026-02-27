#!/usr/bin/env python3
"""
Generate QC panels for Lander presentation.

Creates 6 publication-quality QC panel figures:
(a) Spatial heatmap of cell mapping to single gene
(b) Cell count boxplot (1 Gene vs >1 Gene)
(c) Mapping percentage bar plot
(d) Barcode prefix matching per well
(e) Q-scores across cycles per well
(f) KDE of cells per perturbation

Usage:
    python generate_qc_panels.py [--output-dir figures/] [--sample-reads 0.1]
"""

import argparse
import string
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Import plot_plate_heatmap from generate_heatmaps
from generate_heatmaps import plot_plate_heatmap

# Base path for brieflow output
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
CONFIG_DIR = Path("/mnt/data/blainey/whitney-analysis/analysis/config")

# Plate and well configuration
PLATES = ["P-1", "P-2"]
WELLS = ["A1", "A2", "A3", "B1", "B2", "B3"]

# Color palettes for consistent styling
PALETTE_BC = ["#3498db", "#e74c3c"]  # Blue/Red for panels B/C
PALETTE_DE = {"P-1": "#e74c3c", "P-2": "#3498db"}  # Red/Blue for plates in D/E

# Standard figure size for all panels (slightly wider than tall)
FIGURE_SIZE = (6.5, 5.5)


def setup_style():
    """Set up publication-quality matplotlib style with sans-serif font."""
    plt.style.use("seaborn-v0_8-whitegrid")
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
    })


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_heatmap_data():
    """Load cell mapping heatmap data from both plates (all wells)."""
    dfs = []
    for plate in PLATES:
        path = BRIEFLOW_OUTPUT / "sbs" / "eval" / "mapping" / f"{plate}__cell_mapping_heatmap_one.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df["plate"] = plate
            dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No heatmap data files found")

    combined = pd.concat(dfs, ignore_index=True)
    # Rename column for easier access
    col_name = [c for c in combined.columns if "fraction" in c.lower()][0]
    combined = combined.rename(columns={col_name: "fraction"})

    return combined


def load_cell_density_data():
    """Load cell density heatmap data from phenotype eval (all wells)."""
    dfs = []
    for plate in PLATES:
        path = BRIEFLOW_OUTPUT / "phenotype" / "eval" / "segmentation" / f"{plate}__cell_density_heatmap.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df["plate"] = plate
            dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No cell density heatmap files found")

    return pd.concat(dfs, ignore_index=True)


def load_mapping_overview():
    """Load mapping overview data from both plates."""
    dfs = []
    for plate in PLATES:
        path = BRIEFLOW_OUTPUT / "sbs" / "eval" / "mapping" / f"{plate}__mapping_overview.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df["plate"] = plate
            dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No mapping overview files found")

    return pd.concat(dfs, ignore_index=True)


def load_barcode_library():
    """Load barcode library and extract unique barcodes."""
    path = CONFIG_DIR / "info_CSM_GW_CRISPRko_library_design.csv"
    if not path.exists():
        raise FileNotFoundError(f"Library file not found: {path}")

    df = pd.read_csv(path)
    # Get unique barcodes from iBAR_1 and iBAR_2 columns
    barcodes = set(df["iBAR_1"].dropna().unique()) | set(df["iBAR_2"].dropna().unique())
    return barcodes


def load_reads_data(sample_frac=0.1):
    """
    Load reads data from parquet files.

    Args:
        sample_frac: Fraction of reads to sample (for speed)

    Returns:
        DataFrame with barcodes, Q-scores, and metadata per well
    """
    q_cols = [f"Q_{i}" for i in range(12)] + ["Q_min"]
    all_data = []

    for plate in PLATES:
        for well in WELLS:
            path = BRIEFLOW_OUTPUT / "sbs" / "parquets" / f"{plate}_W-{well}__reads.parquet"
            if path.exists():
                df = pd.read_parquet(path, columns=["barcode"] + q_cols)
                # Sample for speed
                if sample_frac < 1.0:
                    df = df.sample(frac=sample_frac, random_state=42)
                df["plate"] = plate
                df["well"] = well
                all_data.append(df)

    if not all_data:
        raise FileNotFoundError("No reads parquet files found")

    return pd.concat(all_data, ignore_index=True)


def load_perturbation_counts():
    """Load cell counts per perturbation from aggregated TSV."""
    path = BRIEFLOW_OUTPUT / "aggregate" / "tsvs" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__aggregated.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Aggregated file not found: {path}")

    df = pd.read_csv(path, sep="\t", usecols=["gene_symbol_0", "cell_count"])
    return df


# =============================================================================
# Panel A: Spatial Heatmaps (averaged across ALL wells)
# =============================================================================

def compute_tile_average(data, value_col):
    """
    Compute average across all wells for each tile (using whatever wells have data).

    Args:
        data: DataFrame with 'well', 'plate', 'tile', and value column
        value_col: Name of the column to average

    Returns:
        DataFrame with 'tile' and averaged value
    """
    # Create unique well identifier combining plate and well
    data = data.copy()
    data["plate_well"] = data["plate"] + "_" + data["well"]

    # Get all unique wells
    all_wells = data["plate_well"].unique()
    n_wells = len(all_wells)

    # Average each tile across whatever wells have data for it
    avg_data = data.groupby("tile")[value_col].mean().reset_index()

    print(f"  Averaged {len(avg_data)} tiles across {n_wells} wells")

    return avg_data


def load_sbs_tile_positions():
    """
    Load actual SBS tile positions from metadata.
    Returns DataFrame with tile, x_pos, y_pos.
    """
    path = BRIEFLOW_OUTPUT / "preprocess" / "metadata" / "sbs" / "P-1_W-A1__combined_metadata.parquet"
    df = pd.read_parquet(path)
    # Filter to cycle 1 to get unique tile positions
    df_c1 = df[df["cycle"] == 1]
    return df_c1[["tile", "x_pos", "y_pos"]].drop_duplicates().sort_values("tile")


def create_tile_grid_from_positions(tile_positions):
    """
    Create a tile grid from actual physical positions.

    Args:
        tile_positions: DataFrame with 'tile', 'x_pos', 'y_pos'

    Returns:
        2D numpy array where each cell contains the tile number (or NaN)
    """
    # Get unique x and y positions, sorted
    x_vals = sorted(tile_positions["x_pos"].unique())
    y_vals = sorted(tile_positions["y_pos"].unique())

    # Create mapping from position to grid index
    x_to_col = {x: i for i, x in enumerate(x_vals)}
    y_to_row = {y: i for i, y in enumerate(y_vals)}

    # Create grid
    rows, cols = len(y_vals), len(x_vals)
    grid = np.full((rows, cols), np.nan)

    # Fill grid with tile numbers
    for _, row in tile_positions.iterrows():
        r = y_to_row[row["y_pos"]]
        c = x_to_col[row["x_pos"]]
        grid[r, c] = row["tile"]

    return grid


def create_tile_grid(shape, snake_sites=True):
    """
    Create a tile grid mapping for the given shape.

    Args:
        shape: 'squid_sbs' or 'squid_ph'
        snake_sites: Whether to snake the row order

    Returns:
        2D numpy array where each cell contains the tile number (or NaN for empty cells)
    """
    if shape == "squid_sbs":
        # Load actual positions from metadata
        tile_positions = load_sbs_tile_positions()
        return create_tile_grid_from_positions(tile_positions)
    elif shape == "squid_ph":
        # Spatially accurate layout for Squid microscope phenotype tiles (1732 tiles)
        rows = [
            4, 14, 18, 22, 26, 28, 32, 34, 36, 36, 38, 40, 40, 42, 42, 44, 44,
            46, 46, 46, 46, 46, 48, 48, 48, 48, 46, 46, 46, 46, 46, 44, 44,
            42, 42, 40, 40, 38, 36, 36, 34, 32, 28, 26, 22, 18, 14, 4,
        ]

        r, c = len(rows), max(rows)
        grid = np.empty((r, c))
        grid[:] = np.nan

        next_site = 0
        for row, row_sites in enumerate(rows):
            start = int((c - row_sites) / 2)
            grid[row, start : start + row_sites] = range(next_site, next_site + row_sites)
            next_site += row_sites

        if snake_sites:
            grid[1::2] = grid[1::2, ::-1]

        return grid
    else:
        raise ValueError(f"{shape} shape not implemented")


def plot_panel_a_combined(cell_mapping_data, cell_density_data, output_path):
    """
    Panel A: Combined heatmap with cell mapping (left) and cell density (right).
    Horizontal colorbars at the bottom.
    """
    print("  Computing averages across all wells...")

    # Compute averages (using whatever wells have data for each tile)
    avg_mapping = compute_tile_average(cell_mapping_data, "fraction")
    avg_density = compute_tile_average(cell_density_data, "cell count")

    # Create grids
    grid_sbs = create_tile_grid("squid_sbs")
    grid_ph = create_tile_grid("squid_ph")

    # Fill SBS grid with mapping values
    sbs_values = np.full(grid_sbs.shape, np.nan)
    for _, row in avg_mapping.iterrows():
        tile = int(row["tile"])
        positions = np.where(grid_sbs == tile)
        if len(positions[0]) > 0:
            r, c = positions[0][0], positions[1][0]
            sbs_values[r, c] = row["fraction"]

    # Fill phenotype grid with density values
    ph_values = np.full(grid_ph.shape, np.nan)
    for _, row in avg_density.iterrows():
        tile = int(row["tile"])
        positions = np.where(grid_ph == tile)
        if len(positions[0]) > 0:
            r, c = positions[0][0], positions[1][0]
            ph_values[r, c] = row["cell count"]

    # Create masked arrays to hide NaN values properly
    sbs_masked = np.ma.masked_invalid(sbs_values)
    ph_masked = np.ma.masked_invalid(ph_values)

    # Create figure with two subplots - same size as other panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGURE_SIZE)

    # Set background color for masked/NaN areas
    ax1.set_facecolor('#e0e0e0')
    ax2.set_facecolor('#e0e0e0')

    # Left: Cell mapping (viridis colormap)
    vmin1, vmax1 = avg_mapping["fraction"].min(), avg_mapping["fraction"].max()
    im1 = ax1.imshow(sbs_masked, cmap="viridis", vmin=vmin1, vmax=vmax1, aspect="equal")
    ax1.axis("off")
    ax1.set_anchor('C')  # Center vertically

    # Horizontal colorbar for left plot (descriptive label replaces subtitle)
    cbar1 = fig.colorbar(im1, ax=ax1, orientation="horizontal", pad=0.08, shrink=0.7, aspect=20)
    cbar1.ax.tick_params(labelsize=8)
    cbar1.set_label("Fraction of Cells Mapping to 1 Gene (ISS)", fontsize=8)

    # Right: Cell density (plasma colormap)
    vmin2, vmax2 = avg_density["cell count"].min(), avg_density["cell count"].max()
    im2 = ax2.imshow(ph_masked, cmap="plasma", vmin=vmin2, vmax=vmax2, aspect="equal")
    ax2.axis("off")
    ax2.set_anchor('C')  # Center vertically

    # Horizontal colorbar for right plot
    cbar2 = fig.colorbar(im2, ax=ax2, orientation="horizontal", pad=0.08, shrink=0.7, aspect=20)
    cbar2.ax.tick_params(labelsize=8)
    cbar2.set_label("Cells per Tile (Phenotyping)", fontsize=8)

    # Adjust layout first
    fig.subplots_adjust(top=0.95, bottom=0.15, wspace=0.12, left=0.05, right=0.95)

    # Top level title - position relative to figure
    fig.suptitle("Spatial QC (Averaged Across All Wells)", fontsize=12, fontweight="bold", y=0.99)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


# =============================================================================
# Panel B: Cell Count Boxplot
# =============================================================================

def plot_panel_b_boxplot(data, output_path):
    """
    Panel B: Cell count boxplot for all 4 mapping categories (matching Panel C).
    Colors match category, shapes distinguish plates.
    """
    # Define categories and their corresponding columns (same order as Panel C)
    categories = {
        "1+ Barcode": "1_or_more_barcodes__count",
        "1 Barcode": "1_barcode_cells__count",
        "1+ Gene": "1_or_more_genes__count",
        "1 Gene": "1_gene_cells__count"
    }

    # Colors matching Panel C
    colors = {"1+ Barcode": "#3498db", "1 Barcode": "#2ecc71",
              "1+ Gene": "#9b59b6", "1 Gene": "#e74c3c"}
    color_list = list(colors.values())

    # Markers for plates
    markers = {"Plate 1": "o", "Plate 2": "s"}  # circle and square

    # Prepare data
    plot_data = []
    for _, row in data.iterrows():
        for cat_name, col_name in categories.items():
            plot_data.append({
                "Category": cat_name,
                "Cell Count": row[col_name],
                "Plate": row["plate"].replace("P-", "Plate "),
                "Well": row["well"]
            })

    plot_df = pd.DataFrame(plot_data)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Box plot (lower zorder so points appear on top, hide outliers since we plot all points)
    sns.boxplot(
        data=plot_df,
        x="Category",
        y="Cell Count",
        hue="Category",
        ax=ax,
        palette=colors,
        width=0.6,
        order=list(categories.keys()),
        legend=False,
        zorder=1,
        showfliers=False
    )

    # Add individual points - plot each plate separately with different markers
    # Set higher zorder so points appear on top of boxes
    cat_positions = {cat: i for i, cat in enumerate(categories.keys())}

    for plate, marker in markers.items():
        plate_data = plot_df[plot_df["Plate"] == plate]
        for cat_name in categories.keys():
            cat_data = plate_data[plate_data["Category"] == cat_name]
            # Add jitter
            x_pos = cat_positions[cat_name] + np.random.uniform(-0.15, 0.15, len(cat_data))
            ax.scatter(x_pos, cat_data["Cell Count"],
                      marker=marker, s=70, facecolors=colors[cat_name],
                      edgecolors="black", linewidth=1.2, alpha=1.0,
                      zorder=10,
                      label=plate if cat_name == "1+ Barcode" else "")

    # Format y-axis in millions
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}M"))

    ax.set_ylabel("Number of Cells")
    ax.set_xlabel("")
    ax.set_title("Cell Counts by Mapping Category", fontweight="bold")

    # Create custom legend with black markers (shape only, no color)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
               markeredgecolor='black', markersize=8, label='Plate 1'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='black',
               markeredgecolor='black', markersize=8, label='Plate 2')
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Panel C: Mapping Percentage Bar Plot
# =============================================================================

def plot_panel_c_barplot(data, output_path):
    """
    Panel C: Mapping percentage bar plot.
    """
    # Calculate means and stds
    metrics = {
        "1+ Barcode": data["1_or_more_barcodes__percent"],
        "1 Barcode": data["1_barcode_cells__percent"],
        "1+ Gene": data["1_or_more_genes__percent"],
        "1 Gene": data["1_gene_cells__percent"]
    }

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    labels = list(metrics.keys())
    means = [metrics[k].mean() for k in labels]
    stds = [metrics[k].std() for k in labels]

    colors = ["#3498db", "#2ecc71", "#9b59b6", "#e74c3c"]

    bars = ax.bar(labels, means, yerr=stds, capsize=5,
                  color=colors, edgecolor="black", linewidth=1.5, alpha=0.85)

    ax.set_ylim(0, 100)
    ax.set_ylabel("Percent (%)")
    ax.set_xlabel("")
    ax.set_title("Mapping Success Rates", fontweight="bold")

    # Add value labels on bars
    for bar, mean, std in zip(bars, means, stds):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 2,
                f"{mean:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Panel D: Barcode Prefix Matching
# =============================================================================

def compute_prefix_matching(reads_df, library_barcodes):
    """
    Compute barcode prefix matching rates for each prefix length.

    Args:
        reads_df: DataFrame with 'barcode', 'plate', 'well' columns
        library_barcodes: Set of library barcodes

    Returns:
        DataFrame with prefix matching rates per well
    """
    # Build prefix sets for each length
    prefix_sets = {}
    for length in range(1, 13):
        prefix_sets[length] = set(bc[:length] for bc in library_barcodes if len(bc) >= length)

    results = []

    for (plate, well), group in reads_df.groupby(["plate", "well"]):
        barcodes = group["barcode"].dropna().astype(str)
        total = len(barcodes)

        if total == 0:
            continue

        for length in range(1, 13):
            prefixes = barcodes.str[:length]
            matching = prefixes.isin(prefix_sets[length]).sum()
            frac = matching / total

            results.append({
                "plate": plate,
                "well": well,
                "prefix_length": length,
                "fraction_matching": frac
            })

    return pd.DataFrame(results)


def plot_panel_d_prefix_matching(prefix_data, output_path):
    """
    Panel D: Barcode prefix matching plot per well.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Per-well lines (thin, colored by plate)
    for (plate, well), group in prefix_data.groupby(["plate", "well"]):
        color = PALETTE_DE[plate]
        ax.plot(group["prefix_length"], group["fraction_matching"],
                linewidth=1, alpha=0.5, color=color)

    # Per-plate averages
    for plate in PLATES:
        plate_data = prefix_data[prefix_data["plate"] == plate]
        avg = plate_data.groupby("prefix_length")["fraction_matching"].mean()
        ax.plot(avg.index, avg.values, linewidth=2.5, color=PALETTE_DE[plate],
                marker="o", markersize=5, label=f"{plate.replace('P-', 'Plate ')}")

    # Overall average (thick black line)
    overall = prefix_data.groupby("prefix_length")["fraction_matching"].mean()
    ax.plot(overall.index, overall.values, linewidth=3.5, color="black",
            marker="o", markersize=6, label="All Wells Avg")

    # Random expectation line
    lengths = range(1, 13)
    random_exp = [(1/4)**n for n in lengths]
    ax.plot(lengths, random_exp, linewidth=2, color="gray", linestyle="--",
            marker=".", markersize=4, alpha=0.7, label="Random (1/4)ⁿ")

    ax.set_xlabel("Barcode Prefix Length (bases)")
    ax.set_ylabel("Fraction of Reads Matching Library")
    ax.set_title("Barcode Prefix Matching", fontweight="bold")
    ax.set_xlim(1, 12)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(1, 13))
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Panel E: Q-Scores Across Cycles
# =============================================================================

def plot_panel_e_qscores(reads_df, output_path):
    """
    Panel E: Q-scores across cycles per well.
    Style matches Panel D.
    """
    q_cols = [f"Q_{i}" for i in range(12)]
    cycles = range(1, 13)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Per-well lines (thin, colored by plate)
    for (plate, well), group in reads_df.groupby(["plate", "well"]):
        color = PALETTE_DE[plate]
        means = [group[q].mean() for q in q_cols]
        ax.plot(cycles, means, linewidth=1, alpha=0.5, color=color)

    # Per-plate averages
    for plate in PLATES:
        plate_data = reads_df[reads_df["plate"] == plate]
        means = [plate_data[q].mean() for q in q_cols]
        ax.plot(cycles, means, linewidth=2.5, color=PALETTE_DE[plate],
                marker="o", markersize=5, label=f"{plate.replace('P-', 'Plate ')}")

    # Overall average (thick black line)
    overall_means = [reads_df[q].mean() for q in q_cols]
    ax.plot(cycles, overall_means, linewidth=3.5, color="black",
            marker="o", markersize=6, label="All Wells Avg")

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Q-score")
    ax.set_title("Q-Scores Across Cycles", fontweight="bold")
    ax.set_xlim(1, 12)
    ax.set_xticks(cycles)
    ax.legend(loc="upper right", fontsize=9)

    # Set y-axis based on data range
    y_min = min(overall_means) - 0.05
    y_max = max(overall_means) + 0.05
    ax.set_ylim(max(0, y_min), min(1, y_max))

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Panel F: KDE of Cells per Perturbation
# =============================================================================

def plot_panel_f_kde(data, output_path):
    """
    Panel F: KDE of cells per perturbation.
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Plot KDE
    sns.kdeplot(data["cell_count"], ax=ax, fill=True, color="#1abc9c", alpha=0.6, linewidth=2)

    # Add median reference line
    median_val = data["cell_count"].median()
    ax.axvline(median_val, color="#e74c3c", linestyle="--", linewidth=2.5,
               label=f"Median = {median_val:.0f}")

    ax.set_xlabel("Number of Cells per Perturbation")
    ax.set_ylabel("Density")
    ax.set_title("Cells per Perturbation Distribution", fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)

    # Remove y-axis tick labels (keep axis label)
    ax.set_yticklabels([])

    # Set x-axis limit to focus on main distribution
    q99 = data["cell_count"].quantile(0.99)
    ax.set_xlim(0, q99)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Generate QC panels for Lander presentation")
    parser.add_argument("--output-dir", type=str, default="figures",
                       help="Output directory for figures")
    parser.add_argument("--sample-reads", type=float, default=0.1,
                       help="Fraction of reads to sample (0-1, default 0.1)")
    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    setup_style()

    print("=" * 60)
    print("Generating QC Panels")
    print("=" * 60)

    # Load data for Panels A-C (from pre-computed TSVs)
    print("\nLoading data for Panels A-C...")
    print("  Loading cell mapping heatmap data...")
    heatmap_data = load_heatmap_data()

    print("  Loading cell density heatmap data...")
    density_data = load_cell_density_data()

    print("  Loading mapping overview...")
    mapping_overview = load_mapping_overview()

    # Load data for Panels D-E (from parquets)
    print(f"\nLoading reads data (sampling {args.sample_reads*100:.0f}%)...")
    reads_data = load_reads_data(sample_frac=args.sample_reads)

    print("  Loading barcode library...")
    library_barcodes = load_barcode_library()

    # Load data for Panel F
    print("  Loading perturbation counts...")
    perturbation_data = load_perturbation_counts()

    # Generate panels
    print("\n" + "=" * 60)
    print("Generating panels...")
    print("=" * 60)

    # Panel A: Combined heatmap (cell mapping + cell density, averaged across all wells)
    print("\nPanel A: Combined heatmap (cell mapping + cell density)...")
    plot_panel_a_combined(heatmap_data, density_data, output_dir / "qc_panel_a_heatmap.png")

    # Panel B: Cell count boxplot
    print("Panel B: Cell count boxplot...")
    plot_panel_b_boxplot(mapping_overview, output_dir / "qc_panel_b_boxplot.png")

    # Panel C: Mapping percentage bar plot
    print("Panel C: Mapping percentage bar plot...")
    plot_panel_c_barplot(mapping_overview, output_dir / "qc_panel_c_barplot.png")

    # Panel D: Barcode prefix matching
    print("Panel D: Computing barcode prefix matching...")
    prefix_data = compute_prefix_matching(reads_data, library_barcodes)
    plot_panel_d_prefix_matching(prefix_data, output_dir / "qc_panel_d_prefix_matching.png")

    # Panel E: Q-scores
    print("Panel E: Q-scores across cycles...")
    plot_panel_e_qscores(reads_data, output_dir / "qc_panel_e_qscores.png")

    # Panel F: KDE
    print("Panel F: Cells per perturbation KDE...")
    plot_panel_f_kde(perturbation_data, output_dir / "qc_panel_f_kde.png")

    print("\n" + "=" * 60)
    print(f"All panels saved to: {output_dir}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
