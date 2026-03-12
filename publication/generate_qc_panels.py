#!/usr/bin/env python3
"""
Generate QC panels for Lander presentation.

Creates 10 publication-quality QC panel figures:
(a) Spatial heatmap of cell mapping to single gene
(b) Cell count boxplot (1 Gene vs >1 Gene)
(c) Mapping percentage bar plot
(d) Barcode prefix matching per well
(e) Q-scores across cycles per well
(f) KDE of cells per perturbation
(g) Cell retention funnel through deduplication pipeline
(h) Total cell retention barplot (summed across all wells)
(i) Error-corrected prefix matching (raw vs corrected)
(j) Mapping rate and cell retention vs read quality (Q_min)

Usage:
    python generate_qc_panels.py [--output-dir figures/] [--sample-reads 0.1]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Import plot_plate_heatmap from generate_heatmaps

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
    """Set up publication-quality matplotlib style with Myriad Pro font."""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Myriad Pro", "Arial", "DejaVu Sans", "Helvetica"],
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
            "savefig.dpi": 600,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig, output_path, dpi=600):
    """Save figure as PNG (600 dpi), PDF, and SVG."""
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(str(output_path).replace(".png", ".pdf"), bbox_inches="tight")
    fig.savefig(str(output_path).replace(".png", ".svg"), bbox_inches="tight")
    print(f"Saved: {output_path} (.png, .pdf, .svg)")


# =============================================================================
# Data Loading Functions
# =============================================================================


def load_heatmap_data():
    """Load cell mapping heatmap data from both plates (all wells)."""
    dfs = []
    for plate in PLATES:
        path = (
            BRIEFLOW_OUTPUT
            / "sbs"
            / "eval"
            / "mapping"
            / f"{plate}__cell_mapping_heatmap_one.tsv"
        )
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
        path = (
            BRIEFLOW_OUTPUT
            / "phenotype"
            / "eval"
            / "segmentation"
            / f"{plate}__cell_density_heatmap.tsv"
        )
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
        path = (
            BRIEFLOW_OUTPUT
            / "sbs"
            / "eval"
            / "mapping"
            / f"{plate}__mapping_overview.tsv"
        )
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
    barcodes = set(df["iBAR_2"].dropna().unique())
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
            path = (
                BRIEFLOW_OUTPUT
                / "sbs"
                / "parquets"
                / f"{plate}_W-{well}__reads.parquet"
            )
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


def load_retention_funnel_data():
    """
    Load cell counts at each pipeline stage for the retention funnel.

    Stages: Phenotype Cells, SBS Cells, Matched, After SBS Dedup,
    After Phenotype Dedup, 1+ Barcode, 1 Gene.

    Returns DataFrame with columns: plate, well, stage, count.
    """
    rows = []

    for plate in PLATES:
        # Phenotype segmentation overview → final_cells per well
        ph_path = (
            BRIEFLOW_OUTPUT
            / "phenotype"
            / "eval"
            / "segmentation"
            / f"{plate}__segmentation_overview.tsv"
        )
        ph_overview = pd.read_csv(ph_path, sep="\t")

        # SBS segmentation overview → final_cells per well
        sbs_path = (
            BRIEFLOW_OUTPUT
            / "sbs"
            / "eval"
            / "segmentation"
            / f"{plate}__segmentation_overview.tsv"
        )
        sbs_overview = pd.read_csv(sbs_path, sep="\t")

        for well in WELLS:
            ph_row = ph_overview[ph_overview["well"] == well].iloc[0]
            sbs_row = sbs_overview[sbs_overview["well"] == well].iloc[0]

            # Stage 1 & 2: initial segmented cells
            rows.append(
                {"plate": plate, "well": well, "stage": "Phenotype\nCells",
                 "count": ph_row["final_cells"]}
            )
            rows.append(
                {"plate": plate, "well": well, "stage": "SBS\nCells",
                 "count": sbs_row["final_cells"]}
            )

            # Stage 3-5: from deduplication stats
            dedup_path = (
                BRIEFLOW_OUTPUT
                / "merge"
                / "eval"
                / f"{plate}_W-{well}__deduplication_stats.tsv"
            )
            dedup = pd.read_csv(dedup_path, sep="\t")
            stage_map = {
                "Initial": "Matched",
                "After phenotype dedup": "Deduped",
            }
            for _, d_row in dedup.iterrows():
                label = stage_map.get(d_row["stage"])
                if label:
                    rows.append(
                        {"plate": plate, "well": well, "stage": label,
                         "count": d_row["total_cells"]}
                    )

            # Stage 6 & 7: barcode/gene counts from merge_final parquet
            merge_path = (
                BRIEFLOW_OUTPUT
                / "merge"
                / "parquets"
                / f"{plate}_W-{well}__merge_final.parquet"
            )
            merge_df = pd.read_parquet(
                merge_path, columns=["cell_barcode_0", "mapped_single_gene"]
            )
            rows.append(
                {"plate": plate, "well": well, "stage": "1+ Barcode",
                 "count": merge_df["cell_barcode_0"].notna().sum()}
            )
            rows.append(
                {"plate": plate, "well": well, "stage": "1 Gene",
                 "count": merge_df["mapped_single_gene"].sum()}
            )

    return pd.DataFrame(rows)


def load_perturbation_counts():
    """Load cell counts per perturbation from aggregated TSV."""
    path = (
        BRIEFLOW_OUTPUT
        / "aggregate"
        / "tsvs"
        / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__aggregated.tsv"
    )
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
    path = (
        BRIEFLOW_OUTPUT
        / "preprocess"
        / "metadata"
        / "sbs"
        / "P-1_W-A1__combined_metadata.parquet"
    )
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
            4,
            14,
            18,
            22,
            26,
            28,
            32,
            34,
            36,
            36,
            38,
            40,
            40,
            42,
            42,
            44,
            44,
            46,
            46,
            46,
            46,
            46,
            48,
            48,
            48,
            48,
            46,
            46,
            46,
            46,
            46,
            44,
            44,
            42,
            42,
            40,
            40,
            38,
            36,
            36,
            34,
            32,
            28,
            26,
            22,
            18,
            14,
            4,
        ]

        r, c = len(rows), max(rows)
        grid = np.empty((r, c))
        grid[:] = np.nan

        next_site = 0
        for row, row_sites in enumerate(rows):
            start = int((c - row_sites) / 2)
            grid[row, start : start + row_sites] = range(
                next_site, next_site + row_sites
            )
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
    ax1.set_facecolor("#e0e0e0")
    ax2.set_facecolor("#e0e0e0")

    # Left: Cell mapping (viridis colormap)
    vmin1, vmax1 = avg_mapping["fraction"].min(), avg_mapping["fraction"].max()
    im1 = ax1.imshow(sbs_masked, cmap="viridis", vmin=vmin1, vmax=vmax1, aspect="equal")
    ax1.axis("off")
    ax1.set_anchor("C")  # Center vertically

    # Horizontal colorbar for left plot (descriptive label replaces subtitle)
    cbar1 = fig.colorbar(
        im1, ax=ax1, orientation="horizontal", pad=0.08, shrink=0.7, aspect=20
    )
    cbar1.ax.tick_params(labelsize=8)
    cbar1.set_label("Fraction of Cells Mapping to 1 Gene (ISS)", fontsize=8)

    # Right: Cell density (plasma colormap)
    vmin2, vmax2 = avg_density["cell count"].min(), avg_density["cell count"].max()
    im2 = ax2.imshow(ph_masked, cmap="plasma", vmin=vmin2, vmax=vmax2, aspect="equal")
    ax2.axis("off")
    ax2.set_anchor("C")  # Center vertically

    # Horizontal colorbar for right plot
    cbar2 = fig.colorbar(
        im2, ax=ax2, orientation="horizontal", pad=0.08, shrink=0.7, aspect=20
    )
    cbar2.ax.tick_params(labelsize=8)
    cbar2.set_label("Cells per Tile (Phenotyping)", fontsize=8)

    # Adjust layout first
    fig.subplots_adjust(top=0.95, bottom=0.15, wspace=0.12, left=0.05, right=0.95)

    # Top level title - position relative to figure
    fig.suptitle(
        "Spatial QC (Averaged Across All Wells)", fontsize=12, fontweight="bold", y=0.99
    )
    save_figure(fig, output_path)
    plt.close(fig)


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
        "1 Gene": "1_gene_cells__count",
    }

    # Colors matching Panel C
    colors = {
        "1+ Barcode": "#3498db",
        "1 Barcode": "#2ecc71",
        "1+ Gene": "#9b59b6",
        "1 Gene": "#e74c3c",
    }

    # Markers for plates
    markers = {"Plate 1": "o", "Plate 2": "s"}  # circle and square

    # Prepare data
    plot_data = []
    for _, row in data.iterrows():
        for cat_name, col_name in categories.items():
            plot_data.append(
                {
                    "Category": cat_name,
                    "Cell Count": row[col_name],
                    "Plate": row["plate"].replace("P-", "Plate "),
                    "Well": row["well"],
                }
            )

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
        showfliers=False,
    )

    # Add individual points - plot each plate separately with different markers
    # Set higher zorder so points appear on top of boxes
    cat_positions = {cat: i for i, cat in enumerate(categories.keys())}

    for plate, marker in markers.items():
        plate_data = plot_df[plot_df["Plate"] == plate]
        for cat_name in categories.keys():
            cat_data = plate_data[plate_data["Category"] == cat_name]
            # Add jitter
            x_pos = cat_positions[cat_name] + np.random.uniform(
                -0.15, 0.15, len(cat_data)
            )
            ax.scatter(
                x_pos,
                cat_data["Cell Count"],
                marker=marker,
                s=70,
                facecolors=colors[cat_name],
                edgecolors="black",
                linewidth=1.2,
                alpha=1.0,
                zorder=10,
                label=plate if cat_name == "1+ Barcode" else "",
            )

    # Format y-axis in millions
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x / 1e6:.1f}M"))

    ax.set_ylabel("Number of Cells")
    ax.set_xlabel("")
    ax.set_title("Cell Counts by Mapping Category", fontweight="bold")

    # Create custom legend with black markers (shape only, no color)
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=8,
            label="Plate 1",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=8,
            label="Plate 2",
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


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
        "1 Gene": data["1_gene_cells__percent"],
    }

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    labels = list(metrics.keys())
    means = [metrics[k].mean() for k in labels]
    stds = [metrics[k].std() for k in labels]

    colors = ["#3498db", "#2ecc71", "#9b59b6", "#e74c3c"]

    bars = ax.bar(
        labels,
        means,
        yerr=stds,
        capsize=5,
        color=colors,
        edgecolor="black",
        linewidth=1.5,
        alpha=0.85,
    )

    ax.set_ylim(0, 100)
    ax.set_ylabel("Percent (%)")
    ax.set_xlabel("")
    ax.set_title("Mapping Success Rates", fontweight="bold")

    # Add value labels on bars
    for bar, mean, std in zip(bars, means, stds):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + std + 2,
            f"{mean:.1f}%",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


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
        prefix_sets[length] = set(
            bc[:length] for bc in library_barcodes if len(bc) >= length
        )

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

            results.append(
                {
                    "plate": plate,
                    "well": well,
                    "prefix_length": length,
                    "fraction_matching": frac,
                }
            )

    df = pd.DataFrame(results)
    return df, prefix_sets


def plot_panel_d_prefix_matching(prefix_data, prefix_sets, output_path):
    """
    Panel D: Barcode prefix matching plot per well.
    """
    prefix_set_sizes = {length: len(ps) for length, ps in prefix_sets.items()}
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Per-well lines (thin, colored by plate)
    for (plate, well), group in prefix_data.groupby(["plate", "well"]):
        color = PALETTE_DE[plate]
        ax.plot(
            group["prefix_length"],
            group["fraction_matching"],
            linewidth=1,
            alpha=0.5,
            color=color,
        )

    # Per-plate averages
    for plate in PLATES:
        plate_data = prefix_data[prefix_data["plate"] == plate]
        avg = plate_data.groupby("prefix_length")["fraction_matching"].mean()
        ax.plot(
            avg.index,
            avg.values,
            linewidth=2.5,
            color=PALETTE_DE[plate],
            marker="o",
            markersize=5,
            label=f"{plate.replace('P-', 'Plate ')}",
        )

    # Overall average (thick black line)
    overall = prefix_data.groupby("prefix_length")["fraction_matching"].mean()
    ax.plot(
        overall.index,
        overall.values,
        linewidth=3.5,
        color="black",
        marker="o",
        markersize=6,
        label="All Wells Avg",
    )

    # Random expectation line: fraction of possible prefixes covered by library
    lengths = range(1, 13)
    random_exp = [prefix_set_sizes.get(n, 1) / (4**n) for n in lengths]
    ax.plot(
        lengths,
        random_exp,
        linewidth=2,
        color="gray",
        linestyle="--",
        marker=".",
        markersize=4,
        alpha=0.7,
        label="Random Expected",
    )

    ax.set_xlabel("Barcode Prefix Length (bases)")
    ax.set_ylabel("Fraction of Reads Matching Library")
    ax.set_title("Barcode Prefix Matching", fontweight="bold")
    ax.set_xlim(1, 12)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(1, 13))
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


def plot_panel_d_zoomed(prefix_data, prefix_sets, output_path):
    """
    Panel D (zoomed): Same as default but zoomed into cycles 7-12.
    """
    prefix_set_sizes = {length: len(ps) for length, ps in prefix_sets.items()}
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Per-well lines (thin, colored by plate)
    for (plate, well), group in prefix_data.groupby(["plate", "well"]):
        color = PALETTE_DE[plate]
        ax.plot(
            group["prefix_length"],
            group["fraction_matching"],
            linewidth=1,
            alpha=0.5,
            color=color,
        )

    # Per-plate averages
    for plate in PLATES:
        plate_data = prefix_data[prefix_data["plate"] == plate]
        avg = plate_data.groupby("prefix_length")["fraction_matching"].mean()
        ax.plot(
            avg.index,
            avg.values,
            linewidth=2.5,
            color=PALETTE_DE[plate],
            marker="o",
            markersize=5,
            label=f"{plate.replace('P-', 'Plate ')}",
        )

    # Overall average (thick black line)
    overall = prefix_data.groupby("prefix_length")["fraction_matching"].mean()
    ax.plot(
        overall.index,
        overall.values,
        linewidth=3.5,
        color="black",
        marker="o",
        markersize=6,
        label="All Wells Avg",
    )

    # Random expectation line
    lengths = range(1, 13)
    random_exp = [prefix_set_sizes.get(n, 1) / (4**n) for n in lengths]
    ax.plot(
        lengths,
        random_exp,
        linewidth=2,
        color="gray",
        linestyle="--",
        marker=".",
        markersize=4,
        alpha=0.7,
        label="Random Expected",
    )

    ax.set_xlabel("Barcode Prefix Length (bases)")
    ax.set_ylabel("Fraction of Reads Matching Library")
    ax.set_title("Barcode Prefix Matching (Zoomed)", fontweight="bold")
    ax.set_xlim(7, 12)
    ax.set_xticks(range(7, 13))
    ax.legend(loc="upper right", fontsize=9)

    # Auto y-axis to zoom into the separation region
    zoomed = prefix_data[prefix_data["prefix_length"] >= 7]
    y_min = max(0, zoomed["fraction_matching"].min() - 0.05)
    y_max = min(1.05, zoomed["fraction_matching"].max() + 0.05)
    ax.set_ylim(y_min, y_max)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


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
        ax.plot(
            cycles,
            means,
            linewidth=2.5,
            color=PALETTE_DE[plate],
            marker="o",
            markersize=5,
            label=f"{plate.replace('P-', 'Plate ')}",
        )

    # Overall average (thick black line)
    overall_means = [reads_df[q].mean() for q in q_cols]
    ax.plot(
        cycles,
        overall_means,
        linewidth=3.5,
        color="black",
        marker="o",
        markersize=6,
        label="All Wells Avg",
    )

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
    save_figure(fig, output_path)
    plt.close()


# =============================================================================
# Panel F: Cells per Perturbation by Category
# =============================================================================

# Category colors
CATEGORY_COLORS = {
    "Essential genes": "#d67777",
    "Non-essential genes": "#7dd3c0",
    "Intergenic controls": "#999999",
    "Noncutting controls": "#666666",
    "Olfactory receptor controls": "#444444",
}

CATEGORY_ORDER = [
    "Essential genes",
    "Non-essential genes",
    "Intergenic controls",
    "Noncutting controls",
    "Olfactory receptor controls",
]


def categorize_perturbations(data, essential_genes):
    """Categorize perturbations into essential, non-essential, and control types."""
    data = data.copy()
    gs = data["gene_symbol_0"]

    conditions = [
        gs.str.startswith("nontargeting_intergenic", na=False),
        gs.str.startswith("nontargeting_noncutting", na=False),
        gs.str.startswith("nontargeting_or", na=False),
        gs.isin(essential_genes),
    ]
    choices = [
        "Intergenic controls",
        "Noncutting controls",
        "Olfactory receptor controls",
        "Essential genes",
    ]
    data["category"] = np.select(conditions, choices, default="Non-essential genes")
    return data


def plot_panel_f_histogram(data, essential_genes, output_path):
    """
    Panel F: Split histogram of cells per perturbation by category.
    Non-essential plotted first, then essential on top so it's visible.
    Controls plotted last on top.
    """
    data = categorize_perturbations(data, essential_genes)
    bins = np.arange(0, data["cell_count"].quantile(0.99) + 20, 20)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Plot order: non-essential first (back), then essential, then controls on top
    plot_order = [
        "Non-essential genes",
        "Essential genes",
        "Intergenic controls",
        "Noncutting controls",
        "Olfactory receptor controls",
    ]

    for category in plot_order:
        subset = data[data["category"] == category]
        if len(subset) == 0:
            continue

        ax.hist(
            subset["cell_count"],
            bins=bins,
            alpha=0.6,
            color=CATEGORY_COLORS[category],
            label=f"{category} (n={len(subset)})",
            edgecolor="white",
            linewidth=0.5,
        )

        # Mean dashed line
        mean_val = subset["cell_count"].mean()
        ax.axvline(
            mean_val,
            color=CATEGORY_COLORS[category],
            linestyle="--",
            linewidth=1.5,
            alpha=0.8,
        )

    ax.set_xlabel("Cells per Perturbation")
    ax.set_ylabel("Count")
    ax.set_title("Cells per Perturbation by Category", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlim(0, bins[-1])

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


# =============================================================================
# Panel G: Cell Retention Funnel
# =============================================================================


def plot_panel_g_retention_funnel(data, output_path):
    """
    Panel G: Cell retention through the full pipeline.
    Boxplot with individual well points (matching Panel B style).
    Stages: Phenotype, SBS, Matched, After SBS Dedup, After Phenotype Dedup,
    1+ Barcode, 1 Gene.
    """
    from matplotlib.lines import Line2D

    stage_order = [
        "Phenotype\nCells", "SBS\nCells", "Matched",
        "Deduped", "1+ Barcode", "1 Gene",
    ]
    stage_colors = [
        "#3498db", "#2980b9", "#8e44ad",
        "#e74c3c", "#e67e22", "#f39c12",
    ]
    color_map = dict(zip(stage_order, stage_colors))
    markers = {"Plate 1": "o", "Plate 2": "s"}

    # Prepare plot data
    data = data.copy()
    data["Plate"] = data["plate"].str.replace("P-", "Plate ")

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Boxplot
    sns.boxplot(
        data=data,
        x="stage",
        y="count",
        hue="stage",
        ax=ax,
        palette=color_map,
        width=0.6,
        order=stage_order,
        legend=False,
        zorder=1,
        showfliers=False,
    )

    # Individual well points
    stage_positions = {s: i for i, s in enumerate(stage_order)}
    for plate, marker in markers.items():
        plate_data = data[data["Plate"] == plate]
        for stage in stage_order:
            subset = plate_data[plate_data["stage"] == stage]
            x_pos = stage_positions[stage] + np.random.uniform(
                -0.15, 0.15, len(subset)
            )
            ax.scatter(
                x_pos,
                subset["count"],
                marker=marker,
                s=70,
                facecolors=color_map[stage],
                edgecolors="black",
                linewidth=1.2,
                alpha=1.0,
                zorder=10,
            )

    # Format y-axis in millions
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x / 1e6:.1f}M"))

    ax.set_ylabel("Number of Cells")
    ax.set_xlabel("")
    ax.set_title("Cell Retention Through Pipeline (Per Well)", fontweight="bold")

    # Plate shape legend (matching Panel B)
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=8,
            label="Plate 1",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=8,
            label="Plate 2",
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


def plot_panel_h_retention_totals(data, output_path):
    """
    Panel H: Total cell retention barplot (summed across all wells).
    Same stages as Panel G but showing grand totals.
    """
    stage_order = [
        "Phenotype\nCells", "SBS\nCells", "Matched",
        "Deduped", "1+ Barcode", "1 Gene",
    ]
    stage_colors = [
        "#3498db", "#2980b9", "#8e44ad",
        "#e74c3c", "#e67e22", "#f39c12",
    ]

    # Sum across all wells per stage
    totals = data.groupby("stage")["count"].sum()
    totals = totals.reindex(stage_order)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    bars = ax.bar(
        stage_order,
        totals.values,
        color=stage_colors,
        edgecolor="black",
        linewidth=1.5,
        alpha=0.85,
    )

    # Value labels on bars (matching Panel C style)
    for bar, val in zip(bars, totals.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + totals.max() * 0.01,
            f"{val / 1e6:.1f}M",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x / 1e6:.0f}M"))
    ax.set_ylabel("Total Cells")
    ax.set_xlabel("")
    ax.set_title("Cell Retention Through Pipeline (Total)", fontweight="bold")

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


# =============================================================================
# Panel I: Error-Corrected Prefix Matching
# =============================================================================


def error_correct_barcodes(barcodes_series, library_set, max_dist=2):
    """
    Error-correct barcodes using the same logic as the brieflow pipeline.

    Mirrors brieflow's error_correct_reads (sbs/call_cells.py): for each barcode,
    finds the unique closest library match within Hamming distance max_dist.
    Uses neighbor generation + set lookup for speed (equivalent to the pipeline's
    brute-force distance matrix but tractable for large read sets).

    Returns Series of corrected barcodes (uncorrected if no unique match).
    """
    bases = "ACGT"

    # Work on unique barcodes for speed
    unique_bcs = list(barcodes_series.dropna().astype(str).unique())
    correction_map = {}

    for bc in unique_bcs:
        if bc in library_set:
            correction_map[bc] = bc
            continue

        # Distance 1 neighbors
        matches = set()
        for i in range(len(bc)):
            for b in bases:
                if b != bc[i]:
                    variant = bc[:i] + b + bc[i + 1 :]
                    if variant in library_set:
                        matches.add(variant)
        if len(matches) == 1:
            correction_map[bc] = matches.pop()
            continue
        elif len(matches) > 1:
            correction_map[bc] = bc  # ambiguous
            continue

        if max_dist < 2:
            correction_map[bc] = bc
            continue

        # Distance 2 neighbors
        for i in range(len(bc)):
            for j in range(i + 1, len(bc)):
                for bi in bases:
                    for bj in bases:
                        if bi != bc[i] or bj != bc[j]:
                            variant = bc[:i] + bi + bc[i + 1 : j] + bj + bc[j + 1 :]
                            if variant in library_set:
                                matches.add(variant)
        if len(matches) == 1:
            correction_map[bc] = matches.pop()
        else:
            correction_map[bc] = bc  # ambiguous or no match

    return barcodes_series.map(correction_map)


def compute_corrected_prefix_matching(reads_df, library_barcodes):
    """
    Compute prefix matching on error-corrected barcodes.

    Returns DataFrame with prefix matching rates per well, with a 'corrected'
    column indicating whether the data is raw or corrected.
    """
    # Build prefix sets
    prefix_sets = {}
    for length in range(1, 13):
        prefix_sets[length] = set(
            bc[:length] for bc in library_barcodes if len(bc) >= length
        )

    print("  Error-correcting barcodes (this may take a minute)...")
    corrected_barcodes = error_correct_barcodes(
        reads_df["barcode"].dropna().astype(str), library_barcodes
    )

    results = []
    for (plate, well), group in reads_df.groupby(["plate", "well"]):
        raw_barcodes = group["barcode"].dropna().astype(str)
        corr_barcodes = corrected_barcodes.loc[raw_barcodes.index]
        total = len(raw_barcodes)
        if total == 0:
            continue

        for length in range(1, 13):
            ps = prefix_sets[length]
            # Raw
            raw_frac = raw_barcodes.str[:length].isin(ps).sum() / total
            results.append(
                {
                    "plate": plate,
                    "well": well,
                    "prefix_length": length,
                    "fraction_matching": raw_frac,
                    "corrected": False,
                }
            )
            # Corrected
            corr_frac = corr_barcodes.str[:length].isin(ps).sum() / total
            results.append(
                {
                    "plate": plate,
                    "well": well,
                    "prefix_length": length,
                    "fraction_matching": corr_frac,
                    "corrected": True,
                }
            )

    return pd.DataFrame(results), prefix_sets


def plot_panel_i_corrected_prefix(data, prefix_sets, output_path):
    """
    Panel I: Prefix matching comparing raw vs error-corrected barcodes.
    Style matches Panel D.
    """
    prefix_set_sizes = {length: len(ps) for length, ps in prefix_sets.items()}
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Per-well thin lines (both raw and corrected)
    for (plate, well, corrected), group in data.groupby(
        ["plate", "well", "corrected"]
    ):
        color = PALETTE_DE[plate]
        alpha = 0.3 if not corrected else 0.15
        linestyle = "-" if not corrected else "--"
        ax.plot(
            group["prefix_length"],
            group["fraction_matching"],
            linewidth=0.8,
            alpha=alpha,
            color=color,
            linestyle=linestyle,
        )

    # Overall averages: raw vs corrected
    raw_data = data[~data["corrected"]]
    corr_data = data[data["corrected"]]

    raw_avg = raw_data.groupby("prefix_length")["fraction_matching"].mean()
    corr_avg = corr_data.groupby("prefix_length")["fraction_matching"].mean()

    ax.plot(
        raw_avg.index,
        raw_avg.values,
        linewidth=3.5,
        color="black",
        marker="o",
        markersize=6,
        label="Raw (All Wells Avg)",
    )
    ax.plot(
        corr_avg.index,
        corr_avg.values,
        linewidth=3.5,
        color="black",
        marker="s",
        markersize=6,
        linestyle="--",
        label="Error-Corrected (All Wells Avg)",
    )

    # Random expectation
    lengths = range(1, 13)
    random_exp = [prefix_set_sizes.get(n, 1) / (4**n) for n in lengths]
    ax.plot(
        lengths,
        random_exp,
        linewidth=2,
        color="gray",
        linestyle=":",
        marker=".",
        markersize=4,
        alpha=0.7,
        label="Random Expected",
    )

    ax.set_xlabel("Barcode Prefix Length (bases)")
    ax.set_ylabel("Fraction of Reads Matching Library")
    ax.set_title("Prefix Matching: Raw vs Error-Corrected", fontweight="bold")
    ax.set_xlim(1, 12)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(1, 13))
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


# =============================================================================
# Panel J: Mapping Rate & Retention vs Read Quality
# =============================================================================


def compute_mapping_vs_quality(reads_df, library_barcodes):
    """
    Compute mapping rate and cell retention as a function of Q_min threshold.

    Returns DataFrame with columns: plate, well, q_threshold,
    mapping_rate, retention_fraction.
    """
    thresholds = np.arange(0, 0.31, 0.01)
    library_set = set(library_barcodes)

    results = []
    for (plate, well), group in reads_df.groupby(["plate", "well"]):
        barcodes = group["barcode"].dropna().astype(str)
        q_min = group.loc[barcodes.index, "Q_min"]
        total = len(barcodes)
        if total == 0:
            continue

        for thresh in thresholds:
            mask = q_min >= thresh
            n_passing = mask.sum()
            if n_passing == 0:
                continue
            n_mapped = barcodes[mask].isin(library_set).sum()

            results.append(
                {
                    "plate": plate,
                    "well": well,
                    "q_threshold": round(thresh, 3),
                    "mapping_rate": n_mapped / n_passing,
                    "retention_fraction": n_passing / total,
                }
            )

    return pd.DataFrame(results)


def plot_panel_j_mapping_vs_quality(data, output_path):
    """
    Panel J: Mapping rate (left axis) and retention (right axis) vs Q_min.
    Dual-axis plot, style matches Panel D/E.
    """
    fig, ax1 = plt.subplots(figsize=FIGURE_SIZE)
    ax2 = ax1.twinx()

    # Per-well thin lines
    for (plate, well), group in data.groupby(["plate", "well"]):
        color = PALETTE_DE[plate]
        ax1.plot(
            group["q_threshold"],
            group["mapping_rate"],
            linewidth=0.8,
            alpha=0.3,
            color=color,
        )
        ax2.plot(
            group["q_threshold"],
            group["retention_fraction"],
            linewidth=0.8,
            alpha=0.15,
            color=color,
            linestyle="--",
        )

    # Per-plate averages
    for plate in PLATES:
        plate_data = data[data["plate"] == plate]
        avg = plate_data.groupby("q_threshold")[
            ["mapping_rate", "retention_fraction"]
        ].mean()
        ax1.plot(
            avg.index,
            avg["mapping_rate"],
            linewidth=2.5,
            color=PALETTE_DE[plate],
            marker="o",
            markersize=4,
            label=f"{plate.replace('P-', 'Plate ')} Mapping",
        )
        ax2.plot(
            avg.index,
            avg["retention_fraction"],
            linewidth=2.5,
            color=PALETTE_DE[plate],
            marker="s",
            markersize=4,
            linestyle="--",
            label=f"{plate.replace('P-', 'Plate ')} Retention",
        )

    # Overall averages
    overall = data.groupby("q_threshold")[
        ["mapping_rate", "retention_fraction"]
    ].mean()
    ax1.plot(
        overall.index,
        overall["mapping_rate"],
        linewidth=3.5,
        color="black",
        marker="o",
        markersize=6,
        label="Avg Mapping Rate",
    )
    ax2.plot(
        overall.index,
        overall["retention_fraction"],
        linewidth=3.5,
        color="black",
        marker="s",
        markersize=6,
        linestyle="--",
        label="Avg Retention",
    )

    ax1.set_xlabel("Read Quality Threshold (Q_min)")
    ax1.set_ylabel("Mapping Rate (fraction matching library)")
    ax2.set_ylabel("Fraction of Reads Retained")
    ax1.set_title("Mapping Rate & Retention vs Read Quality", fontweight="bold")

    ax1.set_ylim(0, 1.05)
    ax2.set_ylim(0, 1.05)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=8)

    # Keep right spine visible for dual axis
    ax2.spines["right"].set_visible(True)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Generate QC panels for Lander presentation"
    )
    parser.add_argument(
        "--output-dir", type=str, default="figures", help="Output directory for figures"
    )
    parser.add_argument(
        "--sample-reads",
        type=float,
        default=0.1,
        help="Fraction of reads to sample (0-1, default 0.1)",
    )
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
    print(f"\nLoading reads data (sampling {args.sample_reads * 100:.0f}%)...")
    reads_data = load_reads_data(sample_frac=args.sample_reads)

    print("  Loading barcode library...")
    library_barcodes = load_barcode_library()

    # Load data for Panel F
    print("  Loading perturbation counts...")
    perturbation_data = load_perturbation_counts()

    print("  Loading essential genes list...")
    essential_genes_file = Path(__file__).parent / "essential_genes.tsv"
    essential_genes = set(
        pd.read_csv(essential_genes_file, sep="\t")["gene_symbol_0"].dropna().unique()
    )

    # Generate panels
    print("\n" + "=" * 60)
    print("Generating panels...")
    print("=" * 60)

    # Panel A: Combined heatmap (cell mapping + cell density, averaged across all wells)
    print("\nPanel A: Combined heatmap (cell mapping + cell density)...")
    plot_panel_a_combined(
        heatmap_data, density_data, output_dir / "qc_panel_a_heatmap.png"
    )

    # Panel B: Cell count boxplot
    print("Panel B: Cell count boxplot...")
    plot_panel_b_boxplot(mapping_overview, output_dir / "qc_panel_b_boxplot.png")

    # Panel C: Mapping percentage bar plot
    print("Panel C: Mapping percentage bar plot...")
    plot_panel_c_barplot(mapping_overview, output_dir / "qc_panel_c_barplot.png")

    # Panel D: Barcode prefix matching (3 variants)
    print("Panel D: Computing barcode prefix matching...")
    prefix_data, prefix_sets = compute_prefix_matching(reads_data, library_barcodes)
    plot_panel_d_prefix_matching(
        prefix_data, prefix_sets, output_dir / "qc_panel_d_prefix_matching.png"
    )

    # D zoomed
    print("Panel D (zoomed): Zoomed prefix matching...")
    plot_panel_d_zoomed(
        prefix_data, prefix_sets, output_dir / "qc_panel_d_zoomed.png"
    )

    # Panel E: Q-scores
    print("Panel E: Q-scores across cycles...")
    plot_panel_e_qscores(reads_data, output_dir / "qc_panel_e_qscores.png")

    # Panel F: Cells per perturbation by category
    print("Panel F: Cells per perturbation histogram...")
    plot_panel_f_histogram(
        perturbation_data, essential_genes, output_dir / "qc_panel_f_kde.png"
    )

    # Panel G: Cell retention funnel (per well)
    print("Panel G: Cell retention funnel...")
    retention_data = load_retention_funnel_data()
    plot_panel_g_retention_funnel(
        retention_data, output_dir / "qc_panel_g_retention_funnel.png"
    )

    # Panel H: Cell retention totals (summed)
    print("Panel H: Cell retention totals...")
    plot_panel_h_retention_totals(
        retention_data, output_dir / "qc_panel_h_retention_totals.png"
    )

    # Panel I: Error-corrected prefix matching
    print("Panel I: Error-corrected prefix matching...")
    corr_prefix_data, corr_prefix_sets = compute_corrected_prefix_matching(
        reads_data, library_barcodes
    )
    plot_panel_i_corrected_prefix(
        corr_prefix_data, corr_prefix_sets, output_dir / "qc_panel_i_corrected_prefix.png"
    )

    # Panel J: Mapping rate & retention vs read quality
    print("Panel J: Mapping rate & retention vs read quality...")
    quality_data = compute_mapping_vs_quality(reads_data, library_barcodes)
    plot_panel_j_mapping_vs_quality(
        quality_data, output_dir / "qc_panel_j_mapping_vs_quality.png"
    )

    print("\n" + "=" * 60)
    print(f"All panels saved to: {output_dir}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
