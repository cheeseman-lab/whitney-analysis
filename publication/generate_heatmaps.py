#!/usr/bin/env python3
"""Generate cell mapping heatmap PNGs from TSVs using squid SBS layout."""

import string
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_plate_heatmap(
    df, metric=None, shape="square", plate="6W", snake_sites=True, **kwargs
):
    """Plot the heatmap of a summary DataFrame by well and tile in a convenient plate layout.

    Args:
        df (pandas.DataFrame):
            Summary DataFrame of values to plot, expects one row for each (well, tile) combination.
        metric (str, optional):
            Column of `df` to use for plotting the heatmap. If None, attempts to infer based on column names.
        shape (str or list, optional):
            Shape of subplot for each well. Options include 'squid_sbs' for SQUID microscope SBS tiles.
        plate (str):
            Plate type. Options are {'6W', '24W', '96W'}.
        snake_sites (bool, optional):
            If true, plots tiles in a snake order. Defaults to True.
        **kwargs:
            Keyword arguments passed to matplotlib.pyplot.imshow().

    Returns:
        tuple: (fig, cbar) matplotlib figure and colorbar objects.
    """
    tiles = df["tile"].astype(int)
    tiles = max(len(tiles.unique()), tiles.max())

    # Define grid for plotting
    if shape == "square":
        r = c = int(np.ceil(np.sqrt(tiles)))
        grid = np.empty(r * c)
        grid[:] = np.nan
        grid[:tiles] = range(tiles)
        grid = grid.reshape(r, c)
    else:
        if shape == "squid_sbs":
            # Spatially accurate layout for Squid microscope SBS tiles (94 tiles)
            rows = [2, 6, 8, 10, 10, 11, 11, 10, 10, 8, 6, 2]
        elif shape == "squid_ph":
            # Spatially accurate layout for Squid microscope phenotype tiles (1732 tiles)
            rows = [
                4, 14, 18, 22, 26, 28, 32, 34, 36, 36, 38, 40, 40, 42, 42, 44, 44,
                46, 46, 46, 46, 46, 48, 48, 48, 48, 46, 46, 46, 46, 46, 44, 44,
                42, 42, 40, 40, 38, 36, 36, 34, 32, 28, 26, 22, 18, 14, 4,
            ]
        elif isinstance(shape, list):
            rows = shape
        else:
            raise ValueError(f"{shape} shape not implemented")

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

    # Infer metric to plot if necessary
    if not metric:
        metric = [col for col in df.columns if col not in ["plate", "well", "tile"]]
        if len(metric) != 1:
            raise ValueError(
                "Cannot infer metric to plot, can pass metric column name explicitly to metric kwarg"
            )
        metric = metric[0]

    # Define subplots layout
    if df["well"].nunique() == 1:
        wells = df["well"].unique()
        fig, axes = plt.subplots(1, 1, figsize=(10, 10))
        axes = np.array([axes])
    elif plate == "6W":
        wells = [f"{r}{c}" for r in string.ascii_uppercase[:2] for c in range(1, 4)]
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    elif plate == "24W":
        wells = [f"{r}{c}" for r in string.ascii_uppercase[:4] for c in range(1, 7)]
        fig, axes = plt.subplots(4, 6, figsize=(15, 10))
    elif plate == "96W":
        wells = [f"{r}{c}" for r in string.ascii_uppercase[:8] for c in range(1, 13)]
        fig, axes = plt.subplots(8, 12, figsize=(15, 10))
    else:
        wells = sorted(df["well"].unique())
        nr = nc = int(np.ceil(np.sqrt(len(wells))))
        if (nr - 1) * nc >= len(wells):
            nr -= 1
        fig, axes = plt.subplots(nr, nc, figsize=(15, 15))

    # Define colorbar min and max
    cmin, cmax = df[metric].min(), df[metric].max()
    if 0 <= cmin and cmax <= 1:
        cmin, cmax = 0, 1

    # Plot wells
    for ax, well in zip(axes.reshape(-1), wells):
        values = grid.copy()
        df_well = df.query("well==@well")
        if df_well.pipe(len) > 0:
            for tile in range(tiles):
                try:
                    values[grid == tile] = df_well.loc[
                        df_well.tile == tile, metric
                    ].values[0]
                except:
                    values[grid == tile] = np.nan
            plot = ax.imshow(values, vmin=cmin, vmax=cmax, **kwargs)
        ax.set_title(f"Well {well}", fontsize=24)
        ax.axis("off")

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.025, 0.7])
    try:
        cbar = fig.colorbar(plot, cax=cbar_ax)
    except:
        raise ValueError("No data to plot")
    cbar.set_label(metric, fontsize=18)
    cbar_ax.yaxis.set_ticks_position("left")

    return fig, cbar


def generate_heatmap(tsv_path: Path, output_path: Path, shape: str = "squid_sbs", plate: str = "6W"):
    """Generate a heatmap PNG from a TSV file.

    Args:
        tsv_path: Path to the TSV file with well, tile, and metric columns.
        output_path: Path for the output PNG file.
        shape: Layout shape (default: squid_sbs).
        plate: Plate type (default: 6W).
    """
    df = pd.read_csv(tsv_path, sep="\t")
    fig, cbar = plot_plate_heatmap(df, shape=shape, plate=plate, cmap="viridis")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Generated: {output_path}")


def main():
    mapping_dir = Path("analysis/brieflow_output/sbs/eval/mapping")

    # Define TSV to PNG mappings
    heatmaps = [
        ("P-1__cell_mapping_heatmap_one.tsv", "P-1__cell_mapping_heatmap_one.png"),
        ("P-1__cell_mapping_heatmap_any.tsv", "P-1__cell_mapping_heatmap_any.png"),
        ("P-2__cell_mapping_heatmap_one.tsv", "P-2__cell_mapping_heatmap_one.png"),
        ("P-2__cell_mapping_heatmap_any.tsv", "P-2__cell_mapping_heatmap_any.png"),
    ]

    for tsv_name, png_name in heatmaps:
        tsv_path = mapping_dir / tsv_name
        output_path = mapping_dir / png_name

        if tsv_path.exists():
            generate_heatmap(tsv_path, output_path, shape="squid_sbs", plate="6W")
        else:
            print(f"TSV not found: {tsv_path}")


if __name__ == "__main__":
    main()
