"""
Generate matched SBS/Phenotype raw images for FIJI visualization,
plus fast_merge_example alignment plots for each tile-site pair.

Uses validated initial_sites from config.yml (known-good ph_tile/sbs_site pairs),
then picks the pair with the highest SBS mapping rate per well.

SBS output: (12, 5, 4168, 4168) uint16 with TCYX axes
Phenotype output: (4, 2084, 2084) uint16 with CYX axes
"""

import yaml
import pandas as pd
import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from lib.merge.merge_utils import fast_merge_example
from lib.merge.hash import hash_cell_locations, initial_alignment

# Configuration
CONFIG_PATH = Path("analysis/config/config.yml")
BASE_DIR = Path("analysis/brieflow_output")
SBS_IMAGE_DIR = BASE_DIR / "preprocess/images/sbs"
PHENOTYPE_IMAGE_DIR = BASE_DIR / "preprocess/images/phenotype"
MAPPING_DIR = BASE_DIR / "sbs/eval/mapping"
OUTPUT_DIR = Path("publication/matched_images")

PLATES = [1, 2]
NUM_CYCLES = 12
# Only process the best wells (highest mean mapping rate per plate)
SELECTED_WELLS = {1: ["A1", "A2"], 2: ["A1", "A2"]}

THRESHOLD = 10  # Merge matching threshold (from notebook)


def load_initial_sites():
    """Load validated initial_sites from config.yml."""
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    return config["merge"]["initial_sites"]


def load_mapping_rates(plate):
    """Load mapping rate TSV for a plate."""
    path = MAPPING_DIR / f"P-{plate}__cell_mapping_heatmap_one.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"fraction of cells mapping to 1 gene symbols": "mapping_rate"})
    df["plate"] = plate
    return df


def get_best_initial_site(initial_sites, plate, well, mapping_df):
    """From validated initial_sites, pick the pair with highest SBS mapping rate.

    initial_sites format: {plate: {well: [[ph_tile, sbs_site], ...]}}
    Returns (ph_tile, sbs_site, mapping_rate) or None.
    """
    pairs = initial_sites.get(plate, {}).get(well, [])
    if not pairs:
        return None

    sbs_tiles = [p[1] for p in pairs]
    well_mapping = mapping_df[
        (mapping_df["well"] == well) & (mapping_df["tile"].isin(sbs_tiles))
    ]
    if well_mapping.empty:
        return None

    best = well_mapping.loc[well_mapping["mapping_rate"].idxmax()]
    best_sbs = int(best["tile"])
    ph_tile = [p[0] for p in pairs if p[1] == best_sbs][0]
    return ph_tile, best_sbs, best["mapping_rate"]


def load_cell_info(plate, well):
    """Load phenotype_info and sbs_info for a plate-well."""
    ph_path = BASE_DIR / "phenotype/parquets" / f"P-{plate}_W-{well}__phenotype_info.parquet"
    sbs_path = BASE_DIR / "sbs/parquets" / f"P-{plate}_W-{well}__sbs_info.parquet"
    phenotype_info = pd.read_parquet(ph_path)
    sbs_info = pd.read_parquet(sbs_path)
    return phenotype_info, sbs_info


def save_merge_plot(ph_tile, sbs_site, alignment_df, phenotype_info, sbs_info,
                    threshold, output_path):
    """Run fast_merge_example and save the plot to a file."""
    print(f"  Generating merge example plot (ph_tile={ph_tile}, sbs_site={sbs_site})...")

    original_show = plt.show

    def save_show():
        fig = plt.gcf()
        fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
        fig.savefig(output_path.with_suffix(".pdf"), dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    plt.show = save_show
    try:
        success = fast_merge_example(
            ph_tile, sbs_site, alignment_df,
            phenotype_info, sbs_info, threshold
        )
    finally:
        plt.show = original_show

    if success:
        print(f"  Saved merge plot: {output_path}")
    return success


def save_sbs_fiji(plate, well, tile, output_path):
    """Load all 12 cycles for an SBS tile and save as TCYX FIJI TIFF."""
    cycle_images = []
    for cycle in range(1, NUM_CYCLES + 1):
        filename = f"P-{plate}_W-{well}_T-{tile}_C-{cycle}__image.tiff"
        filepath = SBS_IMAGE_DIR / filename
        if not filepath.exists():
            print(f"  WARNING: Missing {filepath}")
            return False
        img = tifffile.imread(filepath)  # (5, 4168, 4168)
        cycle_images.append(img)

    stack = np.stack(cycle_images, axis=0)  # (12, 5, Y, X) = TCYX
    print(f"  SBS stack shape: {stack.shape} dtype: {stack.dtype}")

    tifffile.imwrite(
        output_path, stack, imagej=True,
        metadata={"axes": "TCYX", "ImageJ": "1.53c"},
    )
    return True


def save_phenotype_fiji(plate, well, site, output_path):
    """Load phenotype image and save with CYX FIJI metadata."""
    filename = f"P-{plate}_W-{well}_T-{site}__image.tiff"
    filepath = PHENOTYPE_IMAGE_DIR / filename
    if not filepath.exists():
        print(f"  WARNING: Missing {filepath}")
        return False

    img = tifffile.imread(filepath)  # (4, 2084, 2084)
    print(f"  Phenotype shape: {img.shape} dtype: {img.dtype}")

    tifffile.imwrite(
        output_path, img, imagej=True,
        metadata={"axes": "CYX", "ImageJ": "1.53c"},
    )
    return True


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    initial_sites = load_initial_sites()

    # Load mapping rates for all plates
    all_mapping = pd.concat([load_mapping_rates(p) for p in PLATES], ignore_index=True)
    print(f"Loaded mapping rates: {len(all_mapping)} tile entries across {len(PLATES)} plates")

    summary_rows = []

    for plate in PLATES:
        plate_mapping = all_mapping[all_mapping["plate"] == plate]
        for well in SELECTED_WELLS.get(plate, []):
            result = get_best_initial_site(initial_sites, plate, well, plate_mapping)
            if result is None:
                print(f"\nP-{plate}_W-{well}: No valid initial_sites found, skipping")
                continue

            ph_tile, sbs_tile, mapping_rate = result

            well_dir = OUTPUT_DIR / f"P-{plate}_W-{well}"
            well_dir.mkdir(parents=True, exist_ok=True)

            print(f"\nProcessing P-{plate}_W-{well}: SBS tile={sbs_tile}, "
                  f"PH tile={ph_tile}, mapping_rate={mapping_rate:.3f}")

            # Load cell info and compute alignment for merge plot
            phenotype_info, sbs_info = load_cell_info(plate, well)
            ph_info_hash = hash_cell_locations(phenotype_info)
            sbs_info_hash = hash_cell_locations(sbs_info).rename(columns={"tile": "site"})
            alignment_df = initial_alignment(
                ph_info_hash, sbs_info_hash, initial_sites=[[ph_tile, sbs_tile]]
            )

            # Generate merge example plot
            merge_plot_out = well_dir / f"P-{plate}_W-{well}_T-{sbs_tile}_ph-{ph_tile}_merge_example.png"
            merge_ok = save_merge_plot(
                ph_tile, sbs_tile, alignment_df,
                phenotype_info, sbs_info, THRESHOLD, merge_plot_out
            )

            # Save SBS FIJI image
            sbs_out = well_dir / f"P-{plate}_W-{well}_T-{sbs_tile}_sbs_fiji.tiff"
            sbs_ok = save_sbs_fiji(plate, well, sbs_tile, sbs_out)
            if sbs_ok:
                print(f"  Saved SBS: {sbs_out}")

            # Save matched phenotype FIJI image
            pheno_out = well_dir / f"P-{plate}_W-{well}_ph-{ph_tile}_phenotype_fiji.tiff"
            pheno_ok = save_phenotype_fiji(plate, well, ph_tile, pheno_out)
            if pheno_ok:
                print(f"  Saved phenotype: {pheno_out}")

            summary_rows.append({
                "plate": plate,
                "well": well,
                "sbs_tile": sbs_tile,
                "ph_tile": ph_tile,
                "mapping_rate": mapping_rate,
                "sbs_saved": sbs_ok,
                "phenotype_saved": pheno_ok,
                "merge_plot_saved": merge_ok,
            })

    # Print summary table
    summary = pd.DataFrame(summary_rows)
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(summary.to_string(index=False))
    print(f"\nTotal SBS images saved: {summary['sbs_saved'].sum()}")
    print(f"Total phenotype images saved: {summary['phenotype_saved'].sum()}")
    print(f"Total merge plots saved: {summary['merge_plot_saved'].sum()}")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
