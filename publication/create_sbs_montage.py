"""
Create a 6-panel montage of SBS images across cycles for a high-mapping-rate tile.
Tile 58 from P-1_W-A1 was selected (83.5% mapping rate, 9650 cells).
"""

import tifffile
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
SBS_IMAGE_DIR = Path("analysis/brieflow_output/preprocess/images/sbs")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

# Selected tile with high mapping rate
PLATE = 1
WELL = "A1"
TILE = 58
CYCLES = [1, 3, 5, 7, 9, 11]  # Every other cycle for 6 panels

# Channel names for SBS (typical)
CHANNEL_NAMES = ["DAPI", "G", "T", "A", "C"]


def load_sbs_image(plate, well, tile, cycle):
    """Load a single SBS image."""
    filename = f"P-{plate}_W-{well}_T-{tile}_C-{cycle}__image.tiff"
    filepath = SBS_IMAGE_DIR / filename
    return tifffile.imread(filepath)


def normalize_image(img, percentile_low=1, percentile_high=99.5):
    """Normalize image to 0-1 range using percentile clipping."""
    vmin = np.percentile(img, percentile_low)
    vmax = np.percentile(img, percentile_high)
    img_norm = (img.astype(float) - vmin) / (vmax - vmin + 1e-8)
    return np.clip(img_norm, 0, 1)


def create_composite(img, channel_colors=None):
    """
    Create RGB composite from multi-channel image.
    Default: DAPI=blue, G=green, T=red, A=magenta, C=cyan
    """
    if channel_colors is None:
        # DAPI, G, T, A, C
        channel_colors = [
            [0, 0, 1],      # DAPI - blue
            [0, 1, 0],      # G - green
            [1, 0, 0],      # T - red
            [1, 0, 1],      # A - magenta
            [0, 1, 1],      # C - cyan
        ]

    rgb = np.zeros((*img.shape[1:], 3), dtype=float)
    for ch_idx, color in enumerate(channel_colors):
        ch_norm = normalize_image(img[ch_idx])
        for c_idx, c_val in enumerate(color):
            rgb[..., c_idx] += ch_norm * c_val

    return np.clip(rgb, 0, 1)


def create_sbs_composite(img):
    """
    Create RGB composite showing just the 4 SBS channels (not DAPI).
    G=green, T=red, A=magenta, C=cyan
    """
    # Only use channels 1-4 (SBS bases)
    channel_colors = [
        [0, 1, 0],      # G - green
        [1, 0, 0],      # T - red
        [1, 0, 1],      # A - magenta
        [0, 1, 1],      # C - cyan
    ]

    rgb = np.zeros((*img.shape[1:], 3), dtype=float)
    for ch_idx, color in enumerate(channel_colors):
        ch_norm = normalize_image(img[ch_idx + 1])  # Skip DAPI (channel 0)
        for c_idx, c_val in enumerate(color):
            rgb[..., c_idx] += ch_norm * c_val

    return np.clip(rgb, 0, 1)


def main():
    print(f"Loading SBS images for P-{PLATE}_W-{WELL}_T-{TILE}")
    print(f"Cycles: {CYCLES}")

    # Load images
    images = {}
    for cycle in CYCLES:
        print(f"  Loading cycle {cycle}...")
        images[cycle] = load_sbs_image(PLATE, WELL, TILE, cycle)

    # Create 6-panel figure showing SBS composite for each cycle
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, cycle in enumerate(CYCLES):
        img = images[cycle]
        composite = create_sbs_composite(img)

        # Show a cropped region (center 1000x1000) for better visibility
        h, w = composite.shape[:2]
        crop_size = 1000
        y_start = (h - crop_size) // 2
        x_start = (w - crop_size) // 2
        cropped = composite[y_start:y_start+crop_size, x_start:x_start+crop_size]

        axes[idx].imshow(cropped)
        axes[idx].set_title(f"Cycle {cycle}", fontsize=14, fontweight='bold')
        axes[idx].axis('off')

    plt.suptitle(f"SBS Images Across Cycles\nP-{PLATE}_W-{WELL}_T-{TILE} (83.5% mapping rate)",
                 fontsize=16, fontweight='bold')

    # Add legend for channels
    legend_text = "Channels: G(green) T(red) A(magenta) C(cyan)"
    fig.text(0.5, 0.02, legend_text, ha='center', fontsize=12, style='italic')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    output_path = OUTPUT_DIR / "sbs_cycles_montage.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nSaved montage to: {output_path}")

    # Also create a full-resolution version
    fig2, axes2 = plt.subplots(2, 3, figsize=(24, 16))
    axes2 = axes2.flatten()

    for idx, cycle in enumerate(CYCLES):
        img = images[cycle]
        composite = create_sbs_composite(img)

        axes2[idx].imshow(composite)
        axes2[idx].set_title(f"Cycle {cycle}", fontsize=14, fontweight='bold')
        axes2[idx].axis('off')

    plt.suptitle(f"SBS Images Across Cycles (Full FOV)\nP-{PLATE}_W-{WELL}_T-{TILE} (83.5% mapping rate)",
                 fontsize=16, fontweight='bold')
    fig2.text(0.5, 0.02, legend_text, ha='center', fontsize=12, style='italic')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    output_path_full = OUTPUT_DIR / "sbs_cycles_montage_full.png"
    plt.savefig(output_path_full, dpi=100, bbox_inches='tight', facecolor='white')
    print(f"Saved full FOV montage to: {output_path_full}")

    plt.close('all')


if __name__ == "__main__":
    main()
