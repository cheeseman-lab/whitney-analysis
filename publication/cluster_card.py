"""
Cluster Card Generator - Native SVG output for Illustrator editing.

Generates cards as clean SVGs with named groups, native text elements,
and scatter plots rendered as SVG circles. Fully editable in Illustrator.

Usage:
    python cluster_card.py                         # Default clusters
    python cluster_card.py --clusters 0 221 5      # Specific clusters
    python cluster_card.py --all                    # All high-confidence clusters
    python cluster_card.py --output-dir /path/to/dir
"""

import argparse
import textwrap
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# PATHS
# =============================================================================
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
OUTPUT_DIR = Path("/mnt/data/blainey/whitney-analysis/publication/figures/LLM")

PHATE_FILE = (
    BRIEFLOW_OUTPUT
    / "cluster"
    / "Hoescht_COX4_AGP_ConA"
    / "all"
    / "filtered"
    / "10"
    / "phate_leiden_clustering.tsv"
)
EFFECTS_FILE = (
    BRIEFLOW_OUTPUT
    / "aggregate"
    / "tsvs"
    / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__features_genes.tsv"
)
PVALS_FILE = (
    BRIEFLOW_OUTPUT
    / "aggregate"
    / "bootstrap"
    / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__all_gene_bootstrap_results.tsv"
)

DEFAULT_CLUSTERS = [1, 2, 5, 6, 8, 13]

CLUSTER_COLORS = {
    1: "#2563eb",
    2: "#1d4ed8",
    13: "#3b82f6",
    5: "#dc2626",
    6: "#ea580c",
    8: "#9333ea",
}
FALLBACK_COLORS = [
    "#8b5cf6", "#06b6d4", "#f97316", "#22c55e", "#ef4444",
    "#3b82f6", "#ec4899", "#14b8a6", "#f59e0b", "#6366f1",
    "#84cc16", "#a855f7",
]


def get_cluster_color(cluster_id, pathway):
    if cluster_id in CLUSTER_COLORS:
        return CLUSTER_COLORS[cluster_id]
    return FALLBACK_COLORS[hash(pathway) % len(FALLBACK_COLORS)]


# =============================================================================
# LAYOUT CONSTANTS
# =============================================================================
CARD_RX = 8
MARGIN = 25

HEADER_H = 58

# Square plots — side by side with gap
PLOT_SIZE = 235  # square: width = height
PLOT_GAP = 20
PLOT_TOP = 72
PLOT_L_X = MARGIN
PLOT_R_X = MARGIN + PLOT_SIZE + PLOT_GAP

# Card width derived from plot layout
CARD_W = PLOT_R_X + PLOT_SIZE + MARGIN

# Text area starts below plots
TEXT_START_Y = PLOT_TOP + PLOT_SIZE + 15
TEXT_X = MARGIN
TEXT_RIGHT = CARD_W - MARGIN
TEXT_W = TEXT_RIGHT - TEXT_X

# Character width estimates (px per character at given font size)
# Calibrated for Myriad Pro / Arial; monospace is wider
CHAR_W_REGULAR = 0.48  # regular proportional font
CHAR_W_MONO = 0.55     # monospace font

FONT_FAMILY = "'Myriad Pro', Arial, Helvetica, sans-serif"
FONT_MONO = "'Courier New', Courier, monospace"

COLOR_TEXT_PRIMARY = "#1f2937"
COLOR_TEXT_SECONDARY = "#4b5563"
COLOR_TEXT_MUTED = "#9ca3af"
COLOR_BORDER = "#e5e7eb"
COLOR_BG_POINTS = "#b0b5bd"


# =============================================================================
# TEXT WIDTH HELPERS
# =============================================================================


def _chars_per_line(font_size, mono=False):
    """Estimate characters per line for the text area width."""
    char_w = CHAR_W_MONO if mono else CHAR_W_REGULAR
    return int(TEXT_W / (font_size * char_w))


def _wrap_text(text, font_size, mono=False):
    """Wrap text to fit the card text area."""
    return textwrap.fill(text, width=_chars_per_line(font_size, mono=mono))


# =============================================================================
# SVG BUILDING HELPERS
# =============================================================================


def svg_open(card_h):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {CARD_W} {card_h}" '
        f'width="{CARD_W}" height="{card_h}">\n'
        f"<defs>\n"
        f'  <clipPath id="plot-clip-left">\n'
        f'    <rect x="{PLOT_L_X}" y="{PLOT_TOP}" '
        f'width="{PLOT_SIZE}" height="{PLOT_SIZE}"/>\n'
        f"  </clipPath>\n"
        f'  <clipPath id="plot-clip-right">\n'
        f'    <rect x="{PLOT_R_X}" y="{PLOT_TOP}" '
        f'width="{PLOT_SIZE}" height="{PLOT_SIZE}"/>\n'
        f"  </clipPath>\n"
        f"</defs>\n"
        f"<style>\n"
        f"  text {{ font-family: {FONT_FAMILY}; }}\n"
        f"  .mono {{ font-family: {FONT_MONO}; }}\n"
        f"  .section-title {{ font-size: 10px; font-weight: 600; "
        f"fill: {COLOR_TEXT_PRIMARY}; }}\n"
        f"  .gene-name {{ font-size: 9.5px; font-weight: bold; }}\n"
        f"  .gene-badge {{ font-size: 6.5px; fill: {COLOR_TEXT_MUTED}; "
        f"font-style: italic; font-weight: normal; }}\n"
        f"  .gene-desc {{ font-size: 7px; fill: {COLOR_TEXT_SECONDARY}; }}\n"
        f"  .summary-text {{ font-size: 7.5px; fill: {COLOR_TEXT_SECONDARY}; "
        f"font-style: italic; }}\n"
        f"  .est-gene {{ font-size: 8px; fill: {COLOR_TEXT_PRIMARY}; }}\n"
        f"</style>\n"
    )


def svg_card_bg(card_h, accent_color):
    return (
        f'<g id="card-background">\n'
        f'  <rect x="2" y="2" width="{CARD_W - 4}" height="{card_h - 4}" '
        f'rx="{CARD_RX}" fill="white" stroke="{COLOR_BORDER}" stroke-width="1.2"/>\n'
        f'  <rect x="2" y="{CARD_RX}" width="4" '
        f'height="{card_h - 2 * CARD_RX}" fill="{accent_color}" opacity="0.3"/>\n'
        f"</g>\n"
    )


def svg_header_bar(color, pathway, cluster_id):
    lines = []
    lines.append('<g id="header">')
    lines.append(
        f'  <rect x="2" y="2" width="{CARD_W - 4}" height="{HEADER_H}" '
        f'rx="{CARD_RX}" fill="{color}"/>'
    )
    lines.append(
        f'  <rect x="2" y="{HEADER_H - CARD_RX}" width="{CARD_W - 4}" '
        f'height="{CARD_RX + 2}" fill="{color}"/>'
    )
    lines.append(
        f'  <text x="{CARD_W / 2}" y="{24}" text-anchor="middle" '
        f'font-size="15" font-weight="bold" fill="white">'
        f"{escape(pathway)}</text>"
    )
    lines.append(
        f'  <text x="{CARD_W / 2}" y="{44}" text-anchor="middle" '
        f'font-size="10" fill="white" opacity="0.8">'
        f"Cluster {cluster_id}</text>"
    )
    lines.append("</g>")
    return "\n".join(lines) + "\n"


def _map_coords(values, data_min, data_max, svg_min, svg_max):
    data_range = data_max - data_min
    if data_range == 0:
        return np.full_like(values, (svg_min + svg_max) / 2)
    return svg_min + (values - data_min) / data_range * (svg_max - svg_min)


def svg_scatter_plot(
    group_id, title, x_label, y_label,
    bg_x, bg_y, fg_x, fg_y, fg_color,
    plot_x, plot_y, plot_size,
    clip_id=None,
):
    """SVG scatter plot — square, with axes and clipping."""
    lines = []
    lines.append(f'<g id="{group_id}">')

    pad = 14
    title_h = 16
    label_h = 14
    # Inner plot area (square data region)
    ix = plot_x + pad
    iy = plot_y + title_h
    inner = plot_size - 2 * pad - title_h - label_h

    # Title
    lines.append(
        f'  <text x="{plot_x + plot_size / 2}" y="{plot_y + 12}" '
        f'text-anchor="middle" font-size="8.5" fill="{COLOR_TEXT_SECONDARY}" '
        f'font-weight="500">{escape(title)}</text>'
    )

    # Axes — L-shape
    lines.append(
        f'  <line x1="{ix}" y1="{iy + inner}" x2="{ix + inner}" '
        f'y2="{iy + inner}" stroke="{COLOR_TEXT_MUTED}" stroke-width="0.4"/>'
    )
    lines.append(
        f'  <line x1="{ix}" y1="{iy}" x2="{ix}" y2="{iy + inner}" '
        f'stroke="{COLOR_TEXT_MUTED}" stroke-width="0.4"/>'
    )

    # Axis labels
    lines.append(
        f'  <text x="{ix + inner / 2}" y="{iy + inner + 12}" '
        f'text-anchor="middle" font-size="6.5" fill="{COLOR_TEXT_MUTED}">'
        f"{escape(x_label)}</text>"
    )
    lines.append(
        f'  <text x="{ix - 6}" y="{iy + inner / 2}" text-anchor="middle" '
        f'font-size="6.5" fill="{COLOR_TEXT_MUTED}" '
        f'transform="rotate(-90,{ix - 6},{iy + inner / 2})">'
        f"{escape(y_label)}</text>"
    )

    # Data ranges (filter NaN/inf)
    bg_mask = np.isfinite(bg_x) & np.isfinite(bg_y)
    bg_x_clean, bg_y_clean = bg_x[bg_mask], bg_y[bg_mask]

    fg_mask = np.isfinite(fg_x) & np.isfinite(fg_y)
    fg_x_clean, fg_y_clean = fg_x[fg_mask], fg_y[fg_mask]

    all_x = np.concatenate([bg_x_clean, fg_x_clean]) if len(fg_x_clean) > 0 else bg_x_clean
    all_y = np.concatenate([bg_y_clean, fg_y_clean]) if len(fg_y_clean) > 0 else bg_y_clean

    if len(all_x) == 0:
        lines.append("</g>")
        return "\n".join(lines) + "\n"

    x_min, x_max = all_x.min(), all_x.max()
    y_min, y_max = all_y.min(), all_y.max()
    x_pad = 0.05 * max(x_max - x_min, 0.001)
    y_pad = 0.05 * max(y_max - y_min, 0.001)
    x_min -= x_pad
    x_max += x_pad
    y_min -= y_pad
    y_max += y_pad

    clip_attr = f' clip-path="url(#{clip_id})"' if clip_id else ""

    # Background points
    if len(bg_x_clean) > 0:
        sx = _map_coords(bg_x_clean, x_min, x_max, ix, ix + inner)
        sy = _map_coords(bg_y_clean, y_min, y_max, iy + inner, iy)
        lines.append(
            f'  <g id="{group_id}-background" opacity="0.5"{clip_attr}>'
        )
        for cx, cy in zip(sx, sy):
            lines.append(
                f'    <circle cx="{cx:.1f}" cy="{cy:.1f}" r="0.7" '
                f'fill="{COLOR_BG_POINTS}"/>'
            )
        lines.append("  </g>")

    # Foreground points
    if len(fg_x_clean) > 0:
        sx = _map_coords(fg_x_clean, x_min, x_max, ix, ix + inner)
        sy = _map_coords(fg_y_clean, y_min, y_max, iy + inner, iy)
        lines.append(f'  <g id="{group_id}-cluster"{clip_attr}>')
        for cx, cy in zip(sx, sy):
            lines.append(
                f'    <circle cx="{cx:.1f}" cy="{cy:.1f}" r="2.2" '
                f'fill="{fg_color}" opacity="0.85"/>'
            )
        lines.append("  </g>")

    lines.append("</g>")
    return "\n".join(lines) + "\n"


def svg_section_header(title, y, accent_color):
    svg = (
        f'  <line x1="{TEXT_X}" y1="{y}" x2="{TEXT_RIGHT}" y2="{y}" '
        f'stroke="{COLOR_BORDER}" stroke-width="0.5"/>\n'
        f'  <circle cx="{TEXT_X + 4}" cy="{y + 12}" r="3" '
        f'fill="{accent_color}" opacity="0.7"/>\n'
        f'  <text x="{TEXT_X + 12}" y="{y + 15}" class="section-title">'
        f"{escape(title)}</text>\n"
    )
    return svg, y + 24


def svg_wrapped_text(text, x, y, font_size, css_class, mono=False, max_lines=8):
    """Render wrapped text as SVG <text> elements, properly fitted to card width."""
    wrapped = _wrap_text(text, font_size, mono=mono)
    all_lines = wrapped.split("\n")
    out = []
    line_h = font_size + 2.5
    for line in all_lines[:max_lines]:
        out.append(
            f'  <text x="{x}" y="{y}" class="{css_class}">{escape(line)}</text>'
        )
        y += line_h
    if len(all_lines) > max_lines:
        out.append(
            f'  <text x="{x}" y="{y}" font-size="6.5" fill="{COLOR_TEXT_MUTED}" '
            f'font-style="italic">[truncated]</text>'
        )
        y += 10
    return "\n".join(out) + "\n", y


def _estimate_text_height(text, font_size, mono=False, max_lines=8):
    """Estimate vertical space needed for wrapped text."""
    wrapped = _wrap_text(text, font_size, mono=mono)
    n_lines = min(len(wrapped.split("\n")), max_lines)
    return n_lines * (font_size + 2.5)


# =============================================================================
# CARD GENERATOR
# =============================================================================


def create_card(cluster_id, df_table, df_clustering, df_effects, df_pvals, output_dir):
    cluster_df = df_table[df_table["cluster_id"] == cluster_id].copy()
    if len(cluster_df) == 0:
        print(f"  No data for cluster {cluster_id}, skipping")
        return

    pathway = cluster_df["pathway"].iloc[0]
    feature = cluster_df["best_feature"].iloc[0]
    cluster_summary = cluster_df["cluster_summary"].iloc[0]
    color = get_cluster_color(cluster_id, pathway)

    established = cluster_df[cluster_df["gene_type"] == "established"]
    novel = cluster_df[cluster_df["gene_type"] == "novel"].sort_values("priority")
    uncharacterized = cluster_df[
        cluster_df["gene_type"] == "uncharacterized"
    ].sort_values("priority")
    est_genes = established["gene"].tolist()

    # Volcano data
    pval_col = f"{feature}_pval"
    df_volcano = df_effects[["gene_symbol_0", feature]].merge(
        df_pvals[["gene", pval_col]],
        left_on="gene_symbol_0",
        right_on="gene",
        how="inner",
    )
    df_volcano["neg_log10_pval"] = -np.log10(
        df_volcano[pval_col].clip(lower=1e-10)
    )

    cluster_genes = cluster_df["gene"].tolist()
    cluster_phate = df_clustering[df_clustering["gene_symbol_0"].isin(cluster_genes)]

    # Cluster volcano data — recompute neg_log10_pval to avoid table NaN
    cluster_volcano = cluster_df[cluster_df["effect_size"].notna()].copy()
    cluster_volcano["neg_log10_pval"] = cluster_volcano["pval"].clip(
        lower=1e-10
    ).pipe(lambda s: -np.log10(s))

    feature_name = (
        feature.replace("_", " ").replace("cell ", "")
        .replace("mean", "").strip().title()
    )

    discovery_genes = pd.concat([novel, uncharacterized]).sort_values("priority")
    max_show = 12

    # =========================================================================
    # PRE-COMPUTE CARD HEIGHT
    # =========================================================================
    y_est = TEXT_START_Y

    # Summary
    y_est += 24
    if pd.notna(cluster_summary):
        summary_text = str(cluster_summary)[:600]
        y_est += _estimate_text_height(summary_text, 7.5, max_lines=6)
    y_est += 12

    # Established genes
    y_est += 24
    if est_genes:
        y_est += _estimate_text_height(", ".join(est_genes), 8, mono=True, max_lines=8)
    else:
        y_est += 12
    y_est += 12

    # Discovery genes
    y_est += 24
    for idx, (_, row) in enumerate(discovery_genes.iterrows()):
        if idx >= max_show:
            y_est += 12
            break
        y_est += 13
        desc = row["description"] if pd.notna(row["description"]) else ""
        if desc:
            desc_short = desc[:200] + "..." if len(desc) > 200 else desc
            y_est += _estimate_text_height(desc_short, 7, max_lines=2)
        y_est += 3

    y_est += 35
    card_h = max(y_est, 500)

    # =========================================================================
    # BUILD SVG
    # =========================================================================
    parts = []
    parts.append(svg_open(card_h))
    parts.append(svg_card_bg(card_h, color))
    parts.append(svg_header_bar(color, pathway, cluster_id))

    # PHATE plot (square)
    parts.append(
        svg_scatter_plot(
            group_id="phate-plot",
            title="PHATE embedding",
            x_label="PHATE 1",
            y_label="PHATE 2",
            bg_x=df_clustering["PHATE_0"].values,
            bg_y=df_clustering["PHATE_1"].values,
            fg_x=cluster_phate["PHATE_0"].values,
            fg_y=cluster_phate["PHATE_1"].values,
            fg_color=color,
            plot_x=PLOT_L_X,
            plot_y=PLOT_TOP,
            plot_size=PLOT_SIZE,
            clip_id="plot-clip-left",
        )
    )

    # Volcano plot (square)
    parts.append(
        svg_scatter_plot(
            group_id="volcano-plot",
            title=feature_name,
            x_label="Effect size",
            y_label="-log\u2081\u2080(p)",
            bg_x=df_volcano[feature].values,
            bg_y=df_volcano["neg_log10_pval"].values,
            fg_x=cluster_volcano["effect_size"].values,
            fg_y=cluster_volcano["neg_log10_pval"].values,
            fg_color=color,
            plot_x=PLOT_R_X,
            plot_y=PLOT_TOP,
            plot_size=PLOT_SIZE,
            clip_id="plot-clip-right",
        )
    )

    # =========================================================================
    # TEXT SECTIONS
    # =========================================================================
    text_parts = []
    text_parts.append('<g id="text-content">')
    y = TEXT_START_Y

    # --- Summary ---
    hdr, y = svg_section_header("Summary", y, color)
    text_parts.append(hdr)
    if pd.notna(cluster_summary):
        summary_text = str(cluster_summary)
        if len(summary_text) > 600:
            summary_text = summary_text[:600] + "..."
        block, y = svg_wrapped_text(
            summary_text, TEXT_X, y, font_size=7.5,
            css_class="summary-text", max_lines=6,
        )
        text_parts.append(block)
    y += 12

    # --- Established genes ---
    hdr, y = svg_section_header("Established genes", y, color)
    text_parts.append(hdr)
    if est_genes:
        gene_str = ", ".join(est_genes)
        block, y = svg_wrapped_text(
            gene_str, TEXT_X, y, font_size=8,
            css_class="mono est-gene", mono=True, max_lines=8,
        )
        text_parts.append(block)
    else:
        text_parts.append(
            f'  <text x="{TEXT_X}" y="{y}" font-size="8" '
            f'fill="{COLOR_TEXT_MUTED}">None</text>\n'
        )
        y += 12
    y += 12

    # --- Novel & uncharacterized genes ---
    hdr, y = svg_section_header("Novel & uncharacterized genes", y, color)
    text_parts.append(hdr)

    for idx, (_, row) in enumerate(discovery_genes.iterrows()):
        if idx >= max_show:
            remaining = len(discovery_genes) - idx
            if remaining > 0:
                text_parts.append(
                    f'  <text x="{TEXT_X}" y="{y}" font-size="7" '
                    f'fill="{COLOR_TEXT_MUTED}" font-style="italic">'
                    f"+ {remaining} more</text>\n"
                )
                y += 12
            break

        gene = row["gene"]
        gene_type = row["gene_type"]
        desc = row["description"] if pd.notna(row["description"]) else ""
        badge = "novel" if gene_type == "novel" else "unc."

        text_parts.append(
            f'  <text x="{TEXT_X}" y="{y}" class="mono gene-name" '
            f'fill="{color}">{escape(gene)}'
            f'<tspan class="gene-badge" dx="5">{badge}</tspan>'
            f"</text>\n"
        )
        y += 13

        if desc:
            desc_short = desc[:200] + "..." if len(desc) > 200 else desc
            desc_indent = TEXT_X + 12
            # Use the indented width for wrapping
            indent_w = TEXT_RIGHT - desc_indent
            chars = int(indent_w / (7 * CHAR_W_REGULAR))
            wrapped = textwrap.fill(desc_short, width=chars)
            for line in wrapped.split("\n")[:2]:
                text_parts.append(
                    f'  <text x="{desc_indent}" y="{y}" class="gene-desc">'
                    f"{escape(line)}</text>\n"
                )
                y += 9.5
        y += 3

    # --- Footer ---
    y += 8
    n_total = len(cluster_df)
    footer = (
        f"{n_total} genes: {len(est_genes)} established  \u00b7  "
        f"{len(novel)} novel  \u00b7  {len(uncharacterized)} uncharacterized"
    )
    text_parts.append(
        f'  <text x="{CARD_W / 2}" y="{card_h - 12}" text-anchor="middle" '
        f'font-size="7.5" fill="{COLOR_TEXT_MUTED}">{footer}</text>\n'
    )

    text_parts.append("</g>")
    parts.append("\n".join(text_parts) + "\n")
    parts.append("</svg>\n")

    # Write
    output_dir.mkdir(parents=True, exist_ok=True)
    svg_content = "".join(parts)
    output_path = output_dir / f"card_cluster_{cluster_id}.svg"
    output_path.write_text(svg_content)
    print(f"  Saved: {output_path}")

    # Content TSV
    tsv_path = output_dir / f"card_cluster_{cluster_id}_content.tsv"
    content_rows = []
    for _, row in cluster_df.iterrows():
        content_rows.append({
            "gene": row["gene"],
            "gene_type": row["gene_type"],
            "priority": row["priority"],
            "description": row["description"],
            "effect_size": row["effect_size"],
            "neg_log10_pval": row["neg_log10_pval"],
        })
    pd.DataFrame(content_rows).to_csv(tsv_path, sep="\t", index=False)
    print(f"  Saved: {tsv_path}")


# =============================================================================
# MAIN
# =============================================================================


def main():
    parser = argparse.ArgumentParser(description="Generate publication cluster cards.")
    parser.add_argument(
        "--clusters", type=int, nargs="+", default=None,
        help="Cluster IDs to generate",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Generate cards for all high-confidence clusters",
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help=f"Output directory (default: {OUTPUT_DIR})",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    table_file = output_dir / "cluster_genes_table.tsv"

    print("Loading data...")
    df_table = pd.read_csv(table_file, sep="\t")
    df_clustering = pd.read_csv(PHATE_FILE, sep="\t")
    df_effects = pd.read_csv(EFFECTS_FILE, sep="\t")
    df_pvals = pd.read_csv(PVALS_FILE, sep="\t")

    if args.all:
        cluster_ids = sorted(df_table["cluster_id"].unique())
    elif args.clusters:
        cluster_ids = args.clusters
    else:
        cluster_ids = DEFAULT_CLUSTERS

    print(f"\nGenerating cards for {len(cluster_ids)} cluster(s): {cluster_ids}\n")

    for cluster_id in cluster_ids:
        print(f"Cluster {cluster_id}:")
        create_card(
            cluster_id, df_table, df_clustering, df_effects, df_pvals, output_dir,
        )

    print(f"\nDone! {len(cluster_ids)} card(s) generated in {output_dir}")


if __name__ == "__main__":
    main()
