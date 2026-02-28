"""
Cluster Card Generator - Publication-quality cards from consolidated table.

Usage:
    python cluster_card.py                         # Default: clusters 0, 221
    python cluster_card.py --clusters 0 221 5      # Specific clusters
    python cluster_card.py --all                    # All high-confidence clusters
    python cluster_card.py --output-dir /path/to/dir
"""
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import pandas as pd
import numpy as np
import textwrap
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
OUTPUT_DIR = Path("/mnt/data/blainey/whitney-analysis/publication/figures/LLM")

PHATE_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "phate_leiden_clustering.tsv"
EFFECTS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "tsvs" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__features_genes.tsv"
PVALS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "bootstrap" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__all_gene_bootstrap_results.tsv"

# Selected clusters for the paper (Fig 2g, 2h)
DEFAULT_CLUSTERS = [1, 2, 5, 6, 8, 13]

# Consistent cluster-to-color mapping for paper figures
CLUSTER_COLORS = {
    1: '#2563eb',    # blue — 60S ribosome biogenesis
    2: '#1d4ed8',    # dark blue — 40S ribosome biogenesis
    13: '#3b82f6',   # light blue — SSU processome
    5: '#dc2626',    # red — pre-mRNA splicing
    6: '#ea580c',    # orange — ubiquitin-proteasome
    8: '#9333ea',    # purple — RNA Pol II transcription
}
FALLBACK_COLORS = [
    '#8b5cf6', '#06b6d4', '#f97316', '#22c55e', '#ef4444', '#3b82f6',
    '#ec4899', '#14b8a6', '#f59e0b', '#6366f1', '#84cc16', '#a855f7',
]


def get_cluster_color(cluster_id, pathway):
    """Get color for a cluster — fixed for paper clusters, hash-based otherwise."""
    if cluster_id in CLUSTER_COLORS:
        return CLUSTER_COLORS[cluster_id]
    return FALLBACK_COLORS[hash(pathway) % len(FALLBACK_COLORS)]


# =============================================================================
# STYLE CONFIG
# =============================================================================
FONT_TITLE = 11
FONT_SECTION = 8
FONT_BODY = 6.5
FONT_GENE = 7
FONT_AXIS_LABEL = 6
FONT_AXIS_TICK = 5

COLOR_TEXT_PRIMARY = '#1f2937'
COLOR_TEXT_SECONDARY = '#4b5563'
COLOR_TEXT_MUTED = '#9ca3af'
COLOR_BORDER = '#e5e7eb'
COLOR_BG_POINTS = '#e5e7eb'


# =============================================================================
# CARD GENERATOR
# =============================================================================

def create_card(cluster_id, df_table, df_clustering, df_effects, df_pvals, output_dir):
    """Generate a publication-quality cluster card."""

    cluster_df = df_table[df_table['cluster_id'] == cluster_id].copy()
    if len(cluster_df) == 0:
        print(f"  No data for cluster {cluster_id}, skipping")
        return

    pathway = cluster_df['pathway'].iloc[0]
    feature = cluster_df['best_feature'].iloc[0]
    cluster_summary = cluster_df['cluster_summary'].iloc[0]
    color = get_cluster_color(cluster_id, pathway)

    # Split by gene type
    established = cluster_df[cluster_df['gene_type'] == 'established']
    novel = cluster_df[cluster_df['gene_type'] == 'novel'].sort_values('priority')
    uncharacterized = cluster_df[cluster_df['gene_type'] == 'uncharacterized'].sort_values('priority')

    # Build volcano data
    pval_col = f'{feature}_pval'
    df_volcano = df_effects[['gene_symbol_0', feature]].merge(
        df_pvals[['gene', pval_col]],
        left_on='gene_symbol_0', right_on='gene', how='inner'
    )
    df_volcano['neg_log10_pval'] = -np.log10(df_volcano[pval_col].clip(lower=1e-10))

    # =========================================================================
    # FIGURE LAYOUT
    # =========================================================================
    fig = plt.figure(figsize=(5.5, 7))
    fig.patch.set_facecolor('white')

    # Main coordinate axes (invisible, for text placement)
    ax_main = fig.add_axes([0, 0, 1, 1])
    ax_main.set_xlim(0, 10)
    ax_main.set_ylim(0, 14)
    ax_main.axis('off')

    # Card border
    border = mpatches.FancyBboxPatch(
        (0.15, 0.15), 9.7, 13.7,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor='white', edgecolor=COLOR_BORDER, linewidth=1.2
    )
    ax_main.add_patch(border)

    # Header bar with pathway color
    header = mpatches.FancyBboxPatch(
        (0.15, 13.15), 9.7, 0.7,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor=color, edgecolor='none'
    )
    ax_main.add_patch(header)
    # Cover bottom corners of header (so only top is rounded)
    ax_main.add_patch(mpatches.Rectangle((0.15, 13.15), 9.7, 0.15,
                                          facecolor=color, edgecolor='none'))

    # Title text
    title_text = pathway
    cluster_label = f"Cluster {cluster_id}"
    ax_main.text(5, 13.62, title_text, fontsize=FONT_TITLE, fontweight='bold',
                 color='white', va='center', ha='center',
                 path_effects=[pe.withStroke(linewidth=0.5, foreground='black', alpha=0.15)])
    ax_main.text(5, 13.28, cluster_label, fontsize=8, color='white',
                 va='center', ha='center', alpha=0.85)

    # =========================================================================
    # PHATE PLOT
    # =========================================================================
    ax_phate = fig.add_axes([0.08, 0.76, 0.40, 0.16])

    # Background
    ax_phate.scatter(df_clustering['PHATE_0'], df_clustering['PHATE_1'],
                     c=COLOR_BG_POINTS, s=0.3, alpha=0.4, rasterized=True, linewidths=0)

    # Cluster genes highlighted
    cluster_genes = cluster_df['gene'].tolist()
    cluster_phate = df_clustering[df_clustering['gene_symbol_0'].isin(cluster_genes)]
    ax_phate.scatter(cluster_phate['PHATE_0'], cluster_phate['PHATE_1'],
                     c=color, s=6, alpha=0.9, linewidths=0, zorder=3)

    # Full embedding extent
    x_min, x_max = df_clustering['PHATE_0'].min(), df_clustering['PHATE_0'].max()
    y_min, y_max = df_clustering['PHATE_1'].min(), df_clustering['PHATE_1'].max()
    pad = 0.05 * max(x_max - x_min, y_max - y_min)
    ax_phate.set_xlim(x_min - pad, x_max + pad)
    ax_phate.set_ylim(y_min - pad, y_max + pad)

    _style_subplot(ax_phate, 'PHATE embedding')

    # =========================================================================
    # VOLCANO PLOT
    # =========================================================================
    ax_volcano = fig.add_axes([0.56, 0.76, 0.40, 0.16])

    # Background
    ax_volcano.scatter(df_volcano[feature], df_volcano['neg_log10_pval'],
                       c=COLOR_BG_POINTS, s=0.3, alpha=0.25, rasterized=True, linewidths=0)

    # Cluster genes
    cluster_volcano = cluster_df[cluster_df['effect_size'].notna()]
    ax_volcano.scatter(cluster_volcano['effect_size'], cluster_volcano['neg_log10_pval'],
                       c=color, s=8, alpha=0.9, linewidths=0, zorder=3)

    feature_name = feature.replace('_', ' ').replace('cell ', '').replace('mean', '').strip().title()
    ax_volcano.set_xlabel('Effect size', fontsize=FONT_AXIS_LABEL, color=COLOR_TEXT_SECONDARY, labelpad=2)
    ax_volcano.set_ylabel('$-\\log_{10}(p)$', fontsize=FONT_AXIS_LABEL, color=COLOR_TEXT_SECONDARY, labelpad=2)
    _style_subplot(ax_volcano, feature_name)

    # =========================================================================
    # TEXT SECTIONS
    # =========================================================================
    y = 10.2

    # --- Summary ---
    y = _draw_section_header(ax_main, 'Summary', y, color)
    if pd.notna(cluster_summary):
        summary_text = cluster_summary[:400] + '...' if len(str(cluster_summary)) > 400 else str(cluster_summary)
        wrapped = textwrap.fill(summary_text, width=100)
        for line in wrapped.split('\n')[:5]:
            ax_main.text(0.6, y, line, fontsize=5.5, color=COLOR_TEXT_SECONDARY, style='italic')
            y -= 0.22
    y -= 0.25

    # --- Established genes ---
    y = _draw_section_header(ax_main, 'Established genes', y, color)
    est_genes = established['gene'].tolist()
    if est_genes:
        gene_str = ', '.join(est_genes)
        wrapped = textwrap.fill(gene_str, width=85)
        for line in wrapped.split('\n')[:4]:
            ax_main.text(0.6, y, line, fontsize=FONT_BODY, color=COLOR_TEXT_PRIMARY, family='monospace')
            y -= 0.24
    else:
        ax_main.text(0.6, y, 'None', fontsize=FONT_BODY, color=COLOR_TEXT_MUTED)
        y -= 0.24
    y -= 0.25

    # --- Novel & uncharacterized genes ---
    y = _draw_section_header(ax_main, 'Novel & uncharacterized genes', y, color)

    discovery_genes = pd.concat([novel, uncharacterized]).sort_values('priority')
    max_show = 10

    for idx, (_, row) in enumerate(discovery_genes.iterrows()):
        if idx >= max_show or y < 0.8:
            remaining = len(discovery_genes) - idx
            if remaining > 0:
                ax_main.text(0.6, y, f"+ {remaining} more", fontsize=5.5, color=COLOR_TEXT_MUTED, style='italic')
            break

        gene = row['gene']
        gene_type = row['gene_type']
        desc = row['description'] if pd.notna(row['description']) else ''

        # Gene name with badge
        badge = 'novel' if gene_type == 'novel' else 'unc.'
        ax_main.text(0.6, y, gene, fontsize=FONT_GENE, color=color,
                     family='monospace', fontweight='bold')
        # Badge next to gene name
        gene_width = len(gene) * 0.12 + 0.7
        ax_main.text(gene_width, y + 0.01, badge, fontsize=4.5,
                     color=COLOR_TEXT_MUTED, style='italic')
        y -= 0.24

        # Description
        if desc:
            desc_short = desc[:200] + '...' if len(desc) > 200 else desc
            wrapped = textwrap.fill(desc_short, width=100)
            for line in wrapped.split('\n')[:2]:
                ax_main.text(0.7, y, line, fontsize=5.2, color=COLOR_TEXT_SECONDARY)
                y -= 0.2
        y -= 0.1

    # Gene count footer
    n_est = len(est_genes)
    n_novel = len(novel)
    n_unc = len(uncharacterized)
    footer = f"{n_est} established  ·  {n_novel} novel  ·  {n_unc} uncharacterized"
    ax_main.text(5, 0.4, footer, fontsize=5.5, color=COLOR_TEXT_MUTED,
                 ha='center', va='center')

    # Adjust y limits
    ax_main.set_ylim(max(0, min(y - 0.3, 0.1)), 14)

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'pdf']:
        output_path = output_dir / f'card_cluster_{cluster_id}.{ext}'
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                    facecolor='white', pad_inches=0.05)
        print(f'  Saved: {output_path}')
    plt.close()


def _style_subplot(ax, title):
    """Apply consistent styling to PHATE/volcano subplots."""
    ax.tick_params(axis='both', labelsize=FONT_AXIS_TICK, colors=COLOR_TEXT_MUTED,
                   length=2, width=0.5, pad=1)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color(COLOR_BORDER)
        ax.spines[spine].set_linewidth(0.5)
    ax.set_title(title, fontsize=FONT_AXIS_LABEL + 1, color=COLOR_TEXT_SECONDARY,
                 fontweight='medium', pad=3)


def _draw_section_header(ax, title, y, accent_color):
    """Draw a section header with a subtle accent line."""
    ax.plot([0.45, 9.5], [y + 0.12, y + 0.12], color=COLOR_BORDER, linewidth=0.5)
    ax.text(0.45, y, title, fontsize=FONT_SECTION, color=COLOR_TEXT_PRIMARY, fontweight='bold')
    y -= 0.32
    return y


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Generate publication cluster cards.')
    parser.add_argument('--clusters', type=int, nargs='+', default=None,
                        help='Cluster IDs to generate (default: 0, 221)')
    parser.add_argument('--all', action='store_true',
                        help='Generate cards for all high-confidence clusters')
    parser.add_argument('--output-dir', type=str, default=None,
                        help=f'Output directory (default: {OUTPUT_DIR})')
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    table_file = output_dir / 'cluster_genes_table.tsv'

    print("Loading data...")
    df_table = pd.read_csv(table_file, sep='\t')
    df_clustering = pd.read_csv(PHATE_FILE, sep='\t')
    df_effects = pd.read_csv(EFFECTS_FILE, sep='\t')
    df_pvals = pd.read_csv(PVALS_FILE, sep='\t')

    # Determine which clusters to generate
    if args.all:
        cluster_ids = sorted(df_table['cluster_id'].unique())
    elif args.clusters:
        cluster_ids = args.clusters
    else:
        cluster_ids = DEFAULT_CLUSTERS

    print(f"\nGenerating cards for {len(cluster_ids)} cluster(s): {cluster_ids}\n")

    for cluster_id in cluster_ids:
        print(f"Cluster {cluster_id}:")
        create_card(cluster_id, df_table, df_clustering, df_effects, df_pvals, output_dir)

    print(f'\nDone! {len(cluster_ids)} card(s) generated in {output_dir}')


if __name__ == '__main__':
    main()
