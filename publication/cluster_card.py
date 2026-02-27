"""
Cluster Card Generator - Reads from consolidated table.
Labels all novel/uncharacterized genes on both plots.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import textwrap
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
OUTPUT_DIR = Path("/mnt/data/blainey/whitney-analysis/figures/LLM")

PHATE_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "phate_leiden_clustering.tsv"
EFFECTS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "tsvs" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__features_genes.tsv"
PVALS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "bootstrap" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__all_gene_bootstrap_results.tsv"
TABLE_FILE = OUTPUT_DIR / "cluster_genes_table.tsv"

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...")

# Main table with all gene info
df_table = pd.read_csv(TABLE_FILE, sep='\t')

# Full clustering for background PHATE plot
df_clustering = pd.read_csv(PHATE_FILE, sep='\t')

# Full effects/pvals for background volcano
df_effects = pd.read_csv(EFFECTS_FILE, sep='\t')
df_pvals = pd.read_csv(PVALS_FILE, sep='\t')

# Get unique clusters
cluster_ids = df_table['cluster_id'].unique()

# Colors - vibrant palette for pathway distinction
COLORS = [
    '#8b5cf6', '#06b6d4', '#f97316', '#22c55e', '#ef4444', '#3b82f6',
    '#ec4899', '#14b8a6', '#f59e0b', '#6366f1', '#84cc16', '#a855f7',
    '#0ea5e9', '#f43f5e', '#10b981', '#8b5cf6', '#06b6d4', '#f97316',
    '#22c55e', '#ef4444', '#3b82f6', '#ec4899', '#14b8a6', '#f59e0b',
    '#6366f1', '#84cc16', '#a855f7'
]


def get_pathway_color(pathway):
    """Get consistent color based on pathway name hash."""
    hash_val = hash(pathway) % len(COLORS)
    return COLORS[hash_val]

# =============================================================================
# CARD GENERATOR
# =============================================================================

def create_card(cluster_id):
    """Generate a cluster card from the table."""

    # Get cluster data
    cluster_df = df_table[df_table['cluster_id'] == cluster_id].copy()
    if len(cluster_df) == 0:
        print(f"No data for cluster {cluster_id}")
        return

    pathway = cluster_df['pathway'].iloc[0]
    feature = cluster_df['best_feature'].iloc[0]
    cluster_summary = cluster_df['cluster_summary'].iloc[0]
    color = get_pathway_color(pathway)  # Color based on pathway hash

    # Split by gene type
    established = cluster_df[cluster_df['gene_type'] == 'established']
    novel = cluster_df[cluster_df['gene_type'] == 'novel'].sort_values('priority')
    uncharacterized = cluster_df[cluster_df['gene_type'] == 'uncharacterized'].sort_values('priority')

    # Build volcano data for background
    pval_col = f'{feature}_pval'
    df_volcano = df_effects[['gene_symbol_0', feature]].merge(
        df_pvals[['gene', pval_col]],
        left_on='gene_symbol_0', right_on='gene', how='inner'
    )
    df_volcano['neg_log10_pval'] = -np.log10(df_volcano[pval_col].clip(lower=1e-10))

    # Create figure
    fig = plt.figure(figsize=(5, 6))
    ax_main = fig.add_axes([0, 0, 1, 1])
    ax_main.set_xlim(0, 10)
    ax_main.set_ylim(0, 12)
    ax_main.axis('off')

    # Card background
    bg = mpatches.FancyBboxPatch((0.1, 0.1), 9.8, 11.8,
                                  boxstyle="round,pad=0.02,rounding_size=0.15",
                                  facecolor='white', edgecolor='#d1d5db', linewidth=1)
    ax_main.add_patch(bg)

    # Color accent bar (taller to accommodate two-line title)
    accent = mpatches.Rectangle((0.1, 11.35), 9.8, 0.55, facecolor=color, edgecolor='none')
    ax_main.add_patch(accent)

    # Pathway title - full text, wrap to two lines if needed
    full_title = f"{pathway}  (Cluster {cluster_id})"
    if len(full_title) > 55:
        # Split into two lines
        wrapped_title = textwrap.fill(pathway, width=50)
        title_lines = wrapped_title.split('\n')
        ax_main.text(5, 11.75, title_lines[0], fontsize=9, fontweight='bold',
                     color='white', va='center', ha='center')
        if len(title_lines) > 1:
            second_line = ' '.join(title_lines[1:]) + f"  (Cluster {cluster_id})"
            ax_main.text(5, 11.52, second_line, fontsize=8, fontweight='bold',
                         color='white', va='center', ha='center')
        else:
            ax_main.text(5, 11.52, f"(Cluster {cluster_id})", fontsize=8, fontweight='bold',
                         color='white', va='center', ha='center')
    else:
        ax_main.text(5, 11.62, full_title, fontsize=9, fontweight='bold',
                     color='white', va='center', ha='center')

    # =========================================================================
    # PHATE PLOT
    # =========================================================================
    ax_phate = fig.add_axes([0.06, 0.72, 0.40, 0.18])

    # Background - all genes
    ax_phate.scatter(df_clustering['PHATE_0'], df_clustering['PHATE_1'],
                     c='#e5e7eb', s=0.5, alpha=0.5, rasterized=True)

    # Cluster genes
    cluster_genes = cluster_df['gene'].tolist()
    cluster_phate = df_clustering[df_clustering['gene_symbol_0'].isin(cluster_genes)]
    ax_phate.scatter(cluster_phate['PHATE_0'], cluster_phate['PHATE_1'], c=color, s=5, alpha=0.9)

    # No gene labels on plot - just colored points

    # Show full PHATE embedding (all data) with cluster highlighted
    x_min, x_max = df_clustering['PHATE_0'].min(), df_clustering['PHATE_0'].max()
    y_min, y_max = df_clustering['PHATE_1'].min(), df_clustering['PHATE_1'].max()
    padding = 0.05 * max(x_max - x_min, y_max - y_min)
    ax_phate.set_xlim(x_min - padding, x_max + padding)
    ax_phate.set_ylim(y_min - padding, y_max + padding)
    ax_phate.tick_params(axis='both', labelsize=4, colors='#9ca3af', length=2)
    ax_phate.spines['top'].set_visible(False)
    ax_phate.spines['right'].set_visible(False)
    ax_phate.spines['bottom'].set_color('#d1d5db')
    ax_phate.spines['left'].set_color('#d1d5db')
    ax_phate.set_title('PHATE', fontsize=6, color='#6b7280', pad=2)

    # =========================================================================
    # VOLCANO PLOT
    # =========================================================================
    ax_volcano = fig.add_axes([0.54, 0.72, 0.40, 0.18])

    # Background - all genes
    ax_volcano.scatter(df_volcano[feature], df_volcano['neg_log10_pval'],
                       c='#e5e7eb', s=0.5, alpha=0.3, rasterized=True)

    # Cluster genes
    cluster_volcano = cluster_df[cluster_df['effect_size'].notna()]
    ax_volcano.scatter(cluster_volcano['effect_size'], cluster_volcano['neg_log10_pval'],
                       c=color, s=6, alpha=0.9)

    # No gene labels on plot - just colored points

    ax_volcano.set_xlabel('Effect size', fontsize=5, color='#6b7280')
    ax_volcano.set_ylabel('-log10(p)', fontsize=5, color='#6b7280')
    ax_volcano.tick_params(axis='both', labelsize=4, colors='#9ca3af', length=2)
    ax_volcano.spines['top'].set_visible(False)
    ax_volcano.spines['right'].set_visible(False)
    ax_volcano.spines['bottom'].set_color('#d1d5db')
    ax_volcano.spines['left'].set_color('#d1d5db')
    feature_name = feature.replace('_', ' ').replace('cell ', '').replace('mean', '').strip()
    ax_volcano.set_title(feature_name, fontsize=6, color='#6b7280', pad=2)

    # =========================================================================
    # CLUSTER SUMMARY SECTION
    # =========================================================================
    y = 8.1  # Moved down to give space from plots

    ax_main.text(0.4, y, 'Summary', fontsize=7, color='#374151', fontweight='bold')
    y -= 0.28

    # Display truncated cluster summary
    if pd.notna(cluster_summary):
        summary_short = cluster_summary[:350] + '...' if len(cluster_summary) > 350 else cluster_summary
        wrapped = textwrap.fill(summary_short, width=115)
        lines = wrapped.split('\n')
        for line in lines[:4]:  # Up to 4 lines
            ax_main.text(0.5, y, line, fontsize=5, color='#6b7280', style='italic')
            y -= 0.2
    y -= 0.15

    # =========================================================================
    # ESTABLISHED GENES SECTION
    # =========================================================================
    ax_main.text(0.4, y, 'Established genes', fontsize=7, color='#374151', fontweight='bold')
    y -= 0.28

    est_genes = established['gene'].tolist()
    if est_genes:
        # Wrap long gene lists - balanced wrapping
        gene_str = ', '.join(est_genes)
        wrapped = textwrap.fill(gene_str, width=98)
        lines = wrapped.split('\n')
        for line in lines[:3]:  # Max 3 lines
            ax_main.text(0.5, y, line, fontsize=5.5, color='#374151', family='monospace')
            y -= 0.22
    y -= 0.15

    # =========================================================================
    # NOVEL & UNCHARACTERIZED GENES WITH DESCRIPTIONS
    # =========================================================================
    ax_main.text(0.4, y, 'Novel & Uncharacterized genes', fontsize=7, color='#374151', fontweight='bold')
    y -= 0.32

    # Combine novel and uncharacterized, sort by priority - show fewer genes
    discovery_genes = pd.concat([novel, uncharacterized]).sort_values('priority')
    max_genes_to_show = 8  # Limit to give more breathing room

    for idx, (_, row) in enumerate(discovery_genes.iterrows()):
        if idx >= max_genes_to_show:
            remaining = len(discovery_genes) - max_genes_to_show
            if remaining > 0:
                ax_main.text(0.5, y, f"... and {remaining} more", fontsize=5, color='#9ca3af', style='italic')
            break

        gene = row['gene']
        gene_type = row['gene_type']
        desc = row['description'] if pd.notna(row['description']) else ''

        # Gene name with type indicator
        type_label = '(novel)' if gene_type == 'novel' else '(unc)'
        ax_main.text(0.5, y, f"{gene} {type_label}", fontsize=6, color=color,
                     family='monospace', fontweight='bold')
        y -= 0.22

        # Description (truncated) - balanced wrapping
        if desc:
            desc_short = desc[:180] + '...' if len(desc) > 180 else desc
            wrapped = textwrap.fill(desc_short, width=115)
            lines = wrapped.split('\n')
            for line in lines[:2]:  # Max 2 lines per gene
                ax_main.text(0.6, y, line, fontsize=5, color='#6b7280')
                y -= 0.18
        y -= 0.08

        # Stop if we're running out of space
        if y < 0.5:
            ax_main.text(0.5, y, f"... and {len(discovery_genes) - discovery_genes.index.get_loc(row.name) - 1} more",
                        fontsize=5, color='#9ca3af', style='italic')
            break

    # Adjust y limits
    ax_main.set_ylim(max(0, y - 0.3), 12)

    # Save
    output_path = OUTPUT_DIR / f'card_cluster_{cluster_id}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.02)
    plt.close()
    print(f'Saved: {output_path}')


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print(f"\nGenerating cards for {len(cluster_ids)} clusters...\n")

    for cluster_id in sorted(cluster_ids):
        create_card(cluster_id)

    print(f'\nAll {len(cluster_ids)} cards generated!')
