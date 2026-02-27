"""
Generate consolidated cluster table for card generation.
One row per gene, with all data needed for visualization.
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================
BRIEFLOW_OUTPUT = Path("/mnt/data/blainey/whitney-analysis/analysis/brieflow_output")
OUTPUT_DIR = Path("/mnt/data/blainey/whitney-analysis/figures/LLM")

PHATE_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "phate_leiden_clustering.tsv"
EFFECTS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "tsvs" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__features_genes.tsv"
PVALS_FILE = BRIEFLOW_OUTPUT / "aggregate" / "bootstrap" / "CeCl-all_ChCo-Hoescht_COX4_AGP_ConA__all_gene_bootstrap_results.tsv"
LLM_RESULTS_FILE = BRIEFLOW_OUTPUT / "cluster" / "Hoescht_COX4_AGP_ConA" / "all" / "filtered" / "10" / "mozzarellm" / "claude-sonnet-4-5-20250929_results.json"

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...")

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df_clustering = pd.read_csv(PHATE_FILE, sep='\t')
df_effects = pd.read_csv(EFFECTS_FILE, sep='\t')
df_pvals = pd.read_csv(PVALS_FILE, sep='\t')

with open(LLM_RESULTS_FILE) as f:
    json_data = json.load(f)

# Index clusters by ID
clusters_by_id = {c['cluster_id']: c for c in json_data['clusters']}

# High confidence clusters
high_conf_ids = [c['cluster_id'] for c in json_data['clusters'] if c['pathway_confidence'] == 'High']

# Features to consider
FEATURES = [
    'cell_Hoescht_mean', 'cell_COX4_mean', 'cell_AGP_mean', 'cell_ConA_mean',
    'nucleus_Hoescht_mean', 'nucleus_COX4_mean', 'nucleus_AGP_mean', 'nucleus_ConA_mean',
    'cytoplasm_Hoescht_mean', 'cytoplasm_COX4_mean', 'cytoplasm_AGP_mean', 'cytoplasm_ConA_mean',
]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_best_feature(cluster_genes):
    """Find the feature with most significant p-values for cluster genes."""
    best_feature = None
    best_score = -np.inf

    for feature in FEATURES:
        pval_col = f'{feature}_pval'
        if pval_col not in df_pvals.columns:
            continue
        cluster_pvals = df_pvals[df_pvals['gene'].isin(cluster_genes)][pval_col].dropna()
        if len(cluster_pvals) == 0:
            continue
        score = (-np.log10(cluster_pvals.clip(lower=1e-10))).sum()
        if score > best_score:
            best_score = score
            best_feature = feature

    return best_feature or 'cell_Hoescht_mean'


def get_gene_stats(gene, feature):
    """Get effect size and p-value for a gene."""
    pval_col = f'{feature}_pval'

    eff_row = df_effects[df_effects['gene_symbol_0'] == gene]
    pval_row = df_pvals[df_pvals['gene'] == gene]

    effect_size = eff_row[feature].values[0] if len(eff_row) > 0 and feature in eff_row.columns else np.nan
    pval = pval_row[pval_col].values[0] if len(pval_row) > 0 and pval_col in pval_row.columns else np.nan

    return effect_size, pval


def get_phate_coords(gene):
    """Get PHATE coordinates for a gene."""
    row = df_clustering[df_clustering['gene_symbol_0'] == gene]
    if len(row) > 0:
        return row['PHATE_0'].values[0], row['PHATE_1'].values[0]
    return np.nan, np.nan


# =============================================================================
# BUILD TABLE
# =============================================================================
print("Building table...")

rows = []

for cluster_id_str in high_conf_ids:
    cluster_id = int(cluster_id_str)
    cluster_data = clusters_by_id[cluster_id_str]

    pathway = cluster_data['dominant_process']
    confidence = cluster_data['pathway_confidence']
    cluster_summary = cluster_data.get('summary', '')

    # Get genes from clustering file for this cluster
    cluster_genes = df_clustering[df_clustering['cluster'] == cluster_id]['gene_symbol_0'].tolist()

    # Find best feature for this cluster
    best_feature = get_best_feature(cluster_genes)

    # Process established genes (list of strings)
    for gene in cluster_data.get('established_genes', []):
        effect_size, pval = get_gene_stats(gene, best_feature)
        phate_0, phate_1 = get_phate_coords(gene)
        rows.append({
            'cluster_id': cluster_id,
            'pathway': pathway,
            'confidence': confidence,
            'cluster_summary': cluster_summary,
            'best_feature': best_feature,
            'gene': gene,
            'gene_type': 'established',
            'priority': None,
            'description': '',
            'effect_size': effect_size,
            'pval': pval,
            'neg_log10_pval': -np.log10(pval) if pval and pval > 0 else np.nan,
            'phate_0': phate_0,
            'phate_1': phate_1,
        })

    # Process novel role genes (list of dicts)
    for g in cluster_data.get('novel_role_genes', []):
        gene = g['gene']
        effect_size, pval = get_gene_stats(gene, best_feature)
        phate_0, phate_1 = get_phate_coords(gene)
        rows.append({
            'cluster_id': cluster_id,
            'pathway': pathway,
            'confidence': confidence,
            'cluster_summary': cluster_summary,
            'best_feature': best_feature,
            'gene': gene,
            'gene_type': 'novel',
            'priority': g.get('priority'),
            'description': g.get('rationale', ''),
            'effect_size': effect_size,
            'pval': pval,
            'neg_log10_pval': -np.log10(pval) if pval and pval > 0 else np.nan,
            'phate_0': phate_0,
            'phate_1': phate_1,
        })

    # Process uncharacterized genes (list of dicts)
    for g in cluster_data.get('uncharacterized_genes', []):
        gene = g['gene']
        effect_size, pval = get_gene_stats(gene, best_feature)
        phate_0, phate_1 = get_phate_coords(gene)
        rows.append({
            'cluster_id': cluster_id,
            'pathway': pathway,
            'confidence': confidence,
            'cluster_summary': cluster_summary,
            'best_feature': best_feature,
            'gene': gene,
            'gene_type': 'uncharacterized',
            'priority': g.get('priority'),
            'description': g.get('rationale', ''),
            'effect_size': effect_size,
            'pval': pval,
            'neg_log10_pval': -np.log10(pval) if pval and pval > 0 else np.nan,
            'phate_0': phate_0,
            'phate_1': phate_1,
        })

# Create DataFrame
df = pd.DataFrame(rows)

# Sort by cluster_id, then by gene_type (established first), then by priority
type_order = {'established': 0, 'novel': 1, 'uncharacterized': 2}
df['type_order'] = df['gene_type'].map(type_order)
df = df.sort_values(['cluster_id', 'type_order', 'priority']).drop('type_order', axis=1)

# Save
output_path = OUTPUT_DIR / 'cluster_genes_table.tsv'
df.to_csv(output_path, sep='\t', index=False)
print(f"\nSaved: {output_path}")
print(f"Total rows: {len(df)}")
print(f"Clusters: {df['cluster_id'].nunique()}")
print(f"\nGene types:")
print(df['gene_type'].value_counts())
print(f"\nSample rows:")
print(df[['cluster_id', 'gene', 'gene_type', 'priority', 'best_feature', 'effect_size', 'pval']].head(10))
