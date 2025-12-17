"""Run mozzarellm analysis on OPS clustering data.

This script analyzes gene clusters from an optical pooled screen (OPS).
Gene-wise data is reshaped to cluster format and analyzed using LLMs.

Features:
- Incremental saving: Results saved after each cluster (resumable)
- Progress tracking: Shows which clusters are done
- Resume support: Automatically skips already-analyzed clusters
"""

import argparse
import json
import os
import sys

import pandas as pd
from dotenv import load_dotenv
from tqdm import tqdm

from mozzarellm import ClusterAnalyzer, reshape_to_clusters

# Load environment variables from .env file
load_dotenv()

# Default configuration
DEFAULT_MODEL = "claude-sonnet-4-5-20250929"
DEFAULT_TEMPERATURE = 0.0
DEFAULT_INPUT_FILE = "brieflow_output/cluster/Hoescht_COX4_AGP_ConA/all/filtered/10/phate_leiden_clustering.tsv"
DEFAULT_OUTPUT_DIR = "brieflow_output/cluster/Hoescht_COX4_AGP_ConA/all/filtered/10/mozzarellm"

# Screen context for OPS data
SCREEN_CONTEXT = """
These clusters are from an optical pooled screen (OPS) that measured morphological
phenotypes in human cells. The screen involved perturbing genes using CRISPR knockout
system and imaging the resulting cellular morphology via fluorescence microscopy adapted
from the cell-painting panel, specifically: stains of Hoechst (nucleus), 
COX4 (mitochondria), AGP (actin/golgi/plasma membrane), and ConA (endoplasmic reticulum).
Genes grouped within a cluster tend to exhibit similar phenotypes, suggesting they
may participate in the same biological process or pathway.
"""


def load_existing_results(output_dir, model_name):
    """Load any existing cluster results for resume support."""
    results_dir = os.path.join(output_dir, "clusters")
    completed = {}
    if os.path.exists(results_dir):
        for f in os.listdir(results_dir):
            if f.endswith(".json"):
                try:
                    with open(os.path.join(results_dir, f)) as fh:
                        data = json.load(fh)
                        cluster_id = data.get("cluster_id")
                        if cluster_id is not None:
                            completed[str(cluster_id)] = data  # Ensure string key
                except (json.JSONDecodeError, KeyError):
                    continue
    return completed


def save_cluster_result(output_dir, cluster_id, result_dict):
    """Save a single cluster result to disk."""
    results_dir = os.path.join(output_dir, "clusters")
    os.makedirs(results_dir, exist_ok=True)
    filepath = os.path.join(results_dir, f"cluster_{cluster_id}.json")
    with open(filepath, "w") as f:
        json.dump(result_dict, f, indent=2)


def combine_results(output_dir, model_name):
    """Combine all cluster results into final output files."""
    results_dir = os.path.join(output_dir, "clusters")
    all_results = []

    for f in sorted(os.listdir(results_dir), key=lambda x: int(x.split("_")[1].split(".")[0]) if x.startswith("cluster_") else 0):
        if f.endswith(".json"):
            with open(os.path.join(results_dir, f)) as fh:
                all_results.append(json.load(fh))

    # Save combined JSON
    output_base = os.path.join(output_dir, f"{model_name.replace('/', '_')}_results")
    combined_json = output_base + ".json"
    with open(combined_json, "w") as f:
        json.dump({"clusters": all_results, "model": model_name}, f, indent=2)

    # Save summaries TSV
    summaries = []
    for r in all_results:
        summaries.append({
            "cluster_id": r.get("cluster_id"),
            "dominant_process": r.get("dominant_process"),
            "pathway_confidence": r.get("pathway_confidence"),
            "summary": r.get("summary"),
            "num_established": len(r.get("established_genes", [])),
            "num_novel": len(r.get("novel_role_genes", [])),
            "num_uncharacterized": len(r.get("uncharacterized_genes", [])),
        })
    summaries_df = pd.DataFrame(summaries)
    summaries_tsv = output_base + "_summaries.tsv"
    summaries_df.to_csv(summaries_tsv, sep="\t", index=False)

    return combined_json, summaries_tsv


def main():
    """Run the mozzarellm analysis."""
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Run mozzarellm analysis on OPS clustering data")
    parser.add_argument(
        "--model",
        type=str,
        default=DEFAULT_MODEL,
        help=f"Model to use (default: {DEFAULT_MODEL})",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=DEFAULT_TEMPERATURE,
        help=f"Temperature for model (default: {DEFAULT_TEMPERATURE})",
    )
    parser.add_argument(
        "--input-file",
        type=str,
        default=DEFAULT_INPUT_FILE,
        help=f"Input clustering TSV file (default: {DEFAULT_INPUT_FILE})",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous run (skip completed clusters)",
    )

    args = parser.parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Resolve paths relative to script directory
    input_file = os.path.join(script_dir, args.input_file)
    output_dir = os.path.join(script_dir, args.output_dir)

    print(f"Mozzarellm Analysis - Model: {args.model}")
    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load gene-level data
    print("Loading clustering data...")
    gene_df = pd.read_csv(input_file, sep="\t")

    # Rename gene_symbol_0 to gene_symbol if needed
    if "gene_symbol_0" in gene_df.columns and "gene_symbol" not in gene_df.columns:
        gene_df = gene_df.rename(columns={"gene_symbol_0": "gene_symbol"})

    print(f"Loaded {len(gene_df)} genes across {gene_df['cluster'].nunique()} clusters")

    # Reshape to cluster-level format using mozzarellm's utility
    print("Reshaping data to cluster format...")
    cluster_df, gene_annotations = reshape_to_clusters(
        input_df=gene_df,
        gene_col="gene_symbol",
        cluster_col="cluster",
        uniprot_col="uniprot_function",
        verbose=True,
    )

    print(f"Reshaped to {len(cluster_df)} clusters")

    # Check for existing results (resume support)
    completed = load_existing_results(output_dir, args.model)
    if completed:
        print(f"Found {len(completed)} previously completed clusters")

    # Initialize analyzer (with show_progress=False since we handle our own progress)
    analyzer = ClusterAnalyzer(
        model=args.model, temperature=args.temperature, show_progress=False
    )

    # Process clusters one at a time with incremental saving
    print("Running LLM analysis (with incremental saving)...")
    total_clusters = len(cluster_df)
    skipped = 0

    for i in tqdm(range(total_clusters), desc="Analyzing clusters"):
        row = cluster_df.iloc[i]
        cluster_id = str(row["cluster_id"])  # Ensure string for consistency

        # Skip if already completed
        if cluster_id in completed:
            skipped += 1
            continue

        # Create single-cluster DataFrame
        single_cluster_df = pd.DataFrame([row])

        # Analyze this cluster
        try:
            result = analyzer.analyze(
                single_cluster_df,
                gene_annotations=gene_annotations,
                screen_context=SCREEN_CONTEXT,
            )

            # Extract the cluster result and save immediately
            # result.clusters is a dict[str, ClusterResult]
            if result.clusters:
                # Get the first (and only) cluster result from the dict
                result_cluster_id = list(result.clusters.keys())[0]
                cluster_result = result.clusters[result_cluster_id].model_dump()
                # Use our cluster_id for consistent naming
                cluster_result["cluster_id"] = cluster_id
                save_cluster_result(output_dir, cluster_id, cluster_result)
                completed[cluster_id] = cluster_result
                tqdm.write(f"  Saved cluster {cluster_id}")

        except Exception as e:
            tqdm.write(f"Error analyzing cluster {cluster_id}: {e}")
            import traceback
            traceback.print_exc()
            continue

    if skipped > 0:
        print(f"Skipped {skipped} previously completed clusters")

    # Combine all results into final files
    print("Combining results...")
    combined_json, summaries_tsv = combine_results(output_dir, args.model)

    print(f"\nAnalysis complete!")
    print(f"  - Individual results: {output_dir}/clusters/")
    print(f"  - Combined JSON: {combined_json}")
    print(f"  - Summaries TSV: {summaries_tsv}")


if __name__ == "__main__":
    main()
