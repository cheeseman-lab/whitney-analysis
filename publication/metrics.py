#!/usr/bin/env python3
"""
Standalone metrics script for assessing brieflow output.

This script directly reads from the brieflow output directory structure
and computes metrics without requiring the full brieflow config.
"""

import argparse
import json
import pandas as pd
from pathlib import Path


def get_preprocess_stats(output_dir: Path) -> dict:
    """Get preprocessing statistics including tile counts."""
    preprocess_dir = output_dir / "preprocess"

    # Count TIFF files
    sbs_tiff_dir = preprocess_dir / "images" / "sbs"
    phenotype_tiff_dir = preprocess_dir / "images" / "phenotype"

    sbs_tiff_count = len(list(sbs_tiff_dir.glob("**/*.tiff"))) if sbs_tiff_dir.exists() else 0
    phenotype_tiff_count = len(list(phenotype_tiff_dir.glob("**/*.tiff"))) if phenotype_tiff_dir.exists() else 0

    return {
        "sbs_tiles": sbs_tiff_count,
        "phenotype_tiles": phenotype_tiff_count,
        "total_tiles": sbs_tiff_count + phenotype_tiff_count,
    }


def get_sbs_stats(output_dir: Path) -> dict:
    """Get SBS (Sequencing by Synthesis) statistics."""
    sbs_eval_dir = output_dir / "sbs" / "eval"

    # Load segmentation overview data
    seg_overview_files = list(sbs_eval_dir.glob("**/segmentation/*__segmentation_overview.tsv"))

    total_cells = 0
    plate_stats = {}

    if seg_overview_files:
        for file in seg_overview_files:
            df = pd.read_csv(file, sep="\t")
            plate_name = file.name.split("__")[0]
            plate_cells = df["final_cells"].sum()
            plate_stats[plate_name] = {"segmented_cells": int(plate_cells)}
            total_cells += plate_cells

    # Load mapping overview data
    mapping_overview_files = list(sbs_eval_dir.glob("**/mapping/*__mapping_overview.tsv"))

    mapping_stats = {
        "total_cells_for_mapping": 0,
        "cells_with_1_barcode": 0,
        "cells_with_1_or_more_barcodes": 0,
        "cells_with_1_gene": 0,
        "cells_with_1_or_more_genes": 0,
    }

    if mapping_overview_files:
        for file in mapping_overview_files:
            df = pd.read_csv(file, sep="\t")
            plate_name = file.name.split("__")[0]

            # Sum across all wells in this plate
            mapping_stats["total_cells_for_mapping"] += df["total_cells__count"].sum()
            mapping_stats["cells_with_1_barcode"] += df["1_barcode_cells__count"].sum()
            mapping_stats["cells_with_1_or_more_barcodes"] += df["1_or_more_barcodes__count"].sum()
            mapping_stats["cells_with_1_gene"] += df["1_gene_cells__count"].sum()
            mapping_stats["cells_with_1_or_more_genes"] += df["1_or_more_genes__count"].sum()

            # Store plate-level mapping percentages
            if plate_name in plate_stats:
                plate_stats[plate_name]["avg_1_barcode_pct"] = df["1_barcode_cells__percent"].mean()
                plate_stats[plate_name]["avg_1_gene_pct"] = df["1_gene_cells__percent"].mean()

    # Calculate overall percentages
    total_for_mapping = mapping_stats["total_cells_for_mapping"]
    if total_for_mapping > 0:
        mapping_stats["pct_with_1_barcode"] = 100 * mapping_stats["cells_with_1_barcode"] / total_for_mapping
        mapping_stats["pct_with_1_or_more_barcodes"] = 100 * mapping_stats["cells_with_1_or_more_barcodes"] / total_for_mapping
        mapping_stats["pct_with_1_gene"] = 100 * mapping_stats["cells_with_1_gene"] / total_for_mapping
        mapping_stats["pct_with_1_or_more_genes"] = 100 * mapping_stats["cells_with_1_or_more_genes"] / total_for_mapping

    return {
        "total_segmented_cells": int(total_cells),
        "plate_stats": plate_stats,
        "mapping": mapping_stats,
    }


def get_phenotype_stats(output_dir: Path) -> dict:
    """Get phenotype statistics."""
    phenotype_eval_dir = output_dir / "phenotype" / "eval" / "segmentation"
    phenotype_parquet_dir = output_dir / "phenotype" / "parquets"

    # Load segmentation overview files
    seg_overview_files = list(phenotype_eval_dir.glob("*__segmentation_overview.tsv"))

    total_cells = 0
    plate_stats = {}

    if seg_overview_files:
        for file in seg_overview_files:
            df = pd.read_csv(file, sep="\t")
            plate_name = file.name.split("__")[0]
            plate_cells = df["final_cells"].sum()
            initial_cells = df["initial_cells"].sum() if "initial_cells" in df.columns else 0

            plate_stats[plate_name] = {
                "initial_cells": int(initial_cells),
                "final_cells": int(plate_cells),
                "wells": len(df),
            }
            total_cells += plate_cells

    # Get feature count from parquet schema
    feature_count = 0
    sample_parquet_files = list(phenotype_parquet_dir.glob("**/*__phenotype_cp.parquet"))

    if sample_parquet_files:
        import pyarrow.parquet as pq
        parquet_schema = pq.read_schema(sample_parquet_files[0])
        all_columns = parquet_schema.names
        # Count non-metadata columns (rough estimate)
        metadata_prefixes = ["plate", "well", "tile", "cell", "i", "j", "x", "y", "label"]
        feature_cols = [c for c in all_columns if not any(c.lower().startswith(p) for p in metadata_prefixes)]
        feature_count = len(feature_cols)

    return {
        "total_cells": int(total_cells),
        "feature_count": feature_count,
        "plate_stats": plate_stats,
    }


def get_merge_pipeline_stats(output_dir: Path) -> dict:
    """Get detailed cell counts at each step of the merge pipeline."""
    import pyarrow.parquet as pq

    merge_parquet_dir = output_dir / "merge" / "parquets"
    sbs_parquet_dir = output_dir / "sbs" / "parquets"
    phenotype_parquet_dir = output_dir / "phenotype" / "parquets"
    aggregate_parquet_dir = output_dir / "aggregate" / "parquets"

    pipeline = {}

    # 1. Input: SBS cells (cells with barcodes)
    sbs_cells_files = list(sbs_parquet_dir.glob("*__cells.parquet"))
    sbs_cells_count = 0
    for f in sbs_cells_files:
        try:
            sbs_cells_count += pq.ParquetFile(f).metadata.num_rows
        except:
            pass
    pipeline["sbs_input_cells"] = sbs_cells_count

    # 2. Input: Phenotype cells
    phenotype_min_files = list(phenotype_parquet_dir.glob("*__phenotype_cp_min.parquet"))
    phenotype_cells_count = 0
    for f in phenotype_min_files:
        try:
            phenotype_cells_count += pq.ParquetFile(f).metadata.num_rows
        except:
            pass
    pipeline["phenotype_input_cells"] = phenotype_cells_count

    # 3. Merge pipeline steps
    merge_steps = ["fast_merge", "merge_formatted", "merge_final"]
    for step in merge_steps:
        step_files = list(merge_parquet_dir.glob(f"*__{step}.parquet"))
        step_count = 0
        for f in step_files:
            try:
                step_count += pq.ParquetFile(f).metadata.num_rows
            except:
                pass
        pipeline[step] = step_count

    # 4. Deduplication stats (aggregate from all wells)
    merge_eval_dir = output_dir / "merge" / "eval"
    dedup_files = list(merge_eval_dir.glob("*__deduplication_stats.tsv"))

    dedup_totals = {
        "initial": 0,
        "after_sbs_dedup": 0,
        "after_phenotype_dedup": 0,
        "initial_with_genes": 0,
        "after_sbs_dedup_with_genes": 0,
        "after_phenotype_dedup_with_genes": 0,
    }

    for file in dedup_files:
        try:
            df = pd.read_csv(file, sep="\t")
            for _, row in df.iterrows():
                stage = row["stage"].lower().replace(" ", "_")
                if stage == "initial":
                    dedup_totals["initial"] += row["total_cells"]
                    dedup_totals["initial_with_genes"] += row["mapped_genes"]
                elif stage == "after_sbs_dedup":
                    dedup_totals["after_sbs_dedup"] += row["total_cells"]
                    dedup_totals["after_sbs_dedup_with_genes"] += row["mapped_genes"]
                elif stage == "after_phenotype_dedup":
                    dedup_totals["after_phenotype_dedup"] += row["total_cells"]
                    dedup_totals["after_phenotype_dedup_with_genes"] += row["mapped_genes"]
        except:
            pass

    pipeline["deduplication"] = dedup_totals

    # 5. Aggregate pipeline: merge_data -> filtered
    merge_data_files = list(aggregate_parquet_dir.glob("*__merge_data.parquet"))
    merge_data_count = 0
    for f in merge_data_files:
        try:
            merge_data_count += pq.ParquetFile(f).metadata.num_rows
        except:
            pass
    pipeline["aggregate_merge_data"] = merge_data_count

    filtered_files = list(aggregate_parquet_dir.glob("*__filtered.parquet"))
    filtered_count = 0
    for f in filtered_files:
        try:
            filtered_count += pq.ParquetFile(f).metadata.num_rows
        except:
            pass
    pipeline["aggregate_filtered"] = filtered_count

    # Calculate losses at each step
    pipeline["losses"] = {}
    if pipeline["fast_merge"] > 0 and pipeline["merge_final"] > 0:
        pipeline["losses"]["deduplication"] = pipeline["fast_merge"] - pipeline["merge_final"]
    if pipeline["aggregate_merge_data"] > 0 and pipeline["aggregate_filtered"] > 0:
        pipeline["losses"]["aggregate_filtering"] = pipeline["aggregate_merge_data"] - pipeline["aggregate_filtered"]

    return pipeline


def get_merge_stats(output_dir: Path) -> dict:
    """Get merge statistics."""
    merge_eval_dir = output_dir / "merge" / "eval"

    # Get detailed pipeline stats
    pipeline_stats = get_merge_pipeline_stats(output_dir)

    # Find cell mapping stats files
    mapping_stats_files = list(merge_eval_dir.glob("*__cell_mapping_stats.tsv"))

    if not mapping_stats_files:
        return {"error": "No cell mapping statistics files found", "pipeline": pipeline_stats}

    total_mapped = 0
    total_unmapped = 0
    plate_stats = {}

    for file in mapping_stats_files:
        df = pd.read_csv(file, sep="\t")
        plate_name = file.name.split("__")[0]

        mapped_row = df[df["category"] == "mapped_cells"]
        unmapped_row = df[df["category"] == "unmapped_cells"]

        if not mapped_row.empty and not unmapped_row.empty:
            mapped = mapped_row["count"].iloc[0]
            unmapped = unmapped_row["count"].iloc[0]
            total = mapped + unmapped

            plate_stats[plate_name] = {
                "mapped_cells": int(mapped),
                "unmapped_cells": int(unmapped),
                "total_cells": int(total),
                "mapping_rate": 100 * mapped / total if total > 0 else 0,
            }

            total_mapped += mapped
            total_unmapped += unmapped

    # Load deduplication stats
    dedup_files = list(merge_eval_dir.glob("*__deduplication_stats.tsv"))
    dedup_stats = {}

    for file in dedup_files:
        df = pd.read_csv(file, sep="\t")
        well_name = file.name.replace("__deduplication_stats.tsv", "")
        if not df.empty:
            dedup_stats[well_name] = df.to_dict(orient="records")

    total = total_mapped + total_unmapped

    return {
        "total_cells": int(total),
        "total_mapped_cells": int(total_mapped),
        "total_unmapped_cells": int(total_unmapped),
        "overall_mapping_rate": 100 * total_mapped / total if total > 0 else 0,
        "plate_stats": plate_stats,
        "deduplication_samples": len(dedup_stats),
        "pipeline": pipeline_stats,
    }


def get_aggregate_stats(output_dir: Path) -> dict:
    """Get aggregation statistics."""
    aggregate_dir = output_dir / "aggregate"
    aggregated_tsvs = list((aggregate_dir / "tsvs").glob("*__aggregated.tsv"))

    results = {}

    for tsv_path in aggregated_tsvs:
        # Parse the filename to get cell_class and channel_combo
        # Format: CeCl-{cell_class}_ChCo-{channel_combo}__aggregated.tsv
        name = tsv_path.stem.replace("__aggregated", "")

        try:
            df = pd.read_csv(tsv_path, sep="\t")

            # Basic stats
            distinct_perturbations = len(df["gene_symbol_0"].unique())
            total_cells = df["cell_count"].sum()
            median_cells = df["cell_count"].median()
            min_cells = df["cell_count"].min()
            max_cells = df["cell_count"].max()

            # Count PC features
            pc_cols = [c for c in df.columns if c.startswith("PC_")]

            # Count controls
            controls = df[df["gene_symbol_0"].str.contains("nontargeting", case=False, na=False)]

            results[name] = {
                "distinct_perturbations": int(distinct_perturbations),
                "total_aggregated_cells": int(total_cells),
                "median_cells_per_perturbation": float(median_cells),
                "min_cells_per_perturbation": int(min_cells),
                "max_cells_per_perturbation": int(max_cells),
                "pc_features": len(pc_cols),
                "control_perturbations": len(controls),
            }
        except Exception as e:
            results[name] = {"error": str(e)}

    # Check for NA stats
    na_stats_files = list((aggregate_dir / "eval").glob("*__na_stats.tsv"))
    if na_stats_files:
        for f in na_stats_files:
            name = f.stem.replace("__na_stats", "")
            if name in results:
                try:
                    pd.read_csv(f, sep="\t")  # Verify file is readable
                    results[name]["na_stats_available"] = True
                except Exception:
                    pass

    return results


def get_cluster_stats_for_dir(base_dir: Path) -> list:
    """Get clustering statistics for a specific directory."""
    results = []

    # Find all resolutions (numeric directories)
    for resolution_dir in base_dir.iterdir():
        if not resolution_dir.is_dir():
            continue

        try:
            resolution = float(resolution_dir.name)
        except ValueError:
            continue

        # Read real metrics
        real_metrics_path = resolution_dir / "CB-Real__global_metrics.json"
        shuffled_metrics_path = resolution_dir / "CB-Shuffled__global_metrics.json"
        combined_table_path = resolution_dir / "CB-Real__combined_table.tsv"

        if not real_metrics_path.exists():
            continue

        try:
            with open(real_metrics_path, "r") as f:
                real_metrics = json.load(f)

            # Count clusters
            n_clusters = 0
            if combined_table_path.exists():
                combined_df = pd.read_csv(combined_table_path, sep="\t")
                n_clusters = len(combined_df["cluster"].unique())

            result = {
                "resolution": resolution,
                "n_clusters": n_clusters,
                "real_corum_enriched": real_metrics.get("CORUM", {}).get("num_enriched_clusters", 0),
                "real_corum_proportion": real_metrics.get("CORUM", {}).get("proportion_enriched", 0),
                "real_kegg_enriched": real_metrics.get("KEGG", {}).get("num_enriched_clusters", 0),
                "real_kegg_proportion": real_metrics.get("KEGG", {}).get("proportion_enriched", 0),
            }

            # Add STRING metrics if available
            if "STRING" in real_metrics:
                result["string_precision"] = real_metrics["STRING"].get("precision", 0)
                result["string_recall"] = real_metrics["STRING"].get("recall", 0)
                result["string_f1"] = real_metrics["STRING"].get("f1_score", 0)

            # Read shuffled metrics if available
            if shuffled_metrics_path.exists():
                with open(shuffled_metrics_path, "r") as f:
                    shuffled_metrics = json.load(f)
                result["shuffled_corum_enriched"] = shuffled_metrics.get("CORUM", {}).get("num_enriched_clusters", 0)
                result["shuffled_corum_proportion"] = shuffled_metrics.get("CORUM", {}).get("proportion_enriched", 0)
                result["shuffled_kegg_enriched"] = shuffled_metrics.get("KEGG", {}).get("num_enriched_clusters", 0)
                result["shuffled_kegg_proportion"] = shuffled_metrics.get("KEGG", {}).get("proportion_enriched", 0)

                # Calculate fold enrichment
                if result["shuffled_corum_proportion"] > 0:
                    result["corum_fold_enrichment"] = result["real_corum_proportion"] / result["shuffled_corum_proportion"]
                if result["shuffled_kegg_proportion"] > 0:
                    result["kegg_fold_enrichment"] = result["real_kegg_proportion"] / result["shuffled_kegg_proportion"]

            results.append(result)

        except Exception as e:
            results.append({
                "resolution": resolution,
                "error": str(e),
            })

    # Sort by resolution
    results.sort(key=lambda x: x.get("resolution", 0))
    return results


def get_cluster_stats(output_dir: Path) -> dict:
    """Get clustering statistics for both full and filtered clustering."""
    cluster_dir = output_dir / "cluster"

    if not cluster_dir.exists():
        return {"error": "No cluster directory found"}

    results = {"full": {}, "filtered": {}}

    # Find all channel combos
    for channel_combo_dir in cluster_dir.iterdir():
        if not channel_combo_dir.is_dir():
            continue

        channel_combo = channel_combo_dir.name

        # Find all cell classes
        for cell_class_dir in channel_combo_dir.iterdir():
            if not cell_class_dir.is_dir():
                continue

            cell_class = cell_class_dir.name
            key = f"{channel_combo}/{cell_class}"

            # Full clustering (direct resolutions in cell_class dir)
            full_resolutions = get_cluster_stats_for_dir(cell_class_dir)
            if full_resolutions:
                # Count total genes from aggregate_cleaned.tsv in parent
                aggregate_path = cell_class_dir / "aggregate_cleaned.tsv"
                n_genes = 0
                if aggregate_path.exists():
                    try:
                        n_genes = len(pd.read_csv(aggregate_path, sep="\t"))
                    except Exception:
                        pass

                results["full"][key] = {
                    "resolutions": full_resolutions,
                    "n_genes": n_genes,
                }

            # Filtered clustering (in filtered/ subdirectory)
            filtered_dir = cell_class_dir / "filtered"
            if filtered_dir.exists():
                filtered_resolutions = get_cluster_stats_for_dir(filtered_dir)
                if filtered_resolutions:
                    # Count filtered genes
                    filtered_aggregate = filtered_dir / "aggregate_cleaned.tsv"
                    n_filtered_genes = 0
                    if filtered_aggregate.exists():
                        try:
                            n_filtered_genes = len(pd.read_csv(filtered_aggregate, sep="\t"))
                        except Exception:
                            pass

                    results["filtered"][key] = {
                        "resolutions": filtered_resolutions,
                        "n_genes": n_filtered_genes,
                        "filter_params": {
                            "zscore_threshold": 0.3,
                            "fdr_threshold": 0.5,
                            "filter_mode": "both",
                        },
                    }

                    # Check for MozzareLLM results in each resolution
                    for res_dir in filtered_dir.iterdir():
                        if not res_dir.is_dir():
                            continue
                        mozzarellm_dir = res_dir / "mozzarellm"
                        if mozzarellm_dir.exists():
                            # Find results files
                            summary_files = list(mozzarellm_dir.glob("*_results_summaries.tsv"))
                            if summary_files:
                                try:
                                    summary_df = pd.read_csv(summary_files[0], sep="\t")
                                    model_name = summary_files[0].name.replace("_results_summaries.tsv", "")

                                    # Count by confidence
                                    confidence_counts = summary_df["pathway_confidence"].value_counts().to_dict()

                                    # Top processes
                                    top_processes = summary_df["dominant_process"].value_counts().head(5).to_dict()

                                    results["filtered"][key]["mozzarellm"] = {
                                        "resolution": float(res_dir.name),
                                        "model": model_name,
                                        "total_clusters": len(summary_df),
                                        "high_confidence": confidence_counts.get("High", 0),
                                        "medium_confidence": confidence_counts.get("Medium", 0),
                                        "low_confidence": confidence_counts.get("Low", 0),
                                        "established_genes": int(summary_df["num_established"].sum()),
                                        "novel_genes": int(summary_df["num_novel"].sum()),
                                        "uncharacterized_genes": int(summary_df["num_uncharacterized"].sum()),
                                        "top_processes": top_processes,
                                    }
                                except Exception:
                                    pass

    return results


def format_number(n):
    """Format number with commas."""
    if isinstance(n, float):
        if n == int(n):
            return f"{int(n):,}"
        return f"{n:,.2f}"
    return f"{n:,}"


def generate_markdown_report(stats: dict, output_dir: Path) -> str:
    """Generate a markdown report from the statistics."""
    lines = []
    lines.append("# Brieflow Output Assessment Report")
    lines.append("")
    lines.append(f"**Output Directory:** `{output_dir}`")
    lines.append("")

    # Preprocess
    lines.append("## 1. Preprocessing")
    lines.append("")
    preprocess = stats.get("preprocess", {})
    lines.append(f"- **SBS Tiles**: {format_number(preprocess.get('sbs_tiles', 0))}")
    lines.append(f"- **Phenotype Tiles**: {format_number(preprocess.get('phenotype_tiles', 0))}")
    lines.append(f"- **Total Tiles**: {format_number(preprocess.get('total_tiles', 0))}")
    lines.append("")

    # SBS
    lines.append("## 2. SBS (Sequencing by Synthesis)")
    lines.append("")
    sbs = stats.get("sbs", {})
    lines.append(f"- **Total Segmented Cells**: {format_number(sbs.get('total_segmented_cells', 0))}")
    lines.append("")

    mapping = sbs.get("mapping", {})
    if mapping:
        lines.append("### Barcode Mapping")
        lines.append("")
        lines.append(f"- **Cells with 1 barcode**: {format_number(mapping.get('cells_with_1_barcode', 0))} ({mapping.get('pct_with_1_barcode', 0):.1f}%)")
        lines.append(f"- **Cells with 1+ barcodes**: {format_number(mapping.get('cells_with_1_or_more_barcodes', 0))} ({mapping.get('pct_with_1_or_more_barcodes', 0):.1f}%)")
        lines.append(f"- **Cells with 1 gene**: {format_number(mapping.get('cells_with_1_gene', 0))} ({mapping.get('pct_with_1_gene', 0):.1f}%)")
        lines.append(f"- **Cells with 1+ genes**: {format_number(mapping.get('cells_with_1_or_more_genes', 0))} ({mapping.get('pct_with_1_or_more_genes', 0):.1f}%)")
        lines.append("")

    # Phenotype
    lines.append("## 3. Phenotype")
    lines.append("")
    phenotype = stats.get("phenotype", {})
    lines.append(f"- **Total Cells**: {format_number(phenotype.get('total_cells', 0))}")
    lines.append(f"- **Feature Count**: {format_number(phenotype.get('feature_count', 0))}")
    lines.append("")

    # Merge + Aggregate Pipeline (combined for clarity)
    lines.append("## 4. Merge & Aggregate Pipeline")
    lines.append("")
    merge = stats.get("merge", {})
    aggregate = stats.get("aggregate", {})
    pipeline = merge.get("pipeline", {})

    if pipeline:
        # Calculate key metrics
        sbs_input = pipeline.get('sbs_input_cells', 0)
        ph_input = pipeline.get('phenotype_input_cells', 0)
        fast_merge = pipeline.get('fast_merge', 0)
        merge_final = pipeline.get('merge_final', 0)
        filtered = pipeline.get('aggregate_filtered', 0)
        dedup = pipeline.get("deduplication", {})
        losses = pipeline.get("losses", {})

        lines.append("### Inputs")
        lines.append("")
        lines.append(f"- **SBS cells** (with 1+ barcodes): {format_number(sbs_input)}")
        lines.append(f"- **Phenotype cells**: {format_number(ph_input)}")
        lines.append("")

        lines.append("### Spatial Matching")
        lines.append("")
        lines.append(f"- **Candidate matches** (fast_merge): {format_number(fast_merge)}")
        lines.append("")

        lines.append("### Deduplication")
        lines.append("")
        if dedup.get("initial", 0) > 0:
            lines.append(f"- **Initial matches**: {format_number(dedup.get('initial', 0))} ({format_number(dedup.get('initial_with_genes', 0))} with genes)")
            lines.append(f"- **After SBS dedup** (1 phenotype per SBS cell): {format_number(dedup.get('after_sbs_dedup', 0))} ({format_number(dedup.get('after_sbs_dedup_with_genes', 0))} with genes)")
            lines.append(f"- **After phenotype dedup** (1 SBS per phenotype cell): {format_number(dedup.get('after_phenotype_dedup', 0))} ({format_number(dedup.get('after_phenotype_dedup_with_genes', 0))} with genes)")
        else:
            lines.append(f"- **Final matched cells**: {format_number(merge_final)}")
        if losses.get("deduplication"):
            lines.append(f"- *Removed by deduplication*: {format_number(losses.get('deduplication', 0))}")
        lines.append("")

        lines.append("### Aggregate Filtering")
        lines.append("")
        lines.append(f"- **Input to aggregate**: {format_number(merge_final)}")
        lines.append(f"- **After filtering** (requires gene + quality): {format_number(filtered)}")
        if losses.get("aggregate_filtering"):
            lines.append(f"- *Removed* (no gene / quality): {format_number(losses.get('aggregate_filtering', 0))}")
        lines.append("")

        # Retention summary
        if sbs_input > 0 and filtered > 0:
            overall_retention = 100 * filtered / sbs_input
            lines.append(f"**Overall Retention:** {format_number(filtered)} / {format_number(sbs_input)} = **{overall_retention:.1f}%**")
            lines.append("")

    # Aggregation results
    lines.append("### Aggregation Results")
    lines.append("")
    for name, astats in aggregate.items():
        if "error" in astats:
            lines.append(f"**{name}:** Error - {astats['error']}")
        else:
            lines.append(f"**Dataset:** {name}")
            lines.append("")
            lines.append(f"- **Total Cells**: {format_number(astats.get('total_aggregated_cells', 0))}")
            lines.append(f"- **Distinct Perturbations**: {format_number(astats.get('distinct_perturbations', 0))}")
            lines.append(f"- **Control Perturbations**: {format_number(astats.get('control_perturbations', 0))}")
            lines.append(f"- **Median Cells/Perturbation**: {format_number(astats.get('median_cells_per_perturbation', 0))}")
            lines.append(f"- **Cell Range/Perturbation**: {format_number(astats.get('min_cells_per_perturbation', 0))} - {format_number(astats.get('max_cells_per_perturbation', 0))}")
            lines.append(f"- **PC Features**: {format_number(astats.get('pc_features', 0))}")
        lines.append("")

    # Cluster
    lines.append("## 5. Clustering")
    lines.append("")
    cluster = stats.get("cluster", {})
    if "error" in cluster:
        lines.append(f"**Error:** {cluster['error']}")
    else:
        # Full clustering
        full = cluster.get("full", {})
        if full:
            lines.append("### Full Clustering (all perturbations)")
            lines.append("")
            for key, data in full.items():
                n_genes = data.get("n_genes", 0)
                lines.append(f"**{key}** ({format_number(n_genes)} genes)")
                lines.append("")
                resolutions = data.get("resolutions", [])
                for r in resolutions:
                    if "error" in r:
                        lines.append(f"- **Resolution {r.get('resolution', '?')}**: Error - {r['error']}")
                    else:
                        corum_fold = r.get('corum_fold_enrichment', None)
                        kegg_fold = r.get('kegg_fold_enrichment', None)
                        corum_str = f"{corum_fold:.1f}x" if isinstance(corum_fold, float) else "N/A"
                        kegg_str = f"{kegg_fold:.1f}x" if isinstance(kegg_fold, float) else "N/A"
                        lines.append(f"- **Resolution {r.get('resolution', '?')}**: {r.get('n_clusters', 0)} clusters, CORUM {r.get('real_corum_enriched', 0)} ({corum_str}), KEGG {r.get('real_kegg_enriched', 0)} ({kegg_str})")
                lines.append("")

                # Find best resolution
                valid_results = [r for r in resolutions if "error" not in r]
                if valid_results:
                    best = max(valid_results, key=lambda x: x.get("real_corum_enriched", 0) + x.get("real_kegg_enriched", 0))
                    lines.append(f"**Best:** Resolution {best.get('resolution')} with {best.get('real_corum_enriched', 0) + best.get('real_kegg_enriched', 0)} enriched clusters")
                    lines.append("")

        # Filtered clustering
        filtered = cluster.get("filtered", {})
        if filtered:
            lines.append("### Filtered Clustering (bootstrap-significant perturbations)")
            lines.append("")
            for key, data in filtered.items():
                n_genes = data.get("n_genes", 0)
                params = data.get("filter_params", {})
                lines.append(f"**{key}** ({format_number(n_genes)} genes)")
                lines.append("")
                lines.append(f"*Filter: z-score threshold={params.get('zscore_threshold', 'N/A')}, FDR threshold={params.get('fdr_threshold', 'N/A')}, mode={params.get('filter_mode', 'N/A')}*")
                lines.append("")
                resolutions = data.get("resolutions", [])
                for r in resolutions:
                    if "error" in r:
                        lines.append(f"- **Resolution {r.get('resolution', '?')}**: Error - {r['error']}")
                    else:
                        corum_fold = r.get('corum_fold_enrichment', None)
                        kegg_fold = r.get('kegg_fold_enrichment', None)
                        corum_str = f"{corum_fold:.1f}x" if isinstance(corum_fold, float) else "N/A"
                        kegg_str = f"{kegg_fold:.1f}x" if isinstance(kegg_fold, float) else "N/A"
                        string_prec = r.get('string_precision', 0)
                        lines.append(f"- **Resolution {r.get('resolution', '?')}**: {r.get('n_clusters', 0)} clusters, STRING {string_prec*100:.1f}% precision, CORUM {r.get('real_corum_enriched', 0)} ({corum_str}), KEGG {r.get('real_kegg_enriched', 0)} ({kegg_str})")
                lines.append("")

                # Find best resolution
                valid_results = [r for r in resolutions if "error" not in r]
                if valid_results:
                    best = max(valid_results, key=lambda x: x.get("real_corum_enriched", 0) + x.get("real_kegg_enriched", 0))
                    lines.append(f"**Best:** Resolution {best.get('resolution')} with {best.get('real_corum_enriched', 0) + best.get('real_kegg_enriched', 0)} enriched clusters")
                    lines.append("")

                # MozzareLLM results
                mozzarellm = data.get("mozzarellm", {})
                if mozzarellm:
                    lines.append("#### MozzareLLM Analysis")
                    lines.append("")
                    lines.append(f"*Model: {mozzarellm.get('model', 'unknown')}, Resolution: {mozzarellm.get('resolution', 'N/A')}*")
                    lines.append("")
                    lines.append(f"- **Total clusters analyzed**: {mozzarellm.get('total_clusters', 0)}")
                    lines.append(f"- **High confidence**: {mozzarellm.get('high_confidence', 0)}")
                    lines.append(f"- **Medium confidence**: {mozzarellm.get('medium_confidence', 0)}")
                    lines.append(f"- **Low confidence**: {mozzarellm.get('low_confidence', 0)}")
                    lines.append("")
                    lines.append(f"- **Established genes**: {format_number(mozzarellm.get('established_genes', 0))}")
                    lines.append(f"- **Novel gene roles**: {format_number(mozzarellm.get('novel_genes', 0))}")
                    lines.append(f"- **Uncharacterized genes**: {format_number(mozzarellm.get('uncharacterized_genes', 0))}")
                    lines.append("")

                    top_processes = mozzarellm.get("top_processes", {})
                    if top_processes:
                        lines.append("**Top biological processes:**")
                        lines.append("")
                        for process, count in list(top_processes.items())[:5]:
                            if process != "No coherent biological pathway":
                                lines.append(f"- {process} ({count} clusters)")
                        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Assess brieflow output and generate metrics report")
    parser.add_argument(
        "output_dir",
        type=str,
        help="Path to brieflow output directory",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON instead of markdown",
    )
    parser.add_argument(
        "--save",
        type=str,
        help="Save output to file (markdown or json based on --json flag)",
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    if not output_dir.exists():
        print(f"Error: Output directory does not exist: {output_dir}")
        return 1

    print(f"Assessing brieflow output at: {output_dir}")
    print()

    # Collect all statistics
    stats = {}

    print("Gathering preprocessing statistics...")
    stats["preprocess"] = get_preprocess_stats(output_dir)

    print("Gathering SBS statistics...")
    stats["sbs"] = get_sbs_stats(output_dir)

    print("Gathering phenotype statistics...")
    stats["phenotype"] = get_phenotype_stats(output_dir)

    print("Gathering merge statistics...")
    stats["merge"] = get_merge_stats(output_dir)

    print("Gathering aggregation statistics...")
    stats["aggregate"] = get_aggregate_stats(output_dir)

    print("Gathering clustering statistics...")
    stats["cluster"] = get_cluster_stats(output_dir)

    print()

    if args.json:
        output = json.dumps(stats, indent=2, default=str)
    else:
        output = generate_markdown_report(stats, output_dir)

    if args.save:
        with open(args.save, "w") as f:
            f.write(output)
        print(f"Report saved to: {args.save}")
    else:
        print(output)

    return 0


if __name__ == "__main__":
    exit(main())
