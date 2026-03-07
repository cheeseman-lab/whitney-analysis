#!/bin/bash
# Run all publication figure-generation scripts.
# Usage: bash publication/run_all_figures.sh
#
# Run from the whitney-analysis root directory on the cluster.
# All scripts use hardcoded /mnt/data/blainey/whitney-analysis/ paths.
#
# All matplotlib scripts output .png (600 dpi), .pdf, and .svg.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/publication/figures"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/LLM"

cd "$PROJECT_DIR"

echo "=== 1/6 Metrics report ==="
python3 "$SCRIPT_DIR/metrics.py" \
    "$PROJECT_DIR/analysis/brieflow_output" \
    --save "$SCRIPT_DIR/brieflow_output_report.md"
echo "  Done: brieflow_output_report.md"

echo ""
echo "=== 2/6 QC panels (a–f) ==="
python3 "$SCRIPT_DIR/generate_qc_panels.py" \
    --output-dir "$OUTPUT_DIR" \
    --sample-reads 0.1
echo "  Done: qc_panel_{a_heatmap,b_boxplot,c_barplot,d_prefix_matching,e_qscores,f_kde}"

echo ""
echo "=== 3/6 PHATE plot — all perturbations ==="
python3 "$SCRIPT_DIR/generate_phate_plot.py" \
    --output-dir "$OUTPUT_DIR"
echo "  Done: phate_map_controls"

echo ""
echo "=== 4/6 PHATE plot — significant perturbations ==="
python3 "$SCRIPT_DIR/generate_phate_plot_significant.py" \
    --output-dir "$OUTPUT_DIR"
echo "  Done: phate_map_significant"

echo ""
echo "=== 5/6 Cluster gene table ==="
python3 "$SCRIPT_DIR/generate_cluster_table.py" \
    --output-dir "$OUTPUT_DIR/LLM"
echo "  Done: LLM/cluster_genes_table.tsv"

echo ""
echo "=== 6/6 Cluster cards ==="
python3 "$SCRIPT_DIR/cluster_card.py" \
    --clusters 1 2 5 6 8 13 \
    --output-dir "$OUTPUT_DIR/LLM"
echo "  Done: LLM/card_cluster_{1,2,5,6,8,13}.{png,pdf,svg}"

echo ""
echo "=== All figures generated ==="
echo "Output: $OUTPUT_DIR"
