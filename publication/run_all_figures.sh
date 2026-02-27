#!/bin/bash
# Run all publication figure-generation scripts.
# Usage: bash publication/run_all_figures.sh
#
# Run from the whitney-analysis root directory on the cluster.
# All scripts use hardcoded /mnt/data/blainey/whitney-analysis/ paths.
#
# To switch to PDF output for Illustrator, change .png -> .pdf in the
# output filenames of each script's main(). matplotlib will automatically
# produce vector output when the filename ends in .pdf.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/publication/figures"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/LLM"

cd "$PROJECT_DIR"

echo "=== 1/7 Metrics report ==="
python3 "$SCRIPT_DIR/metrics.py" \
    "$PROJECT_DIR/analysis/brieflow_output" \
    --save "$SCRIPT_DIR/brieflow_output_report.md"
echo "  Done: brieflow_output_report.md"

echo ""
echo "=== 2/7 QC panels (a–f) ==="
python3 "$SCRIPT_DIR/generate_qc_panels.py" \
    --output-dir "$OUTPUT_DIR" \
    --sample-reads 0.1
echo "  Done: qc_panel_{a_heatmap,b_boxplot,c_barplot,d_prefix_matching,e_qscores,f_kde}"

echo ""
echo "=== 3/7 PHATE plot — all perturbations ==="
python3 "$SCRIPT_DIR/generate_phate_plot.py" \
    --output-dir "$OUTPUT_DIR"
echo "  Done: phate_map_controls"

echo ""
echo "=== 4/7 PHATE plot — significant perturbations ==="
python3 "$SCRIPT_DIR/generate_phate_plot_significant.py" \
    --output-dir "$OUTPUT_DIR"
echo "  Done: phate_map_significant"

echo ""
echo "=== 5/7 SBS montage ==="
python3 "$SCRIPT_DIR/create_sbs_montage.py"
echo "  Done: sbs_cycles_montage"

echo ""
echo "=== 6/7 Cluster gene table ==="
python3 "$SCRIPT_DIR/generate_cluster_table.py" \
    --output-dir "$OUTPUT_DIR/LLM"
echo "  Done: LLM/cluster_genes_table.tsv"

echo ""
echo "=== 7/7 Cluster cards ==="
python3 "$SCRIPT_DIR/cluster_card.py" \
    --clusters 10 221 \
    --output-dir "$OUTPUT_DIR/LLM"
echo "  Done: LLM/card_cluster_{0,221}.{png,pdf}"

echo ""
echo "=== All figures generated ==="
echo "Output: $OUTPUT_DIR"
echo ""
echo "NOTE: To produce PDF (vector) output for Illustrator, change"
echo "  .png -> .pdf in each script's savefig filenames."
echo "  matplotlib natively outputs vector when filename ends in .pdf."
