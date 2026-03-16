#!/bin/bash
# Upload only the required visualization data to GCS.
set -euo pipefail

SRC="/mnt/data/blainey/whitney-analysis/analysis"
BUCKET_NAME="${BUCKET_NAME:?Set BUCKET_NAME}"
DEST="gs://${BUCKET_NAME}"

echo "Uploading visualization data to $DEST ..."

# 1. Config files
echo "[1/6] Config files..."
gcloud storage cp "$SRC/screen.yaml" "$DEST/screen.yaml"
gcloud storage cp "$SRC/config/config.yml" "$DEST/config/config.yml"

# 2. QC eval files: */eval/**/*.{png,tsv}  (~12 MB)
echo "[2/6] Eval files..."
for phase_dir in "$SRC/brieflow_output"/*/eval; do
  if [ -d "$phase_dir" ]; then
    phase=$(basename "$(dirname "$phase_dir")")
    # Use cp with glob patterns instead of rsync
    gcloud storage cp -r "$phase_dir" "$DEST/brieflow_output/$phase/" 2>/dev/null || true
  fi
done

# 3. Cluster directory (~578 MB)
echo "[3/6] Cluster directory..."
gcloud storage rsync -r \
  "$SRC/brieflow_output/cluster/" "$DEST/brieflow_output/cluster/"

# 4. Aggregate feature TSVs (~1.8 GB)
echo "[4/6] Aggregate feature TSVs..."
gcloud storage cp "$SRC/brieflow_output/aggregate/tsvs/"*.tsv "$DEST/brieflow_output/aggregate/tsvs/"

# 5. Aggregate montages - PNGs and overlay TIFFs (skip montage_data TSVs)
echo "[5/6] Aggregate montage images..."
gcloud storage rsync -r \
  --exclude=".*montage_data.*" \
  "$SRC/brieflow_output/aggregate/montages/" "$DEST/brieflow_output/aggregate/montages/"

# 6. Aggregate eval files
echo "[6/6] Aggregate eval files..."
gcloud storage rsync -r \
  "$SRC/brieflow_output/aggregate/eval/" "$DEST/brieflow_output/aggregate/eval/"

echo ""
echo "Upload complete."
