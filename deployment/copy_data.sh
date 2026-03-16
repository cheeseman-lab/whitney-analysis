#!/bin/bash
# Copy only the data files required by the Streamlit visualization app
# from the full brieflow_output directory into deployment/data/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC="/mnt/data/blainey/whitney-analysis/analysis"
DEST="$SCRIPT_DIR/data"

echo "Copying visualization data to $DEST ..."
rm -rf "$DEST"
mkdir -p "$DEST/brieflow_output"

# 1. Config files
cp "$SRC/screen.yaml" "$DEST/screen.yaml"
mkdir -p "$DEST/config"
cp "$SRC/config/config.yml" "$DEST/config/config.yml"

# 2. QC eval files: */eval/**/*.{png,tsv}  (~12 MB)
echo "Copying eval files..."
(cd "$SRC/brieflow_output" && find . -path "*/eval/*" \( -name "*.png" -o -name "*.tsv" \) -print0) | \
  while IFS= read -r -d '' f; do
    mkdir -p "$DEST/brieflow_output/$(dirname "$f")"
    cp "$SRC/brieflow_output/$f" "$DEST/brieflow_output/$f"
  done

# 3. Cluster directory (~578 MB) - all of it (TSVs, PNGs, JSONs, mozzarellm)
echo "Copying cluster directory..."
rsync -a "$SRC/brieflow_output/cluster/" "$DEST/brieflow_output/cluster/"

# 4. Aggregate feature TSVs (~1.8 GB)
echo "Copying aggregate feature TSVs..."
mkdir -p "$DEST/brieflow_output/aggregate/tsvs"
cp "$SRC/brieflow_output/aggregate/tsvs/"*.tsv "$DEST/brieflow_output/aggregate/tsvs/"

# 5. Aggregate montages - PNGs and overlay TIFFs (~255 MB)
echo "Copying aggregate montage images..."
rsync -a \
  --include='*/' \
  --include='*.png' \
  --include='overlay_montage.tiff' \
  --exclude='*' \
  "$SRC/brieflow_output/aggregate/montages/" "$DEST/brieflow_output/aggregate/montages/"

# 6. Aggregate eval files
echo "Copying aggregate eval files..."
rsync -a "$SRC/brieflow_output/aggregate/eval/" "$DEST/brieflow_output/aggregate/eval/"

echo ""
echo "Done. Total size:"
du -sh "$DEST"
