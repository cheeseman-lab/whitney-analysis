#!/bin/bash
# Create a merged view: baked-in data + GCS-mounted montages via symlinks.
# The baked-in data has everything except aggregate/montages (too large).
# Montages are served from gcsfuse mount at /gcs_data (lazy-loaded on click).
set -e

# Create merged output directory mirroring baked-in data
mkdir -p /app/merged_data/brieflow_output/aggregate
cp -al /app/data/brieflow_output/* /app/merged_data/brieflow_output/ 2>/dev/null || true

# Link montages from GCS mount if available, otherwise from baked-in data
if [ -d "/gcs_data/brieflow_output/aggregate/montages" ]; then
    ln -sf /gcs_data/brieflow_output/aggregate/montages /app/merged_data/brieflow_output/aggregate/montages
    echo "Montages linked from GCS mount"
else
    echo "WARNING: GCS mount not available, montages will not be shown"
fi

exec streamlit run visualization/Experimental_Overview.py --server.port=8080 --server.address=0.0.0.0
