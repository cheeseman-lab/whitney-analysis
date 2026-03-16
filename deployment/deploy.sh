#!/bin/bash
# Deploy the Streamlit visualization app to Google Cloud Run with GCS-mounted data.
#
# Prerequisites:
#   - gcloud CLI authenticated with a project that has billing enabled
#   - Set GCP_PROJECT to your project ID
#
# Usage:
#   GCP_PROJECT=your-project-id bash deploy.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Configuration
GCP_PROJECT="${GCP_PROJECT:?Set GCP_PROJECT to your Google Cloud project ID}"
REGION="${REGION:-us-central1}"
SERVICE_NAME="${SERVICE_NAME:-whitney-viz}"
export BUCKET_NAME="${BUCKET_NAME:-${GCP_PROJECT}-whitney-viz-data}"
IMAGE_NAME="gcr.io/${GCP_PROJECT}/${SERVICE_NAME}"

echo "=== Configuration ==="
echo "Project:  $GCP_PROJECT"
echo "Region:   $REGION"
echo "Service:  $SERVICE_NAME"
echo "Bucket:   $BUCKET_NAME"
echo "Image:    $IMAGE_NAME"
echo ""

# Step 1: Copy brieflow source code needed by the app
echo "=== Step 1: Preparing code ==="
rm -rf "$SCRIPT_DIR/brieflow"
mkdir -p "$SCRIPT_DIR/brieflow/visualization"
mkdir -p "$SCRIPT_DIR/brieflow/workflow/lib/shared"

cp -r "$REPO_ROOT/brieflow/visualization/"* "$SCRIPT_DIR/brieflow/visualization/"
cp "$REPO_ROOT/brieflow/workflow/lib/shared/file_utils.py" "$SCRIPT_DIR/brieflow/workflow/lib/shared/"
touch "$SCRIPT_DIR/brieflow/workflow/__init__.py"
touch "$SCRIPT_DIR/brieflow/workflow/lib/__init__.py"
touch "$SCRIPT_DIR/brieflow/workflow/lib/shared/__init__.py"

# Step 2: Create GCS bucket and upload data
echo ""
echo "=== Step 2: Creating GCS bucket and uploading data ==="
gcloud storage buckets create "gs://${BUCKET_NAME}" \
  --project="$GCP_PROJECT" \
  --location="$REGION" \
  --uniform-bucket-level-access \
  2>/dev/null || echo "Bucket already exists"

echo "Uploading data to GCS (this may take a while)..."
bash "$SCRIPT_DIR/upload_data.sh"

# Step 3: Build container image with Cloud Build
echo ""
echo "=== Step 3: Building container image with Cloud Build ==="
gcloud builds submit "$SCRIPT_DIR" \
  --project="$GCP_PROJECT" \
  --tag="$IMAGE_NAME"

# Step 4: Deploy to Cloud Run with GCS volume mount
echo ""
echo "=== Step 4: Deploying to Cloud Run ==="
gcloud run deploy "$SERVICE_NAME" \
  --project="$GCP_PROJECT" \
  --region="$REGION" \
  --image="$IMAGE_NAME" \
  --platform=managed \
  --allow-unauthenticated \
  --port=8080 \
  --memory=4Gi \
  --cpu=2 \
  --min-instances=0 \
  --max-instances=3 \
  --execution-environment=gen2 \
  --add-volume=name=data-vol,type=cloud-storage,bucket="$BUCKET_NAME",readonly=true \
  --add-volume-mount=volume=data-vol,mount-path=/data \
  --set-env-vars="BRIEFLOW_OUTPUT_PATH=/data/brieflow_output/,CONFIG_PATH=/data/config/config.yml,SCREEN_PATH=/data/screen.yaml"

echo ""
echo "=== Deployment complete ==="
gcloud run services describe "$SERVICE_NAME" \
  --project="$GCP_PROJECT" \
  --region="$REGION" \
  --format="value(status.url)"
