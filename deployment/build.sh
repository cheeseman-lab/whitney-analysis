#!/bin/bash
# Build and test the Docker image locally (requires Docker installed).
# For deployment, use deploy.sh instead (uses Cloud Build, no local Docker needed).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Copy brieflow source code needed by the app
echo "=== Preparing code ==="
rm -rf "$SCRIPT_DIR/brieflow"
mkdir -p "$SCRIPT_DIR/brieflow/visualization"
mkdir -p "$SCRIPT_DIR/brieflow/workflow/lib/shared"

cp -r "$REPO_ROOT/brieflow/visualization/"* "$SCRIPT_DIR/brieflow/visualization/"
cp "$REPO_ROOT/brieflow/workflow/lib/shared/file_utils.py" "$SCRIPT_DIR/brieflow/workflow/lib/shared/"
touch "$SCRIPT_DIR/brieflow/workflow/__init__.py"
touch "$SCRIPT_DIR/brieflow/workflow/lib/__init__.py"
touch "$SCRIPT_DIR/brieflow/workflow/lib/shared/__init__.py"

# Build Docker image
echo ""
echo "=== Building Docker image ==="
docker build -t whitney-viz "$SCRIPT_DIR"

echo ""
echo "=== Build complete ==="
echo "Test locally with:"
echo "  docker run -p 8080:8080 \\"
echo "    -e BRIEFLOW_OUTPUT_PATH=/data/brieflow_output/ \\"
echo "    -e CONFIG_PATH=/data/config/config.yml \\"
echo "    -e SCREEN_PATH=/data/screen.yaml \\"
echo "    -v /path/to/data:/data \\"
echo "    whitney-viz"
