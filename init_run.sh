#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <run_dir>" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_CONFIG="$SCRIPT_DIR/config/config.yaml"
TEMPLATE_SAMPLES="$SCRIPT_DIR/config/samples.tsv"
RUN_DIR="$1"

mkdir -p "$RUN_DIR"

RUN_CONFIG="$RUN_DIR/config.yaml"
RUN_SAMPLES="$RUN_DIR/samples.tsv"

if [[ -e "$RUN_CONFIG" || -e "$RUN_SAMPLES" ]]; then
    echo "Refusing to overwrite existing run files in $RUN_DIR" >&2
    echo "Expected new files: config.yaml and samples.tsv" >&2
    exit 1
fi

cp "$TEMPLATE_CONFIG" "$RUN_CONFIG"
cp "$TEMPLATE_SAMPLES" "$RUN_SAMPLES"

echo "Initialized run directory: $RUN_DIR"
echo "- Edit: $RUN_CONFIG"
echo "- Edit: $RUN_SAMPLES"
