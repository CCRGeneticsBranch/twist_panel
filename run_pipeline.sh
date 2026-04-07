#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <run_dir> [snakemake args...]" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="$SCRIPT_DIR/workflow/Snakefile"
RUN_DIR="$1"
shift

RUN_DIR_ABS="$(cd "$RUN_DIR" && pwd)"
RUN_CONFIG="$RUN_DIR_ABS/config.yaml"

if [[ ! -f "$RUN_CONFIG" ]]; then
    echo "Missing run config: $RUN_CONFIG" >&2
    echo "Initialize one with: $SCRIPT_DIR/init_run.sh $RUN_DIR_ABS" >&2
    exit 1
fi

cd "$RUN_DIR_ABS"

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$RUN_CONFIG" \
    "$@"
