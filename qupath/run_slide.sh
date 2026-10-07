#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# UM-MAP, QuPath step for one whole-slide image:
#   1. build a single-image QuPath project next to the slide, "<slide>-proj/" (mkproj.groovy);
#   2. run classify_export.groovy on it, which writes <output_dir>/<image name>.tsv
#      (every character of the file name other than A-Z, a-z, 0-9 becomes "_"; columns
#      Class, Centroid X µm, Centroid Y µm);
#   3. delete the project folder.
#
# Usage:
#   qupath/run_slide.sh <slide> <output_dir> [classifier] [tissue_mode] [threshold]
#
#   slide        whole-slide image that QuPath can open (e.g. .svs)
#   output_dir   folder for the per-cell table (created if needed)
#   classifier   object classifier JSON file name in the classifier folder (default Combo-12-08-24.json)
#   tissue_mode  "full" for the whole image, or a pixel classifier JSON file name such as
#                TissueROI.json (default full)
#   threshold    watershed detection threshold (default 0.05)
#
#   Optional arguments that are omitted (or given as "") take the defaults defined in
#   classify_export.groovy, which are the UMICH2 settings. The TCGA settings are
#     qupath/run_slide.sh <slide> <output_dir> Tumor-Fibroblast-Lymphoid-Myeloid-2.json TissueROI.json 0.1
#
# Output: <output_dir>/<image name>.tsv, one row per cell (see step 2).
#
# Environment:
#   QUPATH                 QuPath executable (default: "QuPath" on the PATH)
#   UMMAP_CLASSIFIER_DIR   folder with the classifier JSON files (default: classifiers/ in this repository)
#
# Dependencies: bash; QuPath 0.4.3 (tested version).
#
# Notes:
#   - The folder that holds the slide must be writable, because QuPath creates the project
#     there. An existing "<slide>-proj/" folder is deleted and rebuilt, and an existing output
#     table for the slide is replaced.
#   - Do not run two jobs on the same slide file at the same time; they would share the
#     project folder.
#   - QuPath sizes its Java heap from the machine's total memory. On a shared node with a
#     smaller memory allocation, cap it with e.g. JAVA_TOOL_OPTIONS=-Xmx48g.

set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: qupath/run_slide.sh <slide> <output_dir> [classifier] [tissue_mode] [threshold]
  classifier   object classifier JSON file name (default Combo-12-08-24.json)
  tissue_mode  full, or a pixel classifier JSON file name such as TissueROI.json (default full)
  threshold    watershed detection threshold (default 0.05)
Environment: QUPATH (QuPath executable), UMMAP_CLASSIFIER_DIR (classifier folder)
EOF
}

if [ "$#" -lt 2 ] || [ "$#" -gt 5 ]; then
  usage
  exit 1
fi

SLIDE=$1
OUTPUT_DIR=$2
shift 2
# Remaining arguments (classifier, tissue_mode, threshold) are passed on unchanged.
SETTINGS=("$@")

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MK_PROJECT_SCRIPT="$SCRIPT_DIR/mkproj.groovy"
CLASSIFY_SCRIPT="$SCRIPT_DIR/classify_export.groovy"
QUPATH=${QUPATH:-QuPath}

if ! command -v "$QUPATH" >/dev/null 2>&1; then
  echo "Error: QuPath executable '$QUPATH' not found. Set QUPATH=/path/to/QuPath or put QuPath on the PATH." >&2
  exit 1
fi
if [ ! -f "$SLIDE" ]; then
  echo "Error: slide '$SLIDE' does not exist." >&2
  exit 1
fi

# mkproj.groovy creates the project next to the slide's canonical path (symbolic links resolved).
canonical_path() {
  if command -v realpath >/dev/null 2>&1; then
    realpath "$1"
  elif readlink -f / >/dev/null 2>&1; then
    readlink -f "$1"
  else
    python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
  fi
}
PROJECT_DIR="$(canonical_path "$SLIDE")-proj"
PROJECT_FILE="$PROJECT_DIR/project.qpproj"

# Output table name, as written by classify_export.groovy.
SLIDE_NAME=$(basename "$SLIDE")
TSV="${OUTPUT_DIR%/}/$(printf '%s' "$SLIDE_NAME" | LC_ALL=C sed 's/[^a-zA-Z0-9]/_/g').tsv"

cleanup() {
  if [ -d "$PROJECT_DIR" ]; then
    echo "[run_slide] Removing project folder $PROJECT_DIR"
    rm -rf "$PROJECT_DIR"
  fi
}

echo "[run_slide] $(date '+%Y-%m-%d %H:%M:%S') slide: $SLIDE"
echo "[run_slide] output folder: $OUTPUT_DIR"
if [ "${#SETTINGS[@]}" -gt 0 ]; then
  echo "[run_slide] settings passed on: $(printf "'%s' " "${SETTINGS[@]}")"
else
  echo "[run_slide] settings passed on: none (defaults of classify_export.groovy)"
fi
echo "[run_slide] QuPath: $(command -v "$QUPATH")"

if [ -d "$PROJECT_DIR" ]; then
  echo "[run_slide] Removing old project folder $PROJECT_DIR"
  rm -rf "$PROJECT_DIR"
fi
trap cleanup EXIT
mkdir -p "$OUTPUT_DIR"
# A table left by an earlier run is removed first, so a failed run cannot leave it looking new.
if [ -e "$TSV" ]; then
  echo "[run_slide] Removing earlier output $TSV"
  rm -f "$TSV"
fi

# 1. Single-image project.
echo "[run_slide] Creating QuPath project..."
"$QUPATH" script "$MK_PROJECT_SCRIPT" --args "$SLIDE"
if [ ! -f "$PROJECT_FILE" ]; then
  echo "Error: QuPath did not create $PROJECT_FILE" >&2
  exit 1
fi

# 2. Detection, classification and export.
QUPATH_ARGS=(--args "${OUTPUT_DIR%/}")
for value in "${SETTINGS[@]+"${SETTINGS[@]}"}"; do
  QUPATH_ARGS+=(--args "$value")
done
echo "[run_slide] Running classify_export.groovy..."
"$QUPATH" script -p "$PROJECT_FILE" "$CLASSIFY_SCRIPT" "${QUPATH_ARGS[@]}"

# QuPath 0.4.3 logs an error in a project script but still exits with status 0, so success is
# judged by the output table.
if [ ! -s "$TSV" ]; then
  echo "Error: $TSV was not written; see the QuPath messages above. (If the slide name has" \
       "non-ASCII characters, the table name may differ; check $OUTPUT_DIR.)" >&2
  exit 1
fi

# 3. The project folder is removed by the EXIT trap.
echo "[run_slide] $(date '+%Y-%m-%d %H:%M:%S') wrote $TSV ($(($(wc -l < "$TSV") - 1)) cells) in ${SECONDS} s"
