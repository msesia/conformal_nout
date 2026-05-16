#!/bin/bash
# make_figures_tables.sh
# Read map_scripts.tsv in metadata/ and run each R script in order.
# Outputs go to figures_and_tables/. Per-script R output is redirected to
# figures_and_tables/logs/<script>.log.
set -u

TARGET_DIR="figures_and_tables"
TSV_REL="metadata/map_scripts.tsv"   # relative to the directory this script is run from
TSV_ABS="$(realpath "$TSV_REL" 2>/dev/null || readlink -f "$TSV_REL")"

if [[ ! -d "$TARGET_DIR" ]]; then
  echo "Error: $TARGET_DIR not found from $(pwd)" >&2
  exit 1
fi

if [[ ! -f "$TSV_ABS" ]]; then
  echo "Error: $TSV_REL not found from $(pwd)" >&2
  exit 1
fi

TOTAL=$(($(wc -l < "$TSV_ABS") - 1))
COUNT=0
FAILED=()

echo "Running $TOTAL scripts to make figures/tables from $TSV_REL"
echo "Working directory: $TARGET_DIR"
echo "Per-script logs:   $TARGET_DIR/logs/"
echo "=============================================="

cd "$TARGET_DIR" || exit 1
mkdir -p logs

while IFS=$'\t' read -r SCRIPT REQUIRED PRODUCES; do
  COUNT=$((COUNT + 1))
  LOG="logs/${SCRIPT%.R}.log"

  echo ""
  echo "[$COUNT/$TOTAL] $SCRIPT"
  echo "    requires results from : $REQUIRED"
  echo "    produces output       : $PRODUCES"
  echo "    log:                    $TARGET_DIR/$LOG"
  echo "----------------------------------------------"

  if [[ ! -f "$SCRIPT" ]]; then
    echo "  ! Script not found, skipping."
    FAILED+=("$SCRIPT (missing)")
    continue
  fi

  if Rscript "$SCRIPT" >"$LOG" 2>&1; then
    echo "  ✓ Done: $PRODUCES"
  else
    RC=$?
    echo "  ✗ Failed: $SCRIPT (exit $RC) — see $TARGET_DIR/$LOG"
    FAILED+=("$SCRIPT")
  fi
done < <(tail -n +2 "$TSV_ABS")

echo ""
echo "=============================================="
if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "All $TOTAL scripts completed successfully."
else
  echo "${#FAILED[@]} script(s) failed:"
  printf '  - %s\n' "${FAILED[@]}"
fi

echo ""
echo "Outputs:"
[[ -d figures ]] && echo "  figures: figures/  ($(ls figures 2>/dev/null | wc -l) files)"
[[ -d tables  ]] && echo "  tables:  tables/   ($(ls tables  2>/dev/null | wc -l) files)"
echo "  logs:    logs/"

[[ ${#FAILED[@]} -eq 0 ]] || exit 1
