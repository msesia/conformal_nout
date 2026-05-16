#!/bin/bash

## -------------------------------------------------------------------------------------------------
## submit_exp_C1.sh — Wrapper for the R experiment driver exp_C1_time.R.
##
## This script reproduces the timing comparison reported in Appendix C1.
## Unlike submit_exp_A/B1/B2.sh, this is a single Rscript call rather than a
## parameter sweep; flags are accepted for compatibility with run_all_experiments.sh
## but the script always runs locally on the current machine.
##
##   Result key  Produces (manuscript display items)
##   ----------  ----------------------------------------------------------------
##   figA1       Figure A1
##
## Usage: ./submit_exp_C1.sh [RESULT_KEY] [--cluster|--local|--dry-run] [--force]
##   RESULT_KEY: figA1 (the only key produced; argument optional and ignored).
##   --cluster and --local are accepted but ignored; execution is always local.
##   --dry-run:           print the command without running it.
##   --force:             re-run even if the result file already exists
##                        (default: skip if results/figA1/time_results.csv exists).
##
## By default, if the result file already exists the script exits without
## running anything. Pass --force to re-run and overwrite the existing result.
## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------
## Produces results : figA1
## Output directory : results/<FIGURE>/
## -------------------------------------------------------------------------------------------------

DRY_RUN=0
FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cluster|--local)  ;;     # accepted, ignored
    --dry-run)          DRY_RUN=1 ;;
    --force)            FORCE=1   ;;
    -*) echo "Unknown flag: $1" >&2; exit 1 ;;
    *)  ;;                     # ignore any positional arg (e.g. RESULT_KEY name)
  esac
  shift
done

FIGURE="figA1"
OUT_DIR="results/$FIGURE"
OUT_FILE="$OUT_DIR/time_results.csv"
mkdir -p "$OUT_DIR"

# Counters for parent run_all_experiments.sh summary.
N_SUBMITTED=0
N_SKIPPED=0

if [[ -f "$OUT_FILE" && $FORCE -eq 0 ]]; then
  echo "Found existing result $OUT_FILE, skipping (use --force to re-run)."
  N_SKIPPED=1
else
  # Run the timing experiment that produces Figure A1.
  CMD="Rscript --vanilla exp_C1_time.R"
  echo "$CMD"
  if [[ $DRY_RUN -eq 0 ]]; then
    $CMD
  fi
  N_SUBMITTED=1
fi

# --- Append machine-readable summary for run_all_experiments.sh ---
# When run standalone, RUN_SUMMARY_FILE is unset and this writes to /dev/null.
# When run by run_all_experiments.sh, it writes one tab-separated line per
# invocation to the shared summary file.
RUN_SUMMARY_FILE="${RUN_SUMMARY_FILE:-/dev/null}"
echo -e "${0##*/}\t${FIGURE}\t${N_SUBMITTED}\t${N_SKIPPED}" >> "$RUN_SUMMARY_FILE"
