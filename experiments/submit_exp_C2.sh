#!/bin/bash

## -------------------------------------------------------------------------------------------------
## submit_exp_C2.sh — Wrapper for the R experiment driver exp_C2_simes.R.
##
## This script reproduces the Simes-permutation comparison reported in Appendix C2.
## Unlike submit_exp_A/B1/B2.sh, this is a single Rscript call rather than a
## parameter sweep; flags are accepted for compatibility with run_all_experiments.sh
## but the script always runs locally on the current machine.
##
##   Result key  Produces (manuscript display items)
##   ----------  ----------------------------------------------------------------
##   tabA2       Table A2
##
## Usage: ./submit_exp_C2.sh [RESULT_KEY] [--cluster|--local|--dry-run] [--force]
##   RESULT_KEY: tabA2 (the only key produced; argument optional and ignored).
##   --cluster and --local are accepted but ignored; execution is always local.
##   --dry-run:           print the command without running it.
##   --force:             re-run even if the result file already exists
##                        (default: skip if results/tabA2/simes_perm_results.csv exists).
##
## By default, if the result file already exists the script exits without
## running anything. Pass --force to re-run and overwrite the existing result.
## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------
## Produces results : tabA2
## Output directory : results/<TABLE>/
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

TABLE="tabA2"
OUT_DIR="results/$TABLE"
OUT_FILE="$OUT_DIR/simes_perm_results.csv"
mkdir -p "$OUT_DIR"

# Counters for parent run_all_experiments.sh summary.
N_SUBMITTED=0
N_SKIPPED=0

if [[ -f "$OUT_FILE" && $FORCE -eq 0 ]]; then
  echo "Found existing result $OUT_FILE, skipping (use --force to re-run)."
  N_SKIPPED=1
else
  # Run the Simes-permutation experiment that produces Table A2.
  CMD="Rscript --vanilla exp_C2_simes.R"
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
echo -e "${0##*/}\t${TABLE}\t${N_SUBMITTED}\t${N_SKIPPED}" >> "$RUN_SUMMARY_FILE"
