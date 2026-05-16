#!/bin/bash

## -------------------------------------------------------------------------------------------------
## run_all_experiments.sh — Meta-launcher for the full reproducibility pipeline.
##
## This script reads map_submit.tsv and invokes each submit_exp_*.sh once per result key it
## produces. It is the single command that submits every experiment for every figure and table
## in the manuscript (main paper and appendices). For the result-key → display-item mapping,
## see README.md (also reproduced in the ACC form, Part 3).
##
## In aggregate, this script launches THOUSANDS of jobs. It is intended for cluster use.
## Do not run --local without --dry-run.
##
## Usage: ./run_all_experiments.sh [--cluster|--local] [--dry-run] [--force] [-y|--yes] [-h|--help]
##   Default mode is --cluster. See --help for full details.
## -------------------------------------------------------------------------------------------------

set -u

TSV="metadata/map_submit.tsv"

usage() {
  cat <<EOF
Usage: $0 [--cluster|--local] [--dry-run] [--force] [-y|--yes] [-h|--help]

Reads $TSV and invokes each submit_*.sh once per result it produces.

Options:
  --cluster   Submit jobs via sbatch (default).
  --local     Run jobs sequentially on this machine. NOT recommended:
              there are thousands of jobs. Use only with --dry-run for
              inspection, or if you know what you're doing.
  --dry-run   Echo the commands that would run, but don't execute them.
              Forwarded to each submit_*.sh.
  --force     Force re-running configurations even if their result files
              already exist (default: skip existing results).
              Forwarded to each submit_*.sh.
  -y, --yes   Skip the interactive confirmation prompt.
  -h, --help  Show this message and exit.

Examples:
  $0                          # cluster mode, asks for confirmation
  $0 -y                       # cluster mode, no prompt
  $0 --local --dry-run        # print every command, run nothing
  $0 --cluster --dry-run -y   # print every sbatch line, no submission
  $0 --force -y               # re-run all configurations, even completed ones

Prerequisites:
  - $TSV must exist in the current directory.
EOF
}

if [[ ! -f "$TSV" ]]; then
  echo "Error: $TSV not found in $(pwd). Run the metadata extraction script first." >&2
  echo "" >&2
  usage >&2
  exit 1
fi

# --- Parse arguments ---
MODE_FLAG="--cluster"
DRY_RUN_FLAG=""
FORCE_FLAG=""
YES=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cluster)   MODE_FLAG="--cluster" ;;
    --local)     MODE_FLAG="--local"   ;;
    --dry-run)   DRY_RUN_FLAG="--dry-run" ;;
    --force)     FORCE_FLAG="--force" ;;
    -y|--yes)    YES=1 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown flag: $1" >&2; echo "" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

# --- Warning ---
cat <<EOF
============================================================
  WARNING: run_all_experiments.sh
============================================================
  This script reads $TSV and submits every experiment for
  every result group listed there. The total number of
  individual jobs is in the THOUSANDS.

  This is intended to be run on a SLURM cluster. Running
  --local without --dry-run will execute jobs sequentially 
  on this machine and will take a very long time.

  Mode:    $MODE_FLAG
  Dry run: ${DRY_RUN_FLAG:-no}
  Force:   ${FORCE_FLAG:-no}
============================================================

EOF

if [[ $YES -eq 0 ]]; then
  read -r -p "Proceed? [y/N] " ans
  case "$ans" in
    y|Y|yes|YES) ;;
    *) echo "Aborted."; exit 1 ;;
  esac
fi

# --- Set up shared summary file for children to write to ---
# Each submit_exp_*.sh appends one tab-separated line per invocation:
#   script_name<TAB>result_key<TAB>n_submitted<TAB>n_skipped
export RUN_SUMMARY_FILE="$(mktemp)"
trap 'rm -f "$RUN_SUMMARY_FILE"' EXIT

# --- Count work ---
TOTAL_SCRIPTS=$(($(wc -l < "$TSV") - 1))
SCRIPT_COUNT=0
TOTAL_RESULTS=0
FAILED=()

# --- Main loop ---
while IFS=$'\t' read -r SCRIPT PRODUCES OUTPUT_DIR; do
  SCRIPT_COUNT=$((SCRIPT_COUNT + 1))
  echo ""
  echo "[$SCRIPT_COUNT/$TOTAL_SCRIPTS] $SCRIPT"
  echo "    produces:   $PRODUCES"
  echo "    output dir: $OUTPUT_DIR"
  echo "------------------------------------------------------------"

  if [[ ! -f "$SCRIPT" ]]; then
    echo "  ! Script not found, skipping."
    FAILED+=("$SCRIPT (missing)")
    continue
  fi

  # Split produces list on commas and run script once per result
  IFS=',' read -r -a RESULTS <<< "$PRODUCES"
  for r in "${RESULTS[@]}"; do
    r="${r## }"; r="${r%% }"   # trim surrounding spaces
    [[ -z "$r" ]] && continue
    TOTAL_RESULTS=$((TOTAL_RESULTS + 1))

    CMD="./$SCRIPT $r $MODE_FLAG${DRY_RUN_FLAG:+ $DRY_RUN_FLAG}${FORCE_FLAG:+ $FORCE_FLAG}"
    echo "  -> $CMD"
    if ! eval "$CMD" </dev/null; then
      echo "     ✗ Failed (exit $?)"
      FAILED+=("$SCRIPT $r")
    fi
  done
done < <(tail -n +2 "$TSV")

# --- Aggregate per-child summaries written to RUN_SUMMARY_FILE ---
TOTAL_JOBS_SUBMITTED=0
TOTAL_JOBS_SKIPPED=0
if [[ -s "$RUN_SUMMARY_FILE" ]]; then
  while IFS=$'\t' read -r CHILD FIG SUB SKP; do
    TOTAL_JOBS_SUBMITTED=$((TOTAL_JOBS_SUBMITTED + ${SUB:-0}))
    TOTAL_JOBS_SKIPPED=$((TOTAL_JOBS_SKIPPED + ${SKP:-0}))
  done < "$RUN_SUMMARY_FILE"
fi

# --- Summary ---
echo ""
echo "============================================================"
echo "Run-all summary"
echo "============================================================"
echo "Submit scripts invoked:                     $SCRIPT_COUNT"
echo "Result groups processed:                    $TOTAL_RESULTS"
echo "Individual jobs submitted/launched:         $TOTAL_JOBS_SUBMITTED"
echo "Individual jobs skipped (already complete): $TOTAL_JOBS_SKIPPED"

if [[ -n "$DRY_RUN_FLAG" ]]; then
  echo ""
  echo "(Dry run — no jobs were actually submitted or executed.)"
fi

if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo ""
  echo "${#FAILED[@]} invocation(s) failed:"
  printf '  - %s\n' "${FAILED[@]}"
  exit 1
else
  if [[ -z "$DRY_RUN_FLAG" ]]; then
    if [[ "$MODE_FLAG" == "--cluster" ]]; then
      echo ""
      echo "All invocations completed without error; jobs submitted to SLURM."
    else
      echo ""
      echo "All invocations completed without error."
    fi
  fi
fi
