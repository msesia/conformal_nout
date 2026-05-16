#!/bin/bash

## -------------------------------------------------------------------------------------------------
## submit_exp_B2.sh — Grid launcher for the R experiment driver exp_B2_enumeration.sh.
##
## This script reproduces the enumeration simulation study reported in Appendix B2.
## Each result key below produces the inputs for one manuscript display item.
## Outputs are written to results/<RESULT_KEY>/ and consumed by the corresponding make_*.R
## script in figures_and_tables/. See README.md for the full reproducibility map.
##
##   Result key  Produces (manuscript display items)
##   ----------  ----------------------------------------------------------------
##   figA12      Figure A12
##
## Usage: ./submit_exp_B2.sh RESULT_KEY [--cluster|--local|--dry-run] [--force]
##   RESULT_KEY: one of the keys listed above (e.g. figA12)
##   --cluster (default): submit jobs via sbatch.
##   --local:             run jobs sequentially on this machine.
##   --dry-run:           print the commands without submitting / running.
##   --force:             re-run configurations even if their result files
##                        already exist (default: skip existing results).
##
## By default, the script skips any individual configuration whose output file
## already exists in results/<RESULT_KEY>/. This makes it safe to re-run after
## a partial failure: only missing configurations will be (re)submitted.
## Pass --force to ignore existing results and re-run everything.
## -------------------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------------------
## Produces results : figA12
## Output directory : results/<FIGURE>/
## -------------------------------------------------------------------------------------------------

# Parse arguments
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 RESULT_KEY [--cluster|--local|--dry-run] [--force]" >&2
  exit 1
fi

FIGURE=$1
MODE="cluster"
DRY_RUN=0
FORCE=0
shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cluster)  MODE="cluster" ;;
    --local)    MODE="local"   ;;
    --dry-run)  DRY_RUN=1      ;;
    --force)    FORCE=1        ;;
    *) echo "Unknown flag: $1 (use --cluster, --local, --dry-run, or --force)" >&2; exit 1 ;;
  esac
  shift
done


if [[ $FIGURE == "figA12" ]]; then
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of test sample sizes
  N_TEST_LIST=(200)
  # List of alternative distributions
  ALT_LIST=("uniform" "lehmann_k2" "beta_0.25_0.25" "beta_10_10" "normal_1.5_1" "normal_-1.5_1" "normal_0_0.25" "normal_0_2")
  # List of proportions of outlier
  PROP_OUT_LIST=(0 0.2 0.4 0.6)
  # Sequence of seeds for randomization
  SEED_LIST=$(seq 1 2)
  MEMO=5G

else
  echo "Error: unknown RESULT_KEY '$FIGURE'." >&2
  exit 1
fi


# Slurm parameters
TIME=00-00:20:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/"$FIGURE

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/"$FIGURE

# Counters for end-of-run summary
N_SUBMITTED=0
N_SKIPPED=0

# Loop over configurations
for SEED in $SEED_LIST; do
  for N_CAL in "${N_CAL_LIST[@]}"; do
    for N_TEST in "${N_TEST_LIST[@]}"; do
      for ALT in "${ALT_LIST[@]}"; do
        for PROP_OUT in "${PROP_OUT_LIST[@]}"; do

          # Generate a unique and interpretable file name based on the input parameters
          JOBN="${FIGURE}/n_cal_${N_CAL}_n_test_${N_TEST}_seed_${SEED}_alt_${ALT}_prop_out_${PROP_OUT}"
          OUT_FILE=$OUT_DIR"/"$JOBN".txt"

          # Skip if result already exists, unless --force was specified.
          if [[ -f $OUT_FILE && $FORCE -eq 0 ]]; then
            echo "Found existing result $OUT_FILE, skipping (use --force to re-run)."
            N_SKIPPED=$((N_SKIPPED + 1))
            continue
          fi

          # Script to be run
          SCRIPT="exp_B2_enumeration.sh $FIGURE $N_CAL $N_TEST $SEED $ALT $PROP_OUT"
          # Define job name for this configuration
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"

          if [[ $MODE == "cluster" ]]; then
            # Assemble slurm order for this job
            ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
            # Print order
            echo $ORD
            # Submit order to slurm scheduler (on cluster)
            if [[ $DRY_RUN -eq 0 ]]; then
              $ORD
            fi
          else
            # Run command now
            echo "./$SCRIPT"
            if [[ $DRY_RUN -eq 0 ]]; then
              ./$SCRIPT
            fi
          fi
          N_SUBMITTED=$((N_SUBMITTED + 1))

        done
      done
    done
  done
done

# --- Summary -----------------------------------------------------------------
echo ""
echo "------------------------------------------------------------"
if [[ $FORCE -eq 1 ]]; then
  echo "Force mode: all configurations were (re)submitted regardless of existing results."
fi
echo "Submitted/launched: $N_SUBMITTED job(s)"
echo "Skipped (already complete): $N_SKIPPED job(s)"
if [[ $N_SKIPPED -gt 0 && $FORCE -eq 0 ]]; then
  echo ""
  echo "Note: $N_SKIPPED configuration(s) were skipped because their results already exist."
  echo "      To force re-running them, re-invoke with --force."
fi

# --- Append machine-readable summary for run_all_experiments.sh ---
# When run standalone, RUN_SUMMARY_FILE is unset and this writes to /dev/null.
# When run by run_all_experiments.sh, it writes one tab-separated line per
# invocation to the shared summary file.
RUN_SUMMARY_FILE="${RUN_SUMMARY_FILE:-/dev/null}"
echo -e "${0##*/}\t${FIGURE}\t${N_SUBMITTED}\t${N_SKIPPED}" >> "$RUN_SUMMARY_FILE"
