#!/bin/bash

## -------------------------------------------------------------------------------------------------
## submit_exp_A.sh — Grid launcher for the Python experiment driver exp_A.sh / exp_A.py.
##
## Each result key below produces the inputs for one or more manuscript display items.
## Outputs are written to results/<RESULT_KEY>/ and consumed by the corresponding make_*.R
## script in figures_and_tables/. See README.md for the full reproducibility map.
##
##   Result key  Produces (manuscript display items)
##   ----------  ----------------------------------------------------------------
##   fig1        Figure 1, Figure A13, Table A4
##   fig2        Figure 2, Figure A4, Figure A7
##   fig3        Figure 3, Figure A3
##   fig4        Figure 4
##   figA2       Figure A2
##   figA5       Figure A5, Figure A9
##   figA6       Figure A6, Figure A8
##   figA10      Figure A10
##   figA14      Figure A14
##   figA15      Figures A15, A16, A17, A18, A19, A20
##   figA21      Figures A21, A22
##
## Usage: ./submit_exp_A.sh RESULT_KEY [--cluster|--local|--dry-run] [--force]
##   RESULT_KEY: one of the keys listed above (e.g. fig1, fig2, figA15)
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
## Produces results : fig1, fig2, fig3, fig4, figA2, figA5, figA6, figA10, figA14, figA15, figA21
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


if [[ $FIGURE == "fig1" ]]; then
  DATA_LIST=("lhco")
  N_TRAIN_LIST=(10000) # 100000
  N_CAL_LIST=(2000)
  N_TEST_LIST=(10000)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.0 0.02 0.05 0.1 0.15)
  CLASSIFIER_LIST=("auto" "occ-if" "occ-svm" "occ-lof" "bc-mlp" "bc-rf" "bc-abc")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "fig2" ]]; then
  DATA_LIST=("circles-mixed")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(1000)
  A_LIST=(0.7)
  PURITY_LIST=(0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5)
  CLASSIFIER_LIST=("occ-auto" "bc-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "fig3" ]]; then
  DATA_LIST=("circles-mixed")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(2000)
  N_TEST_LIST=(1000)
  P_LIST=(1000)
  A_LIST=(0.7)
  PURITY_LIST=(0.0 0.2 0.5)
  CLASSIFIER_LIST=("occ-svm") # Note: use fixed model, otherwise the selection may be inconsistent
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "top-100") # Log scale
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "fig4" ]]; then
  DATA_LIST=("adversarial")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(100)
  A_LIST=(3.0)
  PURITY_LIST=(0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5)
  CLASSIFIER_LIST=("occ-auto" "bc-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA2" ]]; then
  DATA_LIST=("circles-mixed")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(1000)
  A_LIST=(0.7)
  PURITY_LIST=(0.1 0.2 0.5)
  CLASSIFIER_LIST=("occ-auto" "bc-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.01 0.02 0.05 0.1 0.2 0.5 0.75 0.99)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA5" ]]; then
  DATA_LIST=("binomial")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(10)
  A_LIST=(6.0)
  PURITY_LIST=(0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5)
  CLASSIFIER_LIST=("occ-auto" "bc-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA6" ]]; then
  DATA_LIST=("binomial")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(2000)
  N_TEST_LIST=(1000)
  P_LIST=(10)
  A_LIST=(6.0)
  PURITY_LIST=(0.0 0.1 0.2 0.5)
  CLASSIFIER_LIST=("bc-mlp") # Note: use fixed model.
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "none") # Log scale
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA10" ]]; then
  DATA_LIST=("mixture-0.0" "mixture-0.125" "mixture-0.25" "mixture-0.375" "mixture-0.5" "mixture-0.625" "mixture-0.75" "mixture-0.875" "mixture-1.0")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(100)
  A_LIST=(3.0)
  PURITY_LIST=(0.5)
  CLASSIFIER_LIST=("bc-auto" "occ-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA14" ]]; then
  DATA_LIST=("lhco")
  N_TRAIN_LIST=(10000) # 100000
  N_CAL_LIST=(2000)
  N_TEST_LIST=(2000)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.1 0.15 0.25)
  CLASSIFIER_LIST=("bc-abc")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "top-100") # Log scale
  BATCH_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $FIGURE == "figA15" ]]; then
  DATA_LIST=("creditcard" "pendigits" "cover" "shuttle" "mammography" "aloi")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(200)
  N_TEST_LIST=(100)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
  CLASSIFIER_LIST=("bc-auto" "occ-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("none")
  BATCH_LIST=$(seq 1 10)
  MEMO=1G

elif [[ $FIGURE == "figA21" ]]; then
  DATA_LIST=("creditcard" "pendigits" "cover" "shuttle" "mammography" "aloi")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(200)
  N_TEST_LIST=(100)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.5)
  CLASSIFIER_LIST=("occ-if") # Note: use fixed model.
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "top-100") # Log scale
  BATCH_LIST=$(seq 1 10)
  MEMO=1G

else
  echo "Error: unknown RESULT_KEY '$FIGURE'." >&2
  exit 1
fi


# Slurm parameters
TIME=00-02:00:00                    # Time required (2h)
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
for BATCH in $BATCH_LIST; do
  for DATA in "${DATA_LIST[@]}"; do
    for SELECTION in "${SELECTION_LIST[@]}"; do
      for TUNE_SIZE in "${TUNE_SIZE_LIST[@]}"; do
        for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
          for N_CAL in "${N_CAL_LIST[@]}"; do
            for N_TEST in "${N_TEST_LIST[@]}"; do
              for P in "${P_LIST[@]}"; do
                for A in "${A_LIST[@]}"; do
                  for PURITY in "${PURITY_LIST[@]}"; do
                    for CLASSIFIER in "${CLASSIFIER_LIST[@]}"; do
                      for ALPHA in "${ALPHA_LIST[@]}"; do

                        JOBN=$FIGURE"/"$DATA"_n"$N_TRAIN"_"$N_CAL"_"$N_TEST"_p"$P"_a"$A"_pt"$PURITY"_"$CLASSIFIER"_ts"$TUNE_SIZE"_alpha"$ALPHA"_"$SELECTION"_s"$BATCH
                        OUT_FILE=$OUT_DIR"/"$JOBN".txt"

                        # Skip if result already exists, unless --force was specified.
                        if [[ -f $OUT_FILE && $FORCE -eq 0 ]]; then
                          echo "Found existing result $OUT_FILE, skipping (use --force to re-run)."
                          N_SKIPPED=$((N_SKIPPED + 1))
                          continue
                        fi

                        # Script to be run
                        SCRIPT="exp_A.sh $FIGURE $DATA $N_TRAIN $N_CAL $N_TEST $P $A $PURITY $CLASSIFIER $TUNE_SIZE $ALPHA $SELECTION $BATCH"
                        # Define job name for this chromosome
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
            done
          done
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
