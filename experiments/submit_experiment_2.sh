#!/bin/bash

# Parameters
SETUP=2

if [[ $SETUP == 2 ]]; then
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of test sample sizes
  N_TEST_LIST=(200)
  # List of alternative distributions
  ALT_LIST=("uniform" "lehmann_k2" "beta_0.5_0.5" "beta_4_4" "normal_0.5_1" "normal_-0.5_1" "normal_0_0.5" "normal_0_1.5")
#  ALT_LIST=("uniform" "lehmann_k2")
  # List of proportions of outlier
#  PROP_OUT_LIST=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7)
  PROP_OUT_LIST=(0 0.2 0.4 0.6)
  # Sequence of seeds for randomization
  SEED_LIST=$(seq 1 1)
  MEMO=5G

fi

# Slurm parameters
TIME=00-00:20:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/setup"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/setup"$SETUP

# Loop over configurations
for SEED in $SEED_LIST; do
  for N_CAL in "${N_CAL_LIST[@]}"; do
    for N_TEST in "${N_TEST_LIST[@]}"; do
      for ALT in "${ALT_LIST[@]}"; do
        for PROP_OUT in "${PROP_OUT_LIST[@]}"; do

          # Generate a unique and interpretable file name based on the input parameters
          JOBN="setup${SETUP}/n_cal_${N_CAL}_n_test_${N_TEST}_seed_${SEED}_alt_${ALT}_prop_out_${PROP_OUT}"
          OUT_FILE=$OUT_DIR"/"$JOBN".txt"
          COMPLETE=0

          if [[ -f $OUT_FILE ]]; then
            COMPLETE=1
          fi

          if [[ $COMPLETE -eq 0 ]]; then
            # Script to be run
            SCRIPT="./exp_2_enumeration.sh $SETUP $N_CAL $N_TEST $SEED $ALT $PROP_OUT"
            # Define job name for this configuration
            OUTF=$LOGS"/"$JOBN".out"
            ERRF=$LOGS"/"$JOBN".err"
            # Assemble slurm order for this job
            ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
            # Print order
            echo $ORD
            # Submit order
            $ORD
            # Run command now
            #./$SCRIPT

          fi

        done
      done
    done
  done
done