#!/bin/bash

# Parameters
SETUP=0

if [[ $SETUP == 0 ]]; then
  DATA_LIST=("adversarial")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(1000)
  N_TEST_LIST=(1000)
  P_LIST=(100)
  A_LIST=(3)
  PURITY_LIST=(0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5)
  CLASSIFIER_LIST=("occ-auto" "bc-auto" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  SEED_LIST=$(seq 1 20)
  MEMO=5G

elif [[ $SETUP == 1 ]]; then
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
  SEED_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $SETUP == 2 ]]; then
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
  SEED_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $SETUP == 3 ]]; then
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
  SEED_LIST=$(seq 1 5)
  MEMO=5G

elif [[ $SETUP == 4 ]]; then
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
  SEED_LIST=$(seq 1 10)
  MEMO=1G

elif [[ $SETUP == 5 ]]; then
  DATA_LIST=("creditcard" "pendigits" "cover")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(10 100 1000)
  N_TEST_LIST=(20 40 60 80 100)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.0 0.1 0.2 0.5)
  CLASSIFIER_LIST=("occ-if" "auto")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SELECTION_LIST=("none")
  SEED_LIST=$(seq 1 5)
  MEMO=1G

elif [[ $SETUP == 6 ]]; then
  DATA_LIST=("creditcard")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(100)
  N_TEST_LIST=(20)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.0)
  CLASSIFIER_LIST=("occ-if")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.25)
  SEED_LIST=$(seq 1 1)
  MEMO=1G


elif [[ $SETUP == 100 ]]; then
  DATA_LIST=("lhco")
  N_TRAIN_LIST=(10000) # 100000
  N_CAL_LIST=(2000)
  N_TEST_LIST=(1000)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.0 0.02 0.05 0.1 0.15)
  CLASSIFIER_LIST=("auto" "bc-abc") # Do not run auto with 
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("none")
  SEED_LIST=$(seq 1 10)
  MEMO=5G


# elif [[ $SETUP == 1001 ]]; then
#   DATA_LIST=("circles-mixed")
#   N_TRAIN_LIST=(1000)
#   N_CAL_LIST=(2000)
#   N_TEST_LIST=(1000)
#   P_LIST=(1000)
#   A_LIST=(0.7)
#   PURITY_LIST=(0.0 0.2 0.5)
#   CLASSIFIER_LIST=("occ-svm") # Note: use fixed model, otherwise the selection may be inconsistent
#   ALPHA_LIST=(0.1)
#   TUNE_SIZE_LIST=(0.5)
#   SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "none") # Log scale
#   SEED_LIST=$(seq 1 10)
#   MEMO=5G

elif [[ $SETUP == 1001 ]]; then
  DATA_LIST=("circles-mixed")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(2000)
  N_TEST_LIST=(1000)
  P_LIST=(1000)
  A_LIST=(0.7)
  PURITY_LIST=(0.2)
  CLASSIFIER_LIST=("occ-svm") # Note: use fixed model, otherwise the selection may be inconsistent
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-50") # Log scale
  SEED_LIST=$(seq 1 1)
  MEMO=5G

elif [[ $SETUP == 1002 ]]; then
  DATA_LIST=("binomial")
  N_TRAIN_LIST=(1000)
  N_CAL_LIST=(2000)
  N_TEST_LIST=(1000)
  P_LIST=(10)
  A_LIST=(6.0)
  PURITY_LIST=(0 0.1 0.2 0.5)
  CLASSIFIER_LIST=("bc-mlp") # Note: use fixed model.
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "none") # Log scale
  SEED_LIST=$(seq 1 10)
  MEMO=5G

elif [[ $SETUP == 1005 ]]; then
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
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "none") # Log scale
  SEED_LIST=$(seq 1 5)
  MEMO=1G

elif [[ $SETUP == 1100 ]]; then
  DATA_LIST=("lhco")
  N_TRAIN_LIST=(10000 100000)
  N_CAL_LIST=(2000)
  N_TEST_LIST=(10000)
  P_LIST=(0)
  A_LIST=(0.0)
  PURITY_LIST=(0.05 0.1 0.15)
  CLASSIFIER_LIST=("bc-abc")
  ALPHA_LIST=(0.1)
  TUNE_SIZE_LIST=(0.5)
  SELECTION_LIST=("top-1" "top-2" "top-5" "top-10" "top-20" "top-50" "none") # Log scale
  SEED_LIST=$(seq 1 1)
  MEMO=5G

fi



# Slurm parameters
TIME=00-01:00:00                    # Time required (1 h)
CORE=1                              # Cores required (1)

# Assemble order prefix
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/setup"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/setup"$SETUP

# Loop over configurations and chromosomes
for SEED in $SEED_LIST; do
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

                        if [[ ( ($SETUP == 100) || ($SETUP == 1100) ) && ($N_TRAIN > 10000) && ($CLASSIFIER == "auto") ]]; then 
                           continue
                        fi

                        JOBN="setup"$SETUP"/"$DATA"_n"$N_TRAIN"_"$N_CAL"_"$N_TEST"_p"$P"_a"$A"_pt"$PURITY"_"$CLASSIFIER"_ts"$TUNE_SIZE"_alpha"$ALPHA"_"$SELECTION"_s"$SEED
                        OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                        COMPLETE=0
                        #ls $OUT_FILE
                        if [[ -f $OUT_FILE ]]; then
                          COMPLETE=1
                        fi

                        if [[ $COMPLETE -eq 0 ]]; then
                          # Script to be run
                          SCRIPT="experiment.sh $SETUP $DATA $N_TRAIN $N_CAL $N_TEST $P $A $PURITY $CLASSIFIER $TUNE_SIZE $ALPHA $SELECTION $SEED"
                          # Define job name for this chromosome
                          OUTF=$LOGS"/"$JOBN".out"
                          ERRF=$LOGS"/"$JOBN".err"
                          # Assemble slurm order for this job
                          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
                          # Print order
                          echo $ORD
                          # Submit order
                          $ORD
                          # Run command now
#                          ./$SCRIPT
                        fi
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
