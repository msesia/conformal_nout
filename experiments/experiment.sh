#!/bin/bash

module purge
module load gcc/11.3.0
module load openblas/0.3.21
eval "$(conda shell.bash hook)"
conda activate default

export OPENBLAS_NUM_THREADS=1

python3 experiment.py $1 $2 $3 $4 $5 $6 $7 $8 $9 "${10}" "${11}" "${12}" "${13}"
