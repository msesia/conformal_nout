#!/bin/bash

module purge
module load gcc/11.3.0
module load openblas/0.3.21
eval "$(conda shell.bash hook)"
conda activate default

Rscript --vanilla exp_1_global_testing.R $1 $2 $3 $4 $5 $6
