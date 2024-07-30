#!/bin/bash

module purge
eval "$(conda shell.bash hook)"
conda activate default

Rscript --vanilla exp_1_global_testing.R $1 $2 $3 $4 $5 $6
