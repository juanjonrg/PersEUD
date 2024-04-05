#!/bin/bash

set -o errexit
set -o xtrace

now=$(date +"%Y%m%d_%H%M%S")
slurm_command=""
slurm_command="srun -p cpu_zen2 --exclusive"

plan=5
plan_folder=~/research/Radiotherapy/plans/$plan

results_folder=results
result_path=$results_folder/x_$plan_$now.txt
mkdir -p $results_folder

fluence_folder=~/research/Radiotherapy/legacy/fluences
fluence_prefix=x_5_start_0001

export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=96
export OMP_NUM_THREADS=96
$slurm_command ./gradient_multicore $plan $plan_folder $result_path $fluence_folder $fluence_prefix
