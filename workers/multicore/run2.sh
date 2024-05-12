#!/bin/bash

set -o errexit
set -o xtrace

now=$(date +"%Y%m%d_%H%M%S")
slurm_command=""
slurm_command="srun -p gpu_volta --exclusive"

plan=5
plan_folder=~/Radiotherapy/plans/$plan

results_folder=results
result_path=$results_folder/x_$plan_$now.txt
mkdir -p $results_folder

fluence_folder=~/Radiotherapy/zeromq-load-balancer/fluences/
fluence_prefix=x_5_start_0001

export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=32
export OMP_NUM_THREADS=32

$slurm_command ./gradient_multicore $plan $plan_folder $result_path $fluence_folder $fluence_prefix
#module load valgrind
#$slurm_command valgrind --leak-check=yes --error-limit=no ./gradient_multicore $plan $plan_folder $result_path $fluence_folder $fluence_prefix
