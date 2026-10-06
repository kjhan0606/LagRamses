#!/bin/bash
#SBATCH --job-name=chimes-cpu-build
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/cpu-build-%j.log
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17
task_root=/gpfs/kjhan/chimes-cvode-gpu-20261001
bash "$task_root/simulation/snrt/tools/chimes_cpu_build_test.sh" \
    "$task_root/.cvode-gpu/cpu-build-$SLURM_JOB_ID"
