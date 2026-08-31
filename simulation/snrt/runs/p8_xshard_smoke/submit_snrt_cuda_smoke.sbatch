#!/bin/bash
#SBATCH --job-name=snrt-p1-smoke
#SBATCH --partition=a100,h100,h200
#SBATCH --gpus=1
#SBATCH --time=00:05:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/runs/p8_xshard_smoke/snrt_cuda_smoke_%j.out
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/runs/p8_xshard_smoke/snrt_cuda_smoke_%j.err
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/13.0.2
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
/home/kjhan/BACKUP/lagRamses/bin/snrt_cuda_smoke
