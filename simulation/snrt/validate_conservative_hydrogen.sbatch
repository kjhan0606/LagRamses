#!/bin/bash
#SBATCH --job-name=validate_h_stromgren
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/validate_h_stromgren_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/validate_h_stromgren_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/validate_conservative_hydrogen_stromgren.py \
  --output data/validate_conservative_hydrogen_stromgren.json \
  --fixed-point-iterations 12
