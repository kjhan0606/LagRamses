#!/bin/bash
#SBATCH --job-name=validate_hhe_ledger
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/validate_hhe_ledger_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/validate_hhe_ledger_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/validate_conservative_primordial_hhe.py \
  --output data/validate_conservative_primordial_hhe.json \
  --fixed-point-iterations 16
