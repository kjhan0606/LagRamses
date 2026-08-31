#!/bin/bash
#SBATCH --job-name=p6_conservative_s8
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p6_conservative_s8_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p6_conservative_s8_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/p6_run_conservative_thermochemical_pilot.py \
  --input data/p4_coeval_static_rt_input.h5 \
  --photon-metadata data/p4_pilot_agn_photon_ledger.json \
  --thermal-atlas data/p4_validation_thermal_atlas.h5 \
  --scale-factor 0.208365238497105 \
  --output data/p6_coeval_conservative_primary_s8_6p37myr.h5 \
  --duration-myr 6.3697581 \
  --sn-order 8 \
  --reduced-light-fraction 0.01 \
  --courant 0.4 \
  --fixed-point-iterations 16 \
  --thermal-implicit-iterations 24
