#!/bin/bash
#SBATCH --job-name=p6_conservative_s8_picard
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=0-2
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p6_conservative_s8_picard_%A_%a.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p6_conservative_s8_picard_%A_%a.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

relaxations=(0.10 0.20 0.35)
relaxation=${relaxations[$SLURM_ARRAY_TASK_ID]}
tag="r${relaxation/./}"

.venv/bin/python tools/p6_run_conservative_thermochemical_pilot.py \
  --input data/p4_coeval_static_rt_input.h5 \
  --photon-metadata data/p4_pilot_agn_photon_ledger.json \
  --thermal-atlas data/p4_validation_thermal_atlas.h5 \
  --scale-factor 0.208365238497105 \
  --output "data/p6_coeval_conservative_primary_s8_${tag}_i128_6p37myr.h5" \
  --duration-myr 6.3697581 \
  --sn-order 8 \
  --reduced-light-fraction 0.01 \
  --courant 0.4 \
  --fixed-point-iterations 128 \
  --fixed-point-relaxation "$relaxation" \
  --fixed-point-tolerance 1.0e-4 \
  --ledger-tolerance 1.0e-3 \
  --thermal-implicit-iterations 24
