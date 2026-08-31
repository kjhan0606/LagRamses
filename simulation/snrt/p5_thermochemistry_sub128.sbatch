#!/bin/bash
#SBATCH --job-name=p5_thermo_sub128
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p5_thermo_sub128_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p5_thermo_sub128_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/p5_run_thermochemical_pilot.py \
  --input data/p4_coeval_static_rt_input.h5 \
  --photon-metadata data/p4_pilot_agn_photon_ledger.json \
  --thermal-atlas data/p4_validation_thermal_atlas.h5 \
  --scale-factor 0.208497764676753 \
  --output data/p5_coeval_s8_sub128_6p37myr.h5 \
  --duration-myr 6.3697581 \
  --sn-order 8 \
  --reduced-light-fraction 0.01 \
  --courant 0.4 \
  --thermal-subcycles 128
