#!/bin/bash
#SBATCH --job-name=p4_hdf5_stage
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p4_hdf5_stage_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/p4_hdf5_stage_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/p4_stage_hdf5_level15.py \
  --snapshot ../pilot_rt_ledger/output_00017/data_00017.h5 \
  --info ../pilot_rt_ledger/output_00017/info_00017.txt \
  --zoom-manifest data/p4_high_density_manifest.json \
  --thermal-atlas data/p4_validation_thermal_atlas.h5 \
  --source-ledger data/p4_pilot_agn_photon_ledger.csv \
  --output data/p4_coeval_static_rt_input.h5 \
  --metadata-output data/p4_coeval_static_rt_input.json \
  --analysis-level 15 \
  --scan-chunk-grids 1000000
