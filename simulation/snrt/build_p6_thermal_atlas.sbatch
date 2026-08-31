#!/bin/bash
#SBATCH --job-name=snrt_p6_grackle_atlas
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/build_p6_thermal_atlas_%j.log
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/build_p6_thermal_atlas_%j.err

set -euo pipefail
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt

.venv/bin/python tools/build_grackle_atlas.py \
  --generator /home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/grackle_cooling_grid \
  --grackle-data /home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5 \
  --scale-factors config/p6_thermal_atlas_scale_factors.txt \
  --work-directory data/p6_grackle_subtables \
  --output data/p6_thermal_atlas_a0p20849776.h5 \
  --threads "${SLURM_CPUS_PER_TASK}"
