#!/bin/bash
#SBATCH --job-name=sidm_zoom_ic_z63
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=02:00:00
#SBATCH --output=music_z63_%j.log
#SBATCH --error=music_z63_%j.err

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

exec /gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11_z63/generate_zoom_ic_z63.sh
