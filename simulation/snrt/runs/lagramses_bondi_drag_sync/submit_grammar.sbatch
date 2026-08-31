#!/bin/bash
#SBATCH --job-name=lagr_bondi_drag
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/runs/lagramses_bondi_drag_sync/slurm_%j.out
#SBATCH --error=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/runs/lagramses_bondi_drag_sync/slurm_%j.err
set -euo pipefail

root=/home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt/runs/lagramses_bondi_drag_sync
run="${root}/job_${SLURM_JOB_ID}"
exe=/home/kjhan/BACKUP/lagRamses/bin/ramses_final3d
checkpoint=/home/kjhan/BACKUP/lagRamses-DE/grow_bondi_validation_20260828/target_drag2_512_baseline_o1/output_00001

mkdir -p "${run}"
cp "${root}/bondi-drag-sync.nml" "${run}/bondi-drag-sync.nml"
ln -s "${checkpoint}" "${run}/output_00001"
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

cd "${run}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 1 -c 1 --cpu-bind=cores \
  "${exe}" bondi-drag-sync.nml > run.log
test -f output_00002/COMPLETE
grep -Eq 'NaN_CHK sink_post_bondi: uold=[[:space:]]*0 f=[[:space:]]*0 d0=[[:space:]]*0' run.log
test -s agn_coarse_state_v1.jsonl
test -n "$(find . -maxdepth 1 -type f -name 'sink_*.dat' -print -quit)"
