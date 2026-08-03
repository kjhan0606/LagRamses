#!/bin/bash
#SBATCH --job-name=bounded_ckpt
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/bounded_checkpoint_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/bounded_checkpoint_%j.err

set -euo pipefail
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_dir=/home/kjhan/BACKUP/FDM/bounded_checkpoint_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp /gpfs/kjhan/Hydro/FDM/lbguard_test2/eng3200000_new/cosmo.nml "$run_dir/cosmo.nml"
ln -s /home/kjhan/BACKUP/FDM/hjm_omp_work/val3_deploy/output_00004 "$run_dir/output_00004"
ln -s /gpfs/kjhan/Hydro/FDM/lbguard_test2/noop_old/yield_table.asc "$run_dir/yield_table.asc"
# This compact checkpoint exercises the same Hilbert remap path under a much
# tighter physical grid capacity than the production run.  Its ownership is
# already nearly balanced, so it validates the remap but does not guarantee
# that the multi-round limiter is activated; production is the trigger case.
sed -i '/^nremap=/a lb_grid_headroom=1.0' "$run_dir/cosmo.nml"
: > "$run_dir/jobcontrol.txt"
cd "$run_dir"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  /home/kjhan/BACKUP/lagRamses/bin/ramses_final3d cosmo.nml | tee run.log
grep -q 'Run completed' run.log
