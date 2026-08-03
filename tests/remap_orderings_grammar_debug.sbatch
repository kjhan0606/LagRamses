#!/bin/bash
#SBATCH --job-name=remap_ord
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --array=0-1
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/remap_ordering_%A_%a.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/remap_ordering_%A_%a.err

set -euo pipefail

# The supplied zoom IC is only geometrically valid for the ksection and
# angular decompositions.  Bisection exits during IC sub-volume validation and
# planar overflows its legacy cost map in the unmodified baseline, before the
# dynamic exchange layer can be assessed.
orderings=(ksection angular)
ordering=${orderings[$SLURM_ARRAY_TASK_ID]}

export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_dir=/home/kjhan/BACKUP/FDM/remap_${ordering}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$run_dir"
cp /home/kjhan/BACKUP/lagRamses/tests/dynamic_particle_cdm.nml "$run_dir/cosmo.nml"
sed -i "s/^ordering=.*/ordering='${ordering}'/" "$run_dir/cosmo.nml"
sed -i 's/^nstepmax=.*/nstepmax=1/' "$run_dir/cosmo.nml"
# Particle redistribution has its own 8-to-16-rank restart test.  Keep this
# matrix focused on AMR ordering/remap, since non-Hilbert ownership of this zoom
# IC is extremely uneven and otherwise turns it into an npartmax stress test.
sed -i 's/^pic=.*/pic=.false./' "$run_dir/cosmo.nml"
# The test isolates AMR remapping.  The base-grid FFTW solver assumes a slab
# decomposition that ksection/planar orderings do not provide and is a separate
# pre-existing limitation, so exclude Poisson from this ordering matrix.
sed -i 's/^poisson=.*/poisson=.false./' "$run_dir/cosmo.nml"
: > "$run_dir/jobcontrol.txt"

echo "Testing ordering=${ordering}"
cd "$run_dir"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  /home/kjhan/BACKUP/lagRamses/bin/ramses_final3d cosmo.nml | tee run.log
grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
