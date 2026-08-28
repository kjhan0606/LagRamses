#!/bin/bash
#SBATCH --job-name=mg_omp_o3
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/mg_omp_o3_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/mg_omp_o3_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
run_dir=/home/kjhan/BACKUP/lagRamses-DE/mg_omp_o3_${SLURM_JOB_ID}
mkdir -p "${run_dir}"

cd "${repo}/bin"
make clean
make -j1 MG_PROFILE=1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_dir}/ramses_final3d.o3"
sha256sum "${run_dir}/ramses_final3d.o3"

ln -s "${paper}/output_00002" "${run_dir}/output_00002"
cp "${paper}/j02_species_boundary256_one.nml" "${run_dir}/cosmo.nml"
sed -i 's/nrestart=0/nrestart=2/; s/nstepmax=100000/nstepmax=66/; s/nremap=10/nremap=1000000/' \
  "${run_dir}/cosmo.nml"
sed -i '/use_fftw=.true./a timer_report_interval=1' "${run_dir}/cosmo.nml"
sed -i "/use_fftw=.true./a exchange_method='auto'" "${run_dir}/cosmo.nml"
sed -i 's/aout=0.07692307692307692/aout=2.0/; s/match_aout=.true./match_aout=.false./' \
  "${run_dir}/cosmo.nml"
sed -i 's/foutput=1000000/foutput=100000000/; s/fbackup=1000000/fbackup=100000000/' \
  "${run_dir}/cosmo.nml"
touch "${run_dir}/jobcontrol.txt"
cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 4 --cpu-bind=cores \
  ./ramses_final3d.o3 cosmo.nml | tee run.log

grep -q 'Run completed' run.log
grep -q 'MG_OMP_PROFILE' run.log
grep -q 'MGPROF coarse smoother' run.log
grep -q 'Dynamic exchange \[MG forward L8 cached\]: sparse-p2p' run.log
if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'MG O3 smoke diagnostic failure' >&2
  exit 2
fi
echo MG_OMP_O3_PASS
