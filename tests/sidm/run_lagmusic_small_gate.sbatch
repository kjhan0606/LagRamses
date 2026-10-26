#!/bin/bash
# Small generation/reader contract only. No production continuation or push.
# Submit an immutable GATE_BUNDLE containing SHA256SUMS and source overrides.
# Alternatively run with LM_MANUAL_RUN_TAG=manual_<tag> on LagEunha only.
# Manual mode needs no fabricated Slurm variables and uses one thread/rank.
#SBATCH --job-name=sidm_lagmusic_gate
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_gate_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_gate_%j.err
set -euo pipefail
bundle=${GATE_BUNDLE:?set frozen input bundle}
if [[ -n ${LM_MANUAL_RUN_TAG:-} ]]; then
  [[ $(hostname -s) == LagEunha ]] || { printf 'manual calculation permitted only on LagEunha\n' >&2; exit 2; }
  [[ ${LM_MANUAL_RUN_TAG} =~ ^manual_[a-zA-Z0-9_]+$ ]] || exit 2
  run_tag=${LM_MANUAL_RUN_TAG}
  thread_count=1
else
  : "${SLURM_JOB_ID:?submit through Slurm or explicitly select LagEunha manual mode}"
  run_tag=${SLURM_JOB_ID}
  thread_count=${SLURM_CPUS_PER_TASK:?}
fi
root=/gpfs/kjhan/Hydro/Sidm/lm_gate/${run_tag}
[[ ! -e ${root} ]] || { printf 'refusing existing root %s\n' "${root}" >&2; exit 2; }
(cd "${bundle}"; sha256sum --check SHA256SUMS)
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=${thread_count}
export OMP_STACKSIZE=128M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_FABRICS=shm
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:/home/kjhan/miniconda3/lib:${LD_LIBRARY_PATH:-}
mkdir -p "${root}/music_source" "${root}/reader_source"
module -t list >"${root}/modules.txt" 2>&1
df -h /gpfs >"${root}/free_space.txt"
git -C /home/kjhan/BACKUP/LagMUSIC/music archive \
  35a4eb978ae0e84673d86a0e4722421118277700 | tar -xf - -C "${root}/music_source"
cp -p "${bundle}/output_grafic2.cc" "${root}/music_source/src/plugins/"
git -C /home/kjhan/BACKUP/lagRamses-SIDM archive \
  eadfc9fdac3e10cac3bf3415c05ab67fb061ae96 | tar -xf - -C "${root}/reader_source"
for relative in patch/lagRamses/init_part.f90 patch/cuRamses/init_part.f90 \
  patch/lagRamses/init_time.kjhan.f90 patch/lagRamses/amr_step.jaehyun.f90 \
  patch/cuRamses/amr_step.jaehyun.f90 patch/lagRamses/newdt_fine.kjhan.f90 \
  patch/cuRamses/newdt_fine.kjhan.f90; do
  cp -p "${bundle}/${relative}" "${root}/reader_source/${relative}"
done
cmake -S "${root}/music_source" -B "${root}/music_build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx \
  -DMUSIC_ENABLE_MPI=ON -DFFTW3_ROOT=/home/kjhan/local \
  -DGSL_INCLUDE_DIR=/home/kjhan/miniconda3/include \
  -DGSL_LIBRARY=/home/kjhan/miniconda3/lib/libgsl.so \
  -DGSL_CBLAS_LIBRARY=/home/kjhan/miniconda3/lib/libgslcblas.so \
  -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=TRUE >"${root}/music_configure.log" 2>&1
cmake --build "${root}/music_build" -j2 >"${root}/music_build.log" 2>&1
make -C "${root}/reader_source/bin" -j1 HDF5=1 USE_FFTW=1 \
  GITHASH=eadfc9f-dmo-reader-lm-test-probe GITBRANCH=isolated-lagmusic-gate \
  GITREPO=git@github.com:kjhan0606/LagRamses.git >"${root}/reader_build.log" 2>&1
music=${root}/music_build/MUSIC
reader=${root}/reader_source/bin/ramses_final3d
sha256sum "${music}" "${reader}" >"${root}/binaries.sha256"
ldd "${music}" >"${root}/music_ldd.txt"
ldd "${reader}" >"${root}/reader_ldd.txt"
if grep -q 'not found' "${root}/music_ldd.txt" "${root}/reader_ldd.txt"; then exit 3; fi
python3 "${bundle}/audit_poisson_binary.py" "${reader}" \
  --require-fftw --json "${root}/binary_gate.json"
python3 "${bundle}/lagmusic_small_gate.py" prepare "${root}"
printf 'purpose=small_zoom_reader_initial_step root=%s noutput=1 aout=1.1 tout=1e100 foutput=1e9 fbackup=1e9 nstepmax=0 expected_dumps=0 expected_dump_bytes=0 estimated_IC_upper_MiB=128\n' "${root}"
for label in z6r1 z6r4 z5r1; do
  ranks=1
  [[ ${label} != z6r4 ]] || ranks=4
  (cd "${root}/${label}"; mpirun -np "${ranks}" "${music}" music.conf >music.log 2>&1)
done
python3 "${bundle}/lagmusic_small_gate.py" ic "${root}"
for label in z6r1 z5r1; do
  (cd "${root}/${label}"; mpirun -np 4 "${reader}" run.nml >reader.log 2>&1)
done
python3 "${bundle}/lagmusic_small_gate.py" reader "${root}"
printf 'LAGMUSIC_SMALL_GATE=passed production_release=false\n'
