#!/bin/bash
# Repeat only the small IC/reader gate using locally built, hashed binaries.
# Never overwrite the failed attempt, and never run directly on syntax.
set -euo pipefail
[[ $(hostname -s) == LagEunha ]] || { echo 'manual replay requires LagEunha' >&2; exit 2; }
bundle=${REPLAY_BUNDLE:?set immutable replay tooling}
reuse=${LM_REUSE_ROOT:?set the inspected build artifact directory}
tag=${LM_MANUAL_RUN_TAG:?set a fresh manual tag}
[[ ${tag} =~ ^manual_[a-zA-Z0-9_]+$ ]] || exit 2
root=/gpfs/kjhan/Hydro/Sidm/lm_gate/${tag}
[[ ! -e ${root} ]] || { echo 'refusing existing root' >&2; exit 2; }
(cd "${bundle}"; sha256sum --check SHA256SUMS)
sha256sum --check "${reuse}/binaries.sha256"
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1 OMP_STACKSIZE=128M OMP_PROC_BIND=close OMP_PLACES=cores
export I_MPI_FABRICS=shm
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:/home/kjhan/miniconda3/lib:${LD_LIBRARY_PATH:-}
mkdir -p "${root}"
cp -p "${reuse}/music_build/MUSIC" "${root}/MUSIC"
cp -p "${reuse}/reader_source/bin/ramses_final3d" "${root}/ramses"
printf 'reused_build_root=%s\n' "${reuse}" >"${root}/build_provenance.txt"
sha256sum "${root}/MUSIC" "${root}/ramses" >"${root}/binaries.sha256"
module -t list >"${root}/modules.txt" 2>&1
df -h /gpfs >"${root}/free_space.txt"
python3 "${bundle}/lagmusic_small_gate.py" prepare "${root}"
printf 'purpose=small_zoom_reader_initial_step root=%s noutput=1 aout=1.1 tout=1e100 foutput=1e9 fbackup=1e9 nstepmax=0 expected_dumps=0 expected_dump_bytes=0 estimated_IC_upper_MiB=128\n' "${root}"
for label in z6r1 z6r4 z5r1; do
  ranks=1
  [[ ${label} != z6r4 ]] || ranks=4
  (cd "${root}/${label}"; mpirun -np "${ranks}" "${root}/MUSIC" music.conf >music.log 2>&1)
done
python3 "${bundle}/lagmusic_small_gate.py" ic "${root}"
for label in z6r1 z5r1; do
  (cd "${root}/${label}"; mpirun -np 4 "${root}/ramses" run.nml >reader.log 2>&1)
done
python3 "${bundle}/lagmusic_small_gate.py" reader "${root}"
printf 'LAGMUSIC_SMALL_GATE=passed production_release=false\n'
