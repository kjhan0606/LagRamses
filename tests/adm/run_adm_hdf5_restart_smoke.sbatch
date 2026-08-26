#!/bin/bash
#SBATCH --job-name=adm_hdf5_restart
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/adm_hdf5_restart/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/adm_hdf5_restart/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/adm_hdf5_restart/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
ic_directory=${root}/ic
continuous=${root}/continuous
restart=${root}/restart
mkdir -p "${snapshot}" "${ic_directory}" "${continuous}" "${restart}"

git -C "${code}" status --short >"${snapshot}/git_status.txt"
git -C "${code}" diff --binary >"${snapshot}/source.diff"
cp -p "${code}/tests/sidm/generate_amr_estimator_ic.py" "${snapshot}/"
cp -p "${code}/tests/adm/adm_hdf5_restart.nml.in" "${snapshot}/"
cp -p "${code}/tests/adm/analyze_adm_hdf5_restart.py" "${snapshot}/"
cp -p "${code}/tests/adm/run_adm_hdf5_restart_smoke.sbatch" "${snapshot}/"

"${code}/bin/build_lageunha.sh"
cp -p "${code}/bin/build_lageunha.log" "${snapshot}/build_lageunha.log"
grep -Fq '=== make exit=0' "${snapshot}/build_lageunha.log"
sha256sum "${code}/bin/ramses_final3d" >"${snapshot}/binary.sha256"

python3 "${snapshot}/generate_amr_estimator_ic.py" "${ic_directory}" \
  --multiplicity 4 --coincident-copies >"${root}/ic_generation.log"

render_namelist() {
  local target=$1
  local nrestart=$2
  local nstepmax=$3
  local informat=$4
  sed -e "s|__IC_DIRECTORY__|${ic_directory}|" \
      -e "s|__NRESTART__|${nrestart}|" \
      -e "s|__NSTEPMAX__|${nstepmax}|" \
      -e "s|__INFORMAT__|${informat}|" \
      "${snapshot}/adm_hdf5_restart.nml.in" >"${target}"
}

render_namelist "${continuous}/adm_continuous.nml" 0 3 original
render_namelist "${restart}/adm_seed.nml" 0 2 original
render_namelist "${restart}/adm_restart.nml" 3 3 hdf5
: >"${continuous}/jobcontrol.txt"
: >"${restart}/jobcontrol.txt"
cp -p "${code}/bin/ramses_final3d" "${continuous}/ramses_final3d"
cp -p "${code}/bin/ramses_final3d" "${restart}/ramses_final3d"

(
  cd "${continuous}"
  I_MPI_FABRICS=shm mpirun -np 1 ./ramses_final3d adm_continuous.nml >run.log 2>&1
)
(
  cd "${restart}"
  I_MPI_FABRICS=shm mpirun -np 1 ./ramses_final3d adm_seed.nml >seed.log 2>&1
  I_MPI_FABRICS=shm mpirun -np 1 ./ramses_final3d adm_restart.nml >restart.log 2>&1
)

python3 "${snapshot}/analyze_adm_hdf5_restart.py" "${root}" \
  | tee "${root}/adm_hdf5_restart_summary.txt"
