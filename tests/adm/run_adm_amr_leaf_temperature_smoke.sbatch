#!/bin/bash
#SBATCH --job-name=adm_amr_temp
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/adm_amr_temperature/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/adm_amr_temperature/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/adm_amr_temperature/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
mkdir -p "${snapshot}"

git -C "${code}" status --short >"${snapshot}/git_status.txt"
git -C "${code}" diff --binary >"${snapshot}/source.diff"
cp -p "${code}/tests/sidm/generate_amr_estimator_ic.py" "${snapshot}/"
cp -p "${code}/tests/adm/amr_leaf_temperature.nml.in" "${snapshot}/"
cp -p "${code}/tests/adm/analyze_amr_leaf_temperature.py" "${snapshot}/"
cp -p "${code}/tests/adm/run_adm_amr_leaf_temperature_smoke.sbatch" "${snapshot}/"

"${code}/bin/build_lageunha.sh"
cp -p "${code}/bin/build_lageunha.log" "${snapshot}/build_lageunha.log"
grep -Fq '=== make exit=0' "${snapshot}/build_lageunha.log"
sha256sum "${code}/bin/ramses_final3d" >"${snapshot}/binary.sha256"

for multiplicity in 1 2 4 8; do
  ic_directory=${root}/ic_${multiplicity}
  case_directory=${root}/multiplicity_${multiplicity}
  mkdir -p "${ic_directory}" "${case_directory}"
  python3 "${snapshot}/generate_amr_estimator_ic.py" "${ic_directory}" \
    --multiplicity "${multiplicity}" --coincident-copies \
    >"${case_directory}/ic_generation.log"
  sed "s|__IC_DIRECTORY__|${ic_directory}|" \
    "${snapshot}/amr_leaf_temperature.nml.in" >"${case_directory}/adm_amr_temperature.nml"
  : >"${case_directory}/jobcontrol.txt"
  cp -p "${code}/bin/ramses_final3d" "${case_directory}/ramses_final3d"
  (
    cd "${case_directory}"
    I_MPI_FABRICS=shm mpirun -np 1 \
      ./ramses_final3d adm_amr_temperature.nml >run.log 2>&1
  )
done

python3 "${snapshot}/analyze_amr_leaf_temperature.py" "${root}" \
  | tee "${root}/adm_amr_temperature_summary.txt"
