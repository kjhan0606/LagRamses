#!/bin/bash
#SBATCH --job-name=adm_hpm_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/adm_hpm_smoke/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/adm_hpm_smoke/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/adm_hpm_smoke/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
ic_directory=${root}/ic
mkdir -p "${snapshot}" "${ic_directory}"

git -C "${code}" status --short >"${snapshot}/git_status.txt"
git -C "${code}" diff --binary >"${snapshot}/source.diff"
for file in \
  tests/adm/generate_adm_hpm_gradient_ic.py \
  tests/adm/adm_hpm_gradient.nml.in \
  tests/adm/analyze_adm_hpm_gradient.py \
  tests/adm/run_adm_hpm_gradient_smoke.sbatch; do
  cp -p "${code}/${file}" "${snapshot}/$(basename "${file}")"
done

"${code}/bin/build_lageunha.sh"
cp -p "${code}/bin/build_lageunha.log" "${snapshot}/build_lageunha.log"
grep -Fq '=== make exit=0' "${snapshot}/build_lageunha.log"
sha256sum "${code}/bin/ramses_final3d" >"${snapshot}/binary.sha256"
python3 "${snapshot}/generate_adm_hpm_gradient_ic.py" "${ic_directory}"

for mode in baseline hpm; do
  case_directory=${root}/${mode}
  mkdir -p "${case_directory}"
  if [[ "${mode}" == hpm ]]; then enabled=.true.; else enabled=.false.; fi
  sed -e "s|__IC_DIRECTORY__|${ic_directory}|" \
      -e "s|__ADM_HPM__|${enabled}|" \
      "${snapshot}/adm_hpm_gradient.nml.in" >"${case_directory}/adm_hpm_gradient.nml"
  : >"${case_directory}/jobcontrol.txt"
  cp -p "${code}/bin/ramses_final3d" "${case_directory}/ramses_final3d"
  (
    cd "${case_directory}"
    I_MPI_FABRICS=shm mpirun -np 1 ./ramses_final3d adm_hpm_gradient.nml </dev/null >run.log 2>&1
  )
done

python3 "${snapshot}/analyze_adm_hpm_gradient.py" "${root}" \
  | tee "${root}/adm_hpm_gradient_summary.txt"
