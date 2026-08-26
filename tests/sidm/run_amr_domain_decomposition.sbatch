#!/bin/bash
# Validate that separated occupied AMR leaves retain the same globally
# reduced SIDM pair accounting with one, two, and four MPI ranks.
#SBATCH --job-name=sidm_amr_domain
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/slurm_amr_domain_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/slurm_amr_domain_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/amr_domain/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
mkdir -p "${snapshot}"

git -C "${code}" status --short >"${snapshot}/git_status.txt"
git -C "${code}" diff --binary >"${snapshot}/source.diff"
cp -p "${code}/tests/sidm/amr_estimator.nml.in" "${snapshot}/"
cp -p "${code}/tests/sidm/generate_amr_estimator_ic.py" "${snapshot}/"
cp -p "${code}/tests/sidm/analyze_amr_domain_decomposition.py" "${snapshot}/"
cp -p "${code}/tests/sidm/run_amr_domain_decomposition.sbatch" "${snapshot}/"

"${code}/bin/build_lageunha.sh"
cp -p "${code}/bin/build_lageunha.log" "${snapshot}/build_lageunha.log"
grep -Fq '=== make exit=0' "${snapshot}/build_lageunha.log"
sha256sum "${code}/bin/ramses_final3d" >"${snapshot}/binary.sha256"

ic_directory=${root}/ic
python3 "${snapshot}/generate_amr_estimator_ic.py" "${ic_directory}" \
  --multiplicity 4 --patches 8 >"${root}/ic_generation.log"

for ranks in 1 2 4; do
  case_directory=${root}/mpi_${ranks}
  mkdir -p "${case_directory}"
  sed -e "s|__IC_DIRECTORY__|${ic_directory}|" -e 's|__LEVELMAX__|7|' \
    "${snapshot}/amr_estimator.nml.in" >"${case_directory}/amr_estimator.nml"
  : >"${case_directory}/jobcontrol.txt"
  cp -p "${code}/bin/ramses_final3d" "${case_directory}/ramses_final3d"
  (
    cd "${case_directory}"
    I_MPI_FABRICS=shm mpirun -np "${ranks}" \
      ./ramses_final3d amr_estimator.nml >run.log 2>&1
  )
done

python3 "${snapshot}/analyze_amr_domain_decomposition.py" "${root}" \
  | tee "${root}/amr_domain_summary.txt"
