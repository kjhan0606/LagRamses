#!/bin/bash
#SBATCH --job-name=sidm_amr_occ
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/amr_occupancy/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/amr_occupancy/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/amr_occupancy/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
mkdir -p "$snapshot"

git -C "$code" status --short >"$snapshot/git_status.txt"
git -C "$code" diff --binary >"$snapshot/source.diff"
cp -p "$code/tests/sidm/amr_estimator.nml.in" "$snapshot/"
cp -p "$code/tests/sidm/generate_amr_estimator_ic.py" "$snapshot/"
cp -p "$code/tests/sidm/analyze_amr_occupancy.py" "$snapshot/"

"$code/bin/build_lageunha.sh"
cp -p "$code/bin/build_lageunha.log" "$snapshot/build_lageunha.log"
grep -Fq '=== make exit=0' "$snapshot/build_lageunha.log"
sha256sum "$code/bin/ramses_final3d" >"$snapshot/binary.sha256"

for multiplicity in 1 2 4 8; do
  ic_directory=${root}/ic_${multiplicity}
  case_directory=${root}/multiplicity_${multiplicity}
  mkdir -p "$ic_directory" "$case_directory"
  python3 "$snapshot/generate_amr_estimator_ic.py" "$ic_directory" \
    --multiplicity "$multiplicity" >"$case_directory/ic_generation.log"
  sed -e "s|__IC_DIRECTORY__|${ic_directory}|" -e 's|__LEVELMAX__|7|' \
    "$snapshot/amr_estimator.nml.in" >"$case_directory/amr_estimator.nml"
  : >"$case_directory/jobcontrol.txt"
  cp -p "$code/bin/ramses_final3d" "$case_directory/ramses_final3d"
  (
    cd "$case_directory"
    I_MPI_FABRICS=shm mpirun -np 1 \
      ./ramses_final3d amr_estimator.nml >run.log 2>&1
  )
done

python3 "$snapshot/analyze_amr_occupancy.py" "$root" \
  | tee "$root/amr_occupancy_summary.txt"
