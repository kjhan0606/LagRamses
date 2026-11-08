#!/bin/bash
#SBATCH --job-name=rar-native-quadratic
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_quadratic_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_quadratic_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
out=$base/analysis/rar_native_quadratic_${SLURM_JOB_ID:?}
run=$base/runs/rar_particle_probe_1279276
test -f "$run/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp /home/kjhan/BackUP/rar/scripts/rar_native_quadratic_study.py "$out/study.py"
cp /home/kjhan/BackUP/rar/scripts/rar_native_particles.py "$out/rar_native_particles.py"
printf '%s\n' 'purpose=native discrete-operator reconstruction study, existing read-only geometry and particle sources' 'no simulation launch or snapshots; output budget=16MiB; fixed mass/position/cloud/AMR geometry' > "$out/execution_policy.txt"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
/home/kjhan/miniconda3/bin/python3 "$out/study.py" "$run" "$out/result.json"
(cd "$out" && sha256sum result.json study.py rar_native_particles.py job.sbatch >artifacts.sha256)
printf '%s\n' 'SPARSE_PROTOTYPE_MEASURED_NOT_COMPILED_IMPLEMENTATION' > "$out/COMPLETE"
