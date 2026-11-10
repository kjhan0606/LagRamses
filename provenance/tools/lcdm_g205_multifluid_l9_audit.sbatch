#!/bin/bash
#SBATCH --job-name=rar-g205-ic-audit
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_twofluid_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_twofluid_audit_%j.err
set -euo pipefail
project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
root=$base/initial_conditions/lcdm_g205_multifluid_l9_1284259
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
trace=$base/analysis/lcdm_z0_candidate_trace_1272106
out=$base/analysis/lcdm_g205_multifluid_l9_audit_${SLURM_JOB_ID:?}
test -f "$root/GENERATION_ONLY" && test -f "$parent/COMPLETE" && test -f "$trace/COMPLETE" && test ! -e "$out"
(cd "$root" && sha256sum -c artifacts.sha256 --status)
(cd "$parent" && sha256sum -c artifacts.sha256 --status)
(cd "$trace" && sha256sum -c artifacts.sha256 --status)
mkdir "$out"
cp "$project/scripts/lcdm_g205_multifluid_l9_audit.py" "$out/audit.py"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
srun env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$root" "$parent" "$trace" "$out" > "$out/audit.log" 2>&1
(cd "$out" && sha256sum audit.py job.sbatch geometry_audit.json audit.log > artifacts.sha256 && sha256sum -c artifacts.sha256)
touch "$out/GEOMETRY_ONLY_COMPLETE"
