#!/bin/bash
#SBATCH --job-name=rar-ddo154-profile
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ddo154_profile_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ddo154_profile_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
obs="$base/observations/ddo154_littlethings_1290927"
sparc=/home/kjhan/BackUP/rar/data/sparc/Rotmod_LTG_official.zip
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/ddo154_hi_radial_probe.py
out="$base/analysis/ddo154_hi_radial_probe_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$obs/COMPLETE"
sha256sum -c "$obs/artifacts.sha256" --status
mkdir "$out"
cp -p "$source" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
sha256sum "$sparc" "$obs/hi_moment0.fits" > "$out/inputs.sha256"
export RAR_HI_IMAGE="$obs/hi_moment0.fits" RAR_SPARC_ZIP="$sparc" RAR_ANALYSIS_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1
sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/inputs.sha256" "$out/summary.json" "$out/radial_probe.csv" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'PRELIMINARY_IMAGE_PROBE_NOT_A_VALIDATED_GAS_MODEL' > "$out/COMPLETE"
