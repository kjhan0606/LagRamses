#!/bin/bash
#SBATCH --job-name=rar-g205-window-spectra
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_windowed_spectra_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_windowed_spectra_%j.err
set -Eeuo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
base=/gpfs/kjhan/RAR_QUMOND
baseline=$base/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084
extended=$base/initial_conditions/lcdm_g205_multifluid_l11_z99_highk_1290431
out=$base/analysis/lcdm_g205_z99_windowed_spectra_${SLURM_JOB_ID:?}
revision=${1:?Pass a full committed revision}
test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test ! -e "$out"
for root in "$baseline" "$extended"; do
    test -f "$root/COMPLETE"
    sha256sum -c "$root/artifacts.sha256" --status
    sha256sum -c "$root/inputs.sha256" --status
done
mkdir "$out"
trap 'touch "$out/FAILED"' ERR
for root in "$baseline" "$extended"; do
    name=$(basename "$root")
    awk '$2 ~ /\/ics\/level_011\/ic_(deltab|vel[bc][xyz]|refmap)$/ {print}' \
        "$root/ic_fields.sha256" > "$out/$name.sha256"
    test "$(wc -l < "$out/$name.sha256")" -eq 8
    sha256sum -c "$out/$name.sha256" --status
done
git -C "$repo" show "$revision:provenance/tools/lcdm_g205_z99_windowed_spectra.py" > "$out/analysis.py"
cp -p "$0" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 \
    "$out/analysis.py" "$baseline" "$extended" "$out" > "$out/analysis.log"
for root in "$baseline" "$extended"; do
    sha256sum -c "$out/$(basename "$root").sha256" --status
done
(cd "$out" && find . -maxdepth 1 -type f ! -name artifacts.sha256 ! -name COMPLETE \
    -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256)
(cd "$out" && sha256sum -c artifacts.sha256 --status)
touch "$out/COMPLETE"
