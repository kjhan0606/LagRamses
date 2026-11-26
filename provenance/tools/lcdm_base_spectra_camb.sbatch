#!/bin/bash
#SBATCH --job-name=rar-base-camb-spectra
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/base_camb_spectra_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/base_camb_spectra_%j.err
set -Eeuo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
base=/gpfs/kjhan/RAR_QUMOND
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
zoom=$base/initial_conditions/lcdm_g205_multifluid_l11_1286049
revision=${1:?Pass a full committed revision}
out=$base/analysis/lcdm_base_camb_spectra_${SLURM_JOB_ID:?}
test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test -f "$parent/COMPLETE"
test -f "$zoom/GEOMETRY_AUDIT_COMPLETE"
test ! -e "$out"
mkdir "$out" "$out/parent" "$out/zoom"
trap 'touch "$out/FAILED"' ERR
awk '$2 ~ /\/ics_lcdm\/level_008\/ic_(deltab|posc[xyz]|vel[bc][xyz])$/ ||
    $2 ~ /\/lcdm_control_transfer_out_z49\.dat$/ {print}' \
    "$parent/artifacts.sha256" > "$out/parent_inputs.sha256"
test "$(wc -l < "$out/parent_inputs.sha256")" -eq 11
sha256sum -c "$out/parent_inputs.sha256" --status
awk '$2 ~ /^\.\/ics\/level_008\/ic_(deltab|posc[xyz]|vel[bc][xyz])$/ ||
    $2 == "./transfer.dat" {print}' "$zoom/artifacts.sha256" > "$out/zoom_inputs.sha256"
test "$(wc -l < "$out/zoom_inputs.sha256")" -eq 11
(cd "$zoom" && sha256sum -c "$out/zoom_inputs.sha256" --status)
git -C "$repo" show "$revision:provenance/tools/lcdm_base_spectra_camb.py" > "$out/analysis.py"
cp -p "$0" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=/home/kjhan/miniconda3/bin/python3
env -u LD_LIBRARY_PATH "$py" "$out/analysis.py" "$parent/ics_lcdm/level_008" \
    "$parent/lcdm_control_transfer_out_z49.dat" "$out/parent" > "$out/parent.log"
env -u LD_LIBRARY_PATH "$py" "$out/analysis.py" "$zoom/ics/level_008" \
    "$zoom/transfer.dat" "$out/zoom" > "$out/zoom.log"
sha256sum -c "$out/parent_inputs.sha256" --status
(cd "$zoom" && sha256sum -c "$out/zoom_inputs.sha256" --status)
(cd "$out" && find . -type f ! -name artifacts.sha256 ! -name COMPLETE \
    -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256)
(cd "$out" && sha256sum -c artifacts.sha256 --status)
touch "$out/COMPLETE"
