#!/bin/bash
#SBATCH --job-name=rar-ic-audit-reference
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/photo_audit_reference_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/photo_audit_reference_%j.err
set -Eeuo pipefail
umask 077

repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
base=/gpfs/kjhan/RAR_QUMOND
old=$base/initial_conditions/lcdm_g205_multifluid_l11_1286049
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
previous=$base/analysis
revision=${1:?Pass the committed full Git revision}
out=$base/analysis/lcdm_photo_audit_reference_${SLURM_JOB_ID:?}
test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test ! -e "$out"
test -f "$old/GEOMETRY_AUDIT_COMPLETE"
test -f "$parent/COMPLETE"
mkdir "$out" "$out/reference" "$out/reference/g205" "$out/reference/200" "$out/analysis"
trap 'touch "$out/FAILED"' ERR
awk '$2 ~ /^\.\/ics\/level_0(08|09|10|11)\//' "$old/artifacts.sha256" > "$out/old_inputs.sha256"
test "$(wc -l < "$out/old_inputs.sha256")" -eq 48
(cd "$old" && sha256sum -c "$out/old_inputs.sha256" --status)
awk -v p="$parent" '$2 == p "/lcdm_n256.conf" ||
    $2 == p "/lcdm_control_transfer_out_z49.dat" ||
    $2 ~ /\/ics_lcdm\/level_008\/ic_(deltab|posc[xyz]|vel[bc][xyz])$/ {print}' \
    "$parent/artifacts.sha256" > "$out/parent_inputs.sha256"
test "$(wc -l < "$out/parent_inputs.sha256")" -eq 12
sha256sum -c "$out/parent_inputs.sha256" --status
cp -p "$old/points.txt" "$out/reference/g205/points.txt"
cp -p "$old/effective.conf" "$out/reference/g205/music.conf"
cp -p "$old/music.log" "$out/reference/music.log"
ln -s "$old/ics" "$out/reference/ics"
ln -s "$old/transfer.dat" "$out/reference/200/transfer_z49.dat"
for tool in lcdm_photometric_l12_geometry lcdm_photometric_base_coherence; do
    git -C "$repo" show "$revision:provenance/tools/$tool.py" > "$out/$tool.py"
done
cp -p "$0" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=/home/kjhan/miniconda3/bin/python3
env -u LD_LIBRARY_PATH "$py" "$out/lcdm_photometric_l12_geometry.py" "$out/reference" 205 4536 > "$out/geometry.log"
env -u LD_LIBRARY_PATH "$py" "$out/lcdm_photometric_base_coherence.py" "$out/reference" "$parent" "$out/analysis" 205 > "$out/analysis.log"
"$py" - "$out" "$previous" <<'PY'
import json
import sys
from pathlib import Path

out, previous = map(Path, sys.argv[1:])
geometry = json.loads((out / "reference/geometry_audit.json").read_text())
old = json.loads((previous / "g205_ic_scalars_1290029/summary.json").read_text())
coherence = json.loads((out / "analysis/summary.json").read_text())
reference = json.loads((previous / "lcdm_g205_l11_basefield_coherence_1292418/summary.json").read_text())
for level, stats in old["levels"].items():
    actual = geometry["levels"][level]["fields"]
    assert actual["ic_deltab"]["min"] == stats["delta_min"]
    assert actual["ic_deltab"]["delta_le_minus_one_count"] == stats["delta_le_minus_one_count"]
    assert actual["ic_refmap"]["one_count"] == stats["refmap_one_count"]
for field, stats in reference["field_results"].items():
    actual = coherence["field_results"][field]
    for key in ("zoom_rms_after_mean_removal", "parent_rms_after_mean_removal",
                "difference_rms_over_parent_rms", "registered_pearson_correlation"):
        assert abs(actual[key] - stats[key]) < 1e-12, (field, key)
assert geometry["levels"]["11"]["fields"]["ic_deltab"]["delta_le_minus_one_count"] == 3
print("Reference geometry, density counts and matched-field statistics agree.")
PY
(cd "$old" && sha256sum -c "$out/old_inputs.sha256" --status)
sha256sum -c "$out/parent_inputs.sha256" --status
(cd "$out" && find . -type f ! -name artifacts.sha256 ! -name COMPLETE \
    -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256)
(cd "$out" && sha256sum -c artifacts.sha256 --status)
touch "$out/COMPLETE"
