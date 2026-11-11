#!/bin/bash
#SBATCH --job-name=rar-g205-highk-ic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_highk_ic_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_highk_ic_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0
module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
base=/gpfs/kjhan/RAR_QUMOND
old="$base/initial_conditions/lcdm_g205_multifluid_l11_1286049"
transfer="$base/transfer_probes/lcdm_g205_z99_highk_1290351_1/transfer_z99.dat"
music=/home/kjhan/BackUP/rar/build/lagmusic_qumond_cpu_20261103/MUSIC
audit=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/g205_ic_scalar_audit.py
root="$base/initial_conditions/lcdm_g205_multifluid_l11_z99_highk_${SLURM_JOB_ID:?}"
test ! -e "$root"
test -f "$old/GENERATION_ONLY"
test -f "$old/GEOMETRY_AUDIT_COMPLETE"
test -f "$base/transfer_probes/lcdm_g205_z99_highk_1290351_1/COMPLETE"
printf '%s  %s\n' 687220cdb68d6ee9fafd4da3a23b1ef751dcb10d3ce621a183ee1fd66f9fe415 "$music" | sha256sum -c --status
sha256sum -c "$base/transfer_probes/lcdm_g205_z99_highk_1290351_1/artifacts.sha256" --status
printf '%s  %s\n' 08c9cb4fcdc49cd6b812d43822e80579aa48c476674ed76cf1ad9cd7611f267e "$old/parent.conf" | sha256sum -c --status
printf '%s  %s\n' fa347d060e30106d1f55ef34e1285bcf35c347af5775d79f2d84911ff03cd689 "$old/points.txt" | sha256sum -c --status
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 42949672960
mkdir "$root" "$root/work"
cp -p "$music" "$root/MUSIC"
cp -p "$transfer" "$root/transfer.dat"
cp -p "$old/parent.conf" "$root/parent.conf"
cp -p "$old/points.txt" "$root/points.txt"
cp -p "$audit" "$root/audit.py"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
sed \
  -e 's/zstart = 49/zstart = 99/' \
  -e 's/levelmax = 8/levelmax = 11\nregion = ellipsoid\nregion_point_level = 8/' \
  -e "/region_point_level = 8/a region_point_file = $root/points.txt" \
  -e 's/seed\[8\] = 20261103/seed[8] = 20261103\nseed[9] = 20261104\nseed[10] = 20261105\nseed[11] = 20261106/' \
  -e "s|^transfer_file = .*|transfer_file = $root/transfer.dat|" \
  -e "s|^filename = .*|filename = $root/ics|" "$root/parent.conf" > "$root/effective.conf"
baseline="$base/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084"
test -f "$baseline/COMPLETE"
sha256sum -c "$baseline/inputs.sha256" --status
sed "s|$baseline|$root|g" "$baseline/effective.conf" > "$root/expected_from_baseline.conf"
cmp "$root/effective.conf" "$root/expected_from_baseline.conf"
for key in 'boxlength = 25' 'zstart = 99' 'levelmin = 8' 'levelmax = 11' \
  'region = ellipsoid' 'region_point_level = 8' "region_point_file = $root/points.txt" \
  'baryons = yes' 'use_2LPT = no' 'Omega_m = 0.315' 'Omega_b = 0.0493017' \
  'seed[8] = 20261103' 'seed[9] = 20261104' 'seed[10] = 20261105' 'seed[11] = 20261106' \
  'format = grafic2' "transfer_file = $root/transfer.dat" "filename = $root/ics"; do
  grep -Fxq "$key" "$root/effective.conf"
done
grep -Fq 'Transfer functions at z = 99' "$root/transfer.dat"
printf 'purpose=z99 LCDM multifluid high-k matched-seed IC diagnostic; no RAMSES run\nexpected_RAMSES_dumps=0\nIC_work_budget_bytes=34359738368\nminimum_free_bytes=42949672960\nfree_bytes_before=%s\ntransfer=GR CAMB high-k sampling-200 at z99; solver accuracy and realized normalization unverified\n' "$free_bytes" > "$root/output_policy.txt"
sha256sum "$root/MUSIC" "$root/transfer.dat" "$root/parent.conf" "$root/points.txt" "$root/audit.py" "$root/effective.conf" "$root/job.sbatch" "$root/output_policy.txt" > "$root/inputs.sha256"
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/miniconda3/lib:${LD_LIBRARY_PATH:-}
(cd "$root/work" && srun --mpi=pmi2 --cpu-bind=cores "$root/MUSIC" "$root/effective.conf" > "$root/music.log" 2>&1)
if grep -Eiq 'ERROR:|FATAL|invalid conversion|conversion of .nan. failed|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$root/music.log"; then exit 1; fi
: > "$root/ic_fields.sha256"
for level in 008 009 010 011; do
  for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
    test -s "$root/ics/level_$level/$field"
    sha256sum "$root/ics/level_$level/$field" >> "$root/ic_fields.sha256"
  done
done
sha256sum -c "$root/ic_fields.sha256" --status
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 - "$root" <<'PY'
import struct, sys
from pathlib import Path
root = Path(sys.argv[1])
lines = []
for level in (8, 9, 10, 11):
    path = root / 'ics' / f'level_{level:03d}' / 'ic_deltab'
    with path.open('rb') as f:
        if struct.unpack('<i', f.read(4))[0] != 48:
            raise ValueError(f'invalid header: {path}')
        header = struct.unpack('<3i9f', f.read(48))
    if abs(header[7] - 0.01) > 1e-6:
        raise ValueError(f'unexpected initial scale factor level {level}: {header[7]}')
    lines.append(f'level={level} astart={header[7]:.10g}')
(root / 'header_epochs.txt').write_text('\n'.join(lines) + '\n')
PY
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$root/audit.py" "$root" "$root/scalar_summary.json"
sha256sum -c "$root/inputs.sha256" --status
sha256sum "$root/scalar_summary.json" "$root/header_epochs.txt" "$root/music.log" "$root/inputs.sha256" "$root/ic_fields.sha256" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'GENERATION_AND_SCALAR_AUDIT_ONLY_NOT_PRODUCTION' > "$root/COMPLETE"
