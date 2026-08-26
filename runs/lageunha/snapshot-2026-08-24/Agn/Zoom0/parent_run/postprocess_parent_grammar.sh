#!/bin/bash
#SBATCH --job-name=sidm_parent_halos
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=postprocess_%j.log
#SBATCH --error=postprocess_%j.err

set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/parent_run
cd "${run_dir}"

printf 'Parent post-processing started on %s at %s\n' "$(hostname)" "$(date)"
if [ "${SKIP_HOP:-0}" != 1 ]; then
    ./find_parent_halos.sh
fi

position_file=$(find . -maxdepth 1 -type f -name 'grp_*.pos' -printf '%T@ %p\n' |
                sort -n | tail -1 | cut -d' ' -f2-)

if [ -z "${position_file}" ]; then
    printf 'HOP position catalogue was not produced\n' >&2
    exit 4
fi

python3 - "${position_file}" > zoom_candidates.txt <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
particle_mass = 1.0792e10
rows = []

for line in path.read_text().splitlines():
    fields = line.split()
    if not fields or fields[0] == "#":
        continue
    if len(fields) < 10:
        continue
    group = int(fields[0])
    npart = int(fields[1])
    x, y, z = map(float, fields[4:7])
    rows.append((npart, group, x, y, z))

rows.sort(reverse=True)
print("# rank group npart approximate_mass_h-1_Msun xc yc zc")
for rank, (npart, group, x, y, z) in enumerate(rows, start=1):
    mass = npart * particle_mass
    print(f"{rank:4d} {group:6d} {npart:9d} {mass:16.7e} "
          f"{x:12.8f} {y:12.8f} {z:12.8f}")
PY

printf 'Candidate summary: %s/zoom_candidates.txt\n' "${run_dir}"
printf 'Parent post-processing finished at %s\n' "$(date)"
