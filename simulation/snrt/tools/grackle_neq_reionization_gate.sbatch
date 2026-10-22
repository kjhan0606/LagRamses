#!/usr/bin/env bash
# Verify the Grackle NEQ bridge honors the configured instantaneous z_reion gate.
# This is a tiny CPU-only library test; A40 is requested as the cluster allocation token.
#SBATCH --job-name=grackle-uv-gate
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-uv-gate-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
run=$root/.grackle-uv-gate-${SLURM_JOB_ID}
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
mkdir -p "$run"
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 OMP_STACKSIZE=64M OMP_PROC_BIND=close

mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_reionization_gate_test.c" \
  "$root/patch/lagRamses/grackle_neq_bridge.c" \
  -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/reionization_gate_test"
sha256sum "$run/reionization_gate_test" "$root/patch/lagRamses/grackle_neq_bridge.c" \
  "$uvdata" > "$run/assets.sha256"

for pair in '0.1 8.5' '0.1 -0.5' '0.1 100' '0.12 8.5' '0.12 -0.5'; do
  read -r a z_reion <<< "$pair"
  srun --ntasks=1 --cpus-per-task=1 "$run/reionization_gate_test" "$uvdata" "$a" "$z_reion" \
    | tee -a "$run/results.txt"
done

/home/kjhan/miniconda3/bin/python3 - "$run/results.txt" "$run/summary.json" <<'PY'
import json
import re
import sys

rows = {}
for line in open(sys.argv[1], encoding="utf-8"):
    m = re.search(r"RESULT a=([^ ]+) z_reion=([^ ]+) u=([^ ]+)((?: s\d+=[^ ]+)+)", line)
    if not m:
        continue
    a, z_reion, u, tail = m.groups()
    species = [float(x) for x in re.findall(r"s\d+=([^ ]+)", tail)]
    rows[(float(a), float(z_reion))] = (float(u), species)

if len(rows) != 5:
    raise SystemExit(f"expected five UV/redshift states, found {list(rows)}")

def difference(a, b):
    ua, sa = a
    ub, sb = b
    de = abs(ua - ub) / max(abs(ub), 1e-300)
    ds = max(abs(x - y) for x, y in zip(sa, sb)) / 1e-4
    return de, ds

pre_off = difference(rows[(0.1, 8.5)], rows[(0.1, -0.5)])
pre_on = difference(rows[(0.1, 100.0)], rows[(0.1, 8.5)])
post_on = difference(rows[(0.12, 8.5)], rows[(0.12, -0.5)])
summary = {
    "z_reion": 8.5,
    "a_before_reion": 0.1,
    "z_before_reion": 9.0,
    "a_after_reion": 0.12,
    "z_after_reion": 1 / 0.12 - 1,
    "pre_reion_config_vs_uv_suppressed_relative_energy": pre_off[0],
    "pre_reion_config_vs_uv_suppressed_max_species_per_rho": pre_off[1],
    "pre_reion_forced_uv_vs_config_relative_energy": pre_on[0],
    "pre_reion_forced_uv_vs_config_max_species_per_rho": pre_on[1],
    "post_reion_config_vs_uv_suppressed_relative_energy": post_on[0],
    "post_reion_config_vs_uv_suppressed_max_species_per_rho": post_on[1],
    "pass": (max(pre_off) < 1e-12 and max(pre_on) > 1e-8 and max(post_on) > 1e-8),
}
with open(sys.argv[2], "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
if not summary["pass"]:
    raise SystemExit("Grackle NEQ reionization gate failed")
PY

printf 'GRACKLE_NEQ_REIONIZATION_GATE_PASS %s\n' "$run"
