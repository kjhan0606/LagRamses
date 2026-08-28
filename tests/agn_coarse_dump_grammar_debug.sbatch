#!/bin/bash
#SBATCH --job-name=agn_dump
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/agn_coarse_dump_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/agn_coarse_dump_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
: "${AGN_COARSE_TEST_EXE:?set AGN_COARSE_TEST_EXE to the linked test_agn_coarse_state executable}"
exe=${AGN_COARSE_TEST_EXE}
root=/home/kjhan/BACKUP/lagRamses-DE/agn_coarse_dump_${SLURM_JOB_ID}

mkdir -p "${root}"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

cd "${root}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 1 -c 1 --cpu-bind=cores \
  "${exe}" "${root}/agn_gate.jsonl" > run.log
grep -q 'AGN_COARSE_STATE_HARNESS_PASS' run.log
if grep -Eiq 'FATAL|SIG(SEGV|ABRT)' run.log; then
  echo 'AGN coarse-dump runtime failure' >&2
  exit 2
fi

python3 - "${root}/agn_gate.jsonl" <<'PY'
import json
import math
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
if not rows:
    raise SystemExit("AGN coarse-state ledger is empty")
if len(rows) != 4:
    raise SystemExit(f"expected two sinks at two coarse steps (4 rows), got {len(rows)}")

required = {
    "schema_version", "record_type", "nstep_coarse", "sink_id",
    "mass_msun", "gas_angular_momentum_code", "bh_spin_direction",
    "bh_spin_magnitude", "bh_angular_momentum_cgs",
    "bh_angular_momentum_norm_cgs", "spin_gas_angle_deg", "eddington_ratio",
    "coarse_eddington_ratio", "coarse_accreted_bh_mass_code",
    "bolometric_luminosity_erg_s", "disk_episode_mass_msun",
    "saved_feedback_energy_code",
}
seen = set()
for index, row in enumerate(rows, 1):
    missing = required.difference(row)
    if missing:
        raise SystemExit(f"row {index}: missing {sorted(missing)}")
    if row["schema_version"] != 1 or row["record_type"] != "agn_coarse_state":
        raise SystemExit(f"row {index}: bad schema or record type")
    key = (row["nstep_coarse"], row["sink_id"])
    if key in seen:
        raise SystemExit(f"row {index}: duplicate coarse-step/sink pair {key}")
    seen.add(key)
    for name in ("mass_msun", "eddington_ratio", "bolometric_luminosity_erg_s"):
        value = row[name]
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise SystemExit(f"row {index}: invalid {name}={value!r}")
    for name in ("gas_angular_momentum_code", "bh_spin_direction", "bh_angular_momentum_cgs"):
        if not isinstance(row[name], list) or len(row[name]) != 3:
            raise SystemExit(f"row {index}: invalid vector {name}")
    if row["bh_angular_momentum_norm_cgs"] <= 0:
        raise SystemExit(f"row {index}: non-positive physical BH angular momentum")

if {row["nstep_coarse"] for row in rows} != {7, 8}:
    raise SystemExit("coarse-step coverage or duplicate suppression failed")
by_sink = {row["sink_id"]: row for row in rows if row["nstep_coarse"] == 7}
if set(by_sink) != {101, 202}:
    raise SystemExit("global sink IDs were not preserved")
if not math.isclose(by_sink[101]["spin_gas_angle_deg"], 90.0, abs_tol=1e-10):
    raise SystemExit("orthogonal spin/gas angle is incorrect")
if not math.isclose(by_sink[202]["spin_gas_angle_deg"], 0.0, abs_tol=1e-10):
    raise SystemExit("aligned spin/gas angle is incorrect")
if by_sink[101]["feedback_mode"] != "THERMAL" or by_sink[202]["feedback_mode"] != "JET":
    raise SystemExit("feedback-mode classification is incorrect")
if not math.isclose(by_sink[101]["coarse_eddington_ratio"], 2.0):
    raise SystemExit("thermal coarse-step Eddington ratio is incorrect")
if not math.isclose(by_sink[202]["coarse_eddington_ratio"], 0.005):
    raise SystemExit("jet coarse-step Eddington ratio is incorrect")

print(f"AGN_COARSE_DUMP_PASS rows={len(rows)} sinks={len({r['sink_id'] for r in rows})}")
PY
