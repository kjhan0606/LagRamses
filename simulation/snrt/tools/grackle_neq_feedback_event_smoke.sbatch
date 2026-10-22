#!/usr/bin/env bash
# Integration control: the previously exercised user-selected Kroupa + KL16 +
# LC18 + effective-SSP SNIa source runs through Grackle NEQ. Declared source
# approximations remain; this checks coupling/inventories, not yield calibration.
#SBATCH --job-name=grackle-fb-smoke
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-feedback-event-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
build=$root/.grackle-neq-production-build-415107
binary=$build/bin/ramses_grackle3d
run=$root/.grackle-neq-feedback-event-${SLURM_JOB_ID}
source_nml=$root/simulation/snrt/config/kl16_lc18_local_smoke.nml
history=/gpfs/kjhan/LRD_JWST/.agb-physical.4LAOTJ/snia-input/history.nml
yield_table=/gpfs/kjhan/LRD_JWST/.agb-physical.4LAOTJ/snia-input/yields.dat
snia_contract=$root/simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_noUVB.h5
test -x "$binary" && test -r "$source_nml" && test -r "$history" && test -r "$yield_table" \
  && test -r "$snia_contract" && test -r "$uvdata"
test ! -e "$run"
mkdir -p "$run"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2 OMP_STACKSIZE=64M OMP_PROC_BIND=close

cp "$source_nml" "$run/effective.nml"
sed -i \
  -e "s|cooling=.false.|cooling=.true., grackle_neq_data_file='${uvdata}', grackle_neq_h2_shielding=0|" \
  -e "s/cooling_method='original'/cooling_method='grackle_neq'/" \
  -e 's/gamma=1.6666666667d0/gamma=1.6666666666666667d0/' \
  -e "s|high_mass_history_path='CHANGE_ME_KL16_INPUT/history.nml'|high_mass_history_path='${history}'|" \
  -e "s/population_model='single_star_ssp'/population_model='binary_ssp'/" \
  -e 's/binary_fraction=0d0/binary_fraction=0.5d0/' \
  -e 's/use_snia=.false./use_snia=.true./' \
  "$run/effective.nml"
sha256sum "$binary" "$source_nml" "$history" "$yield_table" "$snia_contract" "$uvdata" \
  "$root/simulation/snrt/tools/grackle_neq_feedback_event_smoke.sbatch" \
  > "$run/assets.sha256"

cd "$run"
TIMEFORMAT='wall=%3R user=%3U sys=%3S'
{ time timeout 1800s env PHASE0_YIELD_TABLE="$yield_table" \
    PHASE0_SNIA_RUNTIME_CONTRACT="$snia_contract" \
    srun --ntasks=1 --cpus-per-task=2 "$binary" effective.nml; } \
  > run.log 2> wall.txt
grep -q 'Run completed' run.log
grep -q 'Entering thermal_feedback' run.log
if grep -E 'ERROR:|Aborting|Too many errors|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' run.log; then
  exit 1
fi

/home/kjhan/miniconda3/bin/python3 - "$run" <<'PY'
import glob
import json
import os
import sys

import h5py
import numpy as np

root = sys.argv[1]
outputs = sorted(glob.glob(os.path.join(root, "output_*")))
if len(outputs) < 2:
    raise SystemExit(f"expected multiple coarse-step outputs, got {outputs}")

def read_output(path):
    files = glob.glob(os.path.join(path, "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"expected one rank HDF5 file in {path}, got {files}")
    with h5py.File(files[0], "r") as f:
        values = []
        for group in f["hydro"]:
            level = int(group.rsplit("_", 1)[-1])
            amr = f[f"amr/{group}"]
            leaf = amr["son_flag"][:] == 0
            u = np.column_stack([f[f"hydro/{group}/uold_{i}"][:] for i in range(1, 28)])
            values.append(u[leaf])
        gas = np.concatenate(values)
        if "particles" not in f:
            raise SystemExit("HDF5 output lacks the RAMSES particles group")
        p = f["particles"]
        ptype = p["ptypep"][:].astype(np.int64)
        identity = p["identity"][:].astype(np.int64)
        mass = p["mass"][:]
        if ptype.shape != identity.shape or mass.shape != identity.shape:
            raise SystemExit("particle type/identity/mass arrays have inconsistent shapes")
        is_star = ptype == 1  # PTYPE_STAR in pm_commons
        stellar_masses = {int(i): float(m) for i, m in zip(identity[is_star], mass[is_star])}
        return files[0], gas, len(identity), stellar_masses

states = []
particle_masses = []
for out in outputs:
    file, u, particle_count, stellar_masses = read_output(out)
    species = u[:, 18:27]
    h = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
    he = species[:, 3] + species[:, 4] + species[:, 5]
    charge = species[:, 1] - species[:, 2] + species[:, 7]/2 + species[:, 4]/4 + species[:, 5]/2
    rho = u[:, 0]
    eint = u[:, 4] - np.sum(u[:, 1:4]**2, axis=1)/(2*rho)
    metrics = {
        "output": os.path.basename(out), "hdf5": file, "cells": len(u),
        "finite": bool(np.isfinite(u).all()), "min_density": float(rho.min()),
        "min_species": float(species.min()), "min_internal_energy_density": float(eint.min()),
        "max_H_residual": float(np.max(np.abs(h-u[:, 6])/rho)),
        "max_He_residual": float(np.max(np.abs(he-u[:, 7])/rho)),
        "max_charge_residual": float(np.max(np.abs(species[:, 8]-charge)/rho)),
        "particle_count": particle_count,
        "stellar_particle_count": len(stellar_masses),
    }
    if (not metrics["finite"] or metrics["min_density"] <= 0 or metrics["min_species"] < 0
            or metrics["min_internal_energy_density"] <= 0
            or max(metrics["max_H_residual"], metrics["max_He_residual"],
                   metrics["max_charge_residual"]) >= 1e-6):
        raise SystemExit(f"Grackle species/energy budget failed: {metrics}")
    states.append(metrics)
    particle_masses.append(stellar_masses)

mass_drops = []
for before, after in zip(particle_masses, particle_masses[1:]):
    for identity in before.keys() & after.keys():
        drop = before[identity] - after[identity]
        if drop > max(abs(before[identity])*1e-12, 1e-30):
            mass_drops.append({"identity": identity, "mass_drop_code": drop})
if not mass_drops:
        raise SystemExit("no persistent PTYPE_STAR mass return was resolved in output snapshots")

summary = {
    "scope": "user-selected Kroupa + KL16 + LC18 + effective-SSP SNIa / Grackle integration control; explicit source approximations remain",
    "outputs": states,
    "persistent_stellar_particle_mass_drops": len(mass_drops),
    "largest_stellar_particle_mass_drop_code": max(x["mass_drop_code"] for x in mass_drops),
    "feedback_event_and_chemistry_gate": True,
}
with open(os.path.join(root, "feedback_event_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
PY

# Delete only evaluated raw snapshots; retain the inputs, logs and compact evidence.
for output in "$run"/output_*; do
  test -d "$output" && rm -r -- "$output"
done
printf 'GRACKLE_NEQ_FEEDBACK_EVENT_WIRING_PASS %s\n' "$run"
