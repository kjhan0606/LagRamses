#!/usr/bin/env bash
# Extended Grackle NEQ cosmology validation on a scheduled compute node.
# Covers matched UV-on/UV-off short 128^3 GRAFIC trajectories, and a
# cooling-off hydro control for the previously observed intermediate-energy warning.
# Grackle is CPU-only here; the requested GPU is only the cluster partition token.
#SBATCH --job-name=grackle-neq-cosmo-ext
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-cosmo-ext-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
ic=/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007
binary=$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d
run=$root/.grackle-neq-cosmo-extended-${SLURM_JOB_ID}
test -x "$binary"
test ! -e "$run"
mkdir -p "$run/real_ic" "$run/uv_off" "$run/hydro_control"
for name in ic_deltab ic_velcx ic_velcy ic_velcz ic_poscx ic_poscy ic_poscz \
            ic_velbx ic_velby ic_velbz ic_pvar_00001 ic_pvar_00012; do
  test -r "$ic/$name"
done

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close
sha256sum "$binary" > "$run/binary.sha256"
printf 'job=%s\nhost=%s\nthreads=%s\n' "$SLURM_JOB_ID" "$(hostname -s)" "$OMP_NUM_THREADS" \
  > "$run/run_environment.txt"

cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" "$run/real_ic/effective.nml"
sed -i -e 's/nstepmax=1,/nstepmax=6,/' -e 's/foutput=1000000/foutput=6/' \
  "$run/real_ic/effective.nml"
cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" "$run/uv_off/effective.nml"
sed -i -e 's/nstepmax=1,/nstepmax=6,/' -e 's/foutput=1000000/foutput=6/' \
  -e 's/haardt_madau=.true./haardt_madau=.false./' "$run/uv_off/effective.nml"
cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" "$run/hydro_control/effective.nml"
sed -i -e 's/cooling=.true./cooling=.false./' \
  -e "s/cooling_method='grackle_neq'/cooling_method='original'/" \
  -e 's/haardt_madau=.true./haardt_madau=.false./' \
  "$run/hydro_control/effective.nml"

run_case() {
  local name=$1
  ( cd "$run/$name"
    TIMEFORMAT='wall=%3R user=%3U sys=%3S'
    { time timeout 3600s srun --ntasks=1 --cpus-per-task=8 "$binary" effective.nml \
        > run.log 2>&1; } 2> wall.txt
    grep -q 'Run completed' run.log
    if grep -E 'ERROR:|Aborting|Too many errors|failed|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' run.log; then
      exit 1
    fi
  )
}

run_case real_ic
run_case uv_off
run_case hydro_control

/home/kjhan/miniconda3/bin/python3 - "$run" <<'PY'
import glob, json, os, sys
import h5py
import numpy as np

root = sys.argv[1]
summary = {"runs": {}}
raw_states = {}
for label in ("real_ic", "uv_off"):
    outputs = sorted(glob.glob(os.path.join(root, label, "output_*")))
    if len(outputs) != 1:
        raise SystemExit(f"{label}: expected exactly one scheduled final snapshot, got {outputs}")
    files = glob.glob(os.path.join(outputs[0], "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"{label}: expected one HDF5 dump, got {files}")
    with h5py.File(files[0], "r") as f:
        arrays = []
        for group in f["hydro"]:
            level = int(group.split("_")[-1])
            if level < int(f["header"].attrs["levelmin"][0]):
                continue
            leaf = f[f"amr/{group}/son_flag"][:] == 0
            values = np.column_stack([f[f"hydro/{group}/uold_{i}"][:] for i in range(1, 28)])
            arrays.append(values[leaf])
        u = np.concatenate(arrays)
        species = u[:, 18:27]
        hydrogen = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
        helium = species[:, 3] + species[:, 4] + species[:, 5]
        charge = species[:, 1] - species[:, 2] + species[:, 7] / 2 + species[:, 4] / 4 + species[:, 5] / 2
        hres = float(np.max(np.abs(hydrogen-u[:, 6]) / u[:, 0]))
        heres = float(np.max(np.abs(helium-u[:, 7]) / u[:, 0]))
        cres = float(np.max(np.abs(species[:, 8]-charge) / u[:, 0]))
        eint = u[:, 4] - np.sum(u[:, 1:4]**2, axis=1)/(2*u[:, 0])
        result = {"snapshot": files[0], "leaf_cells": len(u),
                  "time": float(f["header"].attrs["time"][0]),
                  "finite": bool(np.isfinite(u).all()),
                  "min_density": float(np.min(u[:, 0])),
                  "min_species": float(np.min(species)),
                  "min_internal_energy_density": float(np.min(eint)),
                  "hydrogen_residual": hres, "helium_residual": heres,
                  "charge_residual": cres}
        if not result["finite"] or result["min_density"] <= 0 or result["min_species"] < 0 or result["min_internal_energy_density"] <= 0 or max(hres, heres, cres) >= 1e-6:
            raise SystemExit(f"{label}: final state failed validity gate: {result}")
        summary["runs"][label] = result
        raw_states[label] = u

reference = raw_states["real_ic"]
control = raw_states["uv_off"]
if reference.shape != control.shape:
    raise SystemExit(f"UV pair mesh/field shapes differ: {reference.shape} vs {control.shape}")
field_scale = np.maximum(np.max(np.abs(reference), axis=0), 1e-30)
summary["uv_on_off_comparison"] = {
    "max_field_scaled_difference": float(np.max(np.abs(reference-control)/field_scale)),
    "comparison_basis": "HDF leaf-cell array order; coordinates were not cross-matched",
    "mean_total_energy_density_relative_difference": float(abs(reference[:,4].mean()-control[:,4].mean()) / max(abs(reference[:,4].mean()), 1e-30)),
    "mean_cell_HII_density_over_gas_density_difference": float(np.mean(reference[:,19]/reference[:,0])-np.mean(control[:,19]/control[:,0]))}
del raw_states, reference, control

for label in ("real_ic", "uv_off", "hydro_control"):
    log = open(os.path.join(root, label, "run.log"), encoding="utf-8", errors="replace").read()
    summary["runs"].setdefault(label, {})["completed"] = "Run completed" in log
    summary["runs"][label]["fine_steps_logged"] = log.count("Fine step=")
    with open(os.path.join(root, label, "wall.txt"), encoding="utf-8") as f:
        summary["runs"][label]["elapsed"] = f.read().strip()
    if not summary["runs"][label]["completed"]:
        raise SystemExit(f"{label}: completion marker missing")
with open(os.path.join(root, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
PY

echo "GRACKLE_NEQ_EXTENDED_COSMO_COMPLETED $run"
