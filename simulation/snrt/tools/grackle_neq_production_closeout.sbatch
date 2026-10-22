#!/usr/bin/env bash
# Production-closeout bundle: stiff chemistry time-step convergence plus a
# pressure-fix-enabled cosmological trajectory through the HM2012 UV epoch.
# Grackle is CPU chemistry; the A40 request is only this cluster's allocation token.
#SBATCH --job-name=grackle-neq-closeout
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-closeout-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
ic=/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007
binary=$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d
bridge=$root/.grackle-neq-build-lageunha-20261006T1943/bin/grackle_neq_bridge.o
run=$root/.grackle-neq-production-closeout-${SLURM_JOB_ID}
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
test -x "$binary"
test -r "$bridge"
test -r "$uvdata"
test ! -e "$run"
mkdir -p "$run/hydro_control" "$run/cosmological_neq"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close

sha256sum "$binary" "$bridge" /home/kjhan/local/lib64/libgrackle.a \
  > "$run/build_assets.sha256"
printf 'job=%s\nhost=%s\nthreads=%s\n' "$SLURM_JOB_ID" "$(hostname -s)" "$OMP_NUM_THREADS" \
  > "$run/run_environment.txt"
printf 'noutput=1 aout=2 tout=1d100; periodic output only at foutput=24; backups suppressed\n' \
  >> "$run/run_environment.txt"
printf 'expected full dumps=2; estimated <=1 GiB each, <=2 GiB total; /gpfs free at submission=231 TiB\n' \
  >> "$run/run_environment.txt"

mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_timestep_convergence_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/grackle_neq_timestep_convergence_test"
sha256sum "$run/grackle_neq_timestep_convergence_test" >> "$run/build_assets.sha256"
(
  TIMEFORMAT='wall=%3R user=%3U sys=%3S'
  { time srun --ntasks=1 --cpus-per-task=1 "$run/grackle_neq_timestep_convergence_test" "$uvdata"; } \
    > "$run/timestep_convergence.log" 2> "$run/timestep_convergence.wall"
)

for name in hydro_control cosmological_neq; do
  cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" \
    "$run/$name/effective.nml"
  sed -i -e 's/nstepmax=1,/nstepmax=24,/' \
    -e 's/foutput=1000000/foutput=24/' \
    -e "s/scheme='muscl', riemann='hllc'/scheme='muscl', riemann='hllc'\n pressure_fix=.true., beta_fix=0d0/" \
    "$run/$name/effective.nml"
done
sed -i -e 's/nstepmax=24,/nstepmax=64,/' -e 's/foutput=24/foutput=64/' \
  "$run/cosmological_neq/effective.nml"
sed -i -e 's/cooling=.true./cooling=.false./' \
  -e "s/cooling_method='grackle_neq'/cooling_method='original'/" \
  -e 's/haardt_madau=.true./haardt_madau=.false./' \
  "$run/hydro_control/effective.nml"

run_case() {
  local name=$1
  (
    cd "$run/$name"
    TIMEFORMAT='wall=%3R user=%3U sys=%3S'
    { time timeout 6900s srun --ntasks=1 --cpus-per-task=8 "$binary" effective.nml; } \
      > run.log 2> wall.txt
  )
  grep -q 'Run completed' "$run/$name/run.log"
  if grep -E 'ERROR:|Aborting|Too many errors|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' "$run/$name/run.log"; then
    return 1
  fi
}
run_case hydro_control
run_case cosmological_neq

/home/kjhan/miniconda3/bin/python3 - "$run" <<'PY'
import glob
import json
import os
import re
import sys

import h5py
import numpy as np

root = sys.argv[1]
summary = {"run_root": root, "runs": {}, "acceptance": {}}

def read_state(path):
    with h5py.File(path, "r") as f:
        rows, coords = [], []
        for group in f["hydro"]:
            level = int(group.rsplit("_", 1)[-1])
            amr = f[f"amr/{group}"]
            centers = np.column_stack([amr[f"xg_{i}"][:] for i in (1, 2, 3)])
            offsets = np.array([[((child >> axis) & 1) - .5 for axis in range(3)]
                                for child in range(8)]) * 2.0 ** (-level)
            xyz = (centers[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
            leaf = amr["son_flag"][:] == 0
            if level < int(f["header"].attrs["levelmin"][0]):
                continue
            coords.append(xyz[leaf])
            rows.append(np.column_stack([f[f"hydro/{group}/uold_{i}"][:]
                                          for i in range(1, 28)])[leaf])
        xyz, u = np.concatenate(coords), np.concatenate(rows)
        order = np.lexsort(xyz.T[::-1])
        attrs = f["header"].attrs
        aexp = float(np.ravel(attrs["aexp"])[0]) if "aexp" in attrs else None
        return xyz[order], u[order], aexp

for name in ("hydro_control", "cosmological_neq"):
    case = os.path.join(root, name)
    outputs = sorted(glob.glob(os.path.join(case, "output_*")))
    if len(outputs) != 1:
        raise SystemExit(f"{name}: expected exactly one scheduled dump, got {outputs}")
    if os.path.basename(outputs[0]) != "output_00001":
        raise SystemExit(f"{name}: unexpected output name; refusing automatic cleanup: {outputs[0]}")
    files = glob.glob(os.path.join(outputs[0], "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"{name}: expected one HDF5 file, got {files}")
    xyz, u, aexp = read_state(files[0])
    rho = u[:, 0]
    species = u[:, 18:27]
    h = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
    he = species[:, 3] + species[:, 4] + species[:, 5]
    charge = species[:, 1] - species[:, 2] + species[:, 7] / 2 + species[:, 4] / 4 + species[:, 5] / 2
    eint = u[:, 4] - np.sum(u[:, 1:4] ** 2, axis=1) / (2 * rho)
    hres = float(np.max(np.abs(h - u[:, 6]) / rho))
    heres = float(np.max(np.abs(he - u[:, 7]) / rho))
    eres = float(np.max(np.abs(species[:, 8] - charge) / rho))
    log = open(os.path.join(case, "run.log"), encoding="utf-8", errors="replace").read()
    negative = [line.strip() for line in log.splitlines()
                if "eint_min=" in line and re.search(r"eint_min=\s*-", line)]
    fine_steps = [float(x) for x in re.findall(r"Fine step=.*?\ba=\s*([0-9.Ee+-]+)", log)]
    summary["runs"][name] = {
        "snapshot": files[0], "leaf_cells": len(u), "aexp_header": aexp,
        "aexp_last_logged": fine_steps[-1] if fine_steps else None,
        "finite": bool(np.isfinite(u).all()), "min_density": float(rho.min()),
        "min_species": float(species.min()), "min_internal_energy_density": float(eint.min()),
        "max_H_residual": hres, "max_He_residual": heres, "max_charge_residual": eres,
        "negative_post_hydro_diagnostics": negative,
        "uv_background_requested": "haardt_madau=.true." in open(
            os.path.join(case, "effective.nml"), encoding="utf-8").read(),
        "wall": open(os.path.join(case, "wall.txt"), encoding="utf-8").read().strip(),
    }

control = read_state(summary["runs"]["hydro_control"]["snapshot"])
cosmo = read_state(summary["runs"]["cosmological_neq"]["snapshot"])
same_grid = np.array_equal(control[0], cosmo[0])
hii_control = float(control[1][:, 19].sum() / control[1][:, 6].sum())
hii_grackle = float(cosmo[1][:, 19].sum() / cosmo[1][:, 6].sum())
summary["comparison"] = {"coordinates_match": bool(same_grid),
                         "mass_weighted_HII_over_H_no_cooling_control": hii_control,
                         "mass_weighted_HII_over_H_Grackle_case_not_a_UV_pair": hii_grackle}
summary["acceptance"] = {
    "one_dump_per_run": True,
    "all_final_hydro_states_finite_positive": all(
        r["finite"] and r["min_density"] > 0 and r["min_internal_energy_density"] > 0
        for r in summary["runs"].values()),
    "Grackle_species_nonnegative_and_budgets": (
        summary["runs"]["cosmological_neq"]["min_species"] >= 0
        and max(summary["runs"]["cosmological_neq"]["max_H_residual"],
                summary["runs"]["cosmological_neq"]["max_He_residual"],
                summary["runs"]["cosmological_neq"]["max_charge_residual"]) < 1e-6),
    "pressure_fix_removed_pre_cooling_negative_energy": all(
        not r["negative_post_hydro_diagnostics"] for r in summary["runs"].values()),
    "reached_post_reionization_epoch_with_UVB_requested": (
        summary["runs"]["cosmological_neq"]["uv_background_requested"]
        and summary["runs"]["cosmological_neq"]["aexp_last_logged"] is not None
        and 1.0/summary["runs"]["cosmological_neq"]["aexp_last_logged"]-1.0 < 8.5),
    "matched_final_mesh": bool(same_grid),
    "Grackle_HII_field_nonzero": hii_grackle > 0.0,
}
summary["acceptance"]["pass"] = all(summary["acceptance"].values())
with open(os.path.join(root, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
if not summary["acceptance"]["pass"]:
    raise SystemExit("Grackle production-closeout integration gates failed; preserve raw output for evaluation")
PY

# The output schedule was asserted to be exactly one dump per fresh test run;
# keep only the compact summary, effective inputs, logs and provenance.
rm -r -- "$run/hydro_control/output_00001" "$run/cosmological_neq/output_00001"
printf 'GRACKLE_NEQ_PRODUCTION_CLOSEOUT_COMPLETED %s\n' "$run"
