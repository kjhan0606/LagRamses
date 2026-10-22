#!/usr/bin/env bash
# Rebuild the Grackle NEQ route with the z_reion gate, then run focused UV and
# timestep evidence plus a matched 128^3 cosmological UV-on/off trajectory.
# A PASS here covers that trajectory only; it is not full production approval.
#SBATCH --job-name=grackle-prod-qual
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
# Job 415107 peaked at 21,003,500 KiB; reserve headroom for production reruns.
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-prod-qual-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
build=$root/.grackle-neq-production-build-${SLURM_JOB_ID}
run=$root/.grackle-neq-production-qualification-${SLURM_JOB_ID}
ic=/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
test -r "$uvdata" && test -d "$ic" && test ! -e "$build" && test ! -e "$run"
mkdir -p "$build" "$run/uv_on" "$run/uv_off" "$build/bin" "$build/run"
for folder in patch hydro mhd pm amr poisson aton utils; do ln -s "$root/$folder" "$build/$folder"; done
ln -s "$root/bin/Makefile" "$build/bin/Makefile"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close

cd "$build/bin"
make -j1 GRACKLE_NEQ=1 HDF5=1 USE_CUDA=0 SNRT=0 CHIMES=0 DUST_LIVE=0 \
  EXEC=ramses_grackle USE_FFTW=0 NVECTOR=32 > "$run/build.log" 2>&1
binary=$build/bin/ramses_grackle3d
bridge=$build/bin/grackle_neq_bridge.o
test -x "$binary" && test -r "$bridge"
sha256sum "$binary" "$bridge" "$uvdata" > "$run/assets.sha256"

# Unit-contract regression: equivalent proper states at a=0.1 and 0.2 must
# produce the same chemistry under Grackle's fixed-time comoving units.
mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_unit_contract_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/unit_contract_test"
sha256sum "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_unit_contract_test.c" \
  "$run/unit_contract_test" >> "$run/assets.sha256"
srun --ntasks=1 --cpus-per-task=1 "$run/unit_contract_test" "$uvdata" \
  | tee "$run/unit_contract.txt"
grep -q 'Grackle cosmological unit contract: PASS' "$run/unit_contract.txt"

# Focused regression: at z=9 (> z_reion=8.5), configured UV must be off;
# at z=7.33 it must be on. The z_reion=-0.5 sentinel is a test-only no-UV case.
mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_reionization_gate_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/reionization_gate_test"
for pair in '0.1 8.5' '0.1 -0.5' '0.1 100' '0.12 8.5' '0.12 -0.5'; do
  read -r a z_reion <<< "$pair"
  srun --ntasks=1 --cpus-per-task=1 "$run/reionization_gate_test" "$uvdata" "$a" "$z_reion" \
    | tee -a "$run/reionization_gate.txt"
done

# Chemistry interval convergence in a physical n_H/T span; evidence is retained
# even if the production-scale accuracy criterion fails.
mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_timestep_convergence_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/timestep_convergence_test"
srun --ntasks=1 --cpus-per-task=1 "$run/timestep_convergence_test" "$uvdata" \
  > "$run/timestep_convergence.log"

for name in uv_on uv_off; do
  cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" "$run/$name/effective.nml"
  sed -i -e 's/nstepmax=1,/nstepmax=80,/' \
    -e 's/aout=2d0/aout=0.12d0/' \
    -e 's/scheme=\x27muscl\x27, riemann=\x27hllc\x27/scheme=\x27muscl\x27, riemann=\x27hllc\x27\n pressure_fix=.true., beta_fix=0d0/' \
    -e 's/haardt_madau=.true./haardt_madau=.true., z_reion=8.5d0/' \
    "$run/$name/effective.nml"
done
sed -i 's/haardt_madau=.true./haardt_madau=.false./' "$run/uv_off/effective.nml"
for name in uv_on uv_off; do
  (
    cd "$run/$name"
    TIMEFORMAT='wall=%3R user=%3U sys=%3S'
    { time timeout 6900s srun --ntasks=1 --cpus-per-task=8 "$binary" effective.nml; } \
      > run.log 2> wall.txt
  )
  grep -q 'Run completed' "$run/$name/run.log"
  if grep -E 'ERROR:|Aborting|Too many errors|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' "$run/$name/run.log"; then
    exit 1
  fi
done

/home/kjhan/miniconda3/bin/python3 - "$run" <<'PY'
import glob
import json
import os
import re
import sys

import h5py
import numpy as np

root = sys.argv[1]
summary = {
    "run_root": root,
    "qualification_scope": "matched Grackle NEQ cosmological UV trajectory only",
    "runs": {},
    "acceptance": {},
}

def state(path):
    with h5py.File(path, "r") as f:
        coords, values = [], []
        for group in f["hydro"]:
            level = int(group.rsplit("_", 1)[-1])
            amr = f[f"amr/{group}"]
            centers = np.column_stack([amr[f"xg_{i}"][:] for i in (1, 2, 3)])
            offsets = np.array([[((child >> axis) & 1) - .5 for axis in range(3)]
                                for child in range(8)]) * 2.0 ** (-level)
            xyz = (centers[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
            leaf = amr["son_flag"][:] == 0
            if level >= int(f["header"].attrs["levelmin"][0]):
                coords.append(xyz[leaf])
                values.append(np.column_stack([f[f"hydro/{group}/uold_{i}"][:]
                                               for i in range(1, 28)])[leaf])
        xyz, u = np.concatenate(coords), np.concatenate(values)
        order = np.lexsort(xyz.T[::-1])
        aexp = float(np.ravel(f["header"].attrs["aexp"])[0])
        return xyz[order], u[order], aexp

for name in ("uv_on", "uv_off"):
    case = os.path.join(root, name)
    outputs = sorted(glob.glob(os.path.join(case, "output_*")))
    if len(outputs) != 1 or os.path.basename(outputs[0]) != "output_00001":
        raise SystemExit(f"{name}: expected exactly one scheduled snapshot, got {outputs}")
    files = glob.glob(os.path.join(outputs[0], "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"{name}: expected one HDF5 snapshot, got {files}")
    xyz, u, aexp = state(files[0])
    rho, species = u[:, 0], u[:, 18:27]
    h = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
    he = species[:, 3] + species[:, 4] + species[:, 5]
    charge = species[:, 1] - species[:, 2] + species[:, 7]/2 + species[:, 4]/4 + species[:, 5]/2
    eint = u[:, 4] - np.sum(u[:, 1:4]**2, axis=1)/(2*rho)
    log = open(os.path.join(case, "run.log"), encoding="utf-8", errors="replace").read()
    negative = [line.strip() for line in log.splitlines()
                if "eint_min=" in line and re.search(r"eint_min=\s*-", line)]
    step_a = [float(v) for v in re.findall(r"Fine step=.*?\ba=\s*([0-9.Ee+-]+)", log)]
    summary["runs"][name] = {
        "snapshot": files[0], "leaf_cells": len(u), "aexp": aexp,
        "aexp_last_logged": step_a[-1] if step_a else None,
        "finite": bool(np.isfinite(u).all()), "min_density": float(rho.min()),
        "min_species": float(species.min()), "min_internal_energy_density": float(eint.min()),
        "max_H_residual": float(np.max(np.abs(h-u[:, 6])/rho)),
        "max_He_residual": float(np.max(np.abs(he-u[:, 7])/rho)),
        "max_charge_residual": float(np.max(np.abs(species[:, 8]-charge)/rho)),
        "negative_post_hydro_diagnostics": negative,
        "wall": open(os.path.join(case, "wall.txt"), encoding="utf-8").read().strip(),
    }

on = state(summary["runs"]["uv_on"]["snapshot"])
off = state(summary["runs"]["uv_off"]["snapshot"])
hii_on = float(on[1][:, 19].sum()/on[1][:, 6].sum())
hii_off = float(off[1][:, 19].sum()/off[1][:, 6].sum())
summary["comparison"] = {
    "coordinates_match": bool(np.array_equal(on[0], off[0])),
    "scale_factor_abs_difference": abs(on[2]-off[2]),
    "mass_weighted_HII_over_H_uv_on": hii_on,
    "mass_weighted_HII_over_H_uv_off": hii_off,
}
summary["acceptance"] = {
    "one_snapshot_per_run": True,
    "reached_post_reionization_target": min(on[2], off[2]) >= 0.12,
    "finite_positive_chemistry_and_budgets": all(
        r["finite"] and r["min_density"] > 0 and r["min_species"] >= 0
        and r["min_internal_energy_density"] > 0
        and max(r["max_H_residual"], r["max_He_residual"], r["max_charge_residual"]) < 1e-6
        for r in summary["runs"].values()),
    "no_negative_post_hydro_energy_diagnostics": all(
        not r["negative_post_hydro_diagnostics"] for r in summary["runs"].values()),
    "matched_fixed_level_mesh": bool(np.array_equal(on[0], off[0])),
    "matched_scale_factor": abs(on[2]-off[2]) < 2e-3,
    "UV_increases_mass_weighted_HII": hii_on > hii_off,
}
summary["acceptance"]["trajectory_gate_pass"] = all(summary["acceptance"].values())
with open(os.path.join(root, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
if not summary["acceptance"]["trajectory_gate_pass"]:
    raise SystemExit("Grackle cosmological trajectory gate failed; preserve raw output for diagnosis")
PY

# Keep compact evidence; remove only the two evaluated raw test output directories.
rm -r -- "$run/uv_on/output_00001" "$run/uv_off/output_00001"
printf 'GRACKLE_NEQ_COSMOLOGICAL_TRAJECTORY_PASS %s\n' "$run"
