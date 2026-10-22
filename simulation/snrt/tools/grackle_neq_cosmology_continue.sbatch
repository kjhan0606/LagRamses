#!/usr/bin/env bash
# Continue the matched pressure-fix 128^3 pair from job 415070's a~0.062
# checkpoint to the HM2012 UV-background epoch, without repeating the prefix.
#SBATCH --job-name=grackle-neq-continue
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-continue-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
source_run=$root/.grackle-neq-production-closeout-415070
binary=$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d
run=$root/.grackle-neq-cosmology-continue-${SLURM_JOB_ID}
test -x "$binary" && test -d "$source_run"
test -r "$source_run/hydro_control/output_00001/data_00001.h5"
test -r "$source_run/cosmological_neq/output_00001/data_00001.h5"
test ! -e "$run"
mkdir -p "$run/hydro_control" "$run/cosmological_neq"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close
sha256sum "$binary" > "$run/binary.sha256"
printf 'continuation_of_job=415070\ninitial_dump=aexp~0.062; final_global_step=64\n' \
  > "$run/run_environment.txt"
printf 'aout=2, noutput=1, tout=1d100; foutput=64; fbackup=1000000; one scheduled final dump/case\n' \
  >> "$run/run_environment.txt"

for name in hydro_control cosmological_neq; do
  cp -a "$source_run/$name/output_00001" "$run/$name/output_00001"
  cp "$source_run/$name/effective.nml" "$run/$name/effective.nml"
  sed -i -e 's/nrestart=0, nstepmax=24,/nrestart=1, nstepmax=64,/' \
    -e 's/foutput=24/foutput=64/' "$run/$name/effective.nml"
done

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

/home/kjhan/miniconda3/bin/python3 - "$run" "$source_run" <<'PY'
import glob
import json
import os
import re
import sys

import h5py
import numpy as np

root, source = sys.argv[1:]
summary = {"continuation_root": root, "prefix_job": 415070, "runs": {}, "acceptance": {}}

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
    if len(outputs) != 2 or not outputs[0].endswith("output_00001"):
        raise SystemExit(f"{name}: expected one restart input and one final output, got {outputs}")
    files = glob.glob(os.path.join(outputs[-1], "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"{name}: expected exactly one final HDF5 file, got {files}")
    xyz, u, aexp = read_state(files[0])
    rho, species = u[:, 0], u[:, 18:27]
    h = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
    he = species[:, 3] + species[:, 4] + species[:, 5]
    charge = species[:, 1] - species[:, 2] + species[:, 7] / 2 + species[:, 4] / 4 + species[:, 5] / 2
    eint = u[:, 4] - np.sum(u[:, 1:4] ** 2, axis=1) / (2 * rho)
    log = open(os.path.join(case, "run.log"), encoding="utf-8", errors="replace").read()
    prior_log = open(os.path.join(source, name, "run.log"), encoding="utf-8", errors="replace").read()
    negative = [line.strip() for text in (prior_log, log) for line in text.splitlines()
                if "eint_min=" in line and re.search(r"eint_min=\s*-", line)]
    fine_steps = [int(x) for x in re.findall(r"Fine step=\s*(\d+)", log)]
    fine_a = [float(x) for x in re.findall(r"Fine step=.*?\ba=\s*([0-9.Ee+-]+)", log)]
    summary["runs"][name] = {
        "snapshot": files[0], "leaf_cells": len(u), "aexp_header": aexp,
        "global_final_step": fine_steps[-1] if fine_steps else None,
        "aexp_last_logged": fine_a[-1] if fine_a else None,
        "finite": bool(np.isfinite(u).all()), "min_density": float(rho.min()),
        "min_species": float(species.min()), "min_internal_energy_density": float(eint.min()),
        "max_H_residual": float(np.max(np.abs(h-u[:, 6])/rho)),
        "max_He_residual": float(np.max(np.abs(he-u[:, 7])/rho)),
        "max_charge_residual": float(np.max(np.abs(species[:, 8]-charge)/rho)),
        "negative_post_hydro_diagnostics": negative,
        "uv_background_requested": "haardt_madau=.true." in open(
            os.path.join(case, "effective.nml"), encoding="utf-8").read(),
        "wall": open(os.path.join(case, "wall.txt"), encoding="utf-8").read().strip(),
    }

c = read_state(summary["runs"]["hydro_control"]["snapshot"])
g = read_state(summary["runs"]["cosmological_neq"]["snapshot"])
hii_control = float(c[1][:, 19].sum() / c[1][:, 6].sum())
hii_grackle = float(g[1][:, 19].sum() / g[1][:, 6].sum())
summary["comparison"] = {"coordinates_match": bool(np.array_equal(c[0], g[0])),
                         "mass_weighted_HII_over_H_no_cooling_control": hii_control,
                         "mass_weighted_HII_over_H_Grackle_case_not_a_UV_pair": hii_grackle}
summary["acceptance"] = {
    "final_global_step_64": all(r["global_final_step"] == 64 for r in summary["runs"].values()),
    "reached_post_reionization_epoch_with_UVB_requested": (
        summary["runs"]["cosmological_neq"]["uv_background_requested"]
        and summary["runs"]["cosmological_neq"]["aexp_last_logged"] is not None
        and 1.0/summary["runs"]["cosmological_neq"]["aexp_last_logged"]-1.0 < 8.5),
    "all_final_hydro_states_finite_positive": all(
        r["finite"] and r["min_density"] > 0 and r["min_internal_energy_density"] > 0
        for r in summary["runs"].values()),
    "Grackle_species_nonnegative_and_budgets": (
        summary["runs"]["cosmological_neq"]["min_species"] >= 0
        and max(summary["runs"]["cosmological_neq"]["max_H_residual"],
                summary["runs"]["cosmological_neq"]["max_He_residual"],
                summary["runs"]["cosmological_neq"]["max_charge_residual"]) < 1e-6),
    "no_negative_post_hydro_energy_diagnostics": all(
        not r["negative_post_hydro_diagnostics"] for r in summary["runs"].values()),
    "matched_final_mesh": bool(np.array_equal(c[0], g[0])),
    "Grackle_HII_field_nonzero": hii_grackle > 0.0,
}
summary["acceptance"]["pass"] = all(summary["acceptance"].values())
with open(os.path.join(root, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
if not summary["acceptance"]["pass"]:
    raise SystemExit("Cosmology/UV continuation gate failed; preserve raw outputs for evaluation")
PY

# Remove only output directories belonging to this completed test and its exact
# 415070 restart prefixes; retain compact inputs, logs, checksums and summaries.
for path in "$run/hydro_control"/output_* "$run/cosmological_neq"/output_* \
            "$source_run/hydro_control/output_00001" "$source_run/cosmological_neq/output_00001"; do
  test -d "$path" || continue
  case "$path" in
    "$run/hydro_control"/output_*|"$run/cosmological_neq"/output_*|\
    "$source_run/hydro_control/output_00001"|"$source_run/cosmological_neq/output_00001") rm -r -- "$path" ;;
    *) echo "Refusing unexpected cleanup target: $path" >&2; exit 2 ;;
  esac
done
printf 'GRACKLE_NEQ_COSMOLOGY_CONTINUATION_COMPLETED %s\n' "$run"
