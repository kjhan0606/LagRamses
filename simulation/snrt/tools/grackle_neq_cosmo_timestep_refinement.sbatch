#!/usr/bin/env bash
# Three cosmological trajectories at Courant factors 0.5, 0.25 and 0.125.
# match_aout lands all three trajectories at the same a=0.1065 epoch, avoiding
# epoch mismatch and output interpolation in the timestep comparison.
#SBATCH --job-name=grackle-dt-conv
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-dt-conv-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
build=$root/.grackle-neq-production-build-415107
binary=$build/bin/ramses_grackle3d
run=$root/.grackle-neq-cosmo-timestep-convergence-${SLURM_JOB_ID}
template=$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml
ic=/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
test -x "$binary" && test -r "$template" && test -d "$ic" && test -r "$uvdata"
test ! -e "$run"
mkdir -p "$run"
variants=(base half_dt quarter_dt)

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close

for variant in "${variants[@]}"; do
  mkdir -p "$run/$variant"
  cp "$template" "$run/$variant/effective.nml"
  sed -i \
    -e 's/nstepmax=1,/nstepmax=220,/' \
    -e 's/noutput=1, tout=1d100, aout=2d0/noutput=1, tout=1d100, aout=0.1065d0, match_aout=.true./' \
    -e 's/haardt_madau=.true./haardt_madau=.true., z_reion=8.5d0/' \
    "$run/$variant/effective.nml"
  sed -i "/scheme='muscl', riemann='hllc'/a\\ pressure_fix=.true., beta_fix=0d0" \
    "$run/$variant/effective.nml"
done
sed -i 's/courant_factor=.5d0/courant_factor=.25d0/' "$run/half_dt/effective.nml"
sed -i 's/courant_factor=.5d0/courant_factor=.125d0/' "$run/quarter_dt/effective.nml"
sha256sum "$binary" "$template" "$uvdata" \
  "$root/patch/cuRamses/cooling_fine.kjhan.f90" \
  "$root/patch/lagRamses/grackle_neq_bridge.c" \
  "$root/simulation/snrt/tools/grackle_neq_cosmo_timestep_refinement.sbatch" \
  > "$run/assets.sha256"
for variant in "${variants[@]}"; do
  (
    cd "$run/$variant"
    TIMEFORMAT='wall=%3R user=%3U sys=%3S'
    { time timeout 4500s srun --ntasks=1 --cpus-per-task=8 "$binary" effective.nml; } \
      > run.log 2> wall.txt
  )
  grep -q 'Run completed' "$run/$variant/run.log"
  if grep -E 'ERROR:|Aborting|Too many errors|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' \
       "$run/$variant/run.log"; then
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
target_aexp = float(os.environ.get("GRACKLE_TIMESTEP_TARGET_AEXP", "0.1065"))
hydrogen_species = {"HI": 0, "HII": 1, "HM": 2, "H2I_H_nuclei": 6, "H2II_H_nuclei": 7}
helium_species = {"HeI": 3, "HeII": 4, "HeIII": 5}
abundance_thresholds = (1e-5, 1e-4, 1e-3, 1e-2)
thermal_thresholds = (0.01, 0.05, 0.10, 0.50)

def snapshot(path):
    files = glob.glob(os.path.join(path, "data_*.h5"))
    if len(files) != 1:
        raise SystemExit(f"{path}: expected one rank HDF5 file, got {files}")
    with h5py.File(files[0], "r") as f:
        coords, values, volumes = [], [], []
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
                volumes.append(np.full(int(leaf.sum()), 2.0 ** (-3*level)))
        xyz, u, vol = np.concatenate(coords), np.concatenate(values), np.concatenate(volumes)
        order = np.lexsort(xyz.T[::-1])
        xyz, u, vol = xyz[order], u[order], vol[order]
        aexp = float(np.ravel(f["header"].attrs["aexp"])[0])
    return xyz, u, vol, aexp

def state(variant):
    case = os.path.join(root, variant)
    outputs = sorted(glob.glob(os.path.join(case, "output_*")))
    if [os.path.basename(x) for x in outputs] != ["output_00001"]:
        raise SystemExit(f"{variant}: expected one exact-epoch snapshot, got {outputs}")
    xyz, u, vol, aexp = snapshot(outputs[0])
    if abs(aexp-target_aexp) > 1e-10:
        raise SystemExit(f"{variant}: match_aout wrote a={aexp:.15g}, requested {target_aexp}")
    species = u[:, 18:27]
    rho = u[:, 0]
    h = species[:, 0] + species[:, 1] + species[:, 2] + species[:, 6] + species[:, 7]
    he = species[:, 3] + species[:, 4] + species[:, 5]
    charge = species[:, 1] - species[:, 2] + species[:, 7]/2 + species[:, 4]/4 + species[:, 5]/2
    eint = u[:, 4] - np.sum(u[:, 1:4]**2, axis=1)/(2*rho)
    abundances = {name: species[:, i]/h for name, i in hydrogen_species.items()}
    abundances.update({name: species[:, i]/he for name, i in helium_species.items()})
    abundances["e_per_H_nucleus"] = species[:, 8]/h
    abundance_closure = max(
        float(np.max(np.abs(sum(abundances[name] for name in hydrogen_species)-1.0))),
        float(np.max(np.abs(sum(abundances[name] for name in helium_species)-1.0))),
    )
    wall = open(os.path.join(case, "wall.txt"), encoding="utf-8").read().strip()
    step_a = [float(v) for v in re.findall(r"Fine step=.*?\ba=\s*([0-9.Ee+-]+)",
                                            open(os.path.join(case, "run.log"), encoding="utf-8").read())]
    return {
        "coordinates": xyz, "fields": u, "volume": vol, "abundances": abundances, "aexp": aexp,
        "finite": bool(np.isfinite(u).all()), "min_density": float(rho.min()),
        "min_H_density": float(h.min()), "min_He_density": float(he.min()),
        "min_species": float(species.min()), "min_internal_energy_density": float(eint.min()),
        "max_H_residual": float(np.max(np.abs(h-u[:, 6])/rho)),
        "max_He_residual": float(np.max(np.abs(he-u[:, 7])/rho)),
        "max_charge_residual": float(np.max(np.abs(species[:, 8]-charge)/rho)),
        "max_elemental_abundance_closure": abundance_closure,
        "mass_weighted_HII_over_H": float(np.sum(species[:, 1]*vol)/np.sum(h*vol)),
        "mass_weighted_HeIII_over_He": float(np.sum(species[:, 5]*vol)/np.sum(he*vol)),
        "wall": wall, "aexp_last_logged": step_a[-1] if step_a else None,
    }

states = {name: state(name) for name in ("base", "half_dt", "quarter_dt")}
valid = all(x["finite"] and x["min_density"] > 0 and x["min_H_density"] > 0
            and x["min_He_density"] > 0 and x["min_species"] >= 0
            and x["min_internal_energy_density"] > 0
            and max(x["max_H_residual"], x["max_He_residual"], x["max_charge_residual"],
                    x["max_elemental_abundance_closure"]) < 1e-6
            for x in states.values())
if not valid:
    raise SystemExit("one timestep-convergence trajectory failed state/chemistry budgets")
if any(not np.array_equal(states["base"]["coordinates"], x["coordinates"])
       for x in (states["half_dt"], states["quarter_dt"])):
    raise SystemExit("Courant matrix did not retain the same fixed leaf mesh")

def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    sorted_values, sorted_weights = values[order], weights[order]
    cumulative = np.cumsum(sorted_weights)
    return float(sorted_values[min(np.searchsorted(cumulative, q*cumulative[-1], side="left"),
                                  len(sorted_values)-1)])

def local_distribution(values, volume, mass, thresholds):
    volume_total, mass_total = float(volume.sum()), float(mass.sum())
    return {
        "p50_volume": float(np.quantile(values, 0.50)),
        "p95_volume": float(np.quantile(values, 0.95)),
        "p99_volume": float(np.quantile(values, 0.99)),
        "max": float(values.max()),
        "p95_mass_weighted": weighted_quantile(values, mass, 0.95),
        "exceedance": {
            f"gt_{threshold:g}": {
                "cell_count": int(np.count_nonzero(values > threshold)),
                "volume_fraction": float(volume[values > threshold].sum()/volume_total),
                "mass_fraction": float(mass[values > threshold].sum()/mass_total),
            } for threshold in thresholds
        },
    }

def compare(coarse, fine):
    delta = np.abs(coarse["fields"] - fine["fields"])
    e0 = coarse["fields"][:, 4] - np.sum(coarse["fields"][:, 1:4]**2, axis=1)/(2*coarse["fields"][:, 0])
    e1 = fine["fields"][:, 4] - np.sum(fine["fields"][:, 1:4]**2, axis=1)/(2*fine["fields"][:, 0])
    relative_energy = np.abs(e0-e1)/np.maximum(np.abs(e1), 1e-30)
    relative_density = np.abs(coarse["fields"][:,0]-fine["fields"][:,0]) / \
                       np.maximum(fine["fields"][:,0], 1e-30)
    volume = 0.5*(coarse["volume"]+fine["volume"])
    # Symmetric pair weight avoids privileging either Courant trajectory when
    # reporting the gas-mass fraction of cells in a local-error tail.
    pair_mass = 0.5*(coarse["fields"][:,0]+fine["fields"][:,0])*volume
    abundance_delta = {
        name: np.abs(coarse["abundances"][name]-fine["abundances"][name])
        for name in coarse["abundances"]
    }
    return {
        "aexp_abs_difference": abs(coarse["aexp"]-fine["aexp"]),
        "max_field_abs_difference_by_index": np.max(delta, axis=0).tolist(),
        "thermal_internal_energy_density_relative_difference": local_distribution(
            relative_energy, volume, pair_mass, thermal_thresholds),
        "gas_density_relative_difference": local_distribution(
            relative_density, volume, pair_mass, thermal_thresholds),
        "elemental_species_abundance_absolute_difference": {
            name: local_distribution(values, volume, pair_mass, abundance_thresholds)
            for name, values in abundance_delta.items()
        },
        "mass_weighted_HII_over_H_abs_difference": abs(coarse["mass_weighted_HII_over_H"]-
                                                        fine["mass_weighted_HII_over_H"]),
        "mass_weighted_HeIII_over_He_abs_difference": abs(coarse["mass_weighted_HeIII_over_He"]-
                                                           fine["mass_weighted_HeIII_over_He"]),
        "local_tail_mass_weighting": "symmetric pair weight: 0.5*(rho_coarse+rho_fine)*cell_volume",
        "thermal_exceedance_bins_are_acceptance_thresholds": False,
        "abundance_exceedance_bins_are_acceptance_thresholds": False,
    }

pairwise = {
    "courant_0.5_vs_0.25": compare(states["base"], states["half_dt"]),
    "courant_0.25_vs_0.125": compare(states["half_dt"], states["quarter_dt"]),
    "courant_0.5_vs_0.125": compare(states["base"], states["quarter_dt"]),
}
denom = pairwise["courant_0.5_vs_0.25"]["thermal_internal_energy_density_relative_difference"]["p95_volume"]
pairwise["convergence_diagnostics"] = {
    "p95_energy_finer_to_coarser_pair_ratio":
        pairwise["courant_0.25_vs_0.125"]["thermal_internal_energy_density_relative_difference"]["p95_volume"] / denom
        if denom > 0 else None,
    "p95_HII_over_H_finer_to_coarser_pair_ratio":
        pairwise["courant_0.25_vs_0.125"]["elemental_species_abundance_absolute_difference"]["HII"]["p95_volume"] /
        pairwise["courant_0.5_vs_0.25"]["elemental_species_abundance_absolute_difference"]["HII"]["p95_volume"]
        if pairwise["courant_0.5_vs_0.25"]["elemental_species_abundance_absolute_difference"]["HII"]["p95_volume"] > 0 else None,
    "interpretation": "descriptive ratios only; no post-hoc physical tolerance or pass threshold",
}
summary = {
    "scope": "three-level live cosmological hydro/Grackle timestep comparison at a directly matched post-reionization scale factor",
    "target_aexp": target_aexp,
    "abundance_definition": "species mass density divided by total mass density of the corresponding element; molecular H species are expressed as hydrogen-nuclei mass fractions",
    "diagnostic_bins": {
        "thermal_internal_energy_density_relative_difference": list(thermal_thresholds),
        "elemental_abundance_absolute_difference": list(abundance_thresholds),
        "meaning": "descriptive distribution bins only, not pass/fail tolerances",
    },
    "variants": {name: {k:v for k,v in state.items()
                         if k not in ("coordinates", "fields", "volume", "abundances")}
                 for name, state in states.items()},
    "pairwise": pairwise,
    "timestep_accuracy_verdict": "three-level numerical trend measured; scientific acceptance remains scoped and explicit",
}
with open(os.path.join(root, "timestep_convergence_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
PY

# Preserve compact evidence and remove only the three raw snapshots after the
# full three-way comparison has succeeded.
for variant in "${variants[@]}"; do
  rm -r -- "$run/$variant/output_00001"
done
"/home/kjhan/miniconda3/bin/python3" - "$run" "$SLURM_JOB_ID" <<'PY'
import json
import os
import sys

run, job_id = sys.argv[1:]
variants = ("base", "half_dt", "quarter_dt")
removed = []
for variant in variants:
    path = os.path.join(run, variant, "output_00001")
    if os.path.exists(path):
        raise SystemExit(f"cleanup incomplete; raw snapshot still exists: {path}")
    removed.append({"variant": variant, "snapshot": path,
                    "removed_after_successful_analysis": True})
manifest = {
    "slurm_job_id": job_id,
    "comparison_summary": os.path.join(run, "timestep_convergence_summary.json"),
    "raw_outputs_removed": removed,
    "retained_per_variant": ["effective.nml", "run.log", "wall.txt"],
    "retained_run_evidence": ["assets.sha256", "timestep_convergence_summary.json"],
}
with open(os.path.join(run, "cleanup_manifest.json"), "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)
PY
printf 'GRACKLE_NEQ_HYDRO_TIMESTEP_CONVERGENCE_EVALUATED %s\n' "$run"
