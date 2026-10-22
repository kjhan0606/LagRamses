"""Small integrated Grackle snapshots: coordinate-aligned state comparisons."""
import json
import sys
import subprocess
from pathlib import Path

import h5py
import numpy as np


def read(path, details=False):
    with h5py.File(path, "r") as f:
        rows = []
        coordinates = []
        volumes = []
        for group in f["hydro"]:
            level = int(group.split("_")[-1])
            amr = f[f"amr/{group}"]
            centers = np.column_stack([amr[f"xg_{i}"][:] for i in (1, 2, 3)])
            offsets = np.array([[((child >> axis) & 1) - .5 for axis in range(3)]
                                for child in range(8)]) * 2.0 ** (-level)
            xyz = (centers[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
            leaf = amr["son_flag"][:] == 0
            if level < int(f["header"].attrs["levelmin"][0]):
                continue
            coordinates.append(xyz[leaf])
            rows.append(np.column_stack([f[f"hydro/{group}/uold_{i}"][:]
                                          for i in range(1, 28)])[leaf])
            volumes.append(np.full(int(leaf.sum()), 2.0**(-3*level)))
        xyz = np.concatenate(coordinates)
        values = np.concatenate(rows)
        order = np.lexsort(xyz.T[::-1])
        if details:
            header = {key: float(np.ravel(f["header"].attrs[key])[0])
                      for key in ("aexp", "omega_m", "h0", "boxlen_ini", "gamma")}
            return xyz[order], values[order], np.concatenate(volumes)[order], header
        return xyz[order], values[order], float(f["header"].attrs["time"][0])


def grackle_temperature(u, header, library, data_file):
    """Mirror VPATH-selected cuRamses/units.f90 and comoving energy units."""
    a, h0 = header["aexp"], header["h0"]
    assert abs(header["gamma"]-5/3) < 1e-14
    density_unit = header["omega_m"]*1.88e-29*(h0/100)**2/a**3
    length_unit = a*header["boxlen_ini"]*3.08e24/(h0/100)
    fixed_time_unit = 3.08e24/(h0*1e5)
    epsilon = (u[:, 4]-np.sum(u[:, 1:4]**2, axis=1)/(2*u[:, 0]))/u[:, 0]
    # RAMSES velocity unit is LU/(a^2 TU); Grackle's is LU/(a TU).
    packed = np.ascontiguousarray(np.column_stack(
        (u[:, 0], epsilon/a**2, u[:, 18:27], u[:, 5])).T, dtype=np.float64)
    process = subprocess.run([str(library), str(len(u)), str(density_unit), str(length_unit),
                              str(fixed_time_unit), str(a), str(data_file)],
                             input=packed.tobytes(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if process.returncode:
        raise RuntimeError(f"Grackle temperature API failed: {process.returncode}, {process.stderr.decode()}")
    temperature = np.frombuffer(process.stdout, dtype=np.float64).copy()
    assert temperature.shape == (len(u),)
    return epsilon, temperature


def weighted_p95(values, weights):
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[min(np.searchsorted(cumulative, .95*cumulative[-1]),
                                 len(values)-1)]])


def distribution(values, volume, mass):
    return {"p95_volume": weighted_p95(values, volume),
            "p95_mass": weighted_p95(values, mass), "max": float(values.max())}


def tail(mask, volume, mass):
    return {"cells": int(mask.sum()), "volume_fraction": float(volume[mask].sum()/volume.sum()),
            "mass_fraction": float(mass[mask].sum()/mass.sum())}


def accuracy_controls(library, data_file):
    u = np.zeros((2, 27))
    u[:, 0] = 1; u[:, 4] = 100
    u[:, 6] = .76; u[:, 7] = .24
    u[:, 18] = .76; u[:, 21] = .24
    u[1, 18] = 0; u[1, 19] = .76; u[1, 26] = .76
    header = {"aexp": .1, "omega_m": .315, "h0": 67.4,
              "boxlen_ini": .001, "gamma": 5/3}
    eps, temp = grackle_temperature(u, header, library, data_file)
    eps_scaled, temp_scaled = grackle_temperature(u*3, header, library, data_file)
    assert np.allclose(eps, eps_scaled, rtol=1e-12, atol=0)
    assert np.allclose(temp, temp_scaled, rtol=1e-12, atol=0)
    assert temp[0] > temp[1] > 1
    # Independent atomic particle-count check: ionizing H adds one electron/H.
    assert abs(temp[1]/temp[0] - (.76+.24/4)/(2*.76+.24/4)) < 1e-12
    header2 = dict(header, aexp=.2)
    # Rescale energy separately: epsilon*a^2 is needed to keep physical epsilon
    # fixed as the RAMSES supercomoving velocity unit changes with epoch.
    u2 = u.copy(); u2[:, 4] *= 4
    _, temp2 = grackle_temperature(u2, header2, library, data_file)
    assert np.allclose(temp, temp2, rtol=1e-12, atol=0)
    volume, mass = np.array([1., 2.]), np.array([3., 1.])
    assert tail(np.array([True, False]) | np.array([False, True]), volume, mass)["cells"] == 2
    assert weighted_p95(np.array([0., 1.]), np.array([100., 1.])) == 0
    print("ACCURACY_ANALYZER_CONTROLS_PASS: density, composition, units, weights, union")


def accuracy_policy(root, library, data_file):
    states, diagnostics = {}, {}
    names = ("candidate", "reference")
    species_names = ("HI", "HII", "HM", "HeI", "HeII", "HeIII", "H2I", "H2II", "e_H")
    for name in names:
        outputs = list((root/name).glob("output_*"))
        assert [p.name for p in outputs] == ["output_00001"], (name, outputs)
        files = list(outputs[0].glob("data_*.h5"))
        assert len(files) == 1
        xyz, u, volume, header = read(files[0], details=True)
        assert abs(header["aexp"]-.1065) < 1e-10
        assert np.isfinite(u).all() and u[:, 0].min() > 0
        s, rho = u[:, 18:27], u[:, 0]
        h = s[:, 0]+s[:, 1]+s[:, 2]+s[:, 6]+s[:, 7]
        he = s[:, 3]+s[:, 4]+s[:, 5]
        assert s.min() >= 0 and h.min() > 0 and he.min() > 0
        eps, temp = grackle_temperature(u, header, library, data_file)
        assert np.isfinite(eps).all() and np.isfinite(temp).all()
        assert eps.min() > 0 and temp.min() > 0
        charge = s[:, 1]-s[:, 2]+s[:, 7]/2+s[:, 4]/4+s[:, 5]/2
        residuals = [float(np.max(np.abs(h-u[:, 6])/rho)),
                     float(np.max(np.abs(he-u[:, 7])/rho)),
                     float(np.max(np.abs(charge-s[:, 8])/rho))]
        diagnostics[name] = {"aexp": header["aexp"], "cells": len(u),
                             "max_H_He_charge_residuals": residuals,
                             "min_epsilon": float(eps.min()), "min_temperature": float(temp.min()),
                             "max_temperature": float(temp.max()),
                             "HII_H": float(np.sum(s[:, 1]*volume)/np.sum(h*volume)),
                             "HeIII_He": float(np.sum(s[:, 5]*volume)/np.sum(he*volume))}
        denominators = np.column_stack((h, h, h, he, he, he, h, h, h))
        states[name] = (xyz, rho, volume, eps, temp, s/denominators)
    c, r = states["candidate"], states["reference"]
    assert np.array_equal(c[0], r[0]) and np.array_equal(c[2], r[2])
    volume, mass = c[2], .5*(c[1]+r[1])*c[2]
    de, dt = np.abs(c[3]-r[3])/r[3], np.abs(c[4]-r[4])/r[4]
    ds = np.abs(c[5]-r[5])
    thermal = {"epsilon": distribution(de, volume, mass), "temperature": distribution(dt, volume, mass)}
    composition = {name: distribution(ds[:, i], volume, mass) for i, name in enumerate(species_names)}
    thermal_tail = tail((de > .10) | (dt > .10), volume, mass)
    composition_tail = tail(np.any(ds > .01, axis=1), volume, mass)
    global_deltas = {key: abs(diagnostics["candidate"][key]-diagnostics["reference"][key])
                     for key in ("HII_H", "HeIII_He")}
    checks = {
        "valid_state_inventory": all(max(d["max_H_He_charge_residuals"]) < 1e-6 for d in diagnostics.values()),
        "thermal_p95": all(d[w] <= .02 for d in thermal.values() for w in ("p95_volume", "p95_mass")),
        "composition_p95": all(d[w] <= .001 for d in composition.values() for w in ("p95_volume", "p95_mass")),
        "global_ionization": all(v <= 1e-4 for v in global_deltas.values()),
        "thermal_union_tail": all(thermal_tail[w] <= .001 for w in ("volume_fraction", "mass_fraction")),
        "composition_union_tail": all(composition_tail[w] <= .0001 for w in ("volume_fraction", "mass_fraction")),
    }
    result = {"scope": "approved 2026-10-07 CPU Grackle NEQ, fixed level-7 a=0.1065; not universal validation",
              "candidate_courant": .125, "reference_courant": .0625, "states": diagnostics,
              "thermal": thermal, "composition": composition, "global_ionization_deltas": global_deltas,
              "thermal_union_tail": thermal_tail, "composition_union_tail": composition_tail,
              "density_diagnostic": distribution(np.abs(c[1]-r[1])/r[1], volume, mass),
              "thermal_density_diagnostic": distribution(np.abs(c[3]*c[1]-r[3]*r[1])/(r[3]*r[1]), volume, mass),
              "checks": checks, "policy_pass": all(checks.values())}
    (root/"accuracy_summary.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    if not result["policy_pass"]:
        raise SystemExit("Approved accuracy criteria failed; retain snapshots and inspect, do not relax limits")


def main():
    root = Path(sys.argv[1])
    reference = read(root / "omp1/output_00004/data_00004.h5")
    result = {"comparisons": {}, "states": {}}
    names = ["omp1", "omp2", "mpi2", "restart"]
    if (root / "amr_fixed1/output_00004/data_00004.h5").exists():
        names += ["amr_fixed1", "amr_fixed2"]
    for name in names:
        if name == "amr_fixed1":
            reference = read(root / "amr_fixed1/output_00004/data_00004.h5")
            assert len(reference[1]) > 512, "AMR test did not refine"
        xyz, u, time = read(root / name / "output_00004/data_00004.h5")
        assert np.array_equal(xyz, reference[0]), name
        relative = float(np.max(np.abs(u-reference[1]) /
                                 np.maximum(np.max(np.abs(reference[1]), axis=0), 1e-30)))
        result["comparisons"][name] = {"max_field_scaled_difference": relative,
                                         "time_difference": time-reference[2]}
        assert relative < 1e-10 and abs(time-reference[2]) < 1e-12, (name, relative)
        s = u[:, 18:27]
        hydrogen = s[:, 0]+s[:, 1]+s[:, 2]+s[:, 6]+s[:, 7]
        helium = s[:, 3]+s[:, 4]+s[:, 5]
        charge = s[:, 1]-s[:, 2]+s[:, 7]/2+s[:, 4]/4+s[:, 5]/2
        hr = float(np.max(np.abs(hydrogen-u[:, 6])/u[:, 0]))
        her = float(np.max(np.abs(helium-u[:, 7])/u[:, 0]))
        er = float(np.max(np.abs(s[:, 8]-charge)/u[:, 0]))
        assert np.isfinite(u).all() and s.min() >= 0 and max(hr, her, er) < 1e-6
        internal = u[:, 4] - np.sum(u[:, 1:4]**2, axis=1)/(2*u[:, 0])
        assert internal.min() > 0
        result["states"][name] = {"cells": len(u), "hydrogen_residual": hr,
                                   "helium_residual": her, "charge_residual": er,
                                   "minimum_internal_energy_density": float(internal.min())}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--accuracy-controls":
        accuracy_controls(Path(sys.argv[2]), Path(sys.argv[3]))
    elif len(sys.argv) > 1 and sys.argv[1] == "--accuracy-policy":
        accuracy_policy(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
    else:
        main()
