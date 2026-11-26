#!/usr/bin/env python3
"""Measure how provisional z=0 M200c candidates depend on halo centering."""
import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree


snapshot, catalog_path, tags_path, m200_path, trace_summary_path, out = map(Path, sys.argv[1:7])
K_DENSITY = 64
CENTROID_NEIGHBORS = (64, 128, 256)
QUERY_RADIUS_MPC = 1.0


def periodic_delta(xyz, center):
    delta = xyz - center
    delta -= np.rint(delta)
    return delta


def m200c_at_center(xyz, center, box_mpc, particle_mass, rho_crit):
    delta = periodic_delta(xyz, center)
    radius_mpc = np.sqrt(np.einsum("ij,ij->i", delta, delta)) * box_mpc
    central_count = np.count_nonzero(radius_mpc == 0)
    radius_mpc = np.sort(radius_mpc[(radius_mpc > 0.0) & (radius_mpc <= QUERY_RADIUS_MPC)])
    if len(radius_mpc) < 2:
        raise ValueError("too few particles inside M200c query radius")
    enclosed_mass = (central_count + np.arange(1, len(radius_mpc) + 1, dtype=np.float64)) * particle_mass
    mean_density = enclosed_mass / ((4.0 / 3.0) * math.pi * radius_mpc**3)
    crossing = np.flatnonzero(mean_density <= 200.0 * rho_crit)
    if not len(crossing):
        raise ValueError("no 200c crossing within query radius")
    j = int(crossing[0])
    if j == 0:
        r200 = float(radius_mpc[0])
    else:
        log_r0, log_r1 = np.log(radius_mpc[j - 1:j + 1])
        log_d0, log_d1 = np.log(mean_density[j - 1:j + 1])
        fraction = (math.log(200.0 * rho_crit) - log_d0) / (log_d1 - log_d0)
        r200 = float(math.exp(log_r0 + fraction * (log_r1 - log_r0)))
    mass = float((4.0 / 3.0) * math.pi * 200.0 * rho_crit * r200**3)
    return {"m200c_msun": mass, "r200c_mpc": r200, "n200c_equiv": mass / particle_mass,
            "crossing_bracketed": bool(j > 0)}


def main():
    if not out.is_dir():
        raise ValueError("Slurm runner must create the output directory")
    trace_summary = json.loads(trace_summary_path.read_text())
    if trace_summary.get("status") != "PRELIMINARY_Z0_LAGRANGIAN_TRACE_NO_TARGET_APPROVAL":
        raise ValueError("unexpected candidate trace status")
    groups = [int(item["group_id"]) for item in trace_summary["candidates"]]
    if not groups or len(set(groups)) != len(groups):
        raise ValueError("empty or duplicate trace candidates")

    with catalog_path.open(newline="") as stream:
        hop_rows = {int(row["group_id"]): row for row in csv.DictReader(stream)}
    with m200_path.open(newline="") as stream:
        m200_rows = {int(row["group_id"]): row for row in csv.DictReader(stream)}
    if not set(groups).issubset(hop_rows) or not set(groups).issubset(m200_rows):
        raise ValueError("candidate missing from HOP/M200 catalogue")
    with h5py.File(snapshot, "r") as handle:
        attrs = handle["header"].attrs
        aexp = float(np.asarray(attrs["aexp"]).item())
        box = float(np.asarray(attrs["boxlen"]).item())
        box_mpc_h = float(np.asarray(attrs["boxlen_ini"]).item())
        h = float(np.asarray(attrs["h0"]).item()) / 100.0
        omega_m = float(np.asarray(attrs["omega_m"]).item())
        p = handle["particles"]
        count = len(p["mass"])
        if (abs(aexp - 1.0) > 1e-8 or abs(box - 1.0) > 1e-12
                or abs(box_mpc_h - 25.0) > 1e-4 or abs(h - 0.6736) > 1e-4
                or abs(omega_m - 0.315) > 1e-6 or count != 256**3):
            raise ValueError("unexpected parent snapshot")
        xyz = np.column_stack([p[f"x_{axis}"][:] for axis in (1, 2, 3)]).astype(np.float64, copy=False)
        if not np.isfinite(xyz).all() or np.any(xyz < 0.0) or np.any(xyz >= 1.0):
            raise ValueError("invalid parent positions")
        mass_values = np.asarray(p["mass"][:])
        if not np.isfinite(mass_values).all() or not np.all(mass_values == 1.0 / count):
            raise ValueError("unexpected one-fluid particle masses")

    # HOP's Fortran tag file is a sequential binary record, read separately.
    from scipy.io import FortranFile
    with FortranFile(tags_path, "r") as stream:
        tag_header = stream.read_ints(np.int32)
        tag = stream.read_ints(np.int32)
    if len(tag_header) != 2 or len(tag) != count or int(tag_header[0]) != count:
        raise ValueError("HOP tag dimensions do not match the snapshot")

    rho_crit = 2.77536627e11 * h**2
    box_mpc = box_mpc_h / h
    particle_mass = omega_m * rho_crit * box_mpc**3 / count
    outputs = []
    for gid in groups:
        hop = hop_rows[gid]
        ref = m200_rows[gid]
        member_mask = tag == gid
        member_xyz = xyz[member_mask]
        if len(member_xyz) != int(hop["particle_count"]) or len(member_xyz) < K_DENSITY:
            raise ValueError(f"wrong member inventory for group {gid}")
        hop_center = np.array([float(hop[k]) for k in ("x_box", "y_box", "z_box")])
        tree = cKDTree(member_xyz, boxsize=1.0, compact_nodes=True, balanced_tree=True)
        neighbor_distances, _ = tree.query(member_xyz, k=K_DENSITY, workers=1)
        r64 = neighbor_distances[:, -1]
        peak_index = int(np.argmin(r64))
        peak_particle = member_xyz[peak_index]
        center_defs = {"hop_circular_mean": hop_center.tolist(), "densest_member_particle_64nn": peak_particle.tolist()}
        for k in CENTROID_NEIGHBORS:
            d, indices = tree.query(peak_particle, k=k, workers=1)
            local = member_xyz[np.asarray(indices, dtype=np.int64)]
            center_defs[f"density_peak_{k}nn_centroid"] = ((peak_particle + periodic_delta(local, peak_particle).mean(axis=0)) % 1.0).tolist()

        # A superset of every recentered 1 Mpc sphere, including periodic images.
        pad = max(np.linalg.norm(periodic_delta(np.array(c), hop_center)) for c in center_defs.values())
        delta = periodic_delta(xyz, hop_center)
        local_xyz = xyz[np.einsum("ij,ij->i", delta, delta) <= (QUERY_RADIUS_MPC/box_mpc + pad + 1e-12)**2]
        del delta
        centers = []
        for name, values in center_defs.items():
            center = np.asarray(values, dtype=np.float64)
            offset_box = float(np.linalg.norm(periodic_delta(center[None, :], hop_center)[0]))
            profile = m200c_at_center(local_xyz, center, box_mpc, particle_mass, rho_crit)
            centers.append({
                "center_method": name,
                "center_box": center.tolist(),
                "offset_from_HOP_center_box": offset_box,
                "offset_in_provisional_R200c": offset_box / (float(ref["r200c_mpc"]) / box_mpc),
                **profile,
                "fractional_M200c_change_from_HOP_center": profile["m200c_msun"] / float(ref["m200c_msun"]) - 1.0,
            })
        hop_recomputed = centers[0]["m200c_msun"]
        if abs(hop_recomputed / float(ref["m200c_msun"]) - 1.0) > 0.01:
            raise ValueError(f"M200c implementation does not reproduce screening result for {gid}")
        outputs.append({
            "group_id": gid,
            "hop_particles": int(hop["particle_count"]),
            "hop_mass_msun": float(hop["hop_mass_msun"]),
            "reference_M200c_msun": float(ref["m200c_msun"]),
            "reference_R200c_mpc": float(ref["r200c_mpc"]),
            "r64_density_peak_box": float(r64[peak_index]),
            "centers": centers,
        })

    result = {
        "status": "CENTER_SENSITIVITY_DIAGNOSTIC_NO_HALO_APPROVAL",
        "snapshot": str(snapshot),
        "aexp": aexp,
        "candidate_groups": groups,
        "density_proxy": "minimum 64th-neighbor radius among HOP member particles; centroids use periodic offsets around that densest member particle",
        "M200c_method": "same all-particle spherical-overdensity profile and 200*critical-density crossing used by the provisional screening catalogue, with a 1 Mpc query radius",
        "candidates": outputs,
        "limits": [
            "A 64-neighbor density peak is not a potential minimum or an unbound halo center.",
            "Center and M200c sensitivity at this parent resolution do not measure concentration or establish a production zoom target.",
            "No gas, galaxy counterpart, high-resolution IC, or contamination rerun is included.",
        ],
    }
    (out / "center_sensitivity.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
