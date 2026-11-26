#!/usr/bin/env python3
"""Trace provisional z=0 LCDM halo candidates to their parent IC cells.

This is a zoom-feasibility diagnostic, not a halo-selection or contamination
approval.  Candidate centers and M200c values inherit the limitations of the
low-resolution HOP/M200c screening catalogue.
"""
import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.io import FortranFile
from scipy.optimize import brentq


snapshot, catalog_path, tags_path, ids_path, m200_path, m200_summary_path, map_path, map_audit_path, out = map(Path, sys.argv[1:10])
N = 256
MAX_CANDIDATES = 12
RADIUS_MULTIPLIERS = (2.0, 3.0)


def scalar(attrs, key):
    return float(np.asarray(attrs[key]).item())


def periodic_extent(values):
    values = np.sort(np.mod(np.asarray(values, dtype=np.float64), 1.0))
    if len(values) == 0:
        raise ValueError("empty periodic coordinate")
    if len(values) == 1:
        return float(values[0]), 0.0
    gaps = np.diff(np.r_[values, values[0] + 1.0])
    largest = int(np.argmax(gaps))
    start = values[(largest + 1) % len(values)]
    width = 1.0 - gaps[largest]
    return float(start), float(width)


def tag_catalog(path, n):
    with FortranFile(path, "r") as stream:
        header = stream.read_ints(np.int32)
        tags = stream.read_ints(np.int32)
    if len(header) != 2 or int(header[0]) != n or len(tags) != n:
        raise ValueError("HOP tag dimensions do not match the parent")
    ngroups = int(header[1])
    if np.any((tags >= ngroups) | (tags < -1)):
        raise ValueError("invalid HOP group tag")
    counts = np.bincount(tags[tags >= 0], minlength=ngroups)
    return tags, counts, ngroups


def photometric_reference(path, parent_h0, distance):
    # Express the same physical DMO NFW density profile at the parent's
    # critical density. This conversion does not change the SHMR calibration.
    with path.open(newline="") as stream:
        row, = [r for r in csv.DictReader(stream)
                if float(r["distance_mpc"]) == distance and
                float(r["stellar_offset_dex"]) == float(r["c_offset_dex"]) == 0]
    mass, radius, c = (float(row[k]) for k in ("dmo_m200c_msun","r200c_kpc","concentration"))
    def profile(x):
        return math.log1p(x) - x / (1+x)
    def ratio(x):
        return profile(c*x)/profile(c)
    wmap_h0 = 70.4
    x = brentq(lambda v: ratio(v) - (parent_h0/wmap_h0)**2 * v**3, 0.5, 2)
    converted = mass * ratio(x)
    g_kpc = 4.30091727003628e-6
    if not math.isclose(converted, 100*(parent_h0/1000)**2*(radius*x)**3/g_kpc, rel_tol=1e-9):
        raise ValueError("converted mass and radius do not satisfy the parent 200c definition")
    return dict(source="NGC 2403 Kroupa photometry and Moster z=0 relation",
                original_m200c_msun=mass, original_r200c_kpc=radius, original_h0=wmap_h0,
                original_concentration=c, parent_h0=parent_h0, r200c_kpc=radius*x,
                comparison_m200c_msun=converted, reference_distance_mpc=distance)


def main():
    if not out.is_dir():
        raise ValueError("Slurm runner must create a new output directory first")
    parent_summary = json.loads(m200_summary_path.read_text())
    map_audit = json.loads(map_audit_path.read_text())
    if parent_summary.get("status") != "PROVISIONAL_Z0_LCDM_200C_SCREENING_CATALOGUE":
        raise ValueError("unexpected M200c screening catalogue")
    if map_audit.get("unique_lagrangian_cells") != N**3 or map_audit.get("duplicate_lagrangian_cells") != 0:
        raise ValueError("upstream ID-to-IC-cell map failed its uniqueness audit")

    hop_rows = list(csv.DictReader(catalog_path.open(newline="")))
    hop = {int(row["group_id"]): row for row in hop_rows}
    if len(hop) != len(hop_rows):
        raise ValueError("duplicate HOP group IDs")
    m200_rows = list(csv.DictReader(m200_path.open(newline="")))
    target_mass = float(parent_summary["ngc2403_reference_m200c_msun"])
    reference = {"source": "conditional analytic NFW rotation-curve fit",
                 "comparison_m200c_msun": target_mass}
    if len(sys.argv) > 10:
        reference = photometric_reference(Path(sys.argv[10]), parent_summary["parent_h0_km_s_mpc"],
                                          float(sys.argv[11]))
        target_mass = reference["comparison_m200c_msun"]
    eligible = []
    for row in m200_rows:
        if row["status"] != "measured_but_provisional" or row["n200c_500_particle_flag"].lower() != "true":
            continue
        mass = float(row["m200c_msun"])
        if not math.isfinite(mass) or mass <= 0:
            continue
        gid = int(row["group_id"])
        if gid not in hop:
            raise ValueError(f"M200c candidate {gid} is missing from HOP catalogue")
        eligible.append((abs(math.log10(mass / target_mass)), gid, row))
    eligible.sort(key=lambda item: (item[0], item[1]))
    candidates = eligible[:MAX_CANDIDATES]
    if len(candidates) < 3:
        raise ValueError("fewer than three resolved mass-screen candidates")

    with h5py.File(snapshot, "r") as handle:
        attrs = handle["header"].attrs
        aexp = scalar(attrs, "aexp")
        box_code = scalar(attrs, "boxlen")
        box_mpc_h = scalar(attrs, "boxlen_ini")
        h = scalar(attrs, "h0") / 100.0
        omega_m = scalar(attrs, "omega_m")
        particles = handle["particles"]
        count = len(particles["identity"])
        if (abs(aexp - 1.0) > 1e-8 or abs(box_code - 1.0) > 1e-12
                or abs(box_mpc_h - 25.0) > 1e-4 or abs(h - 0.6736) > 1e-4
                or abs(omega_m - 0.315) > 1e-6 or count != N**3):
            raise ValueError("unexpected z=0 parent or cosmology")

        row_ids = np.load(ids_path, mmap_mode="r")
        if row_ids.shape != (count,):
            raise ValueError("HOP row-to-ID map has wrong length")
        ids = np.asarray(particles["identity"][:], dtype=np.int64)
        if ids.min() != 1 or ids.max() != count:
            raise ValueError("parent IDs outside the verified map")
        if not np.array_equal(ids, row_ids):
            raise ValueError("HOP tag row order does not match the z=0 snapshot")
        expected_mass = 1.0 / count
        for start in range(0, count, 1_048_576):
            stop = min(count, start + 1_048_576)
            mass = np.asarray(particles["mass"][start:stop])
            if not np.isfinite(mass).all() or not np.all(mass == expected_mass):
                raise ValueError("unexpected one-fluid parent particle masses")
        xyz = np.column_stack([particles[f"x_{axis}"][:] for axis in (1, 2, 3)]).astype(np.float64, copy=False)
        if not np.isfinite(xyz).all() or np.any(xyz < 0.0) or np.any(xyz >= 1.0):
            raise ValueError("invalid periodic parent coordinates")

    tags, tag_counts, ngroups = tag_catalog(tags_path, count)
    cell_by_id = np.load(map_path, mmap_mode="r")
    if cell_by_id.shape != (count, 3) or cell_by_id.dtype != np.uint8:
        raise ValueError("unexpected ID-to-IC-cell map format")
    box_mpc = box_mpc_h / h
    particle_mass_msun = omega_m * (2.77536627e11 * h**2) * box_mpc**3 / count

    # HOP-neighbour distances are only a coarse environment descriptor; HOP
    # member masses are not spherical-overdensity masses or an isolation test.
    hop_centers = np.array([[float(row[k]) for k in ("x_box", "y_box", "z_box")] for row in hop_rows])
    hop_ids = np.array([int(row["group_id"]) for row in hop_rows])
    hop_mass = np.array([float(row["hop_mass_msun"]) for row in hop_rows])
    hop_count = np.array([int(row["particle_count"]) for row in hop_rows])
    saved = {}
    result_candidates = []

    for rank, (mass_delta, gid, mrow) in enumerate(candidates, start=1):
        hrow = hop[gid]
        if int(hrow["particle_count"]) != int(tag_counts[gid]):
            raise ValueError(f"HOP catalogue count mismatch for group {gid}")
        center = np.array([float(hrow[k]) for k in ("x_box", "y_box", "z_box")])
        hop_member = tags == gid
        member_ids = ids[hop_member]
        if len(member_ids) != int(hrow["particle_count"]):
            raise ValueError(f"HOP member identity mismatch for group {gid}")

        neighbor_mask = (hop_ids != gid) & (hop_count >= 500) & (hop_mass >= 0.5 * float(hrow["hop_mass_msun"]))
        neighbor = None
        if np.any(neighbor_mask):
            delta = hop_centers[neighbor_mask] - center
            delta -= np.rint(delta)
            distances = np.sqrt(np.einsum("ij,ij->i", delta, delta))
            j = int(np.argmin(distances))
            neighbor = {
                "group_id": int(hop_ids[neighbor_mask][j]),
                "hop_mass_msun": float(hop_mass[neighbor_mask][j]),
                "center_distance_mpc_h": float(distances[j] * box_mpc_h),
                "definition": "nearest listed HOP center with >=500 members and >=half this candidate's HOP mass; not an isolation criterion",
            }

        delta = xyz - center
        delta -= np.rint(delta)
        distance2 = np.einsum("ij,ij->i", delta, delta)
        r200_mpc = float(mrow["r200c_mpc"])
        r200_box = r200_mpc / box_mpc
        candidate_record = {
            "rank_by_M200c_proximity": rank,
            "group_id": gid,
            "hop_particles": int(hrow["particle_count"]),
            "hop_mass_msun": float(hrow["hop_mass_msun"]),
            "center_box": center.tolist(),
            "center_definition": "HOP periodic circular mean",
            "provisional_M200c_msun": float(mrow["m200c_msun"]),
            "provisional_R200c_mpc": r200_mpc,
            "equivalent_particles_inside_M200c": float(mrow["n200c_equiv"]),
            "log10_mass_offset_from_NGC2403_NFW_reference_dex": float(mass_delta),
            "nearest_HOP_neighbor": neighbor,
            "traces": [],
        }
        for multiple in RADIUS_MULTIPLIERS:
            inside = distance2 <= (multiple * r200_box) ** 2
            selected_ids = np.union1d(ids[inside], member_ids).astype(np.int64, copy=False)
            cells = np.asarray(cell_by_id[selected_ids - 1], dtype=np.uint8)
            if len(np.unique(cells, axis=0)) != len(cells):
                raise ValueError(f"duplicate initial cells in group {gid}, radius {multiple}R200c")
            points = (cells.astype(np.float64) + 0.5) / N
            start_width = [periodic_extent(points[:, axis]) for axis in range(3)]
            starts = [item[0] for item in start_width]
            widths = [item[1] for item in start_width]
            key = f"g{gid}_r{int(multiple)}"
            saved[f"{key}_ids"] = selected_ids
            saved[f"{key}_cell_xyz"] = cells
            candidate_record["traces"].append({
                "radius_multiple_R200c": multiple,
                "sphere_particle_count": int(np.count_nonzero(inside)),
                "HOP_members_outside_sphere": int(np.count_nonzero(hop_member & ~inside)),
                "union_particle_count": int(len(selected_ids)),
                "unique_parent_IC_cells": int(len(cells)),
                "parent_IC_cell_fraction": float(len(cells) / N**3),
                "periodic_minimum_bbox_start_xyz": starts,
                "periodic_minimum_bbox_width_xyz": widths,
                "bbox_volume_fraction": float(np.prod(widths)),
                "saved_array_prefix": key,
            })
        result_candidates.append(candidate_record)

    np.savez_compressed(out / "candidate_traces.npz", **saved)
    summary = {
        "status": "PRELIMINARY_Z0_LAGRANGIAN_TRACE_NO_TARGET_APPROVAL",
        "parent_snapshot": str(snapshot),
        "parent_scale_factor": aexp,
        "parent_IC": "same 256^3 one-fluid GRAFIC input used by the verified ID-to-cell mapping",
        "particle_mass_msun": particle_mass_msun,
        "parent_IC_grid": [N, N, N],
        "reference": reference,
        "eligible_parent_candidates": len(eligible),
        "candidate_selection": f"Top {len(candidates)} resolved provisional M200c candidates by absolute log-mass distance to " + reference["source"],
        "trace_definition": "Union of HOP members and every z=0 parent particle within 2 or 3 provisional R200c from the HOP circular-mean center; IDs are mapped to inferred initial grid cells.",
        "candidates": result_candidates,
        "limitations": [
            "Parent is a 256^3, 25 Mpc/h one-fluid DMO screening run; centers and M200c values are provisional.",
            "HOP circular-mean centers are not potential minima; members are not unbound-cleaned.",
            "The ID-to-cell map is inferred from the same IC's early output, not directly recorded particle provenance.",
            "A 2-3 R200c z=0 trace is a feasibility region only; it does not include a validated dynamical buffer or prove low-resolution contamination control.",
            "No high-resolution IC, rerun, galaxy counterpart, gas, or RAR prediction is produced by this analysis.",
        ],
        "hop_groups_in_tag_file": ngroups,
        "saved_arrays": "candidate_traces.npz stores particle IDs and inferred initial-cell coordinates for each candidate and radius.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
