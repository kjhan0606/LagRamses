#!/usr/bin/env python3
"""Measure H/He group-kernel variation across the active PARSEC Q/E model.

This is an offline source-shape analysis, not a runtime model or admission
gate. Measured PARSEC cumulative-Q rows are reconstructed with the same
Q-constrained Planck prior as the native source builder; the missing late
phase uses the builder's track-L/Teff Planck tail. Results are knot-sample
envelopes, not IMF- or time-weighted SSP averages.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
TOOLS = HERE.parent
sys.path.insert(0, str(TOOLS))
from build_parsec_native_sed import (  # noqa: E402
    EXTRA_PHOTONS,
    EXTRA_TRACKS,
    LOW_PHOTONS,
    LOW_TRACKS,
    PHOTONS,
    PINNED,
    Q_EDGES,
    Spectrum,
    source_name,
)
from build_parsec_pair_feedback import Z_GRIDS  # noqa: E402


EV = 1.602176634e-12
KB = 8.617333262145e-5
LSUN = 3.828e33
SPECIES = ("hydrogen_i", "helium_i", "helium_ii")
FITS = (
    # Verner et al. (1996), Table 1; same values as snrt_core.primordial.
    (13.60, 5.0e4, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0.0, 0.0, 0.0),
    (24.59, 5.0e4, 1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136),
    (54.42, 5.0e4, 1.720, 1.369e4, 3.288e1, 2.963, 0.0, 0.0, 0.0),
)
Q_COLUMNS = (9, 5, 6, 8, 7)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verner(energy: np.ndarray, fit: tuple[float, ...]) -> np.ndarray:
    threshold, maximum, e0, sigma0, ya, power, yw, y0, y1 = fit
    x = energy / e0 - y0
    y = np.sqrt(x * x + y1 * y1)
    profile = ((x - 1.0) ** 2 + yw**2) * y ** (0.5 * power - 5.5)
    profile *= (1.0 + np.sqrt(y / ya)) ** (-power)
    sigma = sigma0 * profile * 1.0e-18
    return np.where((energy >= threshold) & (energy <= maximum), sigma, 0.0)


class KernelIntegrator:
    def __init__(self, order: int, edges: np.ndarray):
        self.spectrum = Spectrum(order)
        self.edges = edges
        self.x, self.w = np.polynomial.legendre.leggauss(order)

    def _interval(self, kt: np.ndarray, low: float, high: float,
                  reference: float) -> tuple[np.ndarray, np.ndarray]:
        energy = low + (self.x[None, :] + 1.0) * (0.5 * (high - low))
        weight = self.w[None, :] * (0.5 * (high - low))
        thermal = kt[:, None]
        planck_photons = energy**2 / (-np.expm1(-energy / thermal))
        scaled_weight = weight * np.exp(-(energy - reference) / thermal)
        number = np.sum(scaled_weight * planck_photons, axis=1)
        sigma = np.empty((len(FITS), kt.size), dtype=np.float64)
        for species, fit in enumerate(FITS):
            sigma[species] = np.sum(
                scaled_weight * planck_photons * verner(energy[0], fit)[None, :],
                axis=1,
            )
        return number, sigma

    def constrained(self, data: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        """Return group Q and sigma*N moments for measured cumulative-Q rows."""
        kt = KB * np.power(10.0, data[:, 3])
        q_cumulative = np.power(10.0, data[:, Q_COLUMNS])
        q_disjoint = np.column_stack((q_cumulative[:, :-1] - q_cumulative[:, 1:],
                                      q_cumulative[:, -1]))
        if not np.isfinite(q_disjoint).all() or np.any(q_disjoint < 0.0):
            raise ValueError("non-monotone or non-finite measured cumulative Q")

        # Builder's group rates supply the exact 0.01--11.2 eV bolometric
        # residual closure and an independent denominator for the opacity.
        rates, _ = self.spectrum.constrained(data, self.edges)
        q_group = rates[:, :9]
        sigma_moment = np.zeros((len(FITS), len(data), 9), dtype=np.float64)
        reconstructed_q = np.zeros_like(q_group)

        for band, (band_low, band_high) in enumerate(zip(Q_EDGES[:-1], Q_EDGES[1:])):
            denominator, _ = self.spectrum.moment(kt, band_low, band_high, band_low)
            if np.any(denominator <= 0.0):
                raise ValueError("unresolved Planck normalization in Q band")
            scale = q_disjoint[:, band] / denominator
            for group, (group_low, group_high) in enumerate(zip(self.edges[:-1], self.edges[1:])):
                low, high = max(band_low, group_low), min(band_high, group_high)
                if high <= low:
                    continue
                number, sigma = self._interval(kt, low, high, band_low)
                reconstructed_q[:, group] += scale * number
                sigma_moment[:, :, group] += scale[None, :] * sigma

        # Above 11.2 eV, group counts must reproduce the generator's native
        # Q/E construction to quadrature precision before ratios are formed.
        active = q_group[:, 3:]
        rebuilt = reconstructed_q[:, 3:]
        relative = np.abs(active - rebuilt) / np.maximum(active, 1.0)
        q_error = float(np.max(relative))
        if q_error > 2.0e-10:
            raise ValueError(f"group-Q reconstruction mismatch {q_error:.3e}")
        return q_group, sigma_moment, q_error

    def planck_tail(self, logl: np.ndarray, logt: np.ndarray,
                    group_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return group Q and sigma*N for actual-track Planck-tail knots."""
        kt = KB * np.power(10.0, logt)
        luminosity_ev_s = np.power(10.0, logl) * LSUN / EV
        scale = luminosity_ev_s / (np.pi**4 * kt**4 / 15.0)
        rates = self.spectrum.planck(logl, logt, group_edges)[:, :9]
        sigma_moment = np.zeros((len(FITS), len(logl), 9), dtype=np.float64)
        for group, (low, high) in enumerate(zip(group_edges[:-1], group_edges[1:])):
            if high <= FITS[0][0]:
                continue
            a = max(low, FITS[0][0])
            number, sigma = self._interval(kt, a, high, a)
            attenuation = np.exp(-a / kt)
            sigma_moment[:, :, group] = (
                scale[None, :] * attenuation[None, :] * sigma
            )
            # `_interval` is scaled at a; q_group from Spectrum.planck uses
            # reference zero, so this is only used for sigma*N, whose scaling
            # above restores the physical Planck normalization.
        return rates, sigma_moment


def read_nodes(path: Path) -> tuple[str, list[dict[str, float]], int]:
    nodes = []
    with path.open() as f:
        if f.readline().strip() != "SNRT_POPULATION_QE_INTERVAL_V2":
            raise ValueError("unexpected active PARSEC node payload version")
        node_count, row_count = map(int, f.readline().split())
        model = f.readline().strip()
        consumed = 0
        for _ in range(node_count):
            fields = f.readline().split()
            if len(fields) != 8:
                raise ValueError("malformed PARSEC node header")
            mass, z, terminal_age_myr = map(float, fields[:3])
            nrow = int(fields[-1])
            nodes.append({"mass": mass, "z": z, "terminal_age_myr": terminal_age_myr,
                          "stored_rows": nrow})
            for _ in range(nrow):
                if not f.readline():
                    raise ValueError("truncated PARSEC node payload")
            consumed += nrow
        if consumed != row_count or f.read().strip():
            raise ValueError("PARSEC node extent/trailing-data mismatch")
    return model, nodes, row_count


def pilot_values(path: Path) -> dict[str, list[float]]:
    result = {}
    lines = path.read_text().splitlines()
    for species, key in (("hydrogen_i", "hi_input"), ("helium_i", "hei_input"),
                         ("helium_ii", "heii_input")):
        start = next((i for i, line in enumerate(lines)
                      if re.match(rf"^\s*{key}\s*=", line)), None)
        if start is None:
            raise ValueError(f"missing {key} in reference-control contract")
        first = re.sub(rf"^\s*{key}\s*=", "", lines[start], count=1).strip()
        tokens = [x.strip() for x in first.split(",") if x.strip()]
        cursor = start + 1
        while len(tokens) < 9 and cursor < len(lines):
            tokens.extend(x.strip() for x in lines[cursor].split(",") if x.strip())
            cursor += 1
        result[species] = [float(x) for x in tokens[:9]]
        if len(result[species]) != 9:
            raise ValueError(f"expected nine pilot kernels for {species}")
    return result


def add_samples(samples: dict, key: tuple[str, str, int], sigma: np.ndarray,
                support: np.ndarray, mass: float, z: float, ages: np.ndarray,
                source_mass: float, q_fraction: np.ndarray) -> None:
    values = np.asarray(sigma[support], dtype=np.float64)
    if values.size == 0:
        return
    indices = np.flatnonzero(support)
    states = [{"mass_msun": mass, "z": z, "age_myr": float(ages[i]),
               "source_mass_msun": source_mass,
               "qgroup_over_qtotal": float(q_fraction[i])} for i in indices]
    bucket = samples[key]
    bucket["values"].append(values)
    local_min, local_max = int(np.argmin(values)), int(np.argmax(values))
    if values[local_min] < bucket["minimum"]:
        bucket["minimum"] = float(values[local_min])
        bucket["minimum_state"] = states[local_min]
    if values[local_max] > bucket["maximum"]:
        bucket["maximum"] = float(values[local_max])
        bucket["maximum_state"] = states[local_max]
    bucket["count"] += int(values.size)


def summarize(samples: dict[tuple[str, str, int], dict], pilot: dict,
              phase: str) -> dict:
    out = {}
    for (kind, species, group), bucket in sorted(samples.items()):
        if kind != phase or not bucket["values"]:
            continue
        values = np.concatenate(bucket["values"])
        ref = pilot[species][group]
        ratios = values / ref if ref > 0.0 else np.full(values.shape, np.nan)
        minimum, maximum = bucket["minimum"], bucket["maximum"]
        out[f"{species}:g{group}"] = {
            "group_ev": [float(GROUP_EDGES[group]), float(GROUP_EDGES[group + 1])],
            "n_supported_knots": bucket["count"],
            "sigma_min_cm2": minimum,
            "sigma_max_cm2": maximum,
            "span_max_over_min": float(maximum / minimum),
            "sigma_percentiles_05_50_95_cm2": np.percentile(values, [5, 50, 95]).tolist(),
            "pilot_reference_cm2": float(ref),
            "ratio_to_pilot_min_max": [float(np.nanmin(ratios)), float(np.nanmax(ratios))],
            "minimum_state": bucket["minimum_state"],
            "maximum_state": bucket["maximum_state"],
        }
    return out


def quadrature_difference(q_ref: np.ndarray, sigma_ref: np.ndarray,
                          q_test: np.ndarray, sigma_test: np.ndarray,
                          edges: np.ndarray) -> tuple[float, dict | None]:
    total = q_ref.sum(axis=1)
    max_error, state = 0.0, None
    for species_index, fit in enumerate(FITS):
        for group in range(9):
            if edges[group + 1] <= fit[0]:
                continue
            support = (q_ref[:, group] / np.maximum(total, np.finfo(float).tiny) > 1.0e-8)
            support &= (q_ref[:, group] > 0.0)
            if not np.any(support):
                continue
            a = sigma_ref[species_index, support, group]
            b = sigma_test[species_index, support, group]
            rel = np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), np.finfo(float).tiny)
            index = int(np.argmax(rel))
            if float(rel[index]) > max_error:
                max_error = float(rel[index])
                original_index = int(np.flatnonzero(support)[index])
                state = {"species": SPECIES[species_index], "group": group,
                         "age_index": original_index}
    return max_error, state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--group-edges", type=Path, default=ROOT / "simulation/snrt/config/p0_photon_group_edges_ev.txt")
    parser.add_argument("--pilot-contract", type=Path, default=ROOT / "simulation/snrt/config/snrt_group_contract_reference_control_v1.nml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    package = args.package.resolve()
    package_manifest = json.loads((package / "manifest.json").read_text())
    model, nodes, stored_rows = read_nodes(package / "nodes.dat")
    if model != "parsec_mixed_lowmass_truncated_v1" or package_manifest.get("model") != model:
        raise ValueError("active runtime package is not the expected mixed PARSEC model")
    if len(nodes) != 979:
        raise ValueError(f"expected 979 active source nodes, found {len(nodes)}")

    edges = np.loadtxt(args.group_edges, dtype=np.float64)
    if edges.shape != (10,) or not np.all(np.diff(edges) > 0.0):
        raise ValueError("expected the exact increasing nine-group edge table")
    global GROUP_EDGES
    GROUP_EDGES = edges

    z_by_float = {float(z): (z, y) for z, y in Z_GRIDS["precision_eleven"]}
    expected = {**PINNED, **EXTRA_TRACKS, **LOW_TRACKS, **PHOTONS, **EXTRA_PHOTONS, **LOW_PHOTONS}
    input_hashes = {}
    for name, expected_sha in expected.items():
        path = args.source_dir / name
        actual = sha256(path)
        if actual != expected_sha:
            raise ValueError(f"PARSEC input checksum mismatch: {name}")
        input_hashes[name] = actual

    pilot = pilot_values(args.pilot_contract)
    primary = KernelIntegrator(128, edges)
    lower_order = KernelIntegrator(64, edges)
    quadrature_targets = {2.0, 8.0, 12.0, 13.0, 14.0, 20.0, 50.0, 80.0,
                          120.0, 200.0, 300.0, 600.0}
    quadrature_masses = {n["mass"] for n in nodes
                         if any(abs(n["mass"] - target) < 1.0e-9 for target in quadrature_targets)}
    quadrature_max = {"measured_q_constrained": {"relative_error": 0.0, "state": None},
                      "track_planck_tail": {"relative_error": 0.0, "state": None}}
    def new_bucket():
        return {"values": [], "count": 0, "minimum": float("inf"),
                "maximum": -float("inf"), "minimum_state": None, "maximum_state": None}

    samples: dict[tuple[str, str, int], dict] = defaultdict(new_bucket)
    qclosure_max = 0.0
    processed = 0

    photon_archives = {}
    track_archives = {}
    for zs, ys in Z_GRIDS["precision_eleven"]:
        z = float(zs)
        photon_archives[z] = zipfile.ZipFile(args.source_dir / source_name(zs, ys, "photons"))
        track_archives[z] = zipfile.ZipFile(args.source_dir / source_name(zs, ys, "tracks"))

    try:
        for node_index, node in enumerate(nodes, 1):
            mass, z = node["mass"], node["z"]
            if z not in z_by_float:
                raise ValueError(f"unexpected metallicity node {z}")
            zs, ys = z_by_float[z]
            source_mass = min(mass, 12.0) if mass < 14.0 else mass
            photon_bytes = photon_archives[z].read(source_name(zs, ys, "photons", source_mass))
            track_bytes = track_archives[z].read(source_name(zs, ys, "tracks", source_mass))
            data = np.loadtxt(io.BytesIO(photon_bytes), skiprows=3)
            track = np.loadtxt(io.BytesIO(track_bytes), skiprows=3)
            data = data[np.r_[np.diff(data[:, 0]) > 0.0, True]]
            if not np.isfinite(data).all() or not np.isfinite(track).all():
                raise ValueError(f"non-finite PARSEC source rows at node {node_index}")

            terminal_myr = node["terminal_age_myr"]
            selected = data[:, 0] <= terminal_myr
            measured = data[selected]
            if len(measured) == 0:
                raise ValueError(f"no measured PARSEC source rows before terminal age at {node}")
            q_group, sigma_moment, q_error = primary.constrained(measured)
            qclosure_max = max(qclosure_max, q_error)
            q_total = q_group.sum(axis=1)
            if mass in quadrature_masses:
                q64, s64, _ = lower_order.constrained(measured)
                error, location = quadrature_difference(q_group, sigma_moment, q64, s64, edges)
                if error > quadrature_max["measured_q_constrained"]["relative_error"]:
                    quadrature_max["measured_q_constrained"] = {
                        "relative_error": error,
                        "state": {"mass_msun": mass, "z": z, **(location or {})},
                    }
            for species_index, species in enumerate(SPECIES):
                sigma = np.divide(
                    sigma_moment[species_index], q_group,
                    out=np.zeros_like(q_group), where=q_group > 0.0,
                )
                threshold = FITS[species_index][0]
                for group in range(9):
                    if edges[group + 1] <= threshold:
                        continue
                    support = (q_group[:, group] > 0.0) & (
                        q_group[:, group] / np.maximum(q_total, np.finfo(float).tiny) > 1.0e-8
                    )
                    add_samples(samples, ("measured_q_constrained", species, group),
                                sigma[:, group], support, mass, z, measured[:, 0],
                                source_mass, q_group[:, group] / np.maximum(q_total, np.finfo(float).tiny))

            # Match the active builder's missing-late-phase closure: only when
            # material outlives the final measured Q row, use full-track L/Teff.
            measured_end = float(data[-1, 0])
            if terminal_myr > measured_end:
                track_age = track[:, 1] / 1.0e6
                tail_age = np.unique(np.r_[measured_end,
                    track_age[(track_age > measured_end) & (track_age < terminal_myr)],
                    terminal_myr])
                logt = np.interp(tail_age, track_age, track[:, 4])
                logl = np.interp(tail_age, track_age, track[:, 3])
                tail_q, tail_sigma_moment = primary.planck_tail(logl, logt, edges)
                tail_total = tail_q.sum(axis=1)
                if mass in quadrature_masses:
                    q64, s64 = lower_order.planck_tail(logl, logt, edges)
                    error, location = quadrature_difference(tail_q, tail_sigma_moment, q64, s64, edges)
                    if error > quadrature_max["track_planck_tail"]["relative_error"]:
                        quadrature_max["track_planck_tail"] = {
                            "relative_error": error,
                            "state": {"mass_msun": mass, "z": z, **(location or {})},
                        }
                for species_index, species in enumerate(SPECIES):
                    tail_sigma = np.divide(
                        tail_sigma_moment[species_index], tail_q,
                        out=np.zeros_like(tail_q), where=tail_q > 0.0,
                    )
                    threshold = FITS[species_index][0]
                    for group in range(9):
                        if edges[group + 1] <= threshold:
                            continue
                        support = (tail_q[:, group] > 0.0) & (
                            tail_q[:, group] / np.maximum(tail_total, np.finfo(float).tiny) > 1.0e-8
                        )
                        add_samples(samples, ("track_planck_tail", species, group),
                                    tail_sigma[:, group], support, mass, z, tail_age,
                                    source_mass, tail_q[:, group] / np.maximum(tail_total, np.finfo(float).tiny))
            processed += 1
            if node_index % 100 == 0 or node_index == len(nodes):
                print(f"processed_nodes={node_index}/{len(nodes)}", flush=True)
    finally:
        for archive in (*photon_archives.values(), *track_archives.values()):
            archive.close()

    result = {
        "schema": "parsec_shared_group_kernel_sensitivity_v1",
        "model": model,
        "model_status": package_manifest.get("status"),
        "package_nodes_sha256": sha256(package / "nodes.dat"),
        "package_manifest_sha256": sha256(package / "manifest.json"),
        "source_inputs_sha256": input_hashes,
        "group_edges_ev": edges.tolist(),
        "group_interval_convention": "left_closed_right_open_except_final_closed",
        "photoionization_fits": "Verner et al. 1996 ground-state H I/He I/He II; constants match snrt_core.primordial",
        "quadrature_order": 128,
        "quadrature_convergence_check": {
            "comparison_orders": [64, 128],
            "source_masses_msun": sorted(quadrature_masses),
            "all_metallicity_nodes_included": True,
            "maximum_relative_difference": quadrature_max,
        },
        "support_cut": "Q_group/Q_total(0.01eV-10keV) > 1e-8",
        "sample_weighting": "unweighted instantaneous source knots; no IMF or time weighting; not an SSP mean",
        "source_shape": {
            "measured_phase": "PARSEC cumulative-Q bands reconstructed with build_parsec_native_sed.Spectrum.constrained Planck prior",
            "late_phase": "actual PARSEC track L/Teff Planck closure used by active source builder; reported separately",
            "low_mass_proxy": "13 Msun radiation uses 12 Msun spectrum; uniform mass scaling cancels in sigma_N",
        },
        "active_node_count": len(nodes),
        "active_knot_count_in_package": stored_rows,
        "processed_node_count": processed,
        "pilot_reference_status": "reference_control_not_science_approved",
        "pilot_reference_contract_sha256": sha256(args.pilot_contract),
        "pilot_cross_sections_cm2": pilot,
        "maximum_Q_reconstruction_relative_error": qclosure_max,
        "analysis_code_sha256": sha256(HERE),
        "source_builder_sha256": sha256(TOOLS / "build_parsec_native_sed.py"),
        "production_verner_implementation_sha256": sha256(ROOT / "simulation/snrt/snrt_core/primordial.py"),
        "measured_q_constrained": summarize(samples, pilot, "measured_q_constrained"),
        "track_planck_tail": summarize(samples, pilot, "track_planck_tail"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote={args.output}", flush=True)
    print(f"processed_nodes={processed}; stored_package_knots={stored_rows}", flush=True)


if __name__ == "__main__":
    main()
