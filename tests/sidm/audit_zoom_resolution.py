#!/usr/bin/env python3
"""Bounded IC-header, unit and cached-particle audit; never reads IC payloads.

The count scenarios are not regenerated ICs or throughput forecasts. Legacy
headers and cached QC cannot release production or certify mass convergence.
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
import io
import json
import math
import re
import struct
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
ZOOM = Path("/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11")
VALIDATION = Path("/gpfs/kjhan/Hydro/Sidm/validation")
QC = PROJECT / "runs/lageunha/snapshot-2026-08-24/Agn/Zoom0/validation/output_00003_particle_qc.json"
FIELDS = ("ic_deltab", "ic_velcx", "ic_velcy", "ic_velcz", "ic_poscx",
          "ic_poscy", "ic_poscz", "ic_refmap", "ic_pvar_00001")
MSUN_G = 1.98847e33
MPC_CM = 3.0856775814913673e24
CODE_MPC_CM = 3.08e24  # patch/cuRamses/units.f90, not the precise physical Mpc


def decode_header(stream, label: str = "header") -> dict:
    marker = stream.read(4)
    if len(marker) != 4:
        raise ValueError(f"{label}: missing header marker")
    candidates = [(endian, struct.unpack(endian + "i", marker)[0])
                  for endian in ("<", ">")]
    supported = [(e, n) for e, n in candidates if n in (44, 48)]
    if len(supported) != 1:
        raise ValueError(f"{label}: unsupported record marker")
    endian, length = supported[0]
    payload, trailer = stream.read(length), stream.read(4)
    if (len(payload) != length or len(trailer) != 4 or
            struct.unpack(endian + "i", trailer)[0] != length):
        raise ValueError(f"{label}: truncated or inconsistent header")
    values = struct.unpack(endian + ("3i8f" if length == 44 else "3i9f"), payload)
    dims = list(values[:3])
    dx, x, y, z, a, om, ol, h0 = values[3:11]
    if (any(n <= 0 for n in dims) or not all(math.isfinite(v) for v in values[3:]) or
            dx <= 0 or not 0 < a <= 1 or om <= 0 or h0 <= 0):
        raise ValueError(f"{label}: invalid header quantities")
    ob = values[11] if length == 48 else None
    if ob is not None and not 0 <= ob <= om:
        raise ValueError(f"{label}: invalid omega_b")
    return {"dimensions": dims, "cell_comoving_mpc": dx,
            "offset_comoving_mpc": [x, y, z], "a_start": a,
            "omega_m": om, "omega_l": ol, "H0": h0, "omega_b": ob,
            "header_record_bytes": length, "bytes_read": length + 8,
            "header_sha256": hashlib.sha256(marker + payload + trailer).hexdigest()}


def particle_scenarios(qc: dict, finest_level: int) -> dict:
    tiers = sorted(qc["mass_tiers"], key=lambda item: item["mass_code"])
    counts = [item["npart"] for item in tiers]
    if any(not isinstance(n, int) or n <= 0 for n in counts):
        raise ValueError("invalid cached particle counts")
    if sum(counts) != qc["npart"] or counts[0] != qc["highres_count"]:
        raise ValueError("cached particle count mismatch")
    expected = 2.0**(-3 * finest_level)
    if qc["highres_mass_code"] != expected:
        raise ValueError("cached highest mass differs from IC-level mass")
    for index, tier in enumerate(tiers):
        if tier["mass_code"] != expected * 8**index:
            raise ValueError("cached masses are not the expected factor-eight hierarchy")
    mass_sum = math.fsum(t["mass_code"] * t["npart"] for t in tiers)
    if not math.isclose(mass_sum, qc["total_mass_code"], rel_tol=1e-12, abs_tol=0):
        raise ValueError("cached tier masses do not sum to total mass")
    if not math.isclose(mass_sum, 1.0, rel_tol=1e-12, abs_tol=0):
        raise ValueError("expected unit total-matter DMO mass")
    total, high, outer = qc["npart"], counts[0], qc["npart"] - counts[0]
    return {
        "cached_total_particles": total, "cached_finest_particles": high,
        "cached_other_particles": outer, "cached_base_particles": counts[-1],
        "finest_particle_fraction": high / total,
        "base_particle_fraction": counts[-1] / total,
        "finest_only_eight_to_one_count_scenario": {
            "new_total_particles": outer + high / 8,
            "old_over_new_particle_count": total / (outer + high / 8),
            "new_total_fraction": (outer + high / 8) / total,
            "assumption": "only finest particles grouped eight-to-one; all other particles unchanged",
            "valid_regenerated_mask_proven": False,
        },
        "fixed_base_absolute_count_reduction_ceiling": total / counts[-1],
        "is_wall_time_or_storage_forecast": False,
    }


def info_fields(path: Path) -> dict:
    text = path.read_text()
    result = {}
    for key in ("aexp", "unit_l", "unit_d", "H0", "omega_m", "omega_b", "boxlen", "levelmin", "levelmax"):
        found = re.findall(rf"^\s*{key}\s*=\s*(\S+)", text, re.M)
        if len(found) != 1:
            raise ValueError(f"{path}: expected one {key}")
        value = float(found[0].replace("D", "E"))
        if not math.isfinite(value):
            raise ValueError(f"{path}: nonfinite {key}")
        result[key] = value
    if (not 0 < result["aexp"] <= 1 or min(result[k] for k in
            ("unit_l", "unit_d", "H0", "omega_m")) <= 0 or result["omega_b"] != 0):
        raise ValueError("expected a positive-unit total-matter DMO snapshot")
    result["path"] = str(path)
    result["sha256"] = hashlib.sha256(text.encode()).hexdigest()
    result["mass_unit_msun"] = result["unit_d"] * result["unit_l"]**3 / MSUN_G
    result["box_comoving_hinv_mpc_code_convention"] = (
        result["unit_l"] / result["aexp"] / CODE_MPC_CM * result["H0"] / 100)
    return result


def saved_namelist_contract(path: Path, ic_root: Path, lower: int, upper: int,
                            force_cap: int) -> dict:
    raw = path.read_bytes()
    text = re.sub(r"!.*", "", raw.decode())
    def number(key):
        found = re.findall(rf"^\s*{key}\s*=\s*([^\s,]+)", text, re.M)
        if len(found) != 1:
            raise ValueError(f"{path}: ambiguous {key}")
        value = float(found[0].replace("d", "e").replace("D", "E"))
        if not math.isfinite(value):
            raise ValueError(f"{path}: nonfinite {key}")
        return value
    if number("levelmin") != lower or number("levelmax") != force_cap:
        raise ValueError("saved namelist IC or force level mismatch")
    entries = re.findall(r"^\s*initfile\((\d+)\)\s*=\s*'([^']+)'", text, re.M)
    expected = [(str(index + 1), str(ic_root / f"level_{level:03d}"))
                for index, level in enumerate(range(lower, upper + 1))]
    if entries != expected:
        raise ValueError("saved namelists do not use the same ordered IC paths")
    cutoff = number("mass_cut_refine")
    high_mass = 2.0**(-3 * upper)
    if not high_mass < cutoff < 8 * high_mass:
        raise ValueError("saved refinement mass cutoff does not isolate finest particles")
    return {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
            "ic_paths": [value for _, value in entries], "force_cap": force_cap,
            "mass_cut_refine": cutoff}


def audit() -> dict:
    config_path = ZOOM / "music_zoom_l11.used.conf"
    raw_config = config_path.read_bytes()
    config = configparser.ConfigParser()
    config.read_file(io.StringIO(raw_config.decode()))
    lower, upper = config.getint("setup", "levelmin"), config.getint("setup", "levelmax")
    if (lower, upper) != (8, 11) or config.getboolean("setup", "baryons"):
        raise ValueError("audit expects the existing 8-to-11 DMO hierarchy")
    # Hash the two small provenance inputs, NOT any multi-GB IC payload.
    point_file = (ZOOM / config.get("setup", "region_point_file")).resolve()
    manifest_path = ZOOM / "zoom_ic_provenance.sha256"
    manifest = dict(line.split(maxsplit=1)[::-1] for line in manifest_path.read_text().splitlines())
    for relative, raw in (("music_zoom_l11.used.conf", raw_config),
                          (config.get("setup", "region_point_file"), point_file.read_bytes())):
        if manifest.get(relative) != hashlib.sha256(raw).hexdigest():
            raise ValueError(f"saved IC provenance mismatch: {relative}")
    headers = []
    for level in range(lower, upper + 1):
        baseline = None
        for field in FIELDS:
            path = ZOOM / "ic_zoom_l11" / f"level_{level:03d}" / field
            with path.open("rb") as stream:
                header = decode_header(stream, str(path))
            comparable = {key: value for key, value in header.items()
                          if key not in ("header_sha256",)}
            if baseline is not None and comparable != baseline:
                raise ValueError(f"field header mismatch at level {level}: {field}")
            baseline = comparable
            header.update({"path": str(path), "level": level, "field": field})
            headers.append(header)
    density_headers = [h for h in headers if h["field"] == "ic_deltab"]
    base = density_headers[0]
    if base["dimensions"] != [2**lower] * 3:
        raise ValueError("base header does not match the IC grid")
    for header in density_headers:
        expected_dx = base["cell_comoving_mpc"] / 2**(header["level"] - lower)
        if not math.isclose(header["cell_comoving_mpc"], expected_dx, rel_tol=2e-6):
            raise ValueError("IC cell widths do not halve across levels")
        for key in ("a_start", "omega_m", "omega_l", "H0", "omega_b"):
            if header[key] != base[key]:
                raise ValueError("IC cosmology differs across levels")
    infos = [info_fields(VALIDATION / directory / label / f"output_{output:05d}" /
                         f"info_{output:05d}.txt")
             for directory, label, output in (("z0_production_l17", "cdm", 5),
                                               ("z0_production_l17", "sidm3", 5),
                                               ("z0_production_l18", "sidm1", 12))]
    ref = infos[0]
    for info in infos:
        if (not math.isclose(info["mass_unit_msun"], ref["mass_unit_msun"], rel_tol=1e-12) or
                info["levelmin"] != lower or
                not math.isclose(info["H0"], base["H0"], rel_tol=2e-13) or
                not math.isclose(info["omega_m"], base["omega_m"], rel_tol=2e-13)):
            raise ValueError("snapshot units/cosmology do not match the common IC")
    box_h = base["cell_comoving_mpc"] * 2**lower * base["H0"] / 100
    if not math.isclose(box_h, config.getfloat("setup", "boxlength"), rel_tol=2e-6):
        raise ValueError("IC header and MUSIC boxlength differ")
    for info in infos:
        if not math.isclose(info["box_comoving_hinv_mpc_code_convention"], box_h, rel_tol=2e-6):
            raise ValueError("snapshot physical-length units differ from IC box")
    qc_raw = QC.read_bytes()
    qc = json.loads(qc_raw)
    counts = particle_scenarios(qc, upper)
    namelists = [saved_namelist_contract(VALIDATION / directory / label / "snapshot" /
                                        f"run_{job}.nml", ZOOM / "ic_zoom_l11", lower, upper, cap)
                 for directory, label, job, cap in
                 (("z0_production_l17", "cdm", 324571, 17),
                  ("z0_production_l17", "sidm3", 324572, 17),
                  ("z0_production_l18", "sidm1", 323699, 18))]
    if len({n["mass_cut_refine"] for n in namelists}) != 1:
        raise ValueError("legacy refinement mass selections differ")
    mass_unit = ref["mass_unit_msun"]
    h = base["H0"] / 100
    candidates = [{"ic_finest_level": level, "equivalent_uniform_nside": 2**level,
                   "particle_mass_code": 2.0**(-3 * level),
                   "particle_mass_msun": mass_unit * 2.0**(-3 * level),
                   "particle_mass_hinv_msun": mass_unit * h * 2.0**(-3 * level),
                   "suggested_cut_between_finest_and_next_mass": 2 * 2.0**(-3 * level),
                   "existing_cut_accepts_candidate_particles":
                       2.0**(-3 * level) < namelists[0]["mass_cut_refine"],
                   "is_realized_ic": level == upper}
                  for level in (9, 10, 11)]
    force = [{"amr_level_cap": level,
              "nominal_cell_comoving_hinv_kpc": box_h * 1000 / 2**level,
              "nominal_cell_comoving_physical_unit_kpc":
                  ref["unit_l"] / ref["aexp"] / MPC_CM * 1000 / 2**level,
              "is_measured_force_softening": False}
             for level in (16, 17, 18)]
    units_source = PROJECT / "patch/cuRamses/units.f90"
    compact_units = "".join(units_source.read_text().split()).lower()
    if "scale_l=aexp*boxlen_ini*3.08d24/(h0/100)" not in compact_units:
        raise ValueError("source length-unit convention changed; audit needs revision")
    return {
        "scope": "IC headers, saved configuration, cached particle QC and snapshot units",
        "audit_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "science_seed_certified": False, "production_launch_authorized": False,
        "config_path": str(config_path), "config_sha256": hashlib.sha256(raw_config).hexdigest(),
        "region_point_sha256": hashlib.sha256(point_file.read_bytes()).hexdigest(),
        "provenance_matches": True,
        "generator_binary_and_source_revision_certified": False,
        "source_sha256": {
            str(path.relative_to(PROJECT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (units_source, PROJECT / "patch/lagRamses/init_part.f90",
                         PROJECT / "patch/lagRamses/init_time.kjhan.f90",
                         PROJECT / "patch/cuRamses/rho_fine.kjhan.f90")},
        "configuration": {"box_hinv_mpc": config.getfloat("setup", "boxlength"),
                          "omega_b_cosmology": config.getfloat("cosmology", "Omega_b"),
                          "transfer_plugin": config.get("cosmology", "transfer"),
                          "seeds": {str(l): config.getint("random", f"seed[{l}]")
                                    for l in range(lower, upper + 1)}},
        "constants": {"Msun_g": MSUN_G, "Mpc_cm": MPC_CM,
                      "RAMSES_code_Mpc_cm": CODE_MPC_CM},
        "headers": headers, "total_header_bytes_read": sum(x["bytes_read"] for x in headers),
        "omega_b_header_contract_satisfied": all(x["omega_b"] is not None for x in headers),
        "total_matter_density_and_velocity_transfer_certified": False,
        "snapshot_units": infos,
        "saved_namelist_contracts": namelists,
        "cached_qc": {"path": str(QC), "sha256": hashlib.sha256(qc_raw).hexdigest(),
                      "snapshot": qc["snapshot"], "recomputed_payload": False,
                      "mass_tiers": qc["mass_tiers"]},
        "particle_count_scenarios": counts, "mass_resolution_candidates": candidates,
        "force_cell_scales": force,
        "limitations": [
            "Header agreement does not certify payloads, masks, phases or transfer-function selection.",
            "Cached particle QC covers the recorded legacy CDM snapshot, not a new SIDM payload audit.",
            "Only existing level-11 ICs are realized; lower levels need consistent regenerated masks.",
            "Mass selection for refinement must scale with particle mass; old cutoff excludes coarser candidates.",
            "Nominal AMR cell spacing is neither occupied-cell distribution nor calibrated force softening.",
            "A particle-count reduction is not a measured wall-time or dump-size reduction.",
            "DMO total-matter transfer and cosmological baryon metadata require a reader/writer contract audit.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an audit")
    payload = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(args.output)
    print(json.dumps(payload["particle_count_scenarios"], indent=2))
    print(json.dumps(payload["mass_resolution_candidates"], indent=2))


if __name__ == "__main__":
    main()
