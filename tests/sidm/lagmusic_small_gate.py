#!/usr/bin/env python3
"""Small IC contract gate only; cannot certify a target halo or release production."""
from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import math
import re
import struct
from pathlib import Path

import numpy as np

from audit_zoom_resolution import decode_header
from analyze_mg_convergence_log import analyze_lines

FIELDS = ("ic_deltab", "ic_velcx", "ic_velcy", "ic_velcz", "ic_poscx",
          "ic_poscy", "ic_poscz", "ic_refmap")
CASES = {"z6r1": (6, 1), "z6r4": (6, 4), "z5r1": (5, 1)}
MAX_CELLS = 2**21


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def prepare(root):
    for label, (upper, ranks) in CASES.items():
        case = root / label
        case.mkdir()
        ic = case / "ic"
        if len(str(ic / f"level_{upper:03d}")) > 79:
            raise ValueError("IC path exceeds RAMSES initfile capacity")
        config = f"""[setup]
boxlength = 128
zstart = 100
levelmin = 4
levelmin_TF = 4
levelmax = {upper}
padding = 2
overlap = 2
ref_center = 0.5, 0.5, 0.5
ref_extent = 0.25, 0.25, 0.25
align_top = no
baryons = no
use_2LPT = yes
dmo_velocity_source = density_2lpt
vfact_scale = 1
use_LLA = no
periodic_TF = yes
kspace_TF = no
[cosmology]
Omega_m = 0.3111
Omega_L = 0.6889
Omega_b = 0.04
H0 = 67.66
sigma_8 = 0.8102
nspec = 0.9665
w0 = -1
wa = 0
transfer = eisenstein
[random]
disk_cached = yes
seed[4] = 20260401
seed[5] = 20260402
seed[6] = 20260403
[output]
format = grafic2
filename = {ic}
ramses_nml = yes
[poisson]
fft_fine = yes
accuracy = 1e-6
pre_smooth = 3
post_smooth = 3
smoother = gs
laplace_order = 6
grad_order = 6
"""
        (case / "music.conf").write_text(config)
        paths = "\n".join(f"initfile({i-3})='{ic}/level_{i:03d}'"
                          for i in range(4, upper+1))
        nml = f"""&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
sidm=.true.
use_adm=.false.
use_fdm=.false.
use_pbh=.false.
sink=.false.
rt=.false.
nrestart=0
nstepmax=0
nremap=0
nsubcycle=3*1
ncontrol=1
ordering='hilbert'
use_fftw=.true.
jobcontrolfile='jobcontrol.txt'
/
&COSMO_PARAMS
omega_b=0.04
/
&OUTPUT_PARAMS
noutput=1
aout=1.1
tout=1.0d100
foutput=1000000000
fbackup=1000000000
/
&INIT_PARAMS
filetype='grafic'
{paths}
/
&AMR_PARAMS
levelmin=4
levelmax={upper}
nexpand=0,1,1
ngridtot=300000
nparttot=400000
/
&REFINE_PARAMS
m_refine=3*8.
ivar_refine=0
mass_cut_refine={2*2.0**(-3*upper):.16e}
interpol_var=1
interpol_type=0
q_refine_holdback=.false.
/
&POISSON_PARAMS
epsilon=1.0d-4
maxiter_fine=80
abort_on_mg_nonconvergence=.true.
restart_phi_warm_start=.false.
/
&SIDM_PARAMS
sidm_cross_section=3.0
sidm_npart_min=2
sidm_type='constant'
sidm_courant=0.1
sidm_angular='isotropic'
/
"""
        (case / "run.nml").write_text(nml)
        (case / "jobcontrol.txt").touch()
    policy = {"purpose": "small_zoom_reader_initial_step_only",
              "highest_ic_levels": [5, 6], "base_level": 4,
              "transfer": "Eisenstein total-matter analytic fit; no CAMB switch",
              "noutput": 1, "aout": [1.1], "requested_redshifts": [-1/11],
              "tout": [1e100], "foutput": 10**9, "fbackup": 10**9,
              "expected_full_dumps": 0, "expected_dump_storage_bytes": 0,
              "estimated_ic_storage_upper_bytes": 128*1024**2,
              "case_ranks": {k: v[1] for k, v in CASES.items()},
              "science_certification": False, "production_release": False}
    write_json(root / "output_policy.json", policy)
    files = [root / "output_policy.json"]
    files += [root / label / name for label in CASES for name in ("music.conf", "run.nml")]
    write_json(root / "inputs.json", {str(p.relative_to(root)): digest(p) for p in files})


def read_field(path):
    with path.open("rb") as stream:
        endian = "<" if struct.unpack("<i", stream.read(4))[0] == 48 else ">"
        stream.seek(0)
        header = decode_header(stream, str(path))
        if header["header_record_bytes"] != 48:
            raise ValueError("new IC requires extended header")
        nx, ny, nz = header["dimensions"]
        if nx*ny*nz > MAX_CELLS:
            raise ValueError("payload exceeds small-gate size cap")
        data = np.empty((nz, ny, nx), dtype=np.float32)
        size = 4*nx*ny
        for k in range(nz):
            marker, payload, trailer = stream.read(4), stream.read(size), stream.read(4)
            if (len(marker) != 4 or len(payload) != size or len(trailer) != 4 or
                    struct.unpack(endian+"i", marker)[0] != size or
                    struct.unpack(endian+"i", trailer)[0] != size):
                raise ValueError(f"bad payload record: {path}, plane {k}")
            data[k] = np.frombuffer(payload, dtype=endian+"f4").reshape(ny, nx)
        if stream.read(1) or not np.isfinite(data).all():
            raise ValueError("trailing bytes or nonfinite IC payload")
    return header, data


def compare_fields(a, b):
    if a.shape != b.shape:
        raise ValueError("rank-replay geometry mismatch")
    delta = a.astype(np.float64)-b
    norm = float(np.linalg.norm(a.astype(np.float64)))
    error = float(np.linalg.norm(delta))/max(norm, 1e-30)
    if error > 2e-5 or float(np.max(np.abs(delta))) > 2e-5*max(float(np.max(np.abs(a))), 1e-20):
        raise ValueError(f"MPI rank replay failed: relative L2 {error}")
    return {"relative_l2": error, "max_abs": float(np.max(np.abs(delta))),
            "bitwise_equal": bool(np.array_equal(a, b))}


def audit_ic(root, diagnostic=False):
    results, loaded = {}, {}
    for label, (upper, _) in CASES.items():
        case = root / label
        config = configparser.ConfigParser()
        config.read(case / "music.conf")
        generated_nml = (case / "ic/ramses.nml").read_text()
        if not re.search(r"omega_b\s*=\s*0\.0*4(?:\D|$)", generated_nml):
            raise ValueError("writer namelist lost physical omega_b")
        rows = []
        for level in range(4, upper+1):
            common = None
            for field in FIELDS:
                path = case / "ic" / f"level_{level:03d}" / field
                header, data = read_field(path)
                if abs(header["omega_b"]-.04) > 1e-7:
                    raise ValueError("physical baryon header mismatch")
                for key, expected in (("omega_m", .3111), ("omega_l", .6889),
                                      ("H0", 67.66), ("a_start", 1/101)):
                    if not math.isclose(header[key], expected, rel_tol=2e-6):
                        raise ValueError(f"IC cosmology mismatch: {key}")
                if not math.isclose(header["cell_comoving_mpc"], 128/.6766/2**level,
                                    rel_tol=2e-6):
                    raise ValueError("IC cell width mismatch")
                if level == 4 and header["dimensions"] != [16]*3:
                    raise ValueError("base grid is not full-box level 4")
                geometry = (header["dimensions"], header["offset_comoving_mpc"])
                if common is not None and geometry != common:
                    raise ValueError("inconsistent field geometry")
                common = geometry
                if field == "ic_refmap":
                    if not np.isin(data, [0, 1]).all() or not np.any(data == 1):
                        raise ValueError("empty or nonbinary refinement mask")
                loaded[label, level, field] = data
                rows.append({"level": level, "field": field, "header": header,
                             "sha256": digest(path), "bytes": path.stat().st_size,
                             "minimum": float(data.min()), "maximum": float(data.max())})
        results[label] = rows
    replay, failed = {}, []
    for level in range(4, 7):
        for field in FIELDS:
            key = f"L{level}/{field}"
            a, b = loaded["z6r1", level, field], loaded["z6r4", level, field]
            try:
                replay[key] = {"pass": True, **compare_fields(a, b)}
            except ValueError as exc:
                failed.append(key)
                replay[key] = {"pass": False, "failure": str(exc),
                               "relative_l2": float(np.linalg.norm(a.astype(float)-b)) /
                                   max(float(np.linalg.norm(a.astype(float))), 1e-30)}
    # Same numeric seeds alone do not certify phase equivalence across generators
    # or across differently truncated hierarchies. Preserve measured common modes.
    a = loaded["z6r1", 4, "ic_deltab"].astype(float).ravel()
    b = loaded["z5r1", 4, "ic_deltab"].astype(float).ravel()
    cross_level = {"base_density_correlation": float(np.corrcoef(a, b)[0, 1]),
                   "certifies_existing_L11_phase": False,
                   "certifies_target_contamination": False}
    report = {"ic_contract_pass": not failed, "header_payload_pass": True,
              "rank_replay_pass": not failed, "failed_rank_fields": failed,
              "rank_replay": replay, "production_release": False,
              "cross_level_diagnostic": cross_level, "fields": results,
              "scope": "small new ICs only; not production mask or science convergence"}
    write_json(root / "ic_gate.json", report)
    if failed and not diagnostic:
        raise ValueError(f"MPI rank replay failed in {len(failed)} fields; see ic_gate.json")


def reader_summary(text, upper, require_gravity=True):
    if ("IC/namelist omega_b consistent:" not in text or
            "legacy grafic" in text or "Run completed" not in text):
        raise ValueError("missing header consistency/completion or legacy fallback")
    hits = re.findall(r"LM_GATE_MASS: mass=\s*(\S+) count=(\d+) bad=(\d+)", text)
    if len(hits) != 1:
        raise ValueError("missing unique actual particle audit")
    mass, count, bad = float(hits[0][0]), int(hits[0][1]), int(hits[0][2])
    if not math.isfinite(mass) or abs(mass-1) > 1e-12 or count <= 0 or bad:
        raise ValueError("non-unit DMO mass or invalid particle state")
    tiers = re.findall(r"LM_GATE_TIER: level=(\d+) mass=\s*(\S+) count=(\d+)", text)
    if [int(l) for l, _, _ in tiers] != list(range(4, upper+1)):
        raise ValueError("missing particle tier audit")
    if any(float(m) != 2.0**(-3*int(l)) or int(n) <= 0 for l, m, n in tiers):
        raise ValueError("wrong or empty particle mass tier")
    if sum(int(n) for _, _, n in tiers) != count:
        raise ValueError("particle count not accounted for by expected tiers")
    if abs(sum(float(m)*int(n) for _, m, n in tiers)-mass) > 1e-12:
        raise ValueError("actual mass disagrees with measured tier accounting")
    gravity = analyze_lines(text.splitlines(), 1e-4, 80)
    if require_gravity and not all(gravity[k] for k in ("convergence_gate_pass", "fftw_solver_verified", "run_completed")):
        raise ValueError("strict FFT/MG gate failed")
    return {"mass_code": mass, "particle_count": count, "invalid_values": bad,
            "mass_tiers": [{"level": int(l), "mass": float(m), "count": int(n)}
                           for l, m, n in tiers], "gravity": gravity}


def audit_readers(root, diagnostic=False):
    results, dumps = {}, []
    for label, upper in (("z6r1", 6), ("z5r1", 5)):
        case = root / label
        found = list(case.glob("output_*"))
        dumps.extend(str(p) for p in found)
        if found and not diagnostic:
            raise ValueError("reader preflight wrote an unexpected full dump")
        results[label] = reader_summary((case / "reader.log").read_text(), upper,
                                        require_gravity=not diagnostic)
    for relative, expected in json.loads((root / "inputs.json").read_text()).items():
        if digest(root / relative) != expected:
            raise ValueError("effective input changed after audit")
    physics = all(case["gravity"][key] for case in results.values()
                  for key in ("convergence_gate_pass", "fftw_solver_verified", "run_completed"))
    write_json(root / "reader_gate.json", {"reader_contract_pass": not dumps and physics,
               "mass_header_checks_pass": True, "reader_physics_pass": physics,
               "no_dump_gate_pass": not dumps,
               "cases": results, "new_full_dumps": len(dumps),
               "unexpected_dump_paths": dumps, "production_release": False})


def diagnostic(root):
    # Record independent successful checks without waiving any failed gate.
    audit_ic(root, diagnostic=True)
    audit_readers(root, diagnostic=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "ic", "reader", "diagnostic"))
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    {"prepare": prepare, "ic": audit_ic, "reader": audit_readers,
     "diagnostic": diagnostic}[args.stage](args.root)
