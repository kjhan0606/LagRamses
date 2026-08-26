#!/usr/bin/env python3
"""Measure DMO density--velocity spectra on a common periodic mesh.

Particle mass and momentum are CIC-deposited.  Their ratio estimates the
cell velocity; empty cells are filled by successive periodic nearest-shell
averages.  This practical estimator is accepted for science only where its
LCDM-relative response is stable to mesh changes.  It is not a DTFE estimator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import struct
import subprocess
import tempfile
import time

import numpy as np
from numba import njit

import measure_dmo_pk


FortranReader = measure_dmo_pk.FortranReader
read_info = measure_dmo_pk.read_info
MEASURE_DMO_PK_PATH = Path(measure_dmo_pk.__file__).resolve()
EXPECTED_MEASURE_DMO_PK_PATH = Path(__file__).with_name("measure_dmo_pk.py").resolve()
if MEASURE_DMO_PK_PATH != EXPECTED_MEASURE_DMO_PK_PATH:
    raise RuntimeError(
        f"measure_dmo_pk import is not source-bound: {MEASURE_DMO_PK_PATH}"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def runtime_pk_preflight(
    output_dir: Path, info: dict, snapshot: dict, expected_boxlen_mpc_h: float,
) -> dict[str, int | float | str | bool]:
    """Validate the exact on-the-fly CIC spectrum used for physical box units."""
    number = output_dir.name.rsplit("_", 1)[-1]
    path = output_dir / f"pk_cic_{number}.dat"
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty source-bound runtime P(k): {path}")
    header = "\n".join(path.read_text(errors="strict").splitlines()[:12])

    def field(pattern: str, label: str) -> str:
        match = re.search(pattern, header, re.MULTILINE)
        if not match:
            raise ValueError(f"missing {label} in {path}")
        return match.group(1).strip()

    aexp = float(field(r"^# Power spectrum at a_exp\s*=\s*([0-9.Ee+-]+)$", "aexp"))
    boxlen = float(field(r"^# boxlen \(Mpc/h\)\s*=\s*([0-9.Ee+-]+)$", "boxlen"))
    nmesh = int(field(r"^# N_mesh\s*=\s*(\d+)$", "N_mesh"))
    npart = int(field(r"^# N_part\s*=\s*(\d+)$", "N_part"))
    assignment = field(r"^# assignment\s*=\s*(\S+)$", "assignment")
    interlaced_text = field(r"^# interlaced\s*=\s*(\S+)$", "interlaced").lower()
    if interlaced_text not in {"true", "false"}:
        raise ValueError(f"invalid interlaced flag in {path}: {interlaced_text}")
    interlaced = interlaced_text == "true"
    values = (aexp, boxlen, expected_boxlen_mpc_h)
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError(f"non-positive or non-finite runtime P(k) metadata in {path}")
    if not math.isclose(aexp, float(info["aexp"]), rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"runtime P(k) aexp mismatch in {path}: {aexp}")
    if not math.isclose(
        boxlen, expected_boxlen_mpc_h, rel_tol=0.0, abs_tol=1e-10,
    ):
        raise ValueError(
            f"runtime P(k) box length {boxlen} != expected {expected_boxlen_mpc_h}"
        )
    if npart != int(snapshot["particle_count"]):
        raise ValueError(f"runtime P(k) particle count mismatch in {path}: {npart}")
    base_mesh = 1 << int(info["levelmin"])
    if nmesh != base_mesh:
        raise ValueError(f"runtime P(k) mesh {nmesh} != base mesh {base_mesh}")
    if assignment != "CIC-deconvolved" or not interlaced:
        raise ValueError(
            f"unsupported runtime P(k) estimator in {path}: "
            f"assignment={assignment}, interlaced={interlaced}"
        )
    return {
        "path": str(path.resolve()), "sha256": sha256(path),
        "aexp": aexp, "boxlen_mpc_h": boxlen, "nmesh": nmesh,
        "particle_count": npart, "assignment": assignment,
        "interlaced": interlaced,
    }


def exact_record(reader: FortranReader, dtype: str, count: int, name: str) -> np.ndarray:
    payload = reader.record()
    expected = np.dtype(dtype).itemsize * count
    if len(payload) != expected:
        raise ValueError(
            f"{name} record has {len(payload)} bytes, expected {expected}: {reader.path}"
        )
    return np.frombuffer(payload, dtype=dtype).copy()


def particle_header(path: Path) -> dict[str, int | float]:
    with FortranReader(path) as reader:
        header = {
            "ncpu": int(exact_record(reader, "=i4", 1, "ncpu")[0]),
            "ndim": int(exact_record(reader, "=i4", 1, "ndim")[0]),
            "npart": int(exact_record(reader, "=i4", 1, "npart")[0]),
        }
        exact_record(reader, "=i4", 4, "localseed")
        header["nstar_total"] = int(
            exact_record(reader, "=i8", 1, "nstar_total")[0]
        )
        header["mstar_total"] = float(
            exact_record(reader, "=f8", 1, "mstar_total")[0]
        )
        header["mstar_lost"] = float(
            exact_record(reader, "=f8", 1, "mstar_lost")[0]
        )
        header["nsink"] = int(exact_record(reader, "=i4", 1, "nsink")[0])
    return header


def validate_particle_header(
    header: dict[str, int | float], path: Path, expected_ncpu: int, expected_ndim: int
) -> None:
    if header["ncpu"] != expected_ncpu or header["ndim"] != expected_ndim:
        raise ValueError(f"RAMSES header mismatch in {path}: {header}")
    if header["npart"] < 0:
        raise ValueError(f"negative particle count in {path}: {header['npart']}")
    if (
        header["nstar_total"] != 0 or header["nsink"] != 0
        or header["mstar_total"] != 0.0 or header["mstar_lost"] != 0.0
    ):
        raise ValueError(f"non-DMO particle layout is unsupported in {path}: {header}")


def skip_exact_record(
    reader: FortranReader, expected: int, name: str, require_zero: bool = False
) -> None:
    marker = np.fromfile(reader.stream, dtype=np.dtype("=i4"), count=1)
    if marker.size != 1 or int(marker[0]) != expected:
        observed = None if marker.size != 1 else int(marker[0])
        raise ValueError(
            f"{name} record marker is {observed}, expected {expected}: {reader.path}"
        )
    if require_zero:
        remaining = expected
        while remaining:
            block_size = min(remaining, 1 << 20)
            block = reader.stream.read(block_size)
            if len(block) != block_size:
                raise EOFError(f"truncated {name} record: {reader.path}")
            if np.frombuffer(block, dtype=np.uint8).any():
                raise ValueError(f"non-DM particle type in {reader.path}")
            remaining -= block_size
    else:
        before = reader.stream.tell()
        reader.stream.seek(expected, os.SEEK_CUR)
        if reader.stream.tell() - before != expected:
            raise EOFError(f"truncated {name} record: {reader.path}")
    trailer = np.fromfile(reader.stream, dtype=np.dtype("=i4"), count=1)
    if trailer.size != 1 or int(trailer[0]) != expected:
        raise ValueError(f"mismatched {name} record trailer: {reader.path}")


def validate_particle_tail(reader: FortranReader, npart: int) -> None:
    skip_exact_record(reader, 8 * npart, "identity")
    skip_exact_record(reader, 4 * npart, "level")
    skip_exact_record(reader, npart, "particle type", require_zero=True)
    skip_exact_record(reader, 8 * npart, "particle potential")
    if reader.stream.read(1):
        raise ValueError(f"unexpected records after particle potential: {reader.path}")


def validate_particle_structure(path: Path, header: dict[str, int | float]) -> None:
    npart = int(header["npart"])
    expected_size = 77 * npart + 208
    if path.stat().st_size != expected_size:
        raise ValueError(
            f"unsupported particle layout size {path.stat().st_size}, "
            f"expected {expected_size}: {path}"
        )
    with FortranReader(path) as reader:
        reader.skip(8)
        for axis in range(3):
            skip_exact_record(reader, 8 * npart, f"position[{axis}]")
        for axis in range(3):
            skip_exact_record(reader, 8 * npart, f"velocity[{axis}]")
        skip_exact_record(reader, 8 * npart, "mass")
        validate_particle_tail(reader, npart)


def particle_phase_space(
    path: Path, ndim: int, expected_ncpu: int,
    expected_header: dict[str, int | float], boxlen_code: float,
):
    with FortranReader(path) as reader:
        ncpu = int(exact_record(reader, "=i4", 1, "ncpu")[0])
        file_ndim = int(exact_record(reader, "=i4", 1, "ndim")[0])
        npart = int(exact_record(reader, "=i4", 1, "npart")[0])
        exact_record(reader, "=i4", 4, "localseed")
        nstar_total = int(exact_record(reader, "=i8", 1, "nstar_total")[0])
        mstar_total = float(exact_record(reader, "=f8", 1, "mstar_total")[0])
        mstar_lost = float(exact_record(reader, "=f8", 1, "mstar_lost")[0])
        nsink = int(exact_record(reader, "=i4", 1, "nsink")[0])
        observed = {
            "ncpu": ncpu, "ndim": file_ndim, "npart": npart,
            "nstar_total": nstar_total, "mstar_total": mstar_total,
            "mstar_lost": mstar_lost, "nsink": nsink,
        }
        validate_particle_header(observed, path, expected_ncpu, ndim)
        if observed != expected_header:
            raise ValueError(f"particle header changed after preflight: {path}")
        positions = np.empty((npart, ndim), dtype=np.float64)
        velocities = np.empty_like(positions)
        for axis in range(ndim):
            positions[:, axis] = exact_record(
                reader, "=f8", npart, f"position[{axis}]"
            )
        for axis in range(ndim):
            velocities[:, axis] = exact_record(
                reader, "=f8", npart, f"velocity[{axis}]"
            )
        masses = exact_record(reader, "=f8", npart, "mass")
        validate_particle_tail(reader, npart)
    tolerance = 8.0 * np.finfo(np.float64).eps * boxlen_code
    if (
        not np.all(np.isfinite(positions)) or np.any(positions < -tolerance)
        or np.any(positions >= boxlen_code + tolerance)
    ):
        raise ValueError(f"particle position outside periodic box: {path}")
    if not np.all(np.isfinite(velocities)):
        raise ValueError(f"non-finite particle velocity: {path}")
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError(f"non-finite or non-positive particle mass: {path}")
    return positions, velocities, masses


def header_particle_total(path: Path) -> int:
    text = path.read_text(encoding="ascii")
    match = re.search(r"Total number of particles\s*\n\s*(\d+)", text)
    if not match:
        raise ValueError(f"missing total particle count in {path}")
    return int(match.group(1))


def amr_expansion_state(path: Path, info: dict) -> dict[str, float | str]:
    digest = hashlib.sha256()

    def record(reader: FortranReader, dtype: str, count: int, name: str):
        payload = reader.record()
        digest.update(struct.pack("=Q", len(payload)))
        digest.update(payload)
        expected = np.dtype(dtype).itemsize * count
        if len(payload) != expected:
            raise ValueError(
                f"AMR {name} record has {len(payload)} bytes, expected {expected}: {path}"
            )
        return np.frombuffer(payload, dtype=dtype).copy()

    with FortranReader(path) as reader:
        ncpu = int(record(reader, "=i4", 1, "ncpu")[0])
        ndim = int(record(reader, "=i4", 1, "ndim")[0])
        record(reader, "=i4", 3, "coarse dimensions")
        nlevelmax = int(record(reader, "=i4", 1, "nlevelmax")[0])
        record(reader, "=i4", 1, "ngridmax")
        record(reader, "=i4", 1, "nboundary")
        record(reader, "=i4", 1, "ngrid_current")
        record(reader, "=f8", 1, "boxlen")
        output_state = record(reader, "=i4", 3, "output state")
        noutput = int(output_state[0])
        if noutput <= 0 or nlevelmax <= 0:
            raise ValueError(f"invalid AMR dimensions in {path}")
        record(reader, "=f8", noutput, "tout")
        record(reader, "=f8", noutput, "aout")
        record(reader, "=f8", 1, "time")
        record(reader, "=f8", nlevelmax, "dtold")
        record(reader, "=f8", nlevelmax, "dtnew")
        record(reader, "=i4", 2, "step counters")
        record(reader, "=f8", 3, "mass state")
        cosmology = record(reader, "=f8", 7, "cosmology")
        expansion = record(reader, "=f8", 5, "expansion state")
    if ncpu != int(info["ncpu"]) or ndim != int(info.get("ndim", 3)):
        raise ValueError(f"AMR header disagrees with info file: {path}")
    if not np.allclose(
        cosmology[:5],
        [info["omega_m"], info["omega_l"], info["omega_k"], info["omega_b"], info["H0"]],
        rtol=2e-7, atol=1e-12,
    ):
        raise ValueError(f"AMR cosmology disagrees with info file: {path}")
    aexp, hexp = float(expansion[0]), float(expansion[1])
    if not math.isclose(aexp, float(info["aexp"]), rel_tol=2e-7, abs_tol=1e-12):
        raise ValueError(f"AMR aexp disagrees with info file: {path}")
    if not math.isfinite(hexp) or hexp <= 0.0:
        raise ValueError(f"invalid serialized RAMSES hexp={hexp}: {path}")
    hubble = float(info["H0"]) * hexp / aexp**2
    return {
        "aexp": aexp, "hexp_supercomoving": hexp,
        "hubble_km_s_mpc": hubble,
        "source": str(path.resolve()),
        "header_records_sha256": digest.hexdigest(),
    }


def completion_provenance(
    output_dir: Path, allow_legacy_completion: bool, run_log: Path | None
) -> dict[str, str]:
    marker = output_dir / "COMPLETE"
    if marker.is_file():
        return {
            "mode": "snapshot_COMPLETE", "path": str(marker.resolve()),
            "sha256": sha256(marker),
        }
    if not allow_legacy_completion:
        raise FileNotFoundError(
            f"missing {marker}; legacy snapshots require "
            "--allow-legacy-completion --run-log PATH"
        )
    if run_log is None or not run_log.is_file() or run_log.stat().st_size == 0:
        raise FileNotFoundError("legacy completion requires a nonempty --run-log")
    model = output_dir.parent.name
    text = run_log.read_text(errors="replace")
    if "Run completed" not in text or not re.search(
        rf"^model={re.escape(model)}\s+end=", text, flags=re.MULTILINE
    ):
        raise ValueError(f"run log does not prove completion of model={model}: {run_log}")
    binary_match = re.search(
        r"^([0-9a-f]{64})\s+(/\S+)\s*$", text, flags=re.MULTILINE
    )
    if not binary_match:
        raise ValueError(f"run log has no binary SHA256 line: {run_log}")
    logged_sha, binary_token = binary_match.groups()
    binary = Path(binary_token)
    if not binary.is_file() or binary.is_symlink():
        raise FileNotFoundError(f"missing or symlinked production binary: {binary}")
    actual_sha = sha256(binary)
    if actual_sha != logged_sha:
        raise ValueError(
            f"production binary SHA mismatch: log={logged_sha}, actual={actual_sha}"
        )
    code_match = re.search(
        rf"^model={re.escape(model)}\s+.*\bcode=([0-9a-f]+)\b", text,
        flags=re.MULTILINE,
    )
    if not code_match:
        raise ValueError(f"run log has no model/code identity: {run_log}")
    return {
        "mode": "legacy_completed_run_log", "path": str(run_log.resolve()),
        "sha256": sha256(run_log),
        "binary_path": str(binary.resolve()), "binary_sha256": actual_sha,
        "logged_code_prefix": code_match.group(1),
    }


def cpl_expansion_crosscheck(info: dict, run_namelist: Path, hubble: float) -> dict:
    text = run_namelist.read_text(errors="strict")
    if not re.search(r"^\s*&CPL_PARAMS\b", text, flags=re.MULTILINE | re.IGNORECASE):
        return {"family": "non_CPL", "performed": False}

    def parameter(name: str) -> float:
        match = re.search(
            rf"^\s*{name}\s*=\s*([-+0-9.eEdD]+)", text,
            flags=re.MULTILINE | re.IGNORECASE,
        )
        if not match:
            raise ValueError(f"missing {name} in {run_namelist}")
        return float(match.group(1).replace("D", "E").replace("d", "e"))

    w0, wa = parameter("w0"), parameter("wa")
    aexp = float(info["aexp"])
    f_de = aexp ** (-3.0 * (1.0 + w0 + wa)) * math.exp(-3.0 * wa * (1.0 - aexp))
    analytic = float(info["H0"]) * math.sqrt(
        float(info["omega_m"]) / aexp**3
        + float(info["omega_l"]) * f_de
        + float(info["omega_k"]) / aexp**2
    )
    relative = abs(hubble / analytic - 1.0)
    if relative >= 1.0e-3:
        raise ValueError(
            f"serialized and analytic CPL H(a) differ by {relative:.3e}: {run_namelist}"
        )
    return {
        "family": "CPL", "performed": True, "w0": w0, "wa": wa,
        "analytic_hubble_km_s_mpc": analytic, "relative_difference": relative,
    }


def namelist_parameter_map(path: Path) -> dict[str, str]:
    parameters = {}
    for raw_line in path.read_text(errors="strict").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        match = re.match(r"([A-Za-z][A-Za-z0-9_]*(?:\(\d+\))?)\s*=\s*(.+)$", line)
        if match:
            parameters[match.group(1).lower()] = match.group(2).rstrip(",").strip()
    if not parameters:
        raise ValueError(f"no namelist assignments found in {path}")
    return parameters


def source_layout_provenance(output_dir: Path, commit: str) -> dict[str, str]:
    repository = Path(__file__).resolve().parents[3]
    makefile = output_dir / "makefile.txt"
    text = makefile.read_text(errors="strict")
    vpath_match = re.search(r"^VPATH\s*=\s*(.+)$", text, flags=re.MULTILINE)
    if not vpath_match or "../patch/cuRamses" not in vpath_match.group(1):
        raise ValueError(f"snapshot makefile does not expose cuRamses VPATH: {makefile}")
    candidates = (
        "patch/lagRamses/output_part.f90", "patch/cuda/output_part.f90",
        "patch/oct_tree/output_part.f90", "patch/cuRamses/output_part.f90",
        "pm/output_part.f90",
    )
    selected = None
    for candidate in candidates:
        result = subprocess.run(
            ["git", "-C", str(repository), "cat-file", "-e", f"{commit}:{candidate}"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
        if result.returncode == 0:
            selected = candidate
            break
    if selected != "patch/cuRamses/output_part.f90":
        raise ValueError(f"unexpected output_part source at {commit}: {selected}")
    source = subprocess.run(
        ["git", "-C", str(repository), "show", f"{commit}:{selected}"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    ).stdout
    source_text = source.decode("utf-8")
    if "! Write particle type" not in source_text or "! Write family" in source_text:
        raise ValueError(f"source does not implement the observed 19-record layout: {selected}")
    blob = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", f"{commit}:{selected}"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    ).stdout.strip()
    return {
        "repository": str(repository), "commit": commit, "selected_source": selected,
        "source_blob": blob, "source_sha256": hashlib.sha256(source).hexdigest(),
        "vpath": vpath_match.group(1).strip(),
    }


def campaign_binary_audit(campaign_manifest: Path, completion: dict) -> dict:
    campaign = json.loads(campaign_manifest.read_text())
    declared_token = campaign.get("ramses")
    if not isinstance(declared_token, str) or not declared_token:
        return {"status": "not_declared"}
    declared = Path(declared_token)
    result = {"declared_path": str(declared), "status": "missing"}
    if declared.is_file():
        result["declared_sha256"] = sha256(declared)
        result["status"] = "available"
    executed = completion.get("binary_path")
    if executed is not None:
        result["executed_path"] = executed
        result["executed_sha256"] = completion["binary_sha256"]
        result["status"] = (
            "match" if Path(executed).resolve() == declared.resolve()
            else "stale_campaign_manifest_overridden_by_hashed_run_log"
        )
    return result


def snapshot_preflight(
    output_dir: Path, info: dict, hash_particles: bool,
    allow_legacy_completion: bool = False, run_log: Path | None = None,
) -> dict:
    number = output_dir.name.rsplit("_", 1)[-1]
    if not re.fullmatch(r"\d{5}", number):
        raise ValueError(f"unsupported output directory name: {output_dir.name}")
    if list(output_dir.glob("*.h5")) or list(output_dir.glob("*.hdf5")):
        raise ValueError(f"HDF5/grouped particle layouts are unsupported: {output_dir}")
    required = [
        output_dir / f"info_{number}.txt", output_dir / f"header_{number}.txt",
        output_dir / "compilation.txt", output_dir / "namelist.txt",
        output_dir / "makefile.txt",
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing or empty snapshot provenance file: {path}")
    completion = completion_provenance(output_dir, allow_legacy_completion, run_log)
    model_dir, campaign_dir = output_dir.parent, output_dir.parent.parent
    campaign_manifest, run_namelist = campaign_dir / "campaign.json", model_dir / "run.nml"
    for path in (campaign_manifest, run_namelist):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing model provenance file: {path}")
    ncpu, ndim = int(info["ncpu"]), int(info.get("ndim", 3))
    headers: list[dict[str, int | float]] = []
    files: list[dict[str, int | str]] = []
    amr_states: list[dict[str, float | str]] = []
    for cpu in range(1, ncpu + 1):
        part = output_dir / f"part_{number}.out{cpu:05d}"
        amr = output_dir / f"amr_{number}.out{cpu:05d}"
        for product in (part, amr):
            if not product.is_file() or product.stat().st_size == 0:
                raise FileNotFoundError(f"missing or empty snapshot file: {product}")
        header = particle_header(part)
        validate_particle_header(header, part, ncpu, ndim)
        validate_particle_structure(part, header)
        headers.append(header)
        item: dict[str, int | str] = {
            "path": str(part.resolve()), "bytes": part.stat().st_size,
        }
        if hash_particles:
            item["sha256"] = sha256(part)
        files.append(item)
        amr_states.append(amr_expansion_state(amr, info))
    particle_count = sum(int(header["npart"]) for header in headers)
    expected = header_particle_total(output_dir / f"header_{number}.txt")
    if particle_count != expected or particle_count <= 0:
        raise ValueError(
            f"particle total mismatch: binary={particle_count}, header={expected}"
        )
    expansion = amr_states[0]
    for state in amr_states[1:]:
        if not (
            math.isclose(float(state["aexp"]), float(expansion["aexp"]), rel_tol=0.0, abs_tol=1e-13)
            and math.isclose(
                float(state["hexp_supercomoving"]),
                float(expansion["hexp_supercomoving"]), rel_tol=0.0, abs_tol=1e-13,
            )
        ):
            raise ValueError("aexp/hexp differs among AMR CPU headers")
    expansion["all_cpu_header_records"] = [
        {"source": state["source"], "sha256": state["header_records_sha256"]}
        for state in amr_states
    ]
    crosscheck = cpl_expansion_crosscheck(
        info, run_namelist, float(expansion["hubble_km_s_mpc"])
    )
    provenance = {
        path.name: {"path": str(path.resolve()), "sha256": sha256(path)}
        for path in required
    }
    compilation_text = (output_dir / "compilation.txt").read_text(errors="strict")
    commit_match = re.search(r"^\s*last commit\s*=\s*([0-9a-f]+)\s*$", compilation_text, re.MULTILINE)
    if not commit_match:
        raise ValueError(f"missing source commit in {output_dir / 'compilation.txt'}")
    commit = commit_match.group(1)
    code_prefix = completion.get("logged_code_prefix")
    if code_prefix is not None and not commit.startswith(code_prefix):
        raise ValueError(
            f"run-log code prefix {code_prefix} disagrees with source commit {commit}"
        )
    provenance["campaign.json"] = {
        "path": str(campaign_manifest.resolve()), "sha256": sha256(campaign_manifest),
    }
    provenance["run.nml"] = {
        "path": str(run_namelist.resolve()), "sha256": sha256(run_namelist),
    }
    return {
        "completion_contract": completion,
        "structural_contract": (
            "all indexed part/amr files nonempty; exact 19-record DMO layout and EOF; "
            "DM type zero; binary particle sum equals header total"
        ),
        "model": model_dir.name, "source_commit": commit,
        "source_layout": source_layout_provenance(output_dir, commit),
        "campaign_binary_audit": campaign_binary_audit(
            campaign_manifest, completion
        ),
        "model_parameters": namelist_parameter_map(run_namelist),
        "particle_count": particle_count, "particle_headers": headers,
        "max_npart_local": max(int(header["npart"]) for header in headers),
        "particle_files": files, "provenance_files": provenance,
        "python_dependencies": {
            "measure_dmo_pk": {
                "path": str(MEASURE_DMO_PK_PATH),
                "sha256": sha256(MEASURE_DMO_PK_PATH),
            },
        },
        "expansion": expansion, "expansion_crosscheck": crosscheck,
    }


def estimated_peak_bytes(nmesh: int, max_npart_local: int = 0) -> int:
    if nmesh <= 0 or nmesh % 2:
        raise ValueError("nmesh must be a positive even integer")
    if max_npart_local < 0:
        raise ValueError("max_npart_local must be non-negative")
    particle_buffer = 56 * max_npart_local
    return int(224 * nmesh**3 + particle_buffer + (1 << 30))


def memory_preflight(
    nmesh: int, memory_limit_gb: float | None, max_npart_local: int = 0
) -> dict[str, float]:
    estimate = estimated_peak_bytes(nmesh, max_npart_local)
    minimum = max(estimate, 40 * 1024**3 if nmesh >= 512 else estimate)
    if nmesh >= 512 and memory_limit_gb is None:
        raise ValueError("nmesh >= 512 requires explicit --memory-limit-gb (48 recommended)")
    if memory_limit_gb is None:
        values = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            values[key] = int(value.split()[0]) * 1024
        limit = int(0.75 * values["MemAvailable"])
        source = "0.75*MemAvailable"
    else:
        if not math.isfinite(memory_limit_gb) or memory_limit_gb <= 0:
            raise ValueError("--memory-limit-gb must be positive")
        limit = int(memory_limit_gb * 1024**3)
        source = "--memory-limit-gb"
    if minimum > limit:
        raise MemoryError(
            f"required memory {minimum / 1024**3:.1f} GiB "
            f"(estimated peak {estimate / 1024**3:.1f} GiB) exceeds "
            f"allowed {limit / 1024**3:.1f} GiB ({source})"
        )
    return {
        "estimated_peak_gib": estimate / 1024**3,
        "minimum_required_gib": minimum / 1024**3,
        "max_npart_local": float(max_npart_local),
        "allowed_gib": limit / 1024**3, "limit_source": source,
    }


@njit(cache=True)
def deposit_mass_momentum_pair(
    positions, velocities, masses, mass_a, mass_b, mom_a, mom_b, nmesh, boxlen_code
):
    scale = nmesh / boxlen_code
    for particle in range(masses.size):
        for shifted in range(2):
            offset = -0.5 * shifted
            ux = positions[particle, 0] * scale + offset
            uy = positions[particle, 1] * scale + offset
            uz = positions[particle, 2] * scale + offset
            i0, j0, k0 = int(np.floor(ux)), int(np.floor(uy)), int(np.floor(uz))
            fx, fy, fz = ux - np.floor(ux), uy - np.floor(uy), uz - np.floor(uz)
            mass_mesh = mass_a if shifted == 0 else mass_b
            momentum = mom_a if shifted == 0 else mom_b
            for di in range(2):
                wx, ii = ((1.0 - fx) if di == 0 else fx), (i0 + di) % nmesh
                for dj in range(2):
                    wy, jj = ((1.0 - fy) if dj == 0 else fy), (j0 + dj) % nmesh
                    for dk in range(2):
                        wz, kk = ((1.0 - fz) if dk == 0 else fz), (k0 + dk) % nmesh
                        weight = masses[particle] * wx * wy * wz
                        mass_mesh[ii, jj, kk] += weight
                        for axis in range(3):
                            momentum[axis, ii, jj, kk] += weight * velocities[particle, axis]


def deposit_output(output_dir: Path, info: dict, nmesh: int, preflight: dict):
    ncpu, ndim = int(info["ncpu"]), int(info.get("ndim", 3))
    if ndim != 3:
        raise ValueError("only three-dimensional DMO snapshots are supported")
    number = output_dir.name.rsplit("_", 1)[-1]
    shape = (nmesh, nmesh, nmesh)
    mass_a, mass_b = np.zeros(shape), np.zeros(shape)
    mom_a, mom_b = np.zeros((3, *shape)), np.zeros((3, *shape))
    particle_count, mass_sum = 0, 0.0
    velocity_factor = float(info["unit_l"]) / float(info["unit_t"]) / 1.0e5
    for cpu in range(1, ncpu + 1):
        path = output_dir / f"part_{number}.out{cpu:05d}"
        positions, velocities, masses = particle_phase_space(
            path, ndim, ncpu, preflight["particle_headers"][cpu - 1],
            float(info["boxlen"]),
        )
        velocities *= velocity_factor
        deposit_mass_momentum_pair(
            positions, velocities, masses, mass_a, mass_b, mom_a, mom_b,
            nmesh, float(info["boxlen"]),
        )
        particle_count += masses.size
        mass_sum += float(np.sum(masses, dtype=np.float64))
    for name, mesh in (("unshifted", mass_a), ("shifted", mass_b)):
        if not np.isclose(np.sum(mesh), mass_sum, rtol=2e-12):
            raise RuntimeError(f"{name} CIC mass conservation failed")
    if particle_count != preflight["particle_count"]:
        raise RuntimeError("particle count changed between preflight and deposition")
    return mass_a, mass_b, mom_a, mom_b, particle_count, mass_sum


def benchmark_deposition(
    output_dir: Path, info: dict, nmesh: int, preflight: dict, cpu_count: int
) -> dict[str, float | int]:
    ncpu, ndim = int(info["ncpu"]), int(info.get("ndim", 3))
    if cpu_count <= 0 or cpu_count > ncpu:
        raise ValueError(f"benchmark cpu count must be in [1,{ncpu}]")
    number = output_dir.name.rsplit("_", 1)[-1]
    shape = (nmesh, nmesh, nmesh)
    mass_a, mass_b = np.zeros(shape), np.zeros(shape)
    mom_a, mom_b = np.zeros((3, *shape)), np.zeros((3, *shape))
    particles = 0
    started = time.perf_counter()
    velocity_factor = float(info["unit_l"]) / float(info["unit_t"]) / 1.0e5
    for cpu in range(1, cpu_count + 1):
        path = output_dir / f"part_{number}.out{cpu:05d}"
        positions, velocities, masses = particle_phase_space(
            path, ndim, ncpu, preflight["particle_headers"][cpu - 1],
            float(info["boxlen"]),
        )
        velocities *= velocity_factor
        deposit_mass_momentum_pair(
            positions, velocities, masses, mass_a, mass_b, mom_a, mom_b,
            nmesh, float(info["boxlen"]),
        )
        particles += masses.size
    elapsed = time.perf_counter() - started
    seconds_per_particle = elapsed / particles
    return {
        "benchmark_cpus": cpu_count, "particles": particles,
        "elapsed_seconds": elapsed, "seconds_per_particle": seconds_per_particle,
        "estimated_full_deposition_seconds": (
            seconds_per_particle * preflight["particle_count"]
        ),
        "peak_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2,
    }


def velocity_mesh(momentum: np.ndarray, mass: np.ndarray):
    occupied = mass > 0
    initial_empty = int(np.size(occupied) - np.count_nonzero(occupied))
    velocity = np.zeros_like(momentum)
    velocity[:, occupied] = momentum[:, occupied] / mass[occupied]
    passes = 0
    while not np.all(occupied):
        count = np.zeros(mass.shape, dtype=np.uint8)
        for axis in range(3):
            count += np.roll(occupied, 1, axis=axis)
            count += np.roll(occupied, -1, axis=axis)
        fill = (~occupied) & (count > 0)
        if not np.any(fill):
            raise RuntimeError("periodic empty-cell filling stalled")
        for component in range(3):
            neighbor_sum = np.zeros(mass.shape)
            for axis in range(3):
                neighbor_sum += np.roll(velocity[component] * occupied, 1, axis=axis)
                neighbor_sum += np.roll(velocity[component] * occupied, -1, axis=axis)
            velocity[component, fill] = neighbor_sum[fill] / count[fill]
        occupied[fill] = True
        passes += 1
    return velocity, initial_empty, passes


def interlaced_fourier(field_a, field_b):
    nmesh = field_a.shape[0]
    fourier_a, fourier_b = np.fft.rfftn(field_a), np.fft.rfftn(field_b)
    phase_x = np.exp(-1j * np.pi * np.fft.fftfreq(nmesh))
    phase_y = phase_x
    phase_z = np.exp(-1j * np.pi * np.fft.rfftfreq(nmesh))
    fourier_b *= phase_x[:, None, None]
    fourier_b *= phase_y[None, :, None]
    fourier_b *= phase_z[None, None, :]
    fourier_a += fourier_b
    fourier_a *= 0.5
    return fourier_a


def spectra_from_meshes(
    mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, nmesh
):
    mean_mass = mass_sum / nmesh**3
    delta = interlaced_fourier(mass_a / mean_mass - 1, mass_b / mean_mass - 1)
    occupied_a, occupied_b = mass_a > 0, mass_b > 0
    velocity_a, empty_a, passes_a = velocity_mesh(mom_a, mass_a)
    velocity_b, empty_b, passes_b = velocity_mesh(mom_b, mass_b)
    del mom_a, mom_b, mass_a, mass_b
    kfund = 2 * np.pi / boxlen
    k_axes = (
        np.fft.fftfreq(nmesh) * nmesh * kfund,
        np.fft.fftfreq(nmesh) * nmesh * kfund,
        np.fft.rfftfreq(nmesh) * nmesh * kfund,
    )
    return delta, velocity_a, velocity_b, occupied_a, occupied_b, k_axes, {
        "empty_cells_unshifted": empty_a, "empty_cells_shifted": empty_b,
        "empty_fraction_unshifted": empty_a / nmesh**3,
        "empty_fraction_shifted": empty_b / nmesh**3,
        "fill_passes_unshifted": passes_a, "fill_passes_shifted": passes_b,
    }


def finish_spectra(
    delta, velocity_a, velocity_b, occupied_a, occupied_b, k_axes,
    boxlen, aexp, hubble, h, kmax,
):
    nmesh = velocity_a.shape[1]
    theta = np.zeros_like(delta)
    theta_zero_fill = np.zeros_like(delta)
    for axis in range(3):
        velocity_fourier = interlaced_fourier(velocity_a[axis], velocity_b[axis])
        shape = [1, 1, 1]
        shape[axis] = k_axes[axis].size
        factor = -1j * h / (aexp * hubble) * k_axes[axis].reshape(shape)
        theta += factor * velocity_fourier
        velocity_a[axis, ~occupied_a] = 0.0
        velocity_b[axis, ~occupied_b] = 0.0
        theta_zero_fill += factor * interlaced_fourier(
            velocity_a[axis], velocity_b[axis]
        )
    del velocity_a, velocity_b, occupied_a, occupied_b
    mode_limit = min(nmesh // 2 - 1, int(np.floor(kmax / (2 * np.pi / boxlen))))
    normalization = boxlen**3 / nmesh**6
    sums: dict[int, list[float]] = {}
    for nx in range(-mode_limit, mode_limit + 1):
        ix = nx % nmesh
        for ny in range(-mode_limit, mode_limit + 1):
            iy = ny % nmesh
            for nz in range(mode_limit + 1):
                n2 = nx * nx + ny * ny + nz * nz
                if n2 == 0 or (2 * np.pi / boxlen) * math.sqrt(n2) > kmax:
                    continue
                multiplicity = 1 if nz == 0 else 2
                window = (
                    np.sinc(nx / nmesh) ** 2 * np.sinc(ny / nmesh) ** 2
                    * np.sinc(nz / nmesh) ** 2
                )
                d = delta[ix, iy, nz] / window
                t = theta[ix, iy, nz] / window
                t_no_deconvolution = theta[ix, iy, nz]
                t_zero_fill = theta_zero_fill[ix, iy, nz] / window
                values = sums.setdefault(
                    n2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                )
                values[0] += multiplicity * abs(d) ** 2 * normalization
                values[1] += multiplicity * float(np.real(d * np.conj(t))) * normalization
                values[2] += multiplicity * abs(t) ** 2 * normalization
                values[3] += multiplicity * float(
                    np.real(d * np.conj(t_no_deconvolution))
                ) * normalization
                values[4] += multiplicity * abs(t_no_deconvolution) ** 2 * normalization
                values[5] += multiplicity * float(
                    np.real(d * np.conj(t_zero_fill))
                ) * normalization
                values[6] += multiplicity * abs(t_zero_fill) ** 2 * normalization
                values[7] += multiplicity
    n2_values = np.asarray(sorted(sums))
    nmodes = np.asarray([sums[n][7] for n in n2_values], dtype=int)
    pdd = np.asarray([sums[n][0] for n in n2_values]) / nmodes
    pdt = np.asarray([sums[n][1] for n in n2_values]) / nmodes
    ptt = np.asarray([sums[n][2] for n in n2_values]) / nmodes
    pdt_no_deconvolution = np.asarray([sums[n][3] for n in n2_values]) / nmodes
    ptt_no_deconvolution = np.asarray([sums[n][4] for n in n2_values]) / nmodes
    pdt_zero_fill = np.asarray([sums[n][5] for n in n2_values]) / nmodes
    ptt_zero_fill = np.asarray([sums[n][6] for n in n2_values]) / nmodes
    k = (2 * np.pi / boxlen) * np.sqrt(n2_values)
    correlation = np.full_like(pdd, np.nan)
    f_cross, f_auto = np.full_like(pdd, np.nan), np.full_like(pdd, np.nan)
    valid_cross = (pdd > 0) & (ptt > 0)
    correlation[valid_cross] = pdt[valid_cross] / np.sqrt(
        pdd[valid_cross] * ptt[valid_cross]
    )
    valid_density = pdd > 0
    f_cross[valid_density] = pdt[valid_density] / pdd[valid_density]
    f_auto[valid_density] = np.sqrt(ptt[valid_density] / pdd[valid_density])
    f_cross_no_deconvolution = np.full_like(pdd, np.nan)
    f_auto_no_deconvolution = np.full_like(pdd, np.nan)
    f_cross_zero_fill = np.full_like(pdd, np.nan)
    f_auto_zero_fill = np.full_like(pdd, np.nan)
    f_cross_no_deconvolution[valid_density] = (
        pdt_no_deconvolution[valid_density] / pdd[valid_density]
    )
    f_auto_no_deconvolution[valid_density] = np.sqrt(
        ptt_no_deconvolution[valid_density] / pdd[valid_density]
    )
    f_cross_zero_fill[valid_density] = pdt_zero_fill[valid_density] / pdd[valid_density]
    f_auto_zero_fill[valid_density] = np.sqrt(
        ptt_zero_fill[valid_density] / pdd[valid_density]
    )
    radius = 8.0
    kr = k * radius
    window8 = 3 * (np.sin(kr) - kr * np.cos(kr)) / kr**3
    sigma8_sq = float(np.sum(nmodes * pdd * window8**2) / boxlen**3)
    cross8 = float(np.sum(nmodes * pdt * window8**2) / boxlen**3)
    theta8_sq = float(np.sum(nmodes * ptt * window8**2) / boxlen**3)
    sigma8 = math.sqrt(max(sigma8_sq, 0.0))
    return {
        "k_h_mpc": k, "p_delta_delta": pdd, "p_delta_theta": pdt,
        "p_theta_theta": ptt, "nmodes": nmodes, "r_delta_theta": correlation,
        "f_cross": f_cross, "f_auto": f_auto,
        "p_delta_theta_velocity_not_deconvolved": pdt_no_deconvolution,
        "p_theta_theta_velocity_not_deconvolved": ptt_no_deconvolution,
        "p_delta_theta_zero_fill": pdt_zero_fill,
        "p_theta_theta_zero_fill": ptt_zero_fill,
        "f_cross_velocity_not_deconvolved": f_cross_no_deconvolution,
        "f_auto_velocity_not_deconvolved": f_auto_no_deconvolution,
        "f_cross_zero_fill": f_cross_zero_fill,
        "f_auto_zero_fill": f_auto_zero_fill,
        "sigma8_box_truncated": sigma8,
        "fsigma8_cross_box_truncated": cross8 / sigma8 if sigma8 > 0 else math.nan,
        "fsigma8_auto_box_truncated": math.sqrt(max(theta8_sq, 0.0)),
    }


def write_products(path: Path, result: dict, metadata: dict):
    if path.suffix != ".npz":
        raise ValueError(f"destination must have .npz suffix: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_target = path.with_suffix(".json")
    manifest_target = path.with_suffix(".manifest.json")
    complete_target = path.with_suffix(".COMPLETE")
    targets = (path, sidecar_target, manifest_target, complete_target)
    if any(target.exists() for target in targets):
        raise FileExistsError(f"refusing to overwrite product set for {path}")
    with tempfile.TemporaryDirectory(prefix=".density-velocity-", dir=path.parent) as name:
        temporary = Path(name)
        npz = temporary / path.name
        arrays = {key: value for key, value in result.items() if isinstance(value, np.ndarray)}
        np.savez(npz, **arrays)
        metadata["integrated_amplitudes"] = {
            key: value for key, value in result.items() if not isinstance(value, np.ndarray)
        }
        sidecar = temporary / sidecar_target.name
        sidecar.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        manifest = temporary / manifest_target.name
        manifest_payload = {
            "status": "complete_marker_required",
            "files": {
                path.name: sha256(npz), sidecar_target.name: sha256(sidecar),
            },
        }
        manifest.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n")
        complete = temporary / complete_target.name
        complete.write_text(
            json.dumps({"manifest_sha256": sha256(manifest)}, sort_keys=True) + "\n"
        )
        os.replace(npz, path)
        os.replace(sidecar, sidecar_target)
        os.replace(manifest, manifest_target)
        os.replace(complete, complete_target)


def write_fortran_record(handle, values: np.ndarray) -> None:
    payload = np.ascontiguousarray(values).tobytes()
    marker = struct.pack("=i", len(payload))
    handle.write(marker); handle.write(payload); handle.write(marker)


def write_synthetic_snapshot(
    root: Path, nmesh: int, mode: tuple[int, int, int],
    *, header_ncpu: int = 1, header_ndim: int = 3, nstar_total: int = 0,
    nsink: int = 0, particle_type: int = 0,
) -> tuple[Path, dict, float, float, float, float]:
    model_dir = root / "synthetic_model"
    output = model_dir / "output_00001"
    output.mkdir(parents=True)
    (root / "campaign.json").write_text('{"model":"synthetic_model"}\n')
    (model_dir / "run.nml").write_text("&RUN_PARAMS\nsynthetic=.true.\n/\n")
    (output / "COMPLETE").touch()
    repository = Path(__file__).resolve().parents[3]
    commit = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        text=True, stdout=subprocess.PIPE, check=True,
    ).stdout.strip()
    (output / "compilation.txt").write_text(
        f" compile date = synthetic\n last commit  = {commit}\n"
    )
    (output / "namelist.txt").write_text("&RUN_PARAMS\n/\n")
    (output / "makefile.txt").write_text(
        "VPATH = ../patch/lagRamses:../patch/cuda:../patch/oct_tree:"
        "../patch/cuRamses:../pm\n"
    )

    aexp, hubble, h, growth = 0.5, 120.0, 0.6766, 0.8
    h0, omega_m, omega_l = 100.0 * h, 0.3111, 0.6889
    level = int(round(math.log2(nmesh)))
    if 1 << level != nmesh:
        raise ValueError("synthetic nmesh must be a power of two")
    info_text = (
        f"ncpu = 1\nndim = 3\nlevelmin = {level}\nlevelmax = {level}\n"
        f"boxlen = 1.0\naexp = {aexp}\nH0 = {h0}\n"
        f"omega_m = {omega_m}\nomega_l = {omega_l}\nomega_k = 0.0\n"
        "omega_b = 0.0\nunit_l = 1.0e5\nunit_d = 1.0\nunit_t = 1.0\n"
    )
    (output / "info_00001.txt").write_text(info_text)

    grid = (np.indices((nmesh, nmesh, nmesh)).reshape(3, -1).T + 0.5) / nmesh
    mode_vector = np.asarray(mode, dtype=float)
    phase = 2.0 * np.pi * (grid @ mode_vector)
    amplitude = 1.0e-3
    masses = 1.0 + amplitude * np.cos(phase)
    boxlen = 160.0
    physical_k = 2.0 * np.pi * mode_vector / boxlen
    velocities = np.tile(np.array([13.0, -7.0, 5.0]), (grid.shape[0], 1))
    velocities += (
        -(aexp * hubble / h) * growth * amplitude
        * np.sin(phase)[:, None] * physical_k[None, :]
        / float(physical_k @ physical_k)
    )
    npart = grid.shape[0]
    (output / "header_00001.txt").write_text(
        f" Total number of particles\n {npart}\n"
        f" Total number of dark matter particles\n {npart}\n"
        " Total number of star particles\n 0\n Total number of sink particles\n 0\n"
    )
    (output / "pk_cic_00001.dat").write_text(
        f"# Power spectrum at a_exp = {aexp:.12E}\n"
        f"# boxlen (Mpc/h) = {boxlen:.12E}\n"
        f"# N_mesh = {nmesh}\n"
        f"# N_part = {npart}\n"
        "# shot_noise (Mpc/h)^3 = 1.0\n"
        "# assignment = CIC-deconvolved\n"
        "# interlaced = true\n"
        "# kmax (h/Mpc) = 1.0\n"
    )
    part = output / "part_00001.out00001"
    with part.open("wb") as handle:
        write_fortran_record(handle, np.asarray([header_ncpu], dtype="=i4"))
        write_fortran_record(handle, np.asarray([header_ndim], dtype="=i4"))
        write_fortran_record(handle, np.asarray([npart], dtype="=i4"))
        write_fortran_record(handle, np.asarray([1, 2, 3, 4], dtype="=i4"))
        write_fortran_record(handle, np.asarray([nstar_total], dtype="=i8"))
        write_fortran_record(handle, np.asarray([0.0], dtype="=f8"))
        write_fortran_record(handle, np.asarray([0.0], dtype="=f8"))
        write_fortran_record(handle, np.asarray([nsink], dtype="=i4"))
        for axis in range(3):
            write_fortran_record(handle, grid[:, axis].astype("=f8"))
        for axis in range(3):
            write_fortran_record(handle, velocities[:, axis].astype("=f8"))
        write_fortran_record(handle, masses.astype("=f8"))
        write_fortran_record(handle, np.arange(1, npart + 1, dtype="=i8"))
        write_fortran_record(handle, np.full(npart, level, dtype="=i4"))
        write_fortran_record(handle, np.full(npart, particle_type, dtype="i1"))
        write_fortran_record(handle, np.zeros(npart, dtype="=f8"))

    hexp = hubble / h0 * aexp**2
    amr = output / "amr_00001.out00001"
    with amr.open("wb") as handle:
        records = (
            np.asarray([1], dtype="=i4"), np.asarray([3], dtype="=i4"),
            np.asarray([1, 1, 1], dtype="=i4"), np.asarray([1], dtype="=i4"),
            np.asarray([16], dtype="=i4"), np.asarray([0], dtype="=i4"),
            np.asarray([1], dtype="=i4"), np.asarray([1.0], dtype="=f8"),
            np.asarray([1, 1, 1], dtype="=i4"), np.asarray([0.0], dtype="=f8"),
            np.asarray([aexp], dtype="=f8"), np.asarray([0.0], dtype="=f8"),
            np.asarray([0.0], dtype="=f8"), np.asarray([0.0], dtype="=f8"),
            np.asarray([1, 1], dtype="=i4"), np.asarray([1.0, 1.0, 1.0], dtype="=f8"),
            np.asarray([omega_m, omega_l, 0.0, 0.0, h0, 0.02, boxlen], dtype="=f8"),
            np.asarray([aexp, hexp, aexp, 0.0, 0.0], dtype="=f8"),
        )
        for values in records:
            write_fortran_record(handle, values)
    return output, read_info(output), boxlen, aexp, hubble, growth


def synthetic_test():
    nmesh, boxlen, aexp, hubble, h, growth = 32, 128.0, 0.5, 120.0, 0.6766, 0.8
    x = np.arange(nmesh)[:, None, None] * boxlen / nmesh
    x_shifted = (np.arange(nmesh)[:, None, None] + 0.5) * boxlen / nmesh
    k = 2 * np.pi / boxlen
    delta_real = 0.05 * np.cos(k * x) * np.ones((1, nmesh, nmesh))
    velocity_x = -(aexp * hubble / h) * growth * 0.05 / k * np.sin(k * x)
    velocity_x_shifted = -(aexp * hubble / h) * growth * 0.05 / k * np.sin(
        k * x_shifted
    )
    velocity_x = velocity_x * np.ones((1, nmesh, nmesh))
    velocity_x_shifted = velocity_x_shifted * np.ones((1, nmesh, nmesh))
    velocity = np.zeros((3, nmesh, nmesh, nmesh)); velocity[0] = velocity_x
    velocity_shifted = np.zeros_like(velocity); velocity_shifted[0] = velocity_x_shifted
    delta = np.fft.rfftn(delta_real)
    axes = (
        np.fft.fftfreq(nmesh) * nmesh * k,
        np.fft.fftfreq(nmesh) * nmesh * k,
        np.fft.rfftfreq(nmesh) * nmesh * k,
    )
    occupied = np.ones((nmesh, nmesh, nmesh), dtype=bool)
    measured = finish_spectra(
        delta, velocity, velocity_shifted, occupied, occupied.copy(), axes,
        boxlen, aexp, hubble, h, 0.1,
    )
    index = int(np.argmin(abs(measured["k_h_mpc"] - k)))
    if abs(measured["f_cross"][index] / growth - 1) > 1e-12:
        raise RuntimeError("analytic density--velocity closure failed")
    print("analytic mesh density--velocity closure passed")


def synthetic_particle_test():
    nmesh, h = 32, 0.6766
    for mode in ((1, 0, 0), (1, 2, 1), (3, 2, 0), (15, 0, 0)):
        with tempfile.TemporaryDirectory(prefix="density-velocity-fixture-") as name:
            output, info, boxlen, aexp, hubble, growth = write_synthetic_snapshot(
                Path(name), nmesh, mode
            )
            preflight = snapshot_preflight(output, info, True)
            runtime = runtime_pk_preflight(output, info, preflight, boxlen)
            if runtime["boxlen_mpc_h"] != boxlen:
                raise RuntimeError("synthetic runtime P(k) box length mismatch")
            try:
                runtime_pk_preflight(output, info, preflight, boxlen + 1.0)
            except ValueError:
                pass
            else:
                raise RuntimeError("runtime P(k) expected-box mismatch was accepted")
            deposited = deposit_output(output, info, nmesh, preflight)
            mass_a, mass_b, mom_a, mom_b, count, mass_sum = deposited
            if count != nmesh**3:
                raise RuntimeError("synthetic parser particle count mismatch")
            products = spectra_from_meshes(
                mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, nmesh
            )
            delta, velocity_a, velocity_b, occupied_a, occupied_b, axes, occupancy = products
            mode_vector = np.asarray(mode, dtype=float)
            grid = (
                np.indices((nmesh, nmesh, nmesh)).reshape(3, -1).T + 0.5
            ) / nmesh
            phase = 2.0 * np.pi * (grid @ mode_vector)
            masses = 1.0 + 1.0e-3 * np.cos(phase)
            expected_density = np.sum(masses * np.exp(-1j * phase))
            index3 = (mode[0] % nmesh, mode[1] % nmesh, mode[2])
            cic_window = np.prod([np.sinc(value / nmesh) ** 2 for value in mode])
            interlaced_assignment = 0.5 * (
                1.0 + np.prod([np.cos(np.pi * value / nmesh) for value in mode])
            )
            expected_deconvolved = (
                expected_density * interlaced_assignment / cic_window
            )
            observed_density = delta[index3] / cic_window
            density_amplitude_error = abs(
                observed_density / expected_deconvolved - 1.0
            )
            theta_fourier = np.zeros_like(delta)
            for axis in range(3):
                axis_shape = [1, 1, 1]
                axis_shape[axis] = axes[axis].size
                theta_fourier += (
                    -1j * h / (aexp * hubble)
                    * axes[axis].reshape(axis_shape)
                    * interlaced_fourier(velocity_a[axis], velocity_b[axis])
                )
            observed_theta = theta_fourier[index3] / cic_window
            theta_amplitude_error = abs(
                observed_theta / (growth * expected_deconvolved) - 1.0
            )
            if density_amplitude_error > 1.0e-4 or theta_amplitude_error > 1.0e-4:
                raise RuntimeError(
                    f"absolute Fourier amplitude failed for mode={mode}: "
                    f"density={density_amplitude_error}, theta={theta_amplitude_error}"
                )
            analysis_kmax = max(0.2, 1.01 * 2.0 * np.pi / boxlen * math.sqrt(
                sum(value * value for value in mode)
            ))
            measured = finish_spectra(
                delta, velocity_a, velocity_b, occupied_a, occupied_b, axes,
                boxlen, aexp, hubble, h, analysis_kmax,
            )
            target_k = 2.0 * np.pi / boxlen * math.sqrt(
                sum(value * value for value in mode)
            )
            index = int(np.argmin(abs(measured["k_h_mpc"] - target_k)))
            cross_error = abs(measured["f_cross"][index] / growth - 1.0)
            auto_power_error = abs(
                measured["p_theta_theta"][index]
                / measured["p_delta_delta"][index] / growth**2 - 1.0
            )
            if cross_error > 1.0e-4 or auto_power_error > 2.0e-4:
                raise RuntimeError(
                    f"particle CIC closure failed for mode={mode}: "
                    f"cross={cross_error}, auto_power={auto_power_error}"
                )
            if measured["r_delta_theta"][index] < 1.0 - 1.0e-8:
                raise RuntimeError(f"particle phase closure failed for mode={mode}")
            if occupancy["empty_cells_unshifted"] or occupancy["empty_cells_shifted"]:
                raise RuntimeError(f"unexpected empty cell in particle fixture mode={mode}")
    print("particle CIC/interlacing axial-oblique-cutoff closure passed")

    with tempfile.TemporaryDirectory(prefix="density-velocity-nyquist-") as name:
        output, info, boxlen, aexp, hubble, _ = write_synthetic_snapshot(
            Path(name), nmesh, (nmesh // 2, 0, 0)
        )
        preflight = snapshot_preflight(output, info, False)
        mass_a, mass_b, mom_a, mom_b, _, mass_sum = deposit_output(
            output, info, nmesh, preflight
        )
        products = spectra_from_meshes(
            mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, nmesh
        )
        measured = finish_spectra(
            *products[:6], boxlen, aexp, hubble, h, 1.0,
        )
        nyquist = 2.0 * np.pi / boxlen * (nmesh // 2)
        if np.any(np.isclose(measured["k_h_mpc"], nyquist, rtol=0.0, atol=1.0e-12)):
            raise RuntimeError("Nyquist shell was not excluded")
    print("near-Nyquist amplitude and Nyquist-exclusion tests passed")


def sparse_constant_velocity_test():
    nmesh, boxlen, aexp, hubble, h = 16, 160.0, 0.5, 120.0, 0.6766
    shape = (nmesh, nmesh, nmesh)
    mass_a, mass_b = np.zeros(shape), np.zeros(shape)
    mass_a[::2, ::2, ::2] = 1.0
    mass_b[1::2, 1::2, 1::2] = 1.0
    bulk = np.asarray([13.0, -7.0, 5.0])
    mom_a = bulk[:, None, None, None] * mass_a
    mom_b = bulk[:, None, None, None] * mass_b
    products = spectra_from_meshes(
        mass_a, mass_b, mom_a, mom_b, float(mass_a.sum()), boxlen, nmesh
    )
    delta, velocity_a, velocity_b, occupied_a, occupied_b, axes, _ = products
    measured = finish_spectra(
        delta, velocity_a, velocity_b, occupied_a, occupied_b, axes,
        boxlen, aexp, hubble, h, 0.2,
    )
    if float(np.max(np.abs(measured["p_theta_theta"]))) > 1.0e-24:
        raise RuntimeError("constant-velocity sparse fill generated spurious divergence")
    print("sparse constant-velocity empty-fill closure passed")


def parser_rejection_test():
    cases = (
        ("wrong_ncpu", {"header_ncpu": 2}, None),
        ("wrong_ndim", {"header_ndim": 2}, None),
        ("nonzero_star", {"nstar_total": 1}, None),
        ("nonzero_sink", {"nsink": 1}, None),
        ("nonzero_type", {"particle_type": 1}, None),
        ("corrupt_marker", {}, "marker"),
        ("truncated", {}, "truncate"),
    )
    for label, options, corruption in cases:
        with tempfile.TemporaryDirectory(prefix=f"density-velocity-{label}-") as name:
            output, info, *_ = write_synthetic_snapshot(
                Path(name), 2, (1, 0, 0), **options
            )
            part = output / "part_00001.out00001"
            if corruption == "marker":
                with part.open("r+b") as handle:
                    handle.seek(8); handle.write(struct.pack("=i", 99))
            elif corruption == "truncate":
                with part.open("r+b") as handle:
                    handle.truncate(part.stat().st_size - 1)
            try:
                snapshot_preflight(output, info, False)
            except (EOFError, FileNotFoundError, ValueError):
                continue
            raise RuntimeError(f"parser failed open for case={label}")
    print("RAMSES parser fail-closed rejection tests passed")


def artifact_publication_test():
    with tempfile.TemporaryDirectory(prefix="density-velocity-artifact-") as name:
        destination = Path(name) / "fixture.npz"
        write_products(destination, {"x": np.asarray([1.0]), "value": 2.0}, {})
        sidecar = destination.with_suffix(".json")
        manifest_path = destination.with_suffix(".manifest.json")
        complete_path = destination.with_suffix(".COMPLETE")
        manifest = json.loads(manifest_path.read_text())
        complete = json.loads(complete_path.read_text())
        if complete["manifest_sha256"] != sha256(manifest_path):
            raise RuntimeError("artifact COMPLETE marker mismatch")
        expected = {
            destination.name: sha256(destination), sidecar.name: sha256(sidecar),
        }
        if manifest["files"] != expected:
            raise RuntimeError("artifact manifest mismatch")
    print("hash-manifest/last-COMPLETE publication test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--nmesh", type=int, default=256)
    parser.add_argument("--kmax", type=float, default=0.2)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--memory-limit-gb", type=float)
    parser.add_argument("--expected-boxlen-mpc-h", type=float)
    parser.add_argument("--particle-hashes", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--allow-legacy-completion", action="store_true")
    parser.add_argument("--run-log", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--benchmark-cpus", type=int)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        synthetic_test(); synthetic_particle_test(); sparse_constant_velocity_test()
        parser_rejection_test(); artifact_publication_test(); return
    no_destination = args.preflight_only or args.benchmark_cpus is not None
    if args.output is None or (args.destination is None and not no_destination):
        parser.error(
            "output and --destination are required unless preflight/benchmark mode is used"
        )
    if args.expected_boxlen_mpc_h is None:
        parser.error("--expected-boxlen-mpc-h is required for snapshot measurements")
    output = args.output.resolve(); info = read_info(output)
    preflight = snapshot_preflight(
        output, info, args.particle_hashes, args.allow_legacy_completion,
        args.run_log.resolve() if args.run_log else None,
    )
    preflight["runtime_pk"] = runtime_pk_preflight(
        output, info, preflight, args.expected_boxlen_mpc_h
    )
    memory = memory_preflight(
        args.nmesh, args.memory_limit_gb, preflight["max_npart_local"]
    )
    if args.preflight_only:
        print(json.dumps({"memory": memory, "snapshot": preflight}, indent=2))
        return
    if args.benchmark_cpus is not None:
        benchmark = benchmark_deposition(
            output, info, args.nmesh, preflight, args.benchmark_cpus
        )
        print(json.dumps({"memory": memory, "benchmark": benchmark}, indent=2))
        return
    boxlen = float(preflight["runtime_pk"]["boxlen_mpc_h"])
    mass_a, mass_b, mom_a, mom_b, count, mass_sum = deposit_output(
        output, info, args.nmesh, preflight
    )
    products = spectra_from_meshes(
        mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, args.nmesh
    )
    delta, velocity_a, velocity_b, occupied_a, occupied_b, axes, occupancy = products
    expansion = preflight["expansion"]
    hubble, h = expansion["hubble_km_s_mpc"], float(info["H0"]) / 100.0
    result = finish_spectra(
        delta, velocity_a, velocity_b, occupied_a, occupied_b, axes,
        boxlen, float(info["aexp"]), hubble, h, args.kmax,
    )
    empty_direct_pass = (
        occupancy["empty_fraction_unshifted"] <= 1.0e-3
        and occupancy["empty_fraction_shifted"] <= 1.0e-3
        and occupancy["fill_passes_unshifted"] <= 2
        and occupancy["fill_passes_shifted"] <= 2
    )
    metadata = {
        "status": "measurement_complete_science_gates_pending",
        "source_output": str(output), "nmesh": args.nmesh, "kmax_h_mpc": args.kmax,
        "boxlen_mpc_h": boxlen, "particle_count": count, "aexp": float(info["aexp"]),
        "hubble_km_s_mpc": hubble, "expansion_source": expansion,
        "occupancy": occupancy, "memory_preflight": memory,
        "peak_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2,
        "snapshot_preflight": preflight,
        "velocity_definition": "CIC momentum divided by CIC mass; periodic nearest-shell fill",
        "velocity_systematics": (
            "primary velocity is CIC-deconvolved; NPZ also stores undeconvolved-velocity "
            "and zero-empty-cell variants"
        ),
        "science_gates": {
            "empty_cell_direct_occupancy_pass": empty_direct_pass,
            "empty_cell_response_sensitivity": "pending model/LCDM comparison",
            "velocity_window_response_sensitivity": "pending model/LCDM comparison",
            "mesh_256_512_response_convergence": "pending paired comparison",
            "science_eligible": False,
        },
        "theta_definition": "-div(v)/(aH), coordinates in Mpc/h",
        "limitations": [
            "not a volume-weighted DTFE estimator",
            "science use requires nmesh convergence of LCDM-relative responses",
            "science use requires <0.5% low-k sensitivity to empty-cell treatment",
            "box-truncated sigma8 amplitudes omit modes below the fundamental and above kmax",
        ],
        "script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).resolve()),
    }
    write_products(args.destination.resolve(), result, metadata)
    print(args.destination)


if __name__ == "__main__":
    main()
