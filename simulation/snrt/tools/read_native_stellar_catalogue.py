#!/usr/bin/env python3
"""Extract an audited stellar metadata catalogue from native cuRAMSES output.

This reader is deliberately downstream of the native particle audit.  It
understands the compact Fortran-record layout emitted by ``output_part.f90``
in the registered Phase 0 build, including the one-byte ``ptypep`` record and
the eight-byte ``nstar_tot`` header.  It extracts only particles tagged
``PTYPE_STAR=1`` and never assigns an RT photon luminosity.

The catalogue is a metadata hand-off for a later stellar-population/SED step,
not a photon ledger.  In particular, ``birth_epoch`` is RAMSES conformal
time, while ``birth_proper_time`` is the proper-time value used by the
feedback code; they are intentionally separate output columns.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import struct
from typing import BinaryIO

import numpy as np


MSUN_G = 1.98847e33
SECONDS_PER_MYR = 1.0e6 * 365.25 * 86400.0
_PART_RE = re.compile(r"part_(\d{5})\.out(\d{5})$")
_INFO_VALUE = re.compile(r"^\s*([A-Za-z_]+)\s*=\s*([^!#]+)")

_HEADER_RECORD_NAMES = (
    "ncpu",
    "ndim",
    "npart",
    "localseed",
    "nstar_tot",
    "mstar_tot",
    "mstar_lost",
    "nsink",
)
_HEADER_RECORD_LENGTHS = (4, 4, 4, 16, 8, 8, 8, 4)
_DATA_SPECS = (
    ("position_x", np.dtype("<f8")),
    ("position_y", np.dtype("<f8")),
    ("position_z", np.dtype("<f8")),
    ("velocity_x", np.dtype("<f8")),
    ("velocity_y", np.dtype("<f8")),
    ("velocity_z", np.dtype("<f8")),
    ("mass", np.dtype("<f8")),
    ("identity", np.dtype("<i8")),
    ("level", np.dtype("<i4")),
    ("ptypep", np.dtype("<i1")),
    ("potential", np.dtype("<f8")),
    ("birth_epoch", np.dtype("<f8")),
    ("metallicity", np.dtype("<f8")),
    ("birth_proper_time", np.dtype("<f8")),
    ("initial_mass", np.dtype("<f8")),
    ("yield_table_index", np.dtype("<f8")),
)


@dataclass(frozen=True)
class Record:
    offset: int
    length: int
    dtype: np.dtype


@dataclass(frozen=True)
class PartFile:
    path: Path
    ncpu: int
    ndim: int
    npart: int
    nstar_tot: int
    nsink: int
    records: dict[str, Record]


def _read_marker(stream: BinaryIO, label: str) -> int:
    payload = stream.read(4)
    if len(payload) != 4:
        raise ValueError(f"unexpected EOF before {label} record marker")
    marker = struct.unpack("<i", payload)[0]
    if marker < 0:
        raise ValueError(f"negative {label} record length {marker}")
    return marker


def _read_header_record(stream: BinaryIO, label: str) -> bytes:
    marker = _read_marker(stream, label)
    payload = stream.read(marker)
    if len(payload) != marker:
        raise ValueError(f"unexpected EOF in {label} record")
    trailer = stream.read(4)
    if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
        raise ValueError(f"{label} record trailer mismatch")
    return payload


def _int_payload(payload: bytes, label: str, allowed_lengths: tuple[int, ...] = (4,)) -> int:
    if len(payload) not in allowed_lengths:
        raise ValueError(f"{label} record has {len(payload)} bytes, expected {allowed_lengths}")
    return int.from_bytes(payload, byteorder="little", signed=True)


def _double_payload(payload: bytes, label: str) -> float:
    if len(payload) != 8:
        raise ValueError(f"{label} record has {len(payload)} bytes, expected 8")
    return struct.unpack("<d", payload)[0]


def _scan_part_file(path: Path) -> PartFile:
    """Index one native part file without reading its large data records."""

    with path.open("rb") as stream:
        headers = [_read_header_record(stream, name) for name in _HEADER_RECORD_NAMES]
        lengths = tuple(len(payload) for payload in headers)
        if lengths != _HEADER_RECORD_LENGTHS:
            raise ValueError(
                f"{path}: header record lengths {lengths} do not match the registered cuRAMSES layout"
            )
        ncpu = _int_payload(headers[0], "ncpu")
        ndim = _int_payload(headers[1], "ndim")
        npart = _int_payload(headers[2], "npart")
        nstar_tot = _int_payload(headers[4], "nstar_tot", allowed_lengths=(4, 8))
        nsink = _int_payload(headers[7], "nsink")
        if ncpu <= 0 or ndim != 3 or npart < 0 or nstar_tot < 0 or nsink < 0:
            raise ValueError(f"{path}: invalid native particle header")

        records: dict[str, Record] = {}
        for name, dtype in _DATA_SPECS:
            marker = _read_marker(stream, name)
            expected = npart * dtype.itemsize
            if marker != expected:
                raise ValueError(f"{path}: {name} record has {marker} bytes, expected {expected}")
            offset = stream.tell()
            stream.seek(marker, 1)
            trailer = stream.read(4)
            if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
                raise ValueError(f"{path}: {name} record trailer mismatch")
            records[name] = Record(offset=offset, length=marker, dtype=dtype)

        if stream.tell() != path.stat().st_size:
            raise ValueError(f"{path}: trailing bytes after native particle records")

    return PartFile(
        path=path,
        ncpu=ncpu,
        ndim=ndim,
        npart=npart,
        nstar_tot=nstar_tot,
        nsink=nsink,
        records=records,
    )


def _read_selected(part: PartFile, name: str, mask: np.ndarray | None = None) -> np.ndarray:
    record = part.records[name]
    with part.path.open("rb") as stream:
        stream.seek(record.offset)
        if mask is None:
            values = np.fromfile(stream, dtype=record.dtype, count=part.npart)
        else:
            if mask.shape != (part.npart,):
                raise ValueError(f"{part.path}: selection mask has the wrong shape for {name}")
            selected: list[np.ndarray] = []
            chunk_size = 4 * 1024 * 1024
            for start in range(0, part.npart, chunk_size):
                stop = min(start + chunk_size, part.npart)
                payload = stream.read((stop - start) * record.dtype.itemsize)
                if len(payload) != (stop - start) * record.dtype.itemsize:
                    raise ValueError(f"{part.path}: unexpected EOF in {name} record")
                chunk = np.frombuffer(payload, dtype=record.dtype)
                selected.append(chunk[mask[start:stop]])
            values = np.concatenate(selected) if selected else np.empty(0, dtype=record.dtype)
        expected = part.npart if mask is None else int(mask.sum())
        if values.size != expected:
            raise ValueError(f"{part.path}: {name} returned {values.size} values, expected {expected}")
        trailer = stream.read(4)
        if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != record.length:
            raise ValueError(f"{part.path}: {name} record trailer mismatch")
    return values


def _read_info(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = _INFO_VALUE.match(line)
        if match is None:
            continue
        try:
            values[match.group(1)] = float(match.group(2).strip().replace("D", "E"))
        except ValueError:
            continue
    required = {"aexp", "unit_l", "unit_d", "unit_t"}
    missing = required - set(values)
    if missing:
        raise ValueError(f"RAMSES info file missing values: {sorted(missing)}")
    if not 0.0 < values["aexp"] <= 1.0:
        raise ValueError("RAMSES aexp must lie in (0, 1]")
    for name in ("unit_l", "unit_d", "unit_t"):
        if values[name] <= 0.0:
            raise ValueError(f"RAMSES {name} must be positive")
    return values


def _dark_energy_factor(a: np.ndarray, w0: float, wa: float) -> np.ndarray:
    if wa == 0.0 and w0 == -1.0:
        return np.ones_like(a)
    return a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))


def _is_cosmological(info: dict[str, float]) -> bool:
    return not (
        info.get("time", 0.0) >= 0.0
        and info.get("H0", 1.0) == 1.0
        and info["aexp"] == 1.0
    )


def _proper_time_code_from_aexp(info: dict[str, float], w0: float, wa: float) -> float:
    """Reproduce RAMSES' negative look-back proper-time convention.

    ``dadt`` in the active RAMSES source is
    ``sqrt(omega_m/a + omega_k + omega_l*a**2*f_de(a))`` in H0 units.
    RAMSES stores the present epoch at proper time zero, so the current value
    is the negative integral from the snapshot scale factor to one.
    """

    if not _is_cosmological(info):
        return info.get("time", 0.0)
    required = {"omega_m", "omega_l", "omega_k"}
    missing = required - set(info)
    if missing:
        raise ValueError(f"cosmological proper-time reconstruction needs {sorted(missing)}")
    omega_m = info["omega_m"]
    omega_l = info["omega_l"]
    omega_k = info["omega_k"]
    if omega_m <= 0.0 or omega_l < 0.0:
        raise ValueError("unsupported cosmological density parameters")
    nodes, weights = np.polynomial.legendre.leggauss(256)
    lo = info["aexp"]
    hi = 1.0
    a = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
    f_de = _dark_energy_factor(a, w0, wa)
    dadt_squared = omega_m / a + omega_k + omega_l * a * a * f_de
    if np.any(dadt_squared <= 0.0):
        raise ValueError("cosmological proper-time integrand is not positive")
    integral = 0.5 * (hi - lo) * np.sum(weights / np.sqrt(dadt_squared))
    return -float(integral)


def _part_files(output_dir: Path) -> list[Path]:
    matches = sorted(output_dir.glob("part_*.out*"))
    result = [path for path in matches if _PART_RE.fullmatch(path.name)]
    if not result:
        raise FileNotFoundError(f"no native part files in {output_dir}")
    return result


def _stats(values: np.ndarray) -> dict[str, float | int]:
    if values.size == 0:
        return {"count": 0}
    numeric = np.asarray(values, dtype=np.float64)
    if not np.isfinite(numeric).all():
        raise ValueError("stellar catalogue contains a non-finite value")
    return {
        "count": int(values.size),
        "min": float(numeric.min()),
        "max": float(numeric.max()),
        "mean": float(numeric.mean()),
    }


def extract_catalogue(
    output_dir: str | Path,
    info_path: str | Path,
    *,
    w0: float = -1.0,
    wa: float = 0.0,
    include_velocity: bool = False,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Extract all ``ptypep == 1`` rows and return rows plus manifest data."""

    output = Path(output_dir).expanduser().resolve()
    info = Path(info_path).expanduser().resolve()
    if not output.is_dir():
        raise FileNotFoundError(output)
    if not info.is_file():
        raise FileNotFoundError(info)
    part_paths = _part_files(output)
    parts = [_scan_part_file(path) for path in part_paths]
    first = parts[0]
    if len(parts) != first.ncpu:
        raise ValueError(f"found {len(parts)} rank files, expected {first.ncpu}")
    if any(part.ncpu != first.ncpu or part.ndim != 3 for part in parts):
        raise ValueError("native particle files disagree on ncpu/ndim")
    if any(part.nstar_tot != first.nstar_tot for part in parts):
        raise ValueError("native particle files disagree on nstar_tot")
    if any(part.nsink != 0 for part in parts):
        raise ValueError("the registered stellar reader expects the stopped checkpoint to have no sinks")

    columns: dict[str, list[np.ndarray]] = {
        "rank": [],
        "local_index": [],
        "source_id": [],
        "level": [],
        "position_x_code": [],
        "position_y_code": [],
        "position_z_code": [],
        "velocity_x_code": [],
        "velocity_y_code": [],
        "velocity_z_code": [],
        "mass_code": [],
        "initial_mass_code": [],
        "birth_epoch_code": [],
        "birth_proper_time_code": [],
        "birth_metallicity_mass_fraction": [],
        "yield_table_index_code": [],
    }
    for rank, part in enumerate(parts, start=1):
        ptype = _read_selected(part, "ptypep")
        mask = ptype == 1
        star_indices = np.flatnonzero(mask).astype(np.int64)
        if star_indices.size == 0:
            continue
        columns["rank"].append(np.full(star_indices.size, rank, dtype=np.int32))
        columns["local_index"].append(star_indices)
        columns["source_id"].append(_read_selected(part, "identity", mask))
        columns["level"].append(_read_selected(part, "level", mask))
        for name in ("position_x", "position_y", "position_z"):
            columns[f"{name}_code"].append(_read_selected(part, name, mask))
        if include_velocity:
            for name in ("velocity_x", "velocity_y", "velocity_z"):
                columns[f"{name}_code"].append(_read_selected(part, name, mask))
        else:
            for name in ("velocity_x_code", "velocity_y_code", "velocity_z_code"):
                columns[name].append(np.full(star_indices.size, np.nan, dtype=np.float64))
        columns["mass_code"].append(_read_selected(part, "mass", mask))
        columns["initial_mass_code"].append(_read_selected(part, "initial_mass", mask))
        columns["birth_epoch_code"].append(_read_selected(part, "birth_epoch", mask))
        columns["birth_proper_time_code"].append(_read_selected(part, "birth_proper_time", mask))
        columns["birth_metallicity_mass_fraction"].append(_read_selected(part, "metallicity", mask))
        columns["yield_table_index_code"].append(_read_selected(part, "yield_table_index", mask))

    rows = {
        name: np.concatenate(values) if values else np.empty(0, dtype=np.float64)
        for name, values in columns.items()
    }
    rows["rank"] = rows["rank"].astype(np.int32)
    rows["local_index"] = rows["local_index"].astype(np.int64)
    rows["source_id"] = rows["source_id"].astype(np.int64)
    rows["level"] = rows["level"].astype(np.int32)
    if len(np.unique(rows["source_id"])) != len(rows["source_id"]):
        raise ValueError("stellar particle identities are not globally unique")
    order = np.argsort(rows["source_id"], kind="stable")
    rows = {name: values[order] for name, values in rows.items()}
    if len(rows["source_id"]) != first.nstar_tot:
        raise ValueError(
            f"stellar type count {len(rows['source_id'])} disagrees with nstar_tot {first.nstar_tot}"
        )

    for name in (
        "position_x_code",
        "position_y_code",
        "position_z_code",
        "mass_code",
        "initial_mass_code",
        "birth_epoch_code",
        "birth_proper_time_code",
        "birth_metallicity_mass_fraction",
        "yield_table_index_code",
    ):
        if not np.isfinite(rows[name]).all():
            raise ValueError(f"stellar field {name} contains non-finite values")
    position = np.stack(
        [rows[f"position_{axis}_code"] for axis in "xyz"], axis=1
    )
    if np.any(position < 0.0) or np.any(position > 1.0):
        raise ValueError("stellar positions lie outside the normalized RAMSES domain")
    if np.any(rows["mass_code"] <= 0.0) or np.any(rows["initial_mass_code"] <= 0.0):
        raise ValueError("stellar masses must be positive")
    negative_metallicity = rows["birth_metallicity_mass_fraction"] < 0.0
    if np.any(rows["birth_metallicity_mass_fraction"] < -1.0e-12):
        raise ValueError("stellar birth metallicities contain a materially negative value")
    negative_metallicity_count = int(np.count_nonzero(negative_metallicity))
    if negative_metallicity_count:
        rows["birth_metallicity_mass_fraction"] = np.maximum(
            rows["birth_metallicity_mass_fraction"], 0.0
        )
    if np.any(rows["yield_table_index_code"] < 0.0):
        raise ValueError("stellar yield-table indices must be non-negative")

    info_values = _read_info(info)
    current_proper_time_code = _proper_time_code_from_aexp(info_values, w0, wa)
    proper_time_unit_s = (
        info_values["unit_t"] / info_values["aexp"] ** 2
        if _is_cosmological(info_values)
        else info_values["unit_t"]
    )
    age_myr = (current_proper_time_code - rows["birth_proper_time_code"]) * proper_time_unit_s / SECONDS_PER_MYR
    if np.any(age_myr < -1.0e-6):
        raise ValueError("a stellar birth proper time lies after the checkpoint")
    rows["age_myr"] = np.maximum(age_myr, 0.0)
    rows["position_x_proper_cm"] = rows["position_x_code"] * info_values["unit_l"]
    rows["position_y_proper_cm"] = rows["position_y_code"] * info_values["unit_l"]
    rows["position_z_proper_cm"] = rows["position_z_code"] * info_values["unit_l"]
    mass_unit_msun = info_values["unit_d"] * info_values["unit_l"] ** 3 / MSUN_G
    rows["mass_msun"] = rows["mass_code"] * mass_unit_msun
    rows["initial_mass_msun"] = rows["initial_mass_code"] * mass_unit_msun
    if include_velocity:
        velocity_unit_cm_s = info_values["unit_l"] / info_values["unit_t"]
        rows["velocity_x_cm_s"] = rows["velocity_x_code"] * velocity_unit_cm_s
        rows["velocity_y_cm_s"] = rows["velocity_y_code"] * velocity_unit_cm_s
        rows["velocity_z_cm_s"] = rows["velocity_z_code"] * velocity_unit_cm_s

    stats_fields = (
        "mass_msun",
        "initial_mass_msun",
        "age_myr",
        "birth_metallicity_mass_fraction",
        "yield_table_index_code",
    )
    manifest: dict[str, object] = {
        "record_type": "native_ramses_stellar_catalogue",
        "schema_version": 1,
        "status": "complete_native_stellar_metadata_extracted",
        "output_dir": str(output),
        "info_path": str(info),
        "rank_files": len(parts),
        "native_format": "cuRamses_particle_binary_v1",
        "header": {
            "nstar_tot": int(first.nstar_tot),
            "nstar_tot_record_bytes": int(_scan_header_nstar_record_length(first.path)),
            "nsink": int(first.nsink),
        },
        "particle_type_selection": {"field": "ptypep", "star_code": 1, "sink_code": 2},
        "source_catalogue": {
            "photon_luminosity_assigned": False,
            "sed_model_assigned": False,
            "source_kind": "star_metadata_only",
        },
        "time_semantics": {
            "birth_epoch_column": "birth_epoch_code",
            "birth_epoch_kind": "RAMSES conformal/supercomoving time value; not a scale factor",
            "birth_proper_time_column": "birth_proper_time_code",
            "current_proper_time_code": current_proper_time_code,
            "cosmological": _is_cosmological(info_values),
            "proper_time_unit_s": proper_time_unit_s,
            "age_formula": "(current_proper_time_code - birth_proper_time_code) * proper_time_unit_s",
            "cosmological_model": {"w0": w0, "wa": wa},
        },
        "unit_conversion": {
            "position_code": "normalized RAMSES user coordinate in [0, 1]",
            "position_proper_cm": float(info_values["unit_l"]),
            "mass_code": "RAMSES particle mass unit",
            "mass_unit_msun": mass_unit_msun,
            "velocity_code": "RAMSES particle velocity unit",
            "velocity_unit_cm_s": float(info_values["unit_l"] / info_values["unit_t"]),
            "birth_metallicity": "native zp value; mass-fraction normalization not converted to solar units",
            "yield_table_index": "native indtab checkpoint/progress value",
        },
        "field_provenance": {
            "positions": "output_part.f90 xp",
            "velocities": "output_part.f90 vp",
            "current_mass": "output_part.f90 mp",
            "identity": "output_part.f90 idp",
            "birth_epoch": "output_part.f90 tp",
            "birth_metallicity": "output_part.f90 zp",
            "birth_proper_time": "output_part.f90 tpp",
            "initial_mass": "output_part.f90 mp0",
            "yield_table_index": "output_part.f90 indtab",
        },
        "star_count": int(len(rows["source_id"])),
        "field_stats": {name: _stats(rows[name]) for name in stats_fields},
        "sanitization": {
            "negative_birth_metallicity_clamp_tolerance": 1.0e-12,
            "negative_birth_metallicity_clamped_count": negative_metallicity_count,
        },
        "scientific_readiness": {
            "stellar_particle_metadata": True,
            "stellar_source_catalogue": False,
            "photon_luminosity_ledger": False,
            "reason": "native star metadata is decoded; SED, escape fraction, and grouped photon luminosity remain explicit downstream inputs",
        },
    }
    return rows, manifest


def _scan_header_nstar_record_length(path: Path) -> int:
    with path.open("rb") as stream:
        for index, name in enumerate(_HEADER_RECORD_NAMES):
            marker = _read_marker(stream, name)
            if index == 4:
                return marker
            stream.seek(marker, 1)
            trailer = stream.read(4)
            if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
                raise ValueError(f"{path}: header {name} trailer mismatch")
    raise RuntimeError("nstar_tot header record was not found")


def _csv_columns(include_velocity: bool) -> tuple[str, ...]:
    columns = (
        "rank",
        "local_index",
        "source_id",
        "source_kind",
        "level",
        "position_x_code",
        "position_y_code",
        "position_z_code",
        "position_x_proper_cm",
        "position_y_proper_cm",
        "position_z_proper_cm",
        "mass_code",
        "mass_msun",
        "initial_mass_code",
        "initial_mass_msun",
        "birth_epoch_code",
        "birth_proper_time_code",
        "age_myr",
        "birth_metallicity_mass_fraction",
        "yield_table_index_code",
    )
    if include_velocity:
        columns += (
            "velocity_x_code",
            "velocity_y_code",
            "velocity_z_code",
            "velocity_x_cm_s",
            "velocity_y_cm_s",
            "velocity_z_cm_s",
        )
    return columns


def _format_csv_value(value: object) -> object:
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return format(float(value), ".17g")
    return value


def write_catalogue(
    output_path: str | Path,
    manifest_path: str | Path,
    rows: dict[str, np.ndarray],
    manifest: dict[str, object],
    *,
    include_velocity: bool,
) -> None:
    output = Path(output_path).expanduser().resolve()
    metadata = Path(manifest_path).expanduser().resolve()
    if output.exists() or metadata.exists():
        raise FileExistsError(f"refusing to overwrite {output if output.exists() else metadata}")
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = _csv_columns(include_velocity)
    with output.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for index in range(len(rows["source_id"])):
            record = {name: _format_csv_value(rows[name][index]) for name in columns if name != "source_kind"}
            record["source_kind"] = "star"
            writer.writerow(record)
    manifest["catalogue_csv"] = str(output)
    manifest["catalogue_csv_sha256"] = hashlib.sha256(output.read_bytes()).hexdigest()
    metadata.parent.mkdir(parents=True, exist_ok=True)
    with metadata.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--info", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--w0", type=float, default=-1.0)
    parser.add_argument("--wa", type=float, default=0.0)
    parser.add_argument(
        "--include-velocity",
        action="store_true",
        help="also decode particle velocities; the default source hand-off omits them",
    )
    args = parser.parse_args()
    rows, manifest = extract_catalogue(
        args.output_dir,
        args.info,
        w0=args.w0,
        wa=args.wa,
        include_velocity=args.include_velocity,
    )
    write_catalogue(
        args.output,
        args.manifest,
        rows,
        manifest,
        include_velocity=args.include_velocity,
    )
    print(
        "NATIVE_STELLAR_CATALOGUE_OK "
        f"stars={len(rows['source_id'])} output={args.output} manifest={args.manifest}"
    )


if __name__ == "__main__":
    main()
