#!/usr/bin/env python3
"""Audit native cuRAMSES particle records without copying the large payload.

The stopped comparison output uses a compact one-byte ``ptypep`` record and
an eight-byte ``nstar_tot`` header field.  This audit streams only that type
record, validates all Fortran record lengths, and does not create a particle
catalogue or infer photon luminosities.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import struct


_PART_RE = re.compile(r"part_(\d{5})\.out\d{5}$")
_DATA_BYTES_PER_PARTICLE = (8, 8, 8, 8, 8, 8, 8, 8, 4, 1, 8, 8, 8, 8, 8, 8)
_DATA_FIELD_NAMES = (
    "position_x",
    "position_y",
    "position_z",
    "velocity_x",
    "velocity_y",
    "velocity_z",
    "mass",
    "identity",
    "level",
    "ptypep",
    "potential",
    "birth_epoch",
    "metallicity",
    "birth_proper_time",
    "initial_mass",
    "yield_table_index",
)


def _read_record(stream, label: str) -> tuple[int, int]:
    marker_bytes = stream.read(4)
    if len(marker_bytes) != 4:
        raise ValueError(f"unexpected EOF before {label} record marker")
    marker = struct.unpack("<i", marker_bytes)[0]
    if marker < 0:
        raise ValueError(f"negative {label} record length {marker}")
    return marker, stream.tell()


def _finish_record(stream, marker: int, payload_offset: int, label: str) -> None:
    stream.seek(payload_offset + marker)
    trailer = stream.read(4)
    if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
        raise ValueError(f"{label} record trailer mismatch")


def _read_payload(stream, marker: int, payload_offset: int, label: str) -> bytes:
    stream.seek(payload_offset)
    payload = stream.read(marker)
    if len(payload) != marker:
        raise ValueError(f"unexpected EOF in {label} record")
    trailer = stream.read(4)
    if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
        raise ValueError(f"{label} record trailer mismatch")
    return payload


def _unpack_int(payload: bytes, label: str, allowed_sizes: tuple[int, ...] = (4,)) -> int:
    if len(payload) not in allowed_sizes:
        raise ValueError(f"{label} record has {len(payload)} bytes, expected {allowed_sizes}")
    return int.from_bytes(payload, "little", signed=True)


def _unpack_double(payload: bytes, label: str) -> float:
    if len(payload) != 8:
        raise ValueError(f"{label} record has {len(payload)} bytes, expected 8")
    return struct.unpack("<d", payload)[0]


def _scan_part_file(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        header_lengths: list[int] = []
        header_payloads: list[bytes] = []
        for index in range(8):
            marker, offset = _read_record(stream, f"header_{index + 1}")
            header_lengths.append(marker)
            header_payloads.append(_read_payload(stream, marker, offset, f"header_{index + 1}"))

        ncpu = _unpack_int(header_payloads[0], "ncpu")
        ndim = _unpack_int(header_payloads[1], "ndim")
        npart = _unpack_int(header_payloads[2], "npart")
        nstar_tot = _unpack_int(header_payloads[4], "nstar_tot", allowed_sizes=(4, 8))
        mstar_tot = _unpack_double(header_payloads[5], "mstar_tot")
        mstar_lost = _unpack_double(header_payloads[6], "mstar_lost")
        nsink = _unpack_int(header_payloads[7], "nsink")
        if npart < 0:
            raise ValueError(f"{path}: negative npart {npart}")

        data_record_lengths: list[int] = []
        ptype_counts: Counter[int] = Counter()
        for field_index, (field_name, bytes_per_particle) in enumerate(
            zip(_DATA_FIELD_NAMES, _DATA_BYTES_PER_PARTICLE, strict=True)
        ):
            marker, offset = _read_record(stream, field_name)
            expected = npart * bytes_per_particle
            if marker != expected:
                raise ValueError(
                    f"{path}: {field_name} record has {marker} bytes, expected {expected}"
                )
            data_record_lengths.append(marker)
            if field_name == "ptypep":
                stream.seek(offset)
                remaining = marker
                while remaining:
                    chunk = stream.read(min(1024 * 1024, remaining))
                    if not chunk:
                        raise ValueError(f"{path}: unexpected EOF in ptypep record")
                    for value in chunk:
                        ptype_counts[value if value < 128 else value - 256] += 1
                    remaining -= len(chunk)
                trailer = stream.read(4)
                if len(trailer) != 4 or struct.unpack("<i", trailer)[0] != marker:
                    raise ValueError(f"{path}: ptypep record trailer mismatch")
            else:
                _finish_record(stream, marker, offset, field_name)

        if stream.tell() != path.stat().st_size:
            raise ValueError(
                f"{path}: trailing bytes after native particle records: "
                f"{path.stat().st_size - stream.tell()}"
            )

    return {
        "path": str(path.resolve()),
        "npart": npart,
        "ncpu": ncpu,
        "ndim": ndim,
        "nstar_tot": nstar_tot,
        "nstar_tot_record_bytes": len(header_payloads[4]),
        "mstar_tot": mstar_tot,
        "mstar_lost": mstar_lost,
        "nsink": nsink,
        "header_record_lengths": header_lengths,
        "data_record_lengths": data_record_lengths,
        "ptype_counts": {str(key): value for key, value in sorted(ptype_counts.items())},
        "file_bytes": path.stat().st_size,
    }


def audit_output(output_dir: Path) -> dict[str, object]:
    output_dir = output_dir.expanduser().resolve()
    if not output_dir.is_dir():
        raise FileNotFoundError(output_dir)
    matches = sorted(output_dir.glob("part_*.out*"))
    part_files = [path for path in matches if _PART_RE.fullmatch(path.name)]
    if not part_files:
        raise FileNotFoundError(f"no native part files in {output_dir}")

    records = [_scan_part_file(path) for path in part_files]
    first = records[0]
    if any(record["ncpu"] != first["ncpu"] for record in records):
        raise ValueError("native particle files disagree on ncpu")
    if any(record["ndim"] != first["ndim"] for record in records):
        raise ValueError("native particle files disagree on ndim")
    if any(record["header_record_lengths"] != first["header_record_lengths"] for record in records):
        raise ValueError("native particle header record widths differ between ranks")
    if any(record["data_record_lengths"][9] != record["npart"] for record in records):
        raise ValueError("native ptypep records do not match local particle counts")

    ptype_counts: Counter[int] = Counter()
    for record in records:
        for key, value in record["ptype_counts"].items():
            ptype_counts[int(key)] += int(value)
    nstar_headers = sorted({int(record["nstar_tot"]) for record in records})
    nsink_headers = sorted({int(record["nsink"]) for record in records})
    npart_total = sum(int(record["npart"]) for record in records)
    expected_ranks = int(first["ncpu"])
    if len(records) != expected_ranks:
        raise ValueError(f"found {len(records)} part files, expected {expected_ranks}")

    output_header = output_dir / f"header_{_PART_RE.fullmatch(part_files[0].name).group(1)}.txt"
    header_text = output_header.read_text(encoding="utf-8", errors="replace") if output_header.is_file() else ""
    expected_star = None
    expected_dm = None
    expected_sink = None
    for label in ("star", "dark matter", "sink"):
        match = re.search(rf"(?mi)^\s*Total number of {label} particles\s*$\s*([0-9]+)", header_text)
        if match:
            value = int(match.group(1))
            if label == "star": expected_star = value
            elif label == "dark matter": expected_dm = value
            else: expected_sink = value

    code_counts = {str(key): int(value) for key, value in sorted(ptype_counts.items())}
    structural_ok = (
        len(records) == expected_ranks
        and npart_total == sum(ptype_counts.values())
        and nstar_headers == [ptype_counts.get(1, 0)]
        and nsink_headers == [ptype_counts.get(2, 0)]
        and (expected_star is None or expected_star == ptype_counts.get(1, 0))
        and (expected_dm is None or expected_dm == ptype_counts.get(0, 0))
        and (expected_sink is None or expected_sink == ptype_counts.get(2, 0))
    )
    return {
        "record_type": "native_ramses_particle_audit",
        "audit_version": 1,
        "status": "complete_native_particle_metadata_audited" if structural_ok else "structural_check_failed",
        "output_dir": str(output_dir),
        "rank_files": len(records),
        "expected_ranks": expected_ranks,
        "native_format": "cuRamses_particle_binary_v1",
        "endianness": "little",
        "header_record_lengths": first["header_record_lengths"],
        "data_fields": list(_DATA_FIELD_NAMES),
        "data_bytes_per_particle": list(_DATA_BYTES_PER_PARTICLE),
        "local_particle_total": npart_total,
        "nstar_tot_headers": nstar_headers,
        "nsink_headers": nsink_headers,
        "ptype_codes": {
            "0": "DM",
            "1": "STAR",
            "2": "SINK",
        },
        "ptype_counts": code_counts,
        "header_totals": {
            "dark_matter": expected_dm,
            "stars": expected_star,
            "sinks": expected_sink,
        },
        "payload_bytes_stat_sum": sum(int(record["file_bytes"]) for record in records),
        "hash_policy": "native particle payload not hashed; headers and ptype record structure audited",
        "yt_frontend_compatibility": {
            "status": "blocked_for_default_handler",
            "reason": "yt 4.4.2 expects a 4-byte nstar_tot header field, while this native LONGINT output stores 8 bytes",
        },
        "scientific_readiness": {
            "particle_type_roster": structural_ok,
            "stellar_source_catalogue": False,
            "photon_luminosity_ledger": False,
            "reason": "type counts are certified; positions, ages, masses, SED, and source normalization remain a separate explicit reader/ledger step",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = audit_output(args.output_dir)
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
