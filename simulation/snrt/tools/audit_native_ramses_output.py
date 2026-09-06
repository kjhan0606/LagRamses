#!/usr/bin/env python3
"""Audit a native RAMSES output without reading the large binary payloads.

This is a provenance/preflight tool, not a converter.  It checks the
completion marker, text metadata, expected per-rank component counts, build
identity, and the physics/source-ledger inventory.  Large ``amr_*``,
``hydro_*``, ``grav_*`` and ``part_*`` files are only stat'ed; their contents
are deliberately not hashed or decoded here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable


_OUTPUT_RE = re.compile(r"(?:info|header|namelist|compilation|makefile)_(\d{5})")
_KEY_VALUE_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _first_output_number(output_dir: Path) -> str:
    candidates: list[str] = []
    for path in sorted(output_dir.iterdir()):
        match = _OUTPUT_RE.search(path.name)
        if match:
            candidates.append(match.group(1))
    if not candidates:
        raise ValueError(f"cannot infer output number from {output_dir}")
    return sorted(set(candidates))[0]


def _parse_scalar_metadata(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text.splitlines():
        match = _KEY_VALUE_RE.match(line)
        if match:
            values[match.group(1)] = match.group(2)
    return values


def _parse_compilation(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    labels = {
        "compile date": "compile_date",
        "patch dir": "patch_dir",
        "remote repo": "remote_repo",
        "local branch": "local_branch",
        "last commit": "last_commit",
    }
    for line in _read_text(path).splitlines():
        if "=" not in line:
            continue
        label, value = (part.strip() for part in line.split("=", 1))
        key = labels.get(label)
        if key:
            values[key] = value
    return values


def _parse_inventory(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in _read_text(path).splitlines():
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key:
            values[key] = value
    return values


def _parse_hydro_descriptor(path: Path) -> dict[str, object]:
    text = _read_text(path)
    nvar_match = re.search(r"(?mi)^\s*nvar\s*=\s*(\d+)", text)
    fields = {
        int(match.group(1)): match.group(2).strip()
        for match in re.finditer(r"(?mi)^\s*variable\s*#\s*(\d+)\s*:\s*(.*?)\s*$", text)
    }
    return {
        "declared_nvar": int(nvar_match.group(1)) if nvar_match else None,
        "field_count": len(fields),
        "fields": {str(index): name for index, name in sorted(fields.items())},
    }


def _find_flag(text: str, name: str) -> str | None:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*(?:\?\s*)?=\s*([^\s#]+)", text)
    return match.group(1) if match else None


def _component_files(output_dir: Path, output_number: str, prefix: str) -> list[Path]:
    pattern = f"{prefix}_{output_number}.out*"
    return sorted(path for path in output_dir.glob(pattern) if path.is_file())


def _relative_sizes(paths: Iterable[Path]) -> tuple[int, int]:
    files = list(paths)
    return len(files), sum(path.stat().st_size for path in files)


def _parse_runtime_log(path: Path | None) -> dict[str, object] | None:
    if path is None:
        return None
    text = _read_text(path)

    def integer_after(label: str) -> int | None:
        match = re.search(rf"(?mi)^\s*{re.escape(label)}\s*=\s*(\d+)", text)
        return int(match.group(1)) if match else None

    def timing_after(label: str) -> float | None:
        match = re.search(rf"(?mi)^\s*{re.escape(label)}\s*:\s*([0-9.eE+-]+)\s*s", text)
        return float(match.group(1)) if match else None

    mode_match = re.search(r"(?mi)^\s*Stellar feedback mode\s*:\s*(\S+)", text)
    return {
        "path": str(path.expanduser().resolve()),
        "sha256": _sha256(path),
        "phase0_enabled_marker": bool(re.search(r"(?mi)^\s*Phase 0 stellar enrichment enabled\s*$", text)),
        "phase0_table_rows": integer_after("table rows"),
        "total_metal_field": integer_after("total-metal field"),
        "first_element_field": integer_after("first element field"),
        "stellar_feedback_mode_line": mode_match.group(1) if mode_match else None,
        "snrt_advance_seconds": timing_after("snrt_advance"),
        "snrt_diagnose_seconds": timing_after("snrt_diagnose"),
    }


def audit_output(output_dir: Path, run_log: Path | None = None) -> dict[str, object]:
    """Return a JSON-serializable audit record for one native output."""

    output_dir = output_dir.expanduser().resolve()
    if not output_dir.is_dir():
        raise FileNotFoundError(f"native output directory does not exist: {output_dir}")

    output_number = _first_output_number(output_dir)
    required = {
        "complete_marker": output_dir / "COMPLETE",
        "info": output_dir / f"info_{output_number}.txt",
        "header": output_dir / f"header_{output_number}.txt",
        "namelist": output_dir / "namelist.txt",
        "compilation": output_dir / "compilation.txt",
        "makefile": output_dir / "makefile.txt",
        "hydro_descriptor": output_dir / "hydro_file_descriptor.txt",
        "physics_inventory": output_dir / f"resolved_physics_inventory_{output_number}.txt",
    }
    missing = [name for name, path in required.items() if not path.is_file()]

    info = _parse_scalar_metadata(_read_text(required["info"])) if required["info"].is_file() else {}
    compilation = (
        _parse_compilation(required["compilation"])
        if required["compilation"].is_file()
        else {}
    )
    inventory = (
        _parse_inventory(required["physics_inventory"])
        if required["physics_inventory"].is_file()
        else {}
    )
    makefile_text = _read_text(required["makefile"]) if required["makefile"].is_file() else ""
    hydro_descriptor = (
        _parse_hydro_descriptor(required["hydro_descriptor"])
        if required["hydro_descriptor"].is_file()
        else {"declared_nvar": None, "field_count": 0, "fields": {}}
    )

    try:
        expected_ranks = int(info["ncpu"])
    except (KeyError, ValueError):
        expected_ranks = None

    components: dict[str, dict[str, int | bool | None]] = {}
    payload_bytes = 0
    for prefix in ("amr", "hydro", "grav", "part", "sink", "rt"):
        files = _component_files(output_dir, output_number, prefix)
        count, size = _relative_sizes(files)
        payload_bytes += size
        components[prefix] = {
            "file_count": count,
            "bytes": size,
            "rank_count_matches_ncpu": expected_ranks is not None and count == expected_ranks,
            "present": bool(files),
        }

    rank_mismatches = [
        prefix
        for prefix, component in components.items()
        if component["present"] and not component["rank_count_matches_ncpu"]
    ]
    required_components = ("amr", "hydro", "grav", "part")
    missing_native_components = [
        prefix
        for prefix in required_components
        if not components[prefix]["present"]
        or not components[prefix]["rank_count_matches_ncpu"]
    ]
    descriptor_failures = []
    if hydro_descriptor["declared_nvar"] is None:
        descriptor_failures.append("hydro descriptor has no nvar")
    elif hydro_descriptor["field_count"] != hydro_descriptor["declared_nvar"]:
        descriptor_failures.append("hydro descriptor field count differs from nvar")
    structural_ok = (
        not missing
        and not rank_mismatches
        and not missing_native_components
        and not descriptor_failures
    )

    record: dict[str, object] = {
        "record_type": "native_ramses_output_audit",
        "audit_version": 1,
        "output_dir": str(output_dir),
        "output_number": output_number,
        "status": "complete_native_metadata_audited" if structural_ok else "structural_check_failed",
        "completion_marker": required["complete_marker"].is_file(),
        "missing_required_metadata": missing,
        "expected_ranks": expected_ranks,
        "components": components,
        "native_payload_bytes_stat_sum": payload_bytes,
        "hash_policy": "text_metadata_only; native binary payloads not hashed",
        "metadata_sha256": {
            name: _sha256(path)
            for name, path in required.items()
            if path.is_file() and name != "complete_marker"
        },
        "info": info,
        "compilation": compilation,
        "makefile_flags": {
            name: _find_flag(makefile_text, name)
            for name in ("NVAR", "SOLVER", "PHASE0_STELLAR_ENRICHMENT", "F90")
        },
        "hydro_descriptor": hydro_descriptor,
        "physics_inventory": inventory,
        "runtime_log_evidence": _parse_runtime_log(run_log),
        "source_ledger_available": inventory.get("force_source_ledger_status") == "available",
        "sink_info_available": inventory.get("sink_info_file", "none") not in {"none", "unavailable"},
        "rt_runtime_evidence": "unverified_from_native_output_metadata",
        "rank_count_mismatches": rank_mismatches,
        "missing_native_components": missing_native_components,
        "hydro_descriptor_failures": descriptor_failures,
        "scientific_readiness": {
            "hydro_native_checkpoint": structural_ok,
            "direct_snrt_input": False,
            "reason": [
                "native binary checkpoint requires explicit RAMSES-to-canonical field mapping",
                "source ledger and sink metadata are unavailable in this output",
                "dust and non-equilibrium chemistry fields are not certified by this audit",
            ],
        },
    }
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="one native RAMSES output_XXXXX directory")
    parser.add_argument("--run-log", type=Path, help="optional parent run.log for runtime markers")
    parser.add_argument("--output", type=Path, help="optional JSON audit path")
    args = parser.parse_args()

    record = audit_output(args.output_dir, args.run_log)
    rendered = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
