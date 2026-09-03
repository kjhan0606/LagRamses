#!/usr/bin/env python3
"""Audit the staged Keegans SNIa supplementary tables without converting them."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

from fp2_provenance import project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "assets" / "review_only" / "fp2_snia" / "keegans2023"
DEFAULT_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_keegans2023_review_v1.json"
EXPECTED_FILES = (
    "decayed_brox_yields_new_latex_table.tex",
    "decayed_dean_yields_new_latex_table.tex",
    "decayed_shen_yields_new_latex_table.tex",
)
EXPECTED_22NE = (0.0, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2, 1.4e-2, 2.0e-2, 5.0e-2, 1.0e-1)
PROJECT_ELEMENTS = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
DATA_ROW_RE = re.compile(r"^\s*([A-Z][a-z]?)\s*&\s*(\d+)\s*&\s*(\d+)\s*&(.+?)\\\\\s*$")


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_manifest(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("manifest must be a JSON object")
    return value


def _header_columns(lines: list[str]) -> tuple[list[float], str | None]:
    for line in lines:
        if "22" not in line or "&" not in line:
            continue
        parts = [part.strip() for part in line.split("&")]
        if len(parts) < 4:
            continue
        first = re.search(r"=\s*([0-9.eE+-]+)", parts[3])
        if first is None:
            continue
        tokens = [first.group(1)]
        tokens.extend(part.split()[0] for part in parts[4:])
        try:
            return [float(token) for token in tokens], line.strip()
        except ValueError:
            return [], line.strip()
    return [], None


def _parse_table(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    columns, header = _header_columns(lines)
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    symbols_to_z: dict[str, int] = {}
    coordinates: set[tuple[str, int, int]] = set()
    element_stats: dict[str, dict[str, int]] = {}
    for line_number, line in enumerate(lines, start=1):
        match = DATA_ROW_RE.match(line)
        if match is None:
            continue
        symbol, z_text, a_text, values_text = match.groups()
        parts = [part.strip() for part in values_text.split("&")]
        try:
            z = int(z_text)
            a = int(a_text)
            values = [float(part) for part in parts]
        except ValueError:
            failures.append(f"line_{line_number}_non_numeric")
            continue
        if len(values) != len(EXPECTED_22NE):
            failures.append(f"line_{line_number}_column_count_{len(values)}")
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            failures.append(f"line_{line_number}_negative_or_nonfinite_yield")
        stats = element_stats.setdefault(
            symbol,
            {"isotope_row_count": 0, "yield_cell_count": 0, "nonzero_yield_cell_count": 0},
        )
        stats["isotope_row_count"] += 1
        stats["yield_cell_count"] += len(values)
        stats["nonzero_yield_cell_count"] += sum(value != 0.0 for value in values)
        if a <= z or z <= 0:
            failures.append(f"line_{line_number}_invalid_nuclide_coordinate")
        if symbol in symbols_to_z and symbols_to_z[symbol] != z:
            failures.append(f"line_{line_number}_symbol_atomic_number_conflict")
        symbols_to_z[symbol] = z
        coordinate = (symbol, z, a)
        if coordinate in coordinates:
            failures.append(f"line_{line_number}_duplicate_isotope")
        coordinates.add(coordinate)
        rows.append({"symbol": symbol, "z": z, "a": a, "line": line_number})
    return {
        "path": project_relative(path),
        "sha256": _sha256(path),
        "header": header,
        "header_column_count": len(columns),
        "header_22ne": columns,
        "header_22ne_matches_expected": len(columns) == len(EXPECTED_22NE)
        and all(math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-15) for left, right in zip(columns, EXPECTED_22NE)),
        "row_count": len(rows),
        "isotope_coordinate_count": len(coordinates),
        "elements": sorted(symbols_to_z),
        "z_by_element": dict(sorted(symbols_to_z.items())),
        "element_stats": {
            element: {
                **stats,
                "all_yield_cells_explicit_zero": stats["nonzero_yield_cell_count"] == 0,
            }
            for element, stats in sorted(element_stats.items())
        },
        "format_failures": failures,
        "format_integrity_passed": not failures
        and len(columns) == len(EXPECTED_22NE)
        and all(math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-15) for left, right in zip(columns, EXPECTED_22NE)),
    }


def audit_format(root: Path = DEFAULT_ROOT, manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    failures: list[str] = []
    if not manifest_path.is_file():
        failures.append("manifest_missing")
        manifest: dict[str, Any] = {}
    else:
        try:
            manifest = _read_manifest(manifest_path)
        except (OSError, ValueError, json.JSONDecodeError):
            failures.append("manifest_invalid")
            manifest = {}
    if manifest.get("status") != "review_only_asset":
        failures.append("manifest_not_review_only")
    entries = manifest.get("files") if manifest else []
    manifest_names = [entry.get("path") for entry in entries if isinstance(entry, dict)]
    if tuple(manifest_names) != EXPECTED_FILES:
        failures.append("manifest_file_order_or_identity_mismatch")

    file_reports = []
    coordinate_sets: list[set[tuple[str, int, int]]] = []
    for name in EXPECTED_FILES:
        path = root / name
        if not path.is_file():
            failures.append(f"missing_file:{name}")
            continue
        report = _parse_table(path)
        file_reports.append(report)
        if not report["format_integrity_passed"]:
            failures.append(f"format_failure:{name}")
        if report["row_count"] != 70:
            failures.append(f"unexpected_row_count:{name}")
        coordinate_sets.append({
            (row["symbol"], row["z"], row["a"])
            for row in _rows_for_coordinate_audit(path)
        })
    if coordinate_sets and any(coords != coordinate_sets[0] for coords in coordinate_sets[1:]):
        failures.append("file_coordinate_sets_differ")
    source_elements = sorted({element for report in file_reports for element in report["elements"]})
    missing_project_elements = [element for element in PROJECT_ELEMENTS if element not in source_elements]
    source_element_stats: dict[str, dict[str, int]] = {}
    for report in file_reports:
        for element, stats in report["element_stats"].items():
            aggregate = source_element_stats.setdefault(
                element,
                {"isotope_row_count": 0, "yield_cell_count": 0, "nonzero_yield_cell_count": 0},
            )
            for key in aggregate:
                aggregate[key] += stats[key]
    project_element_presence: dict[str, dict[str, Any]] = {}
    explicit_zero_project_elements: list[str] = []
    for element in PROJECT_ELEMENTS:
        stats = source_element_stats.get(element)
        if stats is None:
            project_element_presence[element] = {
                "status": "absent_from_source_isotope_rows",
                "isotope_row_count": 0,
                "yield_cell_count": 0,
                "nonzero_yield_cell_count": 0,
            }
        elif stats["nonzero_yield_cell_count"] == 0:
            project_element_presence[element] = {
                "status": "explicit_zero_only",
                **stats,
            }
            explicit_zero_project_elements.append(element)
        else:
            project_element_presence[element] = {
                "status": "present_with_nonzero_or_mixed_data",
                **stats,
            }
    conversion_blockers = [
        "source package is an isotope-valued table, not normalized per-event runtime rows",
        "isotope-to-project-element aggregation and decay horizon are not approved",
        "event energy and momentum are absent from the staged supplementary files",
        "returned mass and normal-SNIa zero-remnant ownership are not closed by these files",
    ]
    if missing_project_elements:
        conversion_blockers.append(
            "project tracked elements absent from the staged package: " + ", ".join(missing_project_elements)
        )
    return {
        "schema": "snrt-fp2-snia-keegans-format-audit",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_source_format_passed" if not failures else "blocked_source_format_integrity",
        "format_integrity_passed": not failures,
        "canonical_conversion_allowed": False,
        "runtime_activation_allowed": False,
        "root": project_relative(root),
        "manifest": project_relative(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "file_reports": file_reports,
        "file_count": len(file_reports),
        "source_element_count": len(source_elements),
        "source_elements": source_elements,
        "project_elements": list(PROJECT_ELEMENTS),
        "missing_project_elements": missing_project_elements,
        "project_element_presence": project_element_presence,
        "missing_project_elements_are_absent_isotope_rows": all(
            project_element_presence[element]["status"] == "absent_from_source_isotope_rows"
            for element in missing_project_elements
        ),
        "explicit_zero_project_elements": explicit_zero_project_elements,
        "inferred_zero_policy": "never infer absent project elements as zero during conversion",
        "conversion_blockers": conversion_blockers,
        "audit_failures": failures,
        "audit_code_sha256": _sha256(TOOL_PATH),
        "interpretation": "The staged bytes and table grammar are intact, but this report explicitly does not authorize isotope conversion or SNIa runtime activation.",
    }


def _rows_for_coordinate_audit(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = DATA_ROW_RE.match(line)
        if match is not None:
            symbol, z_text, a_text, _ = match.groups()
            rows.append({"symbol": symbol, "z": int(z_text), "a": int(a_text)})
    return rows


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    report = audit_format(args.root, args.manifest)
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["audit_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
