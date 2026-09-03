#!/usr/bin/env python3
"""Audit the review-only HESMA F-P2 SNIa source package.

The HESMA record contains source data that are useful for closing the
composition/profile side of an SNIa event, but this audit deliberately does
not derive or approve a runtime event energy, momentum, or population weight.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile

from fp2_provenance import project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92"
DEFAULT_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"
ARCHIVES = {
    "abundances": "ddt_2013_abundances.zip",
    "density": "ddt_2013_density.zip",
    "isotopes": "ddt_2013_isotopes.zip",
}
PROJECT_ELEMENTS = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
SOLAR_MASS_G = 1.98847e33
PROFILE_REVIEW_THRESHOLD = 0.05
GROSS_PROFILE_ANOMALY_THRESHOLD = 1.0
MODEL_IDS = (
    "n1",
    "n3",
    "n5",
    "n10",
    "n20",
    "n40",
    "n100h",
    "n100",
    "n100l",
    "n150",
    "n200",
    "n300c",
    "n1600",
    "n1600c",
    "n100_z0.01",
)
FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
ABUNDANCE_RE = re.compile(rf"^\s*([a-z]{{1,2}})\s+({FLOAT_RE})\s*$")


def _fingerprint(path: Path) -> tuple[int, str] | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def _composite(records: list[dict[str, Any]]) -> str:
    payload = bytearray()
    for record in sorted(records, key=lambda item: item["path"].encode("utf-8")):
        payload.extend(record["path"].encode("utf-8"))
        payload.extend(b"\0")
        payload.extend(str(record["bytes"]).encode("ascii"))
        payload.extend(b"\0")
        payload.extend(record["sha256"].encode("ascii"))
        payload.extend(b"\n")
    return hashlib.sha256(payload).hexdigest()


def _model_file_names(kind: str, names: list[str]) -> dict[str, str]:
    suffix = f"_{kind}.dat"
    result: dict[str, str] = {}
    for name in names:
        base = Path(name).name
        if not base.startswith("ddt_2013_") or not base.endswith(suffix):
            continue
        model = base[len("ddt_2013_") : -len(suffix)]
        result[model] = name
    return result


def _audit_abundance_file(zfile: ZipFile, name: str) -> tuple[dict[str, Any], list[str]]:
    failures: list[str] = []
    text = zfile.read(name).decode("utf-8")
    section = ""
    stable_elements: dict[str, float] = {}
    for line in text.splitlines():
        lowered = line.strip().lower()
        if lowered.startswith("abundances of radioactive isotopes"):
            section = "radioactive_isotopes"
            continue
        if lowered.startswith("abundances of stable isotopes"):
            section = "stable_isotopes"
            continue
        if lowered.startswith("abundances of stable elements"):
            section = "stable_elements"
            continue
        if section != "stable_elements":
            continue
        match = ABUNDANCE_RE.match(line)
        if match is None:
            continue
        element, value_text = match.groups()
        if element in stable_elements:
            failures.append(f"{name}:duplicate_stable_element:{element}")
        value = float(value_text)
        if not math.isfinite(value) or value < 0.0:
            failures.append(f"{name}:invalid_stable_element:{element}")
        stable_elements[element] = value
    missing = [element.lower() for element in PROJECT_ELEMENTS if element.lower() not in stable_elements]
    if missing:
        failures.append(f"{name}:missing_project_stable_elements:{','.join(missing)}")
    return {
        "path": name,
        "stable_element_count": len(stable_elements),
        "project_elements_present": [element for element in PROJECT_ELEMENTS if element.lower() in stable_elements],
        "missing_project_elements": [element for element in PROJECT_ELEMENTS if element.lower() not in stable_elements],
        "stable_element_mass_msun": sum(stable_elements.values()),
        "project_element_mass_msun": sum(stable_elements.get(element.lower(), 0.0) for element in PROJECT_ELEMENTS),
        "decay_section": "stable_elements_after_about_2_gyr",
    }, failures


def _isotope_element(symbol: str) -> str | None:
    if symbol == "p":
        return "H"
    if symbol in {"n"}:
        return None
    match = re.match(r"^([a-z]{1,2})\d+$", symbol)
    if match is None:
        return None
    raw = match.group(1)
    return raw[0].upper() + raw[1:]


def _audit_isotope_file(zfile: ZipFile, name: str) -> tuple[dict[str, Any], list[str]]:
    failures: list[str] = []
    lines = zfile.read(name).decode("utf-8").splitlines()
    nonempty = [line for line in lines if line.strip()]
    if not nonempty:
        return {"path": name, "row_count": 0, "isotope_count": 0}, [f"{name}:empty"]
    header = nonempty[0].split()
    if len(header) < 5 or header[:4] != ["Velocity", "[km/s]", "Density", "[g/cm^3]"]:
        failures.append(f"{name}:header_mismatch")
        isotope_symbols: list[str] = []
    else:
        isotope_symbols = header[4:]
    if len(set(isotope_symbols)) != len(isotope_symbols):
        failures.append(f"{name}:duplicate_isotope_header")
    source_elements = sorted({element for symbol in isotope_symbols if (element := _isotope_element(symbol))})
    isotope_project_presence = {
        element: (
            "present_via_free_proton_column" if element == "H" and "p" in isotope_symbols else "present"
            if element in source_elements else "absent"
        )
        for element in PROJECT_ELEMENTS
    }
    if any(value == "absent" for value in isotope_project_presence.values()):
        failures.append(f"{name}:project_element_missing_from_isotope_profile")
    rows: list[list[float]] = []
    for line_number, line in enumerate(nonempty[1:], start=2):
        try:
            values = [float(token) for token in line.split()]
        except ValueError:
            failures.append(f"{name}:line_{line_number}_non_numeric")
            continue
        if len(values) != len(isotope_symbols) + 2:
            failures.append(f"{name}:line_{line_number}_column_count_{len(values)}")
            continue
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            failures.append(f"{name}:line_{line_number}_negative_or_nonfinite")
        fraction_sum = sum(values[2:])
        if abs(fraction_sum - 1.0) > 1.0e-8:
            failures.append(f"{name}:line_{line_number}_mass_fraction_sum_{fraction_sum:.17g}")
        rows.append(values)
    if len(rows) > 1 and any(rows[index][0] >= rows[index + 1][0] for index in range(len(rows) - 1)):
        failures.append(f"{name}:velocity_not_strictly_increasing")
    return {
        "path": name,
        "row_count": len(rows),
        "isotope_count": len(isotope_symbols),
        "isotope_header_sha256": hashlib.sha256(" ".join(isotope_symbols).encode("ascii")).hexdigest(),
        "isotope_header_first": isotope_symbols[:10],
        "isotope_header_last": isotope_symbols[-10:],
        "source_elements": source_elements,
        "project_element_presence": isotope_project_presence,
        "mass_fraction_sum_tolerance": 1.0e-8,
        "profile_decay_section": "radioactive_at_end_of_simulation_and_stable_after_about_2_gyr",
    }, failures


def _audit_density_file(zfile: ZipFile, name: str, profile_epoch_s: float) -> tuple[dict[str, Any], list[str]]:
    failures: list[str] = []
    rows: list[tuple[float, float]] = []
    for line_number, line in enumerate(zfile.read(name).decode("utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        if line_number == 1 and line.strip().startswith("Velocity [km/s]"):
            continue
        try:
            values = [float(token) for token in line.split()]
        except ValueError:
            failures.append(f"{name}:line_{line_number}_non_numeric")
            continue
        if len(values) != 2:
            failures.append(f"{name}:line_{line_number}_column_count_{len(values)}")
            continue
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            failures.append(f"{name}:line_{line_number}_negative_or_nonfinite")
        rows.append((values[0], values[1]))
    if len(rows) > 1 and any(rows[index][0] >= rows[index + 1][0] for index in range(len(rows) - 1)):
        failures.append(f"{name}:velocity_not_strictly_increasing")
    profile_mass_g = 0.0
    profile_kinetic_energy_erg = 0.0
    if len(rows) > 1:
        velocities_cm_s = [row[0] * 1.0e5 for row in rows]
        edges_cm_s = [0.0]
        edges_cm_s.extend(
            (velocities_cm_s[index] + velocities_cm_s[index + 1]) / 2.0
            for index in range(len(velocities_cm_s) - 1)
        )
        edges_cm_s.append(
            velocities_cm_s[-1]
            + (velocities_cm_s[-1] - velocities_cm_s[-2]) / 2.0
        )
        for index, (_, density) in enumerate(rows):
            shell_volume = 4.0 * math.pi * (profile_epoch_s ** 3) * (
                edges_cm_s[index + 1] ** 3 - edges_cm_s[index] ** 3
            ) / 3.0
            shell_mass = density * shell_volume
            profile_mass_g += shell_mass
            profile_kinetic_energy_erg += 0.5 * shell_mass * velocities_cm_s[index] ** 2
    return {
        "path": name,
        "row_count": len(rows),
        "velocity_min_km_s": rows[0][0] if rows else None,
        "velocity_max_km_s": rows[-1][0] if rows else None,
        "profile_is_spherical_average": True,
        "profile_epoch_s": profile_epoch_s,
        "profile_mass_estimate_msun": profile_mass_g / SOLAR_MASS_G,
        "profile_kinetic_energy_estimate_erg": profile_kinetic_energy_erg,
        "profile_estimator": "piecewise_constant_density_in_velocity_shells; inner edge v=0 and outer edge extrapolated by half a bin",
        "kinetic_energy_status": "review_only_profile_estimate; not an approved scalar event field",
        "signed_momentum_status": "not determined by spherical average",
    }, failures


def _profile_epochs(record: dict[str, Any]) -> dict[str, float]:
    description = str(record.get("metadata", {}).get("description", ""))
    epochs: dict[str, float] = {}
    for model in MODEL_IDS:
        match = re.search(
            rf"ddt_2013_{re.escape(model)}\s*:\s*({FLOAT_RE})\s*s",
            description,
            flags=re.IGNORECASE,
        )
        if match is not None:
            epochs[model] = float(match.group(1))
    return epochs


def audit_source(root: Path = DEFAULT_ROOT, manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    failures: list[str] = []
    manifest: dict[str, Any] = {}
    if not manifest_path.is_file():
        failures.append("manifest_missing")
    else:
        try:
            manifest = _read_json(manifest_path)
        except ValueError:
            failures.append("manifest_invalid")
    if manifest:
        checks = {
            "schema": "snrt-fp2-snia-hesma-asset-manifest",
            "gate": "F-P2",
            "status": "review_only_asset",
            "candidate_id": "hesma_model_archive_snia_profiles",
            "source_record_id": "yysd4-xap92",
            "asset_root": "assets/review_only/fp2_snia/hesma_yysd4_xap92",
        }
        for key, expected in checks.items():
            if manifest.get(key) != expected:
                failures.append(f"manifest_{key}_mismatch")
        source = manifest.get("source", {})
        if source.get("license") != "CC-BY-4.0":
            failures.append("manifest_license_mismatch")
        if source.get("approval_id") is not None:
            failures.append("source_approval_must_remain_null")
        approval = manifest.get("approval", {})
        if approval.get("approval_id") is not None:
            failures.append("approval_id_must_remain_null")
        if approval.get("canonical_conversion_allowed") is not False:
            failures.append("canonical_conversion_must_remain_disabled")
        if approval.get("runtime_activation_allowed") is not False:
            failures.append("runtime_activation_must_remain_disabled")
    entries = manifest.get("files", []) if manifest else []
    if not isinstance(entries, list) or not entries:
        failures.append("manifest_files_missing")
        entries = []
    records: list[dict[str, Any]] = []
    file_reports: list[dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            failures.append("manifest_file_entry_invalid")
            continue
        relative = Path(entry["path"])
        if relative.is_absolute() or ".." in relative.parts:
            failures.append(f"unsafe_manifest_path:{entry['path']}")
            continue
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError:
            failures.append(f"manifest_path_outside_root:{entry['path']}")
            continue
        observed = _fingerprint(path)
        report = {"path": entry["path"], "exists": observed is not None}
        if observed is None:
            failures.append(f"file_missing:{entry['path']}")
            file_reports.append(report)
            continue
        size, digest = observed
        passed = size == entry.get("bytes") and digest == str(entry.get("sha256", "")).lower()
        report.update({"bytes": size, "sha256": digest, "integrity_passed": passed})
        file_reports.append(report)
        if not passed:
            failures.append(f"fingerprint_mismatch:{entry['path']}")
        else:
            records.append({"path": entry["path"], "bytes": size, "sha256": digest})
    observed_package = _composite(records) if len(records) == len(entries) else None
    if observed_package != manifest.get("package_sha256"):
        failures.append("package_fingerprint_mismatch")

    record_path = root / "record.json"
    record: dict[str, Any] = {}
    if record_path.is_file():
        try:
            record = _read_json(record_path)
        except ValueError:
            failures.append("record_invalid")
    if record:
        if record.get("id") != "yysd4-xap92":
            failures.append("record_id_mismatch")
        rights = record.get("metadata", {}).get("rights", [])
        if not any(isinstance(item, dict) and item.get("id") == "cc-by-4.0" for item in rights):
            failures.append("record_cc_by_right_missing")
        if record.get("access", {}).get("record") != "public" or record.get("access", {}).get("files") != "public":
            failures.append("record_not_public")
        description = str(record.get("metadata", {}).get("description", "")).lower()
        for phrase in ("spherical average of density profile", "spherical average of isotopic profile", "stable elements after about 2 gyr"):
            if phrase not in description:
                failures.append(f"record_description_missing:{phrase}")
        remote_entries = record.get("files", {}).get("entries", {})
        for archive_name in ARCHIVES.values():
            remote = remote_entries.get(archive_name, {}) if isinstance(remote_entries, dict) else {}
            local = next((entry for entry in entries if entry.get("path") == archive_name), {})
            if remote.get("size") != local.get("bytes") or str(remote.get("checksum", "")).lower() != f"md5:{local.get('md5')}".lower():
                failures.append(f"record_file_metadata_mismatch:{archive_name}")

    source_models: dict[str, dict[str, str]] = {}
    for kind, archive_name in ARCHIVES.items():
        archive_path = root / archive_name
        if not archive_path.is_file():
            continue
        try:
            with ZipFile(archive_path) as zfile:
                model_files = _model_file_names(kind, zfile.namelist())
                expected_models = set(MODEL_IDS)
                if set(model_files) != expected_models:
                    failures.append(f"{kind}_model_set_mismatch")
                for model in MODEL_IDS:
                    source_models.setdefault(model, {})[kind] = model_files.get(model, "")
        except (BadZipFile, OSError):
            failures.append(f"{kind}_archive_invalid")

    epochs = _profile_epochs(record)
    for model in MODEL_IDS:
        if model not in epochs:
            failures.append(f"record_profile_epoch_missing:{model}")
    model_reports: dict[str, Any] = {}
    physical_warnings: list[dict[str, Any]] = []
    for model in MODEL_IDS:
        report: dict[str, Any] = {}
        for kind, archive_name in ARCHIVES.items():
            name = source_models.get(model, {}).get(kind)
            if not name:
                continue
            try:
                with ZipFile(root / archive_name) as zfile:
                    if kind == "abundances":
                        detail, detail_failures = _audit_abundance_file(zfile, name)
                    elif kind == "isotopes":
                        detail, detail_failures = _audit_isotope_file(zfile, name)
                    else:
                        detail, detail_failures = _audit_density_file(zfile, name, epochs.get(model, 100.0))
                report[kind] = detail
                failures.extend(detail_failures)
            except (BadZipFile, KeyError, OSError, UnicodeDecodeError) as exc:
                failures.append(f"{kind}:{model}:parse_failed:{type(exc).__name__}")
        if report.get("density", {}).get("row_count") != report.get("isotopes", {}).get("row_count"):
            failures.append(f"{model}:density_isotope_row_count_mismatch")
        if "abundances" in report and "density" in report:
            integrated_mass = report["abundances"]["stable_element_mass_msun"]
            profile_mass = report["density"]["profile_mass_estimate_msun"]
            relative_difference = abs(profile_mass - integrated_mass) / max(integrated_mass, 1.0e-30)
            if relative_difference >= GROSS_PROFILE_ANOMALY_THRESHOLD:
                review_classification = "source_data_anomaly_requires_quarantine"
            elif relative_difference > PROFILE_REVIEW_THRESHOLD:
                review_classification = "profile_mass_warning_requires_resolution"
            else:
                review_classification = "profile_consistent_review_candidate"
            report["profile_mass_vs_integrated_abundance"] = {
                "integrated_stable_element_mass_msun": integrated_mass,
                "profile_mass_estimate_msun": profile_mass,
                "relative_difference": relative_difference,
                "status": "within_5_percent" if relative_difference <= PROFILE_REVIEW_THRESHOLD else "warning_profile_mass_mismatch",
                "review_classification": review_classification,
            }
            if relative_difference > PROFILE_REVIEW_THRESHOLD:
                physical_warnings.append({
                    "model": model,
                    "type": "profile_mass_mismatch",
                    "relative_difference": relative_difference,
                    "severity": (
                        "source_data_anomaly"
                        if review_classification == "source_data_anomaly_requires_quarantine"
                        else "profile_mass_warning"
                    ),
                    "requires_quarantine": review_classification == "source_data_anomaly_requires_quarantine",
                    "interpretation": (
                        "source-data anomaly requires quarantine before any population mixture"
                        if review_classification == "source_data_anomaly_requires_quarantine"
                        else "spherical-average density integration is not ready to define returned mass without a source-specific correction or authoritative scalar"
                    ),
                })
        model_reports[model] = report

    return {
        "schema": "snrt-fp2-snia-hesma-source-audit",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_source_format_passed" if not failures else "blocked_source_format_integrity",
        "source_integrity_passed": not failures,
        "canonical_conversion_allowed": False,
        "runtime_activation_allowed": False,
        "root": project_relative(root),
        "manifest": project_relative(manifest_path),
        "record_id": record.get("id"),
        "record_access": record.get("access"),
        "record_rights": record.get("metadata", {}).get("rights", []),
        "model_count": len(MODEL_IDS),
        "model_ids": list(MODEL_IDS),
        "project_elements": list(PROJECT_ELEMENTS),
        "model_reports": model_reports,
        "physical_warnings": physical_warnings,
        "physical_review_status": (
            "review_only_with_physical_warnings" if physical_warnings
            else "review_only_physical_diagnostics_clean"
        ),
        "physical_review_policy": {
            "profile_warning_threshold_relative_difference": PROFILE_REVIEW_THRESHOLD,
            "gross_source_anomaly_threshold_relative_difference": GROSS_PROFILE_ANOMALY_THRESHOLD,
            "thresholds_are_diagnostic_only": True,
            "quarantine_rule": "models at or above the gross threshold cannot enter a production mixture without source resolution",
        },
        "data_semantics": {
            "integrated_abundance_decay_horizon": "stable elements/isotopes after about 2 Gyr; radioactive isotopes at end of simulation",
            "profile_epoch": "about 100 s after explosion",
            "hydro_profile": "spherical-average density and isotope mass fractions",
            "hydrogen_profile_mapping": "free-proton column p is explicitly mapped to H; no implicit missing-element zero",
            "event_energy": "not supplied as a scalar event field; derive only under an approved profile convention",
            "event_momentum": "signed vector is not determined by spherical-average files",
        },
        "conversion_blockers": [
            "select a concrete HESMA model or approved population mixture",
            "pin decay horizon and isotope-to-11-element conversion, including the explicit p-to-H mapping",
            "derive and validate returned mass and kinetic energy from a documented profile convention or an authoritative scalar source",
            "obtain a signed momentum convention/source; spherical average alone is insufficient",
            "bind population weighting, approval id, and deterministic converted-asset checksum",
        ],
        "files": file_reports,
        "package_sha256": observed_package,
        "manifest_sha256": _fingerprint(manifest_path)[1] if _fingerprint(manifest_path) else None,
        "audit_code_sha256": _fingerprint(TOOL_PATH)[1] if _fingerprint(TOOL_PATH) else None,
        "audit_failures": failures,
        "interpretation": "HESMA supplies a public, checksum-bound composition/profile candidate with all 11 integrated project elements, but this report does not approve event conversion or runtime activation.",
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    report = audit_source(args.root, args.manifest)
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["audit_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
