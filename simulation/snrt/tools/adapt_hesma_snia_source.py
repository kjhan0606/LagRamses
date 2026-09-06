#!/usr/bin/env python3
"""Extract one HESMA SNIa model into a review-only source document.

This adapter is deliberately upstream of ``convert_snia_event_yields.py``.
It reads the integrated stable-element section and records the profile
provenance for one explicitly named HESMA model. It does not infer a returned
mass, terminal remnant, event energy, signed momentum, decay horizon, or
population weight, and it never emits a canonical/runtime asset.

The explicit model argument is intentional: the HESMA record is a model
archive, not a unique physical event prescription. A successful extraction
therefore means only that the selected source model is normalized for review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile


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
ELEMENT_ORDER = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
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


class HesmaAdapterError(ValueError):
    """The review-only HESMA extraction cannot be admitted."""


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HesmaAdapterError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HesmaAdapterError(f"JSON object required: {path}")
    return value


def _canonical_model(model_id: str) -> str:
    if not isinstance(model_id, str) or not model_id.strip():
        raise HesmaAdapterError("an explicit non-empty HESMA model id is required")
    candidate = model_id.strip().lower()
    if candidate not in MODEL_IDS:
        raise HesmaAdapterError(
            f"unknown HESMA model {model_id!r}; choose one of {', '.join(MODEL_IDS)}"
        )
    return candidate


def _manifest_entry(manifest: dict[str, Any], path: str) -> dict[str, Any]:
    entries = manifest.get("files")
    if not isinstance(entries, list):
        raise HesmaAdapterError("HESMA manifest files list is missing")
    for entry in entries:
        if isinstance(entry, dict) and entry.get("path") == path:
            return entry
    raise HesmaAdapterError(f"HESMA manifest entry is missing: {path}")


def _member_name(archive_path: Path, model_id: str, kind: str) -> str:
    suffix = f"_{kind}.dat"
    expected_basename = f"ddt_2013_{model_id}{suffix}"
    try:
        with ZipFile(archive_path) as zfile:
            matches = [name for name in zfile.namelist() if Path(name).name == expected_basename]
    except (BadZipFile, OSError) as exc:
        raise HesmaAdapterError(f"cannot open HESMA {kind} archive: {archive_path}") from exc
    if len(matches) != 1:
        raise HesmaAdapterError(
            f"expected exactly one {kind} member for {model_id}, found {len(matches)}"
        )
    return matches[0]


def _stable_element_values(data: bytes, member_name: str) -> dict[str, float]:
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise HesmaAdapterError(f"HESMA abundance member is not UTF-8: {member_name}") from exc
    section = ""
    values: dict[str, float] = {}
    for line in lines:
        lowered = line.strip().lower()
        if lowered.startswith("abundances of stable elements"):
            section = "stable_elements"
            continue
        if lowered.startswith("abundances of "):
            section = "other"
            continue
        if section != "stable_elements":
            continue
        match = ABUNDANCE_RE.match(line)
        if match is None:
            continue
        symbol, value_text = match.groups()
        value = float(value_text)
        if not math.isfinite(value) or value < 0.0:
            raise HesmaAdapterError(f"invalid stable element value {symbol} in {member_name}")
        if symbol in values:
            raise HesmaAdapterError(f"duplicate stable element {symbol} in {member_name}")
        values[symbol] = value
    if not values:
        raise HesmaAdapterError(f"stable element section is missing: {member_name}")
    missing = [element.lower() for element in ELEMENT_ORDER if element.lower() not in values]
    if missing:
        raise HesmaAdapterError(
            f"selected HESMA model lacks project stable elements: {', '.join(missing)}"
        )
    return values


def _read_member(archive_path: Path, member_name: str) -> bytes:
    try:
        with ZipFile(archive_path) as zfile:
            return zfile.read(member_name)
    except (BadZipFile, KeyError, OSError) as exc:
        raise HesmaAdapterError(f"cannot read HESMA member {member_name}") from exc


def _archive_fingerprints(root: Path, manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    fingerprints: dict[str, dict[str, Any]] = {}
    for archive_name in (*ARCHIVES.values(), "record.json"):
        path = root / archive_name
        entry = _manifest_entry(manifest, archive_name)
        if not path.is_file():
            raise HesmaAdapterError(f"HESMA source file is missing: {path}")
        observed_bytes = path.stat().st_size
        observed_sha256 = _sha256(path)
        if (
            observed_bytes != entry.get("bytes")
            or observed_sha256.lower() != str(entry.get("sha256", "")).lower()
        ):
            raise HesmaAdapterError(f"HESMA source fingerprint mismatch: {archive_name}")
        fingerprints[archive_name] = {
            "bytes": observed_bytes,
            "sha256": observed_sha256,
            "md5": entry.get("md5"),
        }
    return fingerprints


def adapt_source(
    model_id: str,
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path = DEFAULT_MANIFEST,
) -> dict[str, Any]:
    """Extract one explicitly selected HESMA model without runtime admission."""

    model = _canonical_model(model_id)
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    manifest = _read_json(manifest_path)
    fingerprints = _archive_fingerprints(root, manifest)

    # Reuse the strict package audit. Physical warnings are deliberately
    # retained in the report; they are not hidden by this adapter.
    try:
        from audit_fp2_snia_hesma import audit_source
    except ImportError as exc:  # pragma: no cover - protects non-repository use
        raise HesmaAdapterError("cannot import the HESMA source audit") from exc
    audit = audit_source(root, manifest_path)
    if audit["audit_failures"]:
        raise HesmaAdapterError(
            "HESMA source audit failed: " + ", ".join(audit["audit_failures"])
        )

    member_paths: dict[str, str] = {}
    member_data: dict[str, bytes] = {}
    for kind, archive_name in ARCHIVES.items():
        member = _member_name(root / archive_name, model, kind)
        member_paths[kind] = member
        member_data[kind] = _read_member(root / archive_name, member)

    stable_elements = _stable_element_values(member_data["abundances"], member_paths["abundances"])
    project_values = [stable_elements[element.lower()] for element in ELEMENT_ORDER]
    source_stable_mass = math.fsum(stable_elements.values())
    project_mass = math.fsum(project_values)
    model_report = audit["model_reports"][model]
    profile = model_report["density"]
    isotopes = model_report["isotopes"]
    abundance = model_report["abundances"]
    review_classification = model_report.get("profile_mass_vs_integrated_abundance", {}).get(
        "review_classification"
    )
    if review_classification == "source_data_anomaly_requires_quarantine":
        raise HesmaAdapterError(
            f"HESMA model {model} is quarantined as a source-data anomaly"
        )
    model_warnings = [
        warning for warning in audit["physical_warnings"]
        if warning.get("model") == model
    ]
    if model_warnings:
        raise HesmaAdapterError(
            f"HESMA model {model} has unresolved physical warnings; source resolution is required"
        )

    member_fingerprints = {
        kind: {
            "archive": archive_name,
            "archive_sha256": fingerprints[archive_name]["sha256"],
            "member": member_paths[kind],
            "member_bytes": len(member_data[kind]),
            "member_sha256": _sha256_bytes(member_data[kind]),
        }
        for kind, archive_name in ARCHIVES.items()
    }

    return {
        "schema": "snrt-fp2-snia-hesma-source-normalized",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_source_normalized_physics_blocked",
        "source": {
            "source_id": f"hesma:{audit['record_id']}",
            "candidate_id": "hesma_model_archive_snia_profiles",
            "record_id": audit["record_id"],
            "record_api_url": manifest["source"]["record_api_url"],
            "citation": manifest["source"]["citation"],
            "source_version": manifest["source"]["version"],
            "license": manifest["source"]["license"],
            "asset_root": manifest["asset_root"],
            "package_sha256": manifest["package_sha256"],
            "manifest_sha256": _sha256(manifest_path),
            "selected_model_id": model,
            "selection_policy": "explicit_model_argument_required; archive default forbidden",
        },
        "source_members": member_fingerprints,
        "decay_and_mapping": {
            "source_integrated_abundance_convention": audit["data_semantics"][
                "integrated_abundance_decay_horizon"
            ],
            "source_profile_convention": audit["data_semantics"]["profile_epoch"],
            "exact_decay_horizon_yr": None,
            "isotope_to_element_conversion": "not performed; source integrated stable-element section copied as review diagnostic",
            "hydrogen_mapping": audit["data_semantics"]["hydrogen_profile_mapping"],
            "decay_approval_status": "pending_physics_approval",
        },
        "composition": {
            "element_order": list(ELEMENT_ORDER),
            "project_stable_element_masses_msun": project_values,
            "project_stable_element_mass_msun": project_mass,
            "all_source_stable_element_masses_msun": dict(sorted(stable_elements.items())),
            "all_source_stable_element_mass_msun": source_stable_mass,
            "non_project_stable_element_mass_msun": max(0.0, source_stable_mass - project_mass),
            "source_mass_semantics": "integrated stable elements after about 2 Gyr; not an approved returned-mass scalar",
            "missing_project_elements": abundance["missing_project_elements"],
        },
        "profiles": {
            "profile_epoch_s": profile["profile_epoch_s"],
            "density": {
                "member": member_paths["density"],
                "row_count": profile["row_count"],
                "velocity_min_km_s": profile["velocity_min_km_s"],
                "velocity_max_km_s": profile["velocity_max_km_s"],
                "spherical_average": profile["profile_is_spherical_average"],
                "mass_estimate_msun": profile["profile_mass_estimate_msun"],
                "kinetic_energy_estimate_erg": profile["profile_kinetic_energy_estimate_erg"],
                "estimator": profile["profile_estimator"],
            },
            "isotopes": {
                "member": member_paths["isotopes"],
                "row_count": isotopes["row_count"],
                "isotope_count": isotopes["isotope_count"],
                "isotope_header_sha256": isotopes["isotope_header_sha256"],
                "spherical_average": True,
            },
            "density_isotope_row_count_match": profile["row_count"] == isotopes["row_count"],
            "profile_mass_vs_integrated_stable_mass": model_report[
                "profile_mass_vs_integrated_abundance"
            ],
            "physical_warnings": [
                warning for warning in audit["physical_warnings"] if warning["model"] == model
            ],
        },
        "event_contract": {
            "returned_mass_msun_per_event": None,
            "terminal_remnant_msun_per_event": None,
            "energy_erg_per_event": None,
            "momentum_g_cm_s_per_event": None,
            "population_weight": None,
            "status": "blocked_missing_approved_event_contract",
            "reasons": [
                "integrated stable-element total is not an authoritative returned-mass field",
                "profile kinetic energy is a review-only estimator and has model-dependent warnings",
                "spherical-average profiles do not determine a signed momentum vector",
                "terminal remnant ownership is not selected by this archive extraction",
                "DTD normalization and population weighting are separate unapproved inputs",
            ],
        },
        "admission": {
            "canonical_conversion_allowed": False,
            "runtime_activation_allowed": False,
            "converter_input_emitted": False,
            "required_next_approvals": [
                "explicit decay horizon and isotope-to-11-element policy",
                "authoritative returned mass and terminal remnant ownership",
                "event energy convention or authoritative scalar",
                "signed momentum convention/source",
                "DTD population weighting and named approval id",
            ],
        },
        "adapter": {
            "code_sha256": _sha256(TOOL_PATH),
            "source_audit_code_sha256": audit["audit_code_sha256"],
            "source_audit_manifest_sha256": audit["manifest_sha256"],
            "normalization_policy": "read one explicit model; preserve source values and provenance; infer no missing physics",
        },
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="explicit HESMA model id, for example n100")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json-out", type=Path, help="optional review-only JSON output")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = adapt_source(args.model, root=args.root, manifest_path=args.manifest)
        payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(payload, encoding="utf-8")
    except (HesmaAdapterError, OSError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
