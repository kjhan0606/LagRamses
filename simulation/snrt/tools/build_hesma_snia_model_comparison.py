#!/usr/bin/env python3
"""Build a review-only comparison matrix for all staged HESMA SNIa models.

The matrix is evidence for choosing a physical event source. It is not a
ranking, population mixture, event-yield table, or runtime input. In
particular, profile kinetic energies remain diagnostics because the HESMA
files are spherical averages and do not provide a signed momentum source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

from fp2_provenance import project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92"
DEFAULT_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"
DEFAULT_AUDIT = SNRT_ROOT / "data" / "fp2_snia_hesma_source_audit.json"
ELEMENT_ORDER = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")


class HesmaComparisonError(ValueError):
    """The HESMA comparison cannot be constructed from the staged source."""


def _sha256(path: Path) -> str:
    if not path.is_file():
        raise HesmaComparisonError(f"file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HesmaComparisonError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HesmaComparisonError(f"JSON object required: {path}")
    return value


def build_comparison(
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path = DEFAULT_MANIFEST,
    audit_path: Path = DEFAULT_AUDIT,
) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    audit_path = Path(audit_path).resolve()
    manifest = _read_json(manifest_path)
    audit = _read_json(audit_path)
    if audit.get("status") != "review_only_source_format_passed":
        raise HesmaComparisonError("HESMA source audit is not clean")
    if audit.get("record_id") != manifest.get("source_record_id"):
        raise HesmaComparisonError("HESMA audit and manifest record ids disagree")
    if audit.get("package_sha256") != manifest.get("package_sha256"):
        raise HesmaComparisonError("HESMA audit and manifest package hashes disagree")

    try:
        from adapt_hesma_snia_source import _member_name, _read_member, _stable_element_values
    except ImportError as exc:  # pragma: no cover - protects non-repository use
        raise HesmaComparisonError("cannot import HESMA source adapter helpers") from exc

    rows: list[dict[str, Any]] = []
    for model_id in audit.get("model_ids", []):
        model_report = audit["model_reports"][model_id]
        abundance_member = model_report["abundances"]["path"]
        abundance_bytes = _read_member(
            root / "ddt_2013_abundances.zip", abundance_member
        )
        all_stable = _stable_element_values(abundance_bytes, abundance_member)
        project_values = [all_stable[element.lower()] for element in ELEMENT_ORDER]
        total_mass = math.fsum(all_stable.values())
        tracked_mass = math.fsum(project_values)
        density = model_report["density"]
        closure = model_report["profile_mass_vs_integrated_abundance"]
        warnings = [
            warning for warning in audit.get("physical_warnings", [])
            if warning.get("model") == model_id
        ]
        rows.append(
            {
                "model_id": model_id,
                "element_order": list(ELEMENT_ORDER),
                "project_stable_element_masses_msun": project_values,
                "project_stable_element_mass_msun": tracked_mass,
                "all_source_stable_element_mass_msun": total_mass,
                "non_project_stable_element_mass_msun": max(0.0, total_mass - tracked_mass),
                "profile_epoch_s": density["profile_epoch_s"],
                "profile_row_count": density["row_count"],
                "profile_mass_estimate_msun": density["profile_mass_estimate_msun"],
                "profile_kinetic_energy_estimate_erg": density[
                    "profile_kinetic_energy_estimate_erg"
                ],
                "profile_mass_closure_status": closure["status"],
                "profile_mass_relative_difference": closure["relative_difference"],
                "profile_review_classification": closure.get("review_classification"),
                "physical_warning_count": len(warnings),
                "physical_warnings": warnings,
                "selection_status": "not_selected",
            }
        )

    if len(rows) != audit.get("model_count"):
        raise HesmaComparisonError("HESMA model count differs between audit and comparison")
    return {
        "schema": "snrt-fp2-snia-hesma-model-comparison",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_no_model_selected",
        "source": {
            "candidate_id": "hesma_model_archive_snia_profiles",
            "record_id": manifest["source_record_id"],
            "manifest_path": project_relative(manifest_path),
            "package_sha256": manifest["package_sha256"],
            "manifest_sha256": _sha256(manifest_path),
            "source_audit_path": project_relative(audit_path),
            "source_audit_sha256": _sha256(audit_path),
        },
        "comparison_basis": {
            "composition": "integrated stable-element section in solar masses after about 2 Gyr",
            "profile": "spherical-average density/isotope profiles at source-reported epoch",
            "profile_mass_estimator": "inherited review-only piecewise constant velocity-shell estimator",
            "energy_semantics": "profile kinetic energy estimate only; not an approved event scalar",
            "momentum_semantics": "signed momentum unavailable from spherical averages",
            "ranking_policy": "no automatic physical ranking or population weighting",
        },
        "model_count": len(rows),
        "models": rows,
        "selection": {
            "selected_model_id": None,
            "population_mixture": None,
            "approval_id": None,
            "selection_status": "pending_physics_review",
        },
        "admission": {
            "canonical_conversion_allowed": False,
            "runtime_activation_allowed": False,
            "canonical_rows_emitted": 0,
        },
        "blockers": [
            "choose an explicit model or approved population mixture",
            "resolve profile-versus-integrated-mass warnings without silent correction",
            "approve decay and isotope-to-project-element policy",
            "approve returned mass, remnant, event energy, and signed momentum semantics",
            "bind DTD normalization, population weighting, and approval id",
        ],
        "tool_sha256": _sha256(TOOL_PATH),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = build_comparison(root=args.root, manifest_path=args.manifest, audit_path=args.audit)
        payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(payload, encoding="utf-8")
    except (HesmaComparisonError, OSError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
