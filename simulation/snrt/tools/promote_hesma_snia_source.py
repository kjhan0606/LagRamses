#!/usr/bin/env python3
"""Promote one reviewed HESMA model into an explicit SNIa event source.

The promotion is intentionally narrow.  It accepts only the already audited
``n100`` review extraction and records the physical conventions that are not
contained in a spherical HESMA profile: the integrated stable-element mass is
the returned mass, the event has no terminal remnant, the documented profile
kinetic-energy estimate is used as the event energy, and the unresolved
isotropic source has zero net vector momentum.  Net yields are explicitly
marked unavailable and are zero only because they are diagnostic, never a gas
mass source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[2]
ELEMENT_ORDER = ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"]
EXPECTED_MODEL = "n100"
EXPECTED_SOURCE = "hesma:yysd4-xap92"


class PromotionError(ValueError):
    """The reviewed source cannot be promoted under the explicit baseline."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PromotionError(f"cannot read normalized source {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PromotionError("normalized source must be a JSON object")
    return value


def promote(
    input_path: Path,
    output_path: Path,
    *,
    approval_id: str,
    source_commit_binding: str,
) -> dict[str, Any]:
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise PromotionError(f"refusing to overwrite existing output: {output_path}")
    if not approval_id.strip():
        raise PromotionError("approval_id is required")
    if len(source_commit_binding) != 40 or any(
        character not in "0123456789abcdefABCDEF" for character in source_commit_binding
    ):
        raise PromotionError("source_commit_binding must be a 40-character hexadecimal revision")

    normalized = _read(input_path)
    source = normalized.get("source", {})
    if normalized.get("schema") != "snrt-fp2-snia-hesma-source-normalized":
        raise PromotionError("normalized source schema mismatch")
    if normalized.get("status") != "review_only_source_normalized_physics_blocked":
        raise PromotionError("promotion input must remain the reviewed normalized source")
    if source.get("source_id") != EXPECTED_SOURCE or source.get("selected_model_id") != EXPECTED_MODEL:
        raise PromotionError("only the explicitly reviewed HESMA yysd4-xap92 n100 model is approved")
    if normalized.get("profiles", {}).get("physical_warnings"):
        raise PromotionError("selected HESMA model has unresolved physical warnings")

    composition = normalized.get("composition", {})
    elements = composition.get("element_order")
    ejecta = composition.get("project_stable_element_masses_msun")
    returned = composition.get("all_source_stable_element_mass_msun")
    if elements != ELEMENT_ORDER or not isinstance(ejecta, list) or len(ejecta) != len(ELEMENT_ORDER):
        raise PromotionError("review source does not provide the complete project element vector")
    if not isinstance(returned, (int, float)) or returned <= 0.0:
        raise PromotionError("review source returned mass is invalid")
    if sum(ejecta) > returned + 1.0e-12:
        raise PromotionError("tracked ejecta exceed returned mass")

    profile = normalized.get("profiles", {}).get("density", {})
    energy = profile.get("kinetic_energy_estimate_erg")
    if not isinstance(energy, (int, float)) or energy <= 0.0:
        raise PromotionError("review source profile energy is invalid")

    event_source = {
        "schema": "snrt-fp2-snia-hesma-event-source",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "approved_physical_baseline_runtime_gated",
        "source": {
            "source_id": f"{EXPECTED_SOURCE}:{EXPECTED_MODEL}",
            "candidate_id": source.get("candidate_id"),
            "record_id": source.get("record_id"),
            "model_id": EXPECTED_MODEL,
            "citation": source.get("citation"),
            "license": source.get("license"),
            "source_version": source.get("source_version"),
            "manifest_sha256": source.get("manifest_sha256"),
            "package_sha256": source.get("package_sha256"),
            "normalized_source_path": "simulation/snrt/data/fp2_snia_hesma_n100_review_normalized.json",
            "normalized_source_sha256": _sha256(input_path),
            "source_commit_binding": source_commit_binding.lower(),
        },
        "event": {
            "quantity_basis": "per_event",
            "decay_convention": "HESMA stable-element section after approximately 2 Gyr; adopted horizon is 2 Gyr",
            "decay_horizon_yr": 2.0e9,
            "isotope_to_project_element_policy": "use HESMA integrated stable-element section; free-proton column is H",
            "returned_mass_msun_per_event": returned,
            "terminal_remnant_msun_per_event": 0.0,
            "wd_debit_msun_per_event": returned,
            "energy_erg_per_event": energy,
            "energy_basis": "inner_zero_outer_half_bin spherical-average profile estimator",
            "momentum_policy": "isotropic_zero_vector",
            "momentum_g_cm_s_per_event": [0.0, 0.0, 0.0],
            "ejected_mass_msun_per_event": ejecta,
            "untracked_ejecta_msun_per_event": returned - sum(ejecta),
            "net_yield_msun_per_event": [0.0] * len(ELEMENT_ORDER),
            "net_yield_status": "not supplied by the integrated stable-element source; zero is diagnostic-only and never a gas-mass source",
            "population_weight": 1.0,
        },
        "element_order": ELEMENT_ORDER,
        "approval": {
            "approval_id": approval_id,
            "canonical_conversion_allowed": True,
            "runtime_activation_allowed": False,
            "production_ready": True,
            "publication_ready": False,
            "scope": "physical SNIa baseline source; runtime AMR caller and full net-yield/metallicity sensitivity remain separately gated",
        },
        "conversion": {
            "conversion_code_path": "simulation/snrt/tools/promote_hesma_snia_source.py",
            "conversion_code_sha256": _sha256(TOOL_PATH),
            "policy": "deterministic promotion from one explicit reviewed HESMA model; no profile correction or missing-element inference",
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(event_source, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return event_source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--approval-id", required=True)
    parser.add_argument("--source-commit-binding", required=True)
    args = parser.parse_args()
    try:
        result = promote(
            args.input,
            args.output,
            approval_id=args.approval_id,
            source_commit_binding=args.source_commit_binding,
        )
    except (PromotionError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": result["status"], "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
