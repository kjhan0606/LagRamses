#!/usr/bin/env python3
"""Convert normalized SNIa event yields into a separate, fail-closed asset.

This converter is intentionally not a per-star canonical-yield converter.  A
SNIa event is driven by the independently integrated DTD and has per-event
ejecta, energy, and momentum.  The input must therefore already be normalized
by a source-specific adapter; this step only validates, sorts, and records the
immutable provenance needed by the runtime admission gate.

Input JSON::

    {
      "source": {
        "source_id": "...",
        "citation": "...",
        "source_version": "...",
        "source_path": "/absolute/path/to/source",
        "source_sha256": "64 hex characters",
        "license_status": "approved",
        "provenance_status": "approved",
        "approval_id": "...",
        "decay_convention": "...",
        "decay_horizon_yr": 1.0e9,
        "metallicity_definition": "...",
        "population_model": "...",
        "model_selection": "...",
        "element_order": ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"]
      },
      "rows": [{
        "model_id": "...",
        "metallicity_mass_fraction": 0.001,
        "returned_mass_msun_per_event": 1.4,
        "remnant_mass_msun_per_event": 0.0,
        "energy_erg_per_event": 1.0e51,
        "momentum_g_cm_s_per_event": [0.0, 0.0, 0.0],
        "ejecta_msun_per_event": [0.0] * 11,
        "net_yield_msun_per_event": [0.0] * 11
      }]
    }

The output JSON and its sidecar are separate from the ordinary stellar yield
table.  No runtime activation is implied by a successful conversion; the
F-P2 contract must still name this asset and an approved DTD/population.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
ELEMENT_ORDER = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
ELEMENT_COUNT = len(ELEMENT_ORDER)
REQUIRED_SOURCE_FIELDS = (
    "source_id",
    "citation",
    "source_version",
    "source_path",
    "source_sha256",
    "license_status",
    "provenance_status",
    "approval_id",
    "decay_convention",
    "decay_horizon_yr",
    "metallicity_definition",
    "population_model",
    "model_selection",
    "element_order",
)
ROW_FIELDS = (
    "model_id",
    "metallicity_mass_fraction",
    "returned_mass_msun_per_event",
    "remnant_mass_msun_per_event",
    "energy_erg_per_event",
    "momentum_g_cm_s_per_event",
    "ejecta_msun_per_event",
    "net_yield_msun_per_event",
)


class ConversionError(ValueError):
    """The normalized SNIa event source violates the admission contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_input(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConversionError(f"cannot read input JSON {path}: {exc}") from exc
    if not isinstance(value, dict) or not isinstance(value.get("source"), dict):
        raise ConversionError("input must contain an object named source")
    if not isinstance(value.get("rows"), list) or not value["rows"]:
        raise ConversionError("input must contain a non-empty rows list")
    return value


def _finite(value: Any, field: str) -> float:
    if isinstance(value, bool):
        raise ConversionError(f"{field} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ConversionError(f"{field} must be a finite number") from exc
    if not math.isfinite(number):
        raise ConversionError(f"{field} must be finite")
    return number


def _source_metadata(source: dict[str, Any]) -> dict[str, Any]:
    missing = []
    for field in REQUIRED_SOURCE_FIELDS:
        value = source.get(field)
        if field not in source or value is None or (isinstance(value, str) and not value.strip()):
            missing.append(field)
    if missing:
        raise ConversionError(f"missing source metadata: {', '.join(missing)}")
    if source["license_status"] != "approved":
        raise ConversionError("license_status must be exactly 'approved'")
    if source["provenance_status"] != "approved":
        raise ConversionError("provenance_status must be exactly 'approved'")
    if list(source["element_order"]) != list(ELEMENT_ORDER):
        raise ConversionError("element_order does not match the Phase-0 11-element order")
    decay_horizon = _finite(source["decay_horizon_yr"], "decay_horizon_yr")
    if decay_horizon < 0.0:
        raise ConversionError("decay_horizon_yr must be non-negative")
    source_path = Path(source["source_path"])
    if not source_path.is_absolute() or not source_path.is_file():
        raise ConversionError("source_path must be an existing absolute file")
    source_hash = source["source_sha256"]
    if (
        not isinstance(source_hash, str)
        or len(source_hash) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in source_hash)
    ):
        raise ConversionError("source_sha256 must be a 64-character hexadecimal hash")
    observed_hash = _sha256(source_path)
    if observed_hash.lower() != source_hash.lower():
        raise ConversionError("source_sha256 does not match source_path")
    normalized = dict(source)
    normalized["decay_horizon_yr"] = decay_horizon
    normalized["source_sha256"] = observed_hash
    normalized["element_order"] = list(ELEMENT_ORDER)
    return normalized


def _normalize_rows(rows: list[Any]) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    seen: set[tuple[str, float]] = set()
    for row_number, raw in enumerate(rows, start=1):
        if not isinstance(raw, dict):
            raise ConversionError(f"row {row_number} must be an object")
        missing = [field for field in ROW_FIELDS if field not in raw]
        if missing:
            raise ConversionError(f"row {row_number} missing fields: {', '.join(missing)}")
        model_id = raw["model_id"]
        if not isinstance(model_id, str) or not model_id:
            raise ConversionError(f"row {row_number} model_id must be non-empty")
        metallicity = _finite(raw["metallicity_mass_fraction"], f"row {row_number} metallicity")
        returned = _finite(raw["returned_mass_msun_per_event"], f"row {row_number} returned mass")
        remnant = _finite(raw["remnant_mass_msun_per_event"], f"row {row_number} remnant mass")
        energy = _finite(raw["energy_erg_per_event"], f"row {row_number} energy")
        momentum = [_finite(value, f"row {row_number} momentum") for value in raw["momentum_g_cm_s_per_event"]]
        ejecta = [_finite(value, f"row {row_number} ejecta") for value in raw["ejecta_msun_per_event"]]
        net = [_finite(value, f"row {row_number} net yield") for value in raw["net_yield_msun_per_event"]]
        if len(momentum) != 3 or len(ejecta) != ELEMENT_COUNT or len(net) != ELEMENT_COUNT:
            raise ConversionError(f"row {row_number} requires momentum=3, ejecta=11, net_yield=11")
        if metallicity < 0.0 or returned < 0.0 or remnant < 0.0 or energy < 0.0:
            raise ConversionError(f"row {row_number} has a negative physical value")
        if remnant != 0.0:
            raise ConversionError("normal SNIa event source must have zero terminal remnant")
        if any(value < 0.0 for value in ejecta):
            raise ConversionError(f"row {row_number} has negative actual ejecta")
        tracked_ejecta = sum(ejecta)
        tolerance = 1.0e-12 + 1.0e-8 * max(abs(tracked_ejecta), abs(returned), 1.0)
        if tracked_ejecta > returned + tolerance:
            raise ConversionError(f"row {row_number} has tracked ejecta exceeding returned mass")
        coordinate = (model_id, metallicity)
        if coordinate in seen:
            raise ConversionError(f"duplicate model/metallicity coordinate at row {row_number}")
        seen.add(coordinate)
        normalized.append(
            {
                "model_id": model_id,
                "metallicity_mass_fraction": metallicity,
                "returned_mass_msun_per_event": returned,
                "remnant_mass_msun_per_event": 0.0,
                "energy_erg_per_event": energy,
                "momentum_g_cm_s_per_event": momentum,
                "ejecta_msun_per_event": ejecta,
                "net_yield_msun_per_event": net,
                "untracked_ejecta_msun_per_event": max(0.0, returned - tracked_ejecta),
            }
        )
    return sorted(normalized, key=lambda row: (row["model_id"], row["metallicity_mass_fraction"]))


def convert(input_path: Path, output_path: Path, sidecar_path: Path) -> dict[str, Any]:
    document = _read_input(input_path)
    source = _source_metadata(document["source"])
    rows = _normalize_rows(document["rows"])
    if output_path.exists() or sidecar_path.exists():
        raise ConversionError("refusing to overwrite an existing event asset or sidecar")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    event_asset = {
        "schema": "snrt-fp2-snia-event-yield-asset",
        "schema_version": 1,
        "status": "converted_awaiting_fp2_runtime_admission",
        "source": source,
        "element_order": list(ELEMENT_ORDER),
        "event_semantics": {
            "quantity_basis": "per_event",
            "actual_ejecta_is_gas_source": True,
            "net_yield_is_diagnostic": True,
            "terminal_remnant_msun_per_event": 0.0,
            "untracked_ejecta_policy": "returned_mass_minus_sum_tracked_ejecta",
        },
        "rows": rows,
        "conversion": {
            "input_json_sha256": _sha256(input_path),
            "conversion_code_sha256": _sha256(TOOL_PATH),
            "normalization_policy": "validate_and_sort_only; source-specific isotope conversion must be upstream",
        },
    }
    output_path.write_text(json.dumps(event_asset, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    asset_hash = _sha256(output_path)
    residuals = [row["untracked_ejecta_msun_per_event"] for row in rows]
    sidecar = {
        "schema": "snrt-fp2-snia-event-yield-sidecar-v1",
        "status": "converted_awaiting_fp2_runtime_admission",
        "source_id": source["source_id"],
        "citation": source["citation"],
        "source_version": source["source_version"],
        "source_path": source["source_path"],
        "source_sha256": source["source_sha256"],
        "license_status": source["license_status"],
        "provenance_status": source["provenance_status"],
        "approval_id": source["approval_id"],
        "decay_convention": source["decay_convention"],
        "decay_horizon_yr": source["decay_horizon_yr"],
        "conversion_code_sha256": _sha256(TOOL_PATH),
        "input_json_sha256": _sha256(input_path),
        "asset_sha256": asset_hash,
        "asset_bytes": output_path.stat().st_size,
        "row_count": len(rows),
        "element_order": list(ELEMENT_ORDER),
        "maximum_untracked_ejecta_msun_per_event": max(residuals, default=0.0),
        "conversion_policy": "event source only; never an ordinary IMF/per-star canonical yield table",
    }
    sidecar_path.write_text(json.dumps(sidecar, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return sidecar


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="source-normalized event-row JSON")
    parser.add_argument("--output", type=Path, required=True, help="event-yield JSON asset")
    parser.add_argument("--sidecar", type=Path, required=True, help="event-yield provenance sidecar")
    args = parser.parse_args()
    try:
        sidecar = convert(args.input, args.output, args.sidecar)
    except (ConversionError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": sidecar["status"], "rows": sidecar["row_count"], "asset_sha256": sidecar["asset_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
