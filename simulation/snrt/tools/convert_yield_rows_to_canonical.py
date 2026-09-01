#!/usr/bin/env python3
"""Convert explicitly normalized yield rows to the canonical Phase-0 table.

This is a deterministic normalization step, not a source-paper adapter.  It
refuses rate tables because silently integrating a rate a second time is a
scientific error.  A source-specific converter must first produce the JSON
row contract below and record its own conversion method and source hash.

Input JSON::

    {
      "source": {
        "citation": "...",
        "source_version": "...",
        "source_sha256": "64 hex characters",
        "release_history_semantics": "cumulative",
        "license_status": "approved",
        "provenance_status": "approved",
        "approval_id": "...",
        "units": "...",
        "IMF": "...",
        "population_model": "...",
        "channel_boundaries": {},
        "metallicity_definition": "...",
        "solar_abundance_set": "...",
        "remnant_model": "...",
        "untracked_ejecta_policy":
          "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
      },
      "rows": [{
        "channel": 1,
        "initial_mass_msun_per_star": 1.0,
        "birth_metallicity_mass_fraction": 0.001,
        "age_yr": 0.0,
        "returned_mass_msun_per_star": 0.0,
        "remnant_mass_msun_per_star": 0.0,
        "energy_erg_per_star": 0.0,
        "momentum_g_cm_s_per_star": [0.0, 0.0, 0.0],
        "ejecta_msun_per_star": [0.0] * 11,
        "net_yield_msun_per_star": [0.0] * 11
      }]
    }

The converter emits the 32-field ASCII table and a sidecar containing both
the canonical asset hash and the original source hash.  Missing approval or
source metadata stays missing and is subsequently rejected by the production
auditor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
ELEMENT_COUNT = 11
UNTRACKED_EJECTA_POLICY = (
    "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
)
REQUIRED_SOURCE_FIELDS = (
    "citation",
    "source_version",
    "source_sha256",
    "license_status",
    "provenance_status",
    "units",
    "IMF",
    "population_model",
    "channel_boundaries",
    "metallicity_definition",
    "solar_abundance_set",
    "remnant_model",
    "untracked_ejecta_policy",
)
ROW_FIELDS = (
    "channel",
    "initial_mass_msun_per_star",
    "birth_metallicity_mass_fraction",
    "age_yr",
    "returned_mass_msun_per_star",
    "remnant_mass_msun_per_star",
    "energy_erg_per_star",
    "momentum_g_cm_s_per_star",
    "ejecta_msun_per_star",
    "net_yield_msun_per_star",
)


class ConversionError(ValueError):
    """The normalized source rows violate the converter contract."""


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


def _finite_float(value: Any, field: str) -> float:
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
    missing = [field for field in REQUIRED_SOURCE_FIELDS if not source.get(field)]
    if missing:
        raise ConversionError(f"missing source metadata: {', '.join(missing)}")
    if source.get("release_history_semantics") != "cumulative":
        raise ConversionError(
            "release_history_semantics must be exactly 'cumulative'; "
            "rate tables require a separately reviewed source-specific converter"
        )
    if source.get("untracked_ejecta_policy") != UNTRACKED_EJECTA_POLICY:
        raise ConversionError(
            "untracked_ejecta_policy must be exactly "
            f"'{UNTRACKED_EJECTA_POLICY}'"
        )
    source_hash = source["source_sha256"]
    if (
        not isinstance(source_hash, str)
        or len(source_hash) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in source_hash)
    ):
        raise ConversionError("source_sha256 must be a 64-character hexadecimal hash")
    return dict(source)


def _normalize_rows(rows: list[Any]) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    seen: set[tuple[int, float, float, float]] = set()
    for row_number, raw in enumerate(rows, start=1):
        if not isinstance(raw, dict):
            raise ConversionError(f"row {row_number} must be an object")
        missing = [field for field in ROW_FIELDS if field not in raw]
        if missing:
            raise ConversionError(f"row {row_number} missing fields: {', '.join(missing)}")
        channel_float = _finite_float(raw["channel"], f"row {row_number} channel")
        if not channel_float.is_integer() or not 1 <= channel_float <= 5:
            raise ConversionError(f"row {row_number} channel must be an integer in 1..5")
        channel = int(channel_float)
        scalar_fields = ROW_FIELDS[1:7]
        values = {
            field: _finite_float(raw[field], f"row {row_number} {field}")
            for field in scalar_fields
        }
        momentum = [
            _finite_float(value, f"row {row_number} momentum")
            for value in raw["momentum_g_cm_s_per_star"]
        ]
        ejecta = [
            _finite_float(value, f"row {row_number} ejecta")
            for value in raw["ejecta_msun_per_star"]
        ]
        net = [
            _finite_float(value, f"row {row_number} net yield")
            for value in raw["net_yield_msun_per_star"]
        ]
        if len(momentum) != 3 or len(ejecta) != ELEMENT_COUNT or len(net) != ELEMENT_COUNT:
            raise ConversionError(
                f"row {row_number} requires momentum=3, ejecta=11, net_yield=11"
            )
        if (
            values["initial_mass_msun_per_star"] <= 0.0
            or values["birth_metallicity_mass_fraction"] < 0.0
            or values["age_yr"] < 0.0
            or values["returned_mass_msun_per_star"] < 0.0
            or values["remnant_mass_msun_per_star"] < 0.0
            or values["energy_erg_per_star"] < 0.0
            or any(value < 0.0 for value in ejecta)
        ):
            raise ConversionError(f"row {row_number} has a negative physical value")
        returned = values["returned_mass_msun_per_star"]
        tracked_ejecta = sum(ejecta)
        tolerance = 1.0e-12 + 1.0e-8 * max(abs(tracked_ejecta), abs(returned), 1.0)
        if tracked_ejecta > returned + tolerance:
            raise ConversionError(
                f"row {row_number} has tracked ejecta exceeding returned_mass"
            )
        if returned + values["remnant_mass_msun_per_star"] > values["initial_mass_msun_per_star"] + 1.0e-12:
            raise ConversionError(f"row {row_number} exceeds the initial-mass budget")
        coordinate = (
            channel,
            values["initial_mass_msun_per_star"],
            values["birth_metallicity_mass_fraction"],
            values["age_yr"],
        )
        if coordinate in seen:
            raise ConversionError(f"duplicate coordinate tuple at row {row_number}")
        seen.add(coordinate)
        normalized.append(
            {
                "channel": channel,
                **values,
                "momentum_g_cm_s_per_star": momentum,
                "ejecta_msun_per_star": ejecta,
                "net_yield_msun_per_star": net,
                "untracked_ejecta_msun_per_star": max(0.0, returned - tracked_ejecta),
            }
        )
    return sorted(
        normalized,
        key=lambda row: (
            row["channel"],
            row["initial_mass_msun_per_star"],
            row["birth_metallicity_mass_fraction"],
            row["age_yr"],
        ),
    )


def _format_table(rows: list[dict[str, Any]]) -> str:
    lines = [
        "# canonical_phase0_ascii; cumulative actual ejecta; age_yr on disk",
        "# channel initial_mass birth_metallicity age_yr returned remnant energy",
        "# momentum_x momentum_y momentum_z ejecta[H..Fe] net_yield[H..Fe]",
    ]
    for row in rows:
        values: list[float | int] = [
            row["channel"],
            row["initial_mass_msun_per_star"],
            row["birth_metallicity_mass_fraction"],
            row["age_yr"],
            row["returned_mass_msun_per_star"],
            row["remnant_mass_msun_per_star"],
            row["energy_erg_per_star"],
            *row["momentum_g_cm_s_per_star"],
            *row["ejecta_msun_per_star"],
            *row["net_yield_msun_per_star"],
        ]
        lines.append(" ".join(f"{value:.17g}" if isinstance(value, float) else str(value) for value in values))
    return "\n".join(lines) + "\n"


def convert(input_path: Path, output_path: Path, sidecar_path: Path) -> dict[str, Any]:
    document = _read_input(input_path)
    source = _source_metadata(document["source"])
    rows = _normalize_rows(document["rows"])
    table_text = _format_table(rows)
    if output_path.exists() or sidecar_path.exists():
        raise ConversionError("refusing to overwrite an existing table or sidecar")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(table_text, encoding="utf-8")
    asset_hash = _sha256(output_path)
    untracked = [row["untracked_ejecta_msun_per_star"] for row in rows]
    untracked_fractions = [
        residual / row["returned_mass_msun_per_star"]
        for row, residual in zip(rows, untracked)
        if row["returned_mass_msun_per_star"] > 0.0
    ]
    metadata = {
        "schema": "phase0_stellar_yield_asset_sidecar_v1",
        "approval_id": source.get("approval_id"),
        "citation": source["citation"],
        "source_version": source["source_version"],
        "source_sha256": source["source_sha256"],
        "license_status": source["license_status"],
        "provenance_status": source["provenance_status"],
        "units": source["units"],
        "IMF": source["IMF"],
        "population_model": source["population_model"],
        "channel_boundaries": source["channel_boundaries"],
        "metallicity_definition": source["metallicity_definition"],
        "solar_abundance_set": source["solar_abundance_set"],
        "remnant_model": source["remnant_model"],
        "untracked_ejecta_policy": source["untracked_ejecta_policy"],
        "release_history_semantics": "cumulative",
        "conversion_code_sha256": _sha256(TOOL_PATH),
        "asset_sha256": asset_hash,
        "sha256": asset_hash,
        "asset_bytes": output_path.stat().st_size,
        "row_count": len(rows),
        "canonical_field_count": 32,
        "rows_with_untracked_ejecta": sum(value > 1.0e-12 for value in untracked),
        "maximum_untracked_ejecta_msun_per_star": max(untracked, default=0.0),
        "maximum_untracked_ejecta_fraction_of_returned_mass": max(
            untracked_fractions, default=0.0
        ),
        "conversion_policy": "deterministic normalization only; no rate integration",
    }
    sidecar_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="normalized source-row JSON")
    parser.add_argument("--output", type=Path, required=True, help="canonical ASCII table")
    parser.add_argument("--sidecar", type=Path, required=True, help="canonical provenance sidecar")
    args = parser.parse_args()
    try:
        metadata = convert(args.input, args.output, args.sidecar)
    except (ConversionError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": "converted", "rows": metadata["row_count"], "asset_sha256": metadata["asset_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
