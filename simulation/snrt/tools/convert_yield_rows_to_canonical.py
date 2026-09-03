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
import os
from pathlib import Path
from typing import Any

from audit_fp1_source_node_contract import (
    SourceNodeContractError,
    audit_source_node_contract,
)
from audit_fp1_physical_package_admission import (
    PhysicalPackageAdmissionError,
    audit_physical_package_admission,
)
from fp1_source_node_projection import (
    SourceNodeProjectionError,
    validate_canonical_row_against_source_node,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_SOURCE_NODE_CONTRACT = SNRT_ROOT / "config" / "fp1_source_node_contract_v1.json"
DEFAULT_PHYSICAL_PACKAGE_CONTRACT = (
    SNRT_ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
)
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
    "approval_id",
    "units",
    "IMF",
    "population_model",
    "channel_boundaries",
    "metallicity_definition",
    "solar_abundance_set",
    "remnant_model",
    "untracked_ejecta_policy",
    "axis_reduction_policy",
    "energy_semantics",
    "momentum_deposition_contract",
)
ALLOWED_ENERGY_SEMANTICS = {
    "cumulative_physical_erg_per_initial_star",
    "cumulative_injected_erg_per_initial_star",
}
ALLOWED_MOMENTUM_DEPOSITION_CONTRACTS = {
    "source_frame_vector_only_no_scalar_radial_deposition",
    "approved_fp1_terminal_deposition_contract",
}
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
    for hash_field in ("source_sha256",):
        source_hash = source[hash_field]
        if (
            not isinstance(source_hash, str)
            or len(source_hash) != 64
            or any(character not in "0123456789abcdefABCDEF" for character in source_hash)
        ):
            raise ConversionError(f"{hash_field} must be a 64-character hexadecimal hash")
    axis_reduction = source["axis_reduction_policy"]
    if not isinstance(axis_reduction, dict) or axis_reduction.get("mode") not in {
        "none",
        "explicit_frozen_axis",
        "approved_population_marginalization",
    }:
        raise ConversionError("axis_reduction_policy must declare an approved mode")
    if source["energy_semantics"] not in ALLOWED_ENERGY_SEMANTICS:
        raise ConversionError("energy_semantics is not an allowed canonical convention")
    if source["momentum_deposition_contract"] not in ALLOWED_MOMENTUM_DEPOSITION_CONTRACTS:
        raise ConversionError("momentum_deposition_contract is not an allowed convention")
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
        source_node_id = raw.get("source_node_id")
        if not isinstance(source_node_id, str) or not source_node_id:
            raise ConversionError(f"row {row_number} requires a non-empty source_node_id")
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
                "source_node_id": source_node_id,
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


def convert(
    input_path: Path,
    output_path: Path,
    sidecar_path: Path,
    node_mapping_path: Path | None = None,
    node_contract_path: Path = DEFAULT_SOURCE_NODE_CONTRACT,
) -> dict[str, Any]:
    document = _read_input(input_path)
    source = _source_metadata(document["source"])
    rows = _normalize_rows(document["rows"])
    table_text = _format_table(rows)
    if node_mapping_path is None:
        raise ConversionError("a source-node mapping output path is required")
    node_mapping_path = Path(node_mapping_path)
    node_contract_path = Path(node_contract_path).resolve()
    if node_contract_path != DEFAULT_SOURCE_NODE_CONTRACT.resolve():
        raise ConversionError("converter must bind the repository F-P1 source-node contract")
    try:
        node_contract = json.loads(node_contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConversionError(f"cannot read source-node contract {node_contract_path}: {exc}") from exc
    if (
        not isinstance(node_contract, dict)
        or node_contract.get("schema") != "snrt-fp1-source-node-contract"
        or node_contract.get("schema_version") != 1
    ):
        raise ConversionError("unsupported source-node contract")
    try:
        node_contract_audit = audit_source_node_contract(
            node_contract_path=node_contract_path
        )
    except SourceNodeContractError as exc:
        raise ConversionError(f"source-node contract audit failed: {exc}") from exc
    approval = node_contract.get("approval")
    if (
        node_contract_audit.get("status") != "approved_physical_nodes"
        or node_contract_audit.get("production_ready") is not True
        or not isinstance(approval, dict)
        or approval.get("canonical_conversion_allowed") is not True
        or approval.get("approval_id") != source["approval_id"]
    ):
        raise ConversionError(
            "canonical conversion requires the approved repository source-node contract "
            "and matching approval identity"
        )
    physical_nodes = node_contract.get("physical_nodes")
    assert isinstance(physical_nodes, list)
    nodes_by_id = {node["source_node_id"]: node for node in physical_nodes}
    for row_number, row in enumerate(rows, start=1):
        node = nodes_by_id.get(row["source_node_id"])
        if node is None:
            raise ConversionError(
                f"row {row_number} source_node_id is absent from the approved contract"
            )
        if not math.isclose(
            row["initial_mass_msun_per_star"],
            float(node["zams_mass_msun"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ) or not math.isclose(
            row["birth_metallicity_mass_fraction"],
            float(node["birth_metallicity_value"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise ConversionError(
                f"row {row_number} mass/metallicity disagrees with its source node"
            )
        try:
            validate_canonical_row_against_source_node(
                row, node, source["energy_semantics"]
            )
        except SourceNodeProjectionError as exc:
            raise ConversionError(
                f"row {row_number} violates its source-node projection: {exc}"
            ) from exc
    node_contract_hash = _sha256(node_contract_path)
    physical_package_path = DEFAULT_PHYSICAL_PACKAGE_CONTRACT.resolve()
    try:
        physical_package_contract = json.loads(
            physical_package_path.read_text(encoding="utf-8")
        )
        physical_package_audit = audit_physical_package_admission(
            contract_path=physical_package_path
        )
    except (
        OSError,
        json.JSONDecodeError,
        PhysicalPackageAdmissionError,
    ) as exc:
        raise ConversionError(f"physical-package admission audit failed: {exc}") from exc
    package_node_evidence = physical_package_contract.get("evidence_artifacts", {}).get(
        "source_node_contract", {}
    )
    if (
        package_node_evidence.get("path") != "config/fp1_source_node_contract_v1.json"
        or package_node_evidence.get("sha256") != node_contract_hash
    ):
        raise ConversionError(
            "admitted physical package does not bind the converter source-node contract"
        )
    selection = physical_package_contract.get("selection")
    if (
        physical_package_audit.get("status") != "admitted_physical_package"
        or physical_package_audit.get("production_ready") is not True
        or physical_package_audit.get("canonical_conversion_allowed") is not True
        or not isinstance(selection, dict)
        or selection.get("approval_id") != source["approval_id"]
        or selection.get("selected_package_sha256") != source["source_sha256"]
    ):
        raise ConversionError(
            "canonical conversion requires an admitted F-P1 physical package "
            "with matching approval identity"
        )
    physical_package_hash = _sha256(physical_package_path)
    asset_hash = hashlib.sha256(table_text.encode("utf-8")).hexdigest()
    mapping = {
        "schema": "snrt-fp1-source-node-row-mapping",
        "schema_version": 1,
        "source_node_contract_sha256": node_contract_hash,
        "source_node_contract_approval_id": approval["approval_id"],
        "physical_package_approval_id": selection["approval_id"],
        "canonical_asset_sha256": asset_hash,
        "canonical_row_count": len(rows),
        "rows": [
            {
                "canonical_coordinate": [
                    row["channel"],
                    row["initial_mass_msun_per_star"],
                    row["birth_metallicity_mass_fraction"],
                    row["age_yr"],
                ],
                "source_node_id": row["source_node_id"],
            }
            for row in rows
        ],
    }
    mapping_text = json.dumps(mapping, indent=2, sort_keys=True) + "\n"
    node_mapping_hash = hashlib.sha256(mapping_text.encode("utf-8")).hexdigest()
    if selection.get("source_node_mapping_sha256") != node_mapping_hash:
        raise ConversionError(
            "canonical source-node mapping disagrees with the admitted package selection"
        )
    if output_path.exists() or sidecar_path.exists() or node_mapping_path.exists():
        raise ConversionError("refusing to overwrite an existing table, sidecar, or node mapping")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    node_mapping_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(table_text, encoding="utf-8")
    node_mapping_path.write_text(mapping_text, encoding="utf-8")
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
        "source_node_contract_path": os.path.relpath(
            node_contract_path, start=sidecar_path.parent.resolve()
        ),
        "source_node_contract_sha256": node_contract_hash,
        "physical_package_contract_path": os.path.relpath(
            physical_package_path, start=sidecar_path.parent.resolve()
        ),
        "physical_package_contract_sha256": physical_package_hash,
        "source_node_mapping_path": os.path.relpath(
            node_mapping_path.resolve(), start=sidecar_path.parent.resolve()
        ),
        "source_node_mapping_sha256": node_mapping_hash,
        "axis_reduction_policy": source["axis_reduction_policy"],
        "energy_semantics": source["energy_semantics"],
        "momentum_deposition_contract": source["momentum_deposition_contract"],
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
    parser.add_argument(
        "--node-mapping", type=Path, required=True,
        help="write canonical-row to source-node mapping JSON",
    )
    parser.add_argument(
        "--node-contract", type=Path, default=DEFAULT_SOURCE_NODE_CONTRACT,
        help=f"source-node contract JSON (default: {DEFAULT_SOURCE_NODE_CONTRACT})",
    )
    args = parser.parse_args()
    try:
        metadata = convert(
            args.input,
            args.output,
            args.sidecar,
            args.node_mapping,
            args.node_contract,
        )
    except (ConversionError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": "converted", "rows": metadata["row_count"], "asset_sha256": metadata["asset_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
