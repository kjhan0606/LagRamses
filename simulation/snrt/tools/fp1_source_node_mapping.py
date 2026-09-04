#!/usr/bin/env python3
"""Canonical, code-owned serialization for an admitted source-node mapping."""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any


MAPPING_SCHEMA = "snrt-fp1-source-node-row-mapping"
MAPPING_SCHEMA_VERSION = 1
MAPPING_REQUIRED_KEYS = {
    "schema",
    "schema_version",
    "source_node_contract_sha256",
    "source_node_contract_approval_id",
    "physical_package_approval_id",
    "physical_package_sha256",
    "canonical_asset_sha256",
    "canonical_row_count",
    "rows",
}
MAPPING_ROW_KEYS = {"canonical_coordinate", "source_node_id"}


class SourceNodeMappingError(ValueError):
    """The admitted source-node mapping is malformed or non-canonical."""


def _sha256(value: Any, field: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in value)
    ):
        raise SourceNodeMappingError(f"{field} must be a 64-character hexadecimal hash")
    return value.lower()


def _finite_number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SourceNodeMappingError(f"{field} must be a finite JSON number")
    number = float(value)
    if not math.isfinite(number):
        raise SourceNodeMappingError(f"{field} must be a finite JSON number")
    # JSON distinguishes -0.0 from 0.0 at the byte level.  They are the same
    # physical coordinate, so the code-owned representation uses positive zero.
    return 0.0 if number == 0.0 else number


def _coordinate(value: Any, row_number: int) -> list[int | float]:
    if not isinstance(value, list) or len(value) != 4:
        raise SourceNodeMappingError(
            f"mapping row {row_number} canonical_coordinate must have four values"
        )
    channel_number = _finite_number(value[0], f"mapping row {row_number} channel")
    if not channel_number.is_integer() or not 1 <= channel_number <= 5:
        raise SourceNodeMappingError(
            f"mapping row {row_number} channel must be an integer in 1..5"
        )
    return [
        int(channel_number),
        _finite_number(value[1], f"mapping row {row_number} mass"),
        _finite_number(value[2], f"mapping row {row_number} metallicity"),
        _finite_number(value[3], f"mapping row {row_number} age"),
    ]


def normalize_mapping_document(mapping: Any) -> dict[str, Any]:
    """Validate and normalize one admitted mapping before hashing it.

    Numeric coordinates are converted to a single representation, rows are
    sorted by canonical coordinate and source-node id, and duplicate canonical
    coordinates are rejected.  The result is safe for ``allow_nan=False`` JSON
    serialization and contains no negative zero.
    """

    if not isinstance(mapping, dict) or set(mapping) != MAPPING_REQUIRED_KEYS:
        raise SourceNodeMappingError("source-node mapping field set is not exact")
    if mapping.get("schema") != MAPPING_SCHEMA:
        raise SourceNodeMappingError("source-node mapping schema is unsupported")
    if mapping.get("schema_version") != MAPPING_SCHEMA_VERSION:
        raise SourceNodeMappingError("source-node mapping schema version is unsupported")
    normalized: dict[str, Any] = {
        "schema": MAPPING_SCHEMA,
        "schema_version": MAPPING_SCHEMA_VERSION,
        "source_node_contract_sha256": _sha256(
            mapping.get("source_node_contract_sha256"),
            "source_node_contract_sha256",
        ),
        "source_node_contract_approval_id": mapping.get(
            "source_node_contract_approval_id"
        ),
        "physical_package_approval_id": mapping.get("physical_package_approval_id"),
        "physical_package_sha256": _sha256(
            mapping.get("physical_package_sha256"), "physical_package_sha256"
        ),
        "canonical_asset_sha256": _sha256(
            mapping.get("canonical_asset_sha256"), "canonical_asset_sha256"
        ),
    }
    for field in (
        "source_node_contract_approval_id",
        "physical_package_approval_id",
    ):
        if not isinstance(normalized[field], str) or not normalized[field]:
            raise SourceNodeMappingError(f"{field} must be a non-empty string")

    rows = mapping.get("rows")
    row_count = mapping.get("canonical_row_count")
    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count <= 0:
        raise SourceNodeMappingError("canonical_row_count must be a positive integer")
    if not isinstance(rows, list) or not rows or row_count != len(rows):
        raise SourceNodeMappingError("canonical_row_count disagrees with mapping rows")

    normalized_rows: list[dict[str, Any]] = []
    coordinate_keys: set[tuple[int, float, float, float]] = set()
    for row_number, row in enumerate(rows, start=1):
        if not isinstance(row, dict) or set(row) != MAPPING_ROW_KEYS:
            raise SourceNodeMappingError(
                f"source-node mapping row {row_number} field set is not exact"
            )
        source_node_id = row.get("source_node_id")
        if not isinstance(source_node_id, str) or not source_node_id:
            raise SourceNodeMappingError(
                f"source-node mapping row {row_number} source_node_id is invalid"
            )
        coordinate = _coordinate(row.get("canonical_coordinate"), row_number)
        key = tuple(coordinate)  # type: ignore[arg-type]
        if key in coordinate_keys:
            raise SourceNodeMappingError(
                f"duplicate canonical coordinate in source-node mapping: {key}"
            )
        coordinate_keys.add(key)
        normalized_rows.append(
            {
                "canonical_coordinate": coordinate,
                "source_node_id": source_node_id,
            }
        )

    normalized_rows.sort(
        key=lambda row: (
            *row["canonical_coordinate"],
            row["source_node_id"],
        )
    )
    normalized["canonical_row_count"] = row_count
    normalized["rows"] = normalized_rows
    return normalized


def canonical_mapping_text(mapping: Any) -> str:
    """Return the only byte representation admitted for a mapping document."""

    normalized = normalize_mapping_document(mapping)
    return json.dumps(
        normalized,
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def mapping_sha256(mapping: Any) -> str:
    return hashlib.sha256(canonical_mapping_text(mapping).encode("utf-8")).hexdigest()
