#!/usr/bin/env python3
"""Reject stale FS2010 interpolation, wiring, or table provenance."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "furlanetto_stoever_validation.json"
TABLE_DIRECTORY = ROOT / "data" / "furlanetto_stoever_2010"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_finite(value: object) -> None:
    if isinstance(value, dict):
        for child in value.values():
            assert_finite(child)
    elif isinstance(value, list):
        for child in value:
            assert_finite(child)
    elif isinstance(value, float):
        assert math.isfinite(value)


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert_finite(payload)
    assert payload["schema"] == "snrt_furlanetto_stoever_2010_validation_v1"
    assert payload["passed"] is True
    assert payload["table"]["energy_count"] == 258
    assert payload["table"]["energy_minimum_ev"] == 10.0
    assert payload["table"]["energy_maximum_ev"] == 9937.21
    assert len(payload["table"]["ionized_fraction_grid"]) == 14
    assert payload["continuity"][
        "maximum_absolute_channel_delta_99p9_to_100p1_ev"
    ] < payload["continuity"]["acceptance_threshold"]
    assert payload["table_floor_continuity"][
        "maximum_absolute_channel_delta_9p999_to_10p001_ev"
    ] < payload["table_floor_continuity"]["acceptance_threshold"]
    assert payload["pinned_21cmfast_reference"][
        "production_maximum_absolute_error"
    ] < 2.0e-15
    assert payload["independent_host_interpolation_maximum_absolute_error"] < 2.0e-12
    assert payload["energy_closure_maximum_absolute_error"] < 2.0e-15
    assert payload["conservative_solver"][
        "photoelectron_energy_ledger_relative_error"
    ] < 2.0e-14
    assert payload["multiphysics_solver"][
        "photoelectron_energy_ledger_relative_error"
    ] < 2.0e-14
    assert payload["multiphysics_solver"][
        "zero_helium_secondary_channels_are_zero"
    ] is True

    provenance = payload["provenance"]
    for key, relative in (
        ("test_sha256", "tests/secondary_furlanetto_stoever.py"),
        ("secondary_sha256", "snrt_core/secondary.py"),
        ("multiphysics_sha256", "snrt_core/multiphysics.py"),
        ("conservative_primordial_sha256", "snrt_core/conservative_primordial.py"),
        ("implicit_sha256", "snrt_core/implicit.py"),
    ):
        assert provenance[key] == sha256(ROOT / relative)
    assert provenance["table_manifest_sha256"] == sha256(
        TABLE_DIRECTORY / "TABLE_MANIFEST.json"
    )
    assert provenance["table_sha256"] == {
        path.name: sha256(path) for path in sorted(TABLE_DIRECTORY.glob("*.dat"))
    }
    print(
        "SECONDARY_FS2010_ARTIFACT_OK "
        "continuity_delta="
        f"{payload['continuity']['maximum_absolute_channel_delta_99p9_to_100p1_ev']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
