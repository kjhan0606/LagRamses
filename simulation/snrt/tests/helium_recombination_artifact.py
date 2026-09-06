#!/usr/bin/env python3
"""Fail closed when the canonical helium-recombination artifact is stale."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "helium_case_b_recombination_validation.json"


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


def main() -> None:
    report = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert report["schema"] == "snrt_helium_case_b_recombination_validation_v2"
    assert report["passed"] is True
    assert_finite(report)
    assert report["temperature_k"] == [1.0e4, 2.0e4, 4.0e4, 1.0e5]

    one_zone = report["one_zone"]
    fixed_time = report["fixed_elapsed_time_one_zone"]
    assert one_zone["helium_ii_maximum_relative_error"] < one_zone["acceptance_threshold"]
    assert one_zone["helium_iii_maximum_relative_error"] < one_zone["acceptance_threshold"]
    assert fixed_time["helium_ii_maximum_relative_error"] < fixed_time["acceptance_threshold"]
    assert fixed_time["helium_iii_maximum_relative_error"] < fixed_time["acceptance_threshold"]
    assert np.ptp(fixed_time["helium_ii_final_fraction"]) > 0.1
    assert np.ptp(fixed_time["helium_iii_final_fraction"]) > 0.1

    provenance = report["provenance"]
    expected_hashes = {
        "test_sha256": ROOT / "tests" / "helium_recombination.py",
        "primordial_sha256": ROOT / "snrt_core" / "primordial.py",
        "implicit_sha256": ROOT / "snrt_core" / "implicit.py",
        "primordial_cooling_sha256": ROOT / "snrt_core" / "primordial_cooling.py",
        "b1_thermal_coupling_test_sha256": ROOT / "tests" / "b1_thermal_coupling.py",
    }
    for key, path in expected_hashes.items():
        assert provenance[key] == sha256(path), f"stale {key}: {path}"

    print(
        "HELIUM_RECOMBINATION_ARTIFACT_OK "
        f"heii_fixed_error={fixed_time['helium_ii_maximum_relative_error']:.6g} "
        f"heiii_fixed_error={fixed_time['helium_iii_maximum_relative_error']:.6g}"
    )


if __name__ == "__main__":
    main()
