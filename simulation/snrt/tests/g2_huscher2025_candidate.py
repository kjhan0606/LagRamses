#!/usr/bin/env python3
"""Integration checks for the review-only Huscher et al. (2025) adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_huscher2025_candidate import audit_huscher2025_candidate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_huscher2025_candidate()

    assert report["status"] == (
        "candidate_acquired_license_verified_population_normalization_blocked"
    )
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["license"] == "cc-by-4.0"
    assert report["source_identity"]["license_verified"] is True
    assert report["source_identity"]["archive_sha256"] == (
        "dc559ee272d602bcfe95ab0050cb388eed670986e3e62234b4bd9126d0128199"
    )

    grid = report["single_star_grid"]
    assert grid["model_count"] == 120
    assert grid["mass_msun"] == [
        0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0
    ]
    assert len(grid["metallicity_mass_fraction"]) == 10
    assert len(grid["isotopes"]) == 16
    assert grid["tracked_elements_absent"] == ["S", "Ca", "Fe"]
    assert grid["negative_source_species_residual_count"] == 38
    assert grid["negative_residual_outside_printed_quantization_count"] == 0
    assert grid["mass_coverage_gap_msun"] == [7.0, 8.0]

    population = report["population_tables"]
    assert population["row_count_per_table"] == 153
    assert population["normalization_semantics_pass"] is False
    assert population["normalization_inference_applied"] is False
    assert population["metallicity_columns_exceeding_unit_return"] == 10
    assert 1326.0 < population["integrated_return_minimum"] < 1327.0
    assert 2895.0 < population["integrated_return_maximum"] < 2896.0
    assert report["semantic_firewalls"]["second_imf_convolution_forbidden"] is True

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_HUSCHER2025_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
