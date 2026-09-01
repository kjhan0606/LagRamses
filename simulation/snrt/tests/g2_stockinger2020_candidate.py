#!/usr/bin/env python3
"""Integration checks for the review-only Stockinger et al. (2020) adapter."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_stockinger2020_candidate import audit_stockinger2020_candidate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_stockinger2020_candidate()

    assert report["status"] == "candidate_acquired_energy_yields_audited_license_unresolved"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["file_count"] == 12
    assert report["source_identity"]["redistribution_license_verified"] is False

    grid = report["model_grid"]
    assert grid["model_count"] == 3
    assert grid["zams_mass_msun"] == [8.8, 9.0, 9.6]
    assert grid["cross_model_interpolation_allowed"] is False
    for model in grid["models"].values():
        assert model["tracked_elements_absent"] == ["N"]

    closure = report["yield_mass_closure"]
    assert closure["pass"] is True
    assert closure["maximum_absolute_species_sum_residual_msun"] < 0.0011
    energy = report["diagnostic_explosion_energy"]
    assert energy["vsh_metadata_semantics_pass"] is False
    assert energy["vsh_quarantined"] is True
    expected = {"e8.8": 0.9539027937090561, "z9.6": 0.8602597812628994, "s9.0": 0.5250178516256706}
    for model, value in expected.items():
        observed = energy["models"][model]["last_finite_diagnostic_energy_1e50_erg"]
        assert math.isclose(observed, value, rel_tol=0.0, abs_tol=1.0e-15)
        assert energy["models"][model]["vsh_used"] is False
    assert report["quality_findings"]["radioactive_decay_projection_complete"] is False

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_STOCKINGER2020_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
