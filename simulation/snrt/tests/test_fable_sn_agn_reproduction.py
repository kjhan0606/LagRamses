#!/usr/bin/env python3
"""Regression assertions for the independent Fable-finding reproduction."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from reproduce_fable_sn_agn_findings import reproduce  # noqa: E402


EXPECTED_REPRODUCED = {
    "F1", "F2", "F5", "F6", "F9", "F10", "F11", "F12", "F15", "F16",
    "F17",
}
EXPECTED_PARTIAL = {"F13", "F14"}
EXPECTED_NOT_REPRODUCED = {"F3", "F4", "F7", "F8"}


def main() -> int:
    payload = reproduce()
    assert payload["independent_checks"]["production_makefile_selects_patch_tree"]
    assert payload["independent_checks"]["g1_runner_selects_separate_native_mirror"]
    assert payload["independent_checks"]["g1_runner_excludes_ramses_runtime"]
    assert payload["independent_checks"]["compiled_runtime_uses_gyr"]
    assert payload["independent_checks"]["table_axis_declares_years"]
    assert not payload["independent_checks"]["compiled_interpolator_clamps"]
    assert payload["independent_checks"]["production_converts_age_once"]
    assert payload["independent_checks"]["current_interval_contract"]
    assert payload["independent_checks"]["production_requires_external_table"]
    assert payload["independent_checks"]["mirror_converts_year_axis"]
    assert payload["independent_checks"]["mirror_rejects_domain"]
    assert payload["independent_checks"]["production_nvar18"]
    assert payload["independent_checks"]["production_nener0"]
    assert payload["independent_checks"]["hydro_index_formula_present"]
    assert payload["independent_checks"]["mirror_uses_inener"]
    assert payload["independent_checks"]["compiled_literal_energy_field"]

    assert set(payload["summary"]["reproduced"]) == EXPECTED_REPRODUCED
    assert set(payload["summary"]["partially_reproduced"]) == EXPECTED_PARTIAL
    assert set(payload["summary"]["not_reproduced"]) == EXPECTED_NOT_REPRODUCED

    interval = payload["numerical_reproductions"]["forward_cumulative_interval"]
    assert interval["implementation_total"] == 64.0
    assert interval["correct_telescoping_total"] == 36.0
    assert interval["reproduced"] is True

    units = payload["numerical_reproductions"]["year_gyr_coordinate_mismatch"]
    assert units["compiled_query_if_untagged"] == 1.0
    assert units["intended_year_axis_query"] == 1.0e9
    assert units["coordinate_ratio"] == 1.0e-9
    assert units["reproduced"] is True

    print(
        "FABLE_SN_AGN_INDEPENDENT_REPRODUCTION_OK "
        f"reproduced={len(EXPECTED_REPRODUCED)} "
        f"partial={len(EXPECTED_PARTIAL)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
