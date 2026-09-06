#!/usr/bin/env python3
"""Checks for the G2 reduced-chemistry mass-scope audit."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from adapt_g2_candidate_sources import LIMONGI_ID, NUGRID_ID  # noqa: E402
from audit_g2_reduced_chemistry_scope import audit_reduced_chemistry_scope  # noqa: E402


def main() -> int:
    report = audit_reduced_chemistry_scope(root=ROOT.parents[1] / "external" / "g2_candidates")
    assert report["status"] == "blocked_decay_horizon_and_source_approval_required"
    assert report["canonical_ejecta_sum_equals_returned_mass_possible_with_only_tracked_elements"] is False
    assert report["radioactive_decay_applied"] is False
    assert report["untracked_ejecta_residual_contract_implemented"] is True
    assert report["maximum_observed_omitted_mass_fraction"] > 0.0

    limongi = report["candidate_component_summaries"][LIMONGI_ID]
    assert limongi["source_supported_wind"]["record_count"] == 108
    assert limongi["source_supported_terminal_set_R"]["record_count"] == 48
    assert "Ni" in limongi["source_supported_terminal_set_R"]["omitted_elements_by_unweighted_source_grid_sum"]

    nugrid = report["candidate_component_summaries"][NUGRID_ID]
    assert nugrid["integrated_agb_ejecta_candidate"]["record_count"] == 41
    assert nugrid["massive_star_wind"]["record_count"] == 20
    assert nugrid["delayed_explosion_terminal_ejecta"]["record_count"] == 20
    assert all(
        component["maximum_absolute_tracked_plus_omitted_closure_residual"] < 1.0e-12
        for component in nugrid.values()
    )

    print("G2_REDUCED_CHEMISTRY_SCOPE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
