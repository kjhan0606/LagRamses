#!/usr/bin/env python3
"""Integration checks for lossless, fail-closed G2 source adapters."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from adapt_g2_candidate_sources import (  # noqa: E402
    LIMONGI_ID,
    NUGRID_ID,
    SourceAdapterError,
    adapt_candidate,
    require_canonical_promotion_allowed,
)
from convert_yield_rows_to_canonical import ConversionError, convert  # noqa: E402


def _assert_review_only(report: dict) -> None:
    assert report["status"] == "review_only"
    assert report["canonical_conversion_allowed"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["blockers"]
    availability = report["source_availability"]
    assert all(value is False for value in availability.values())
    try:
        require_canonical_promotion_allowed(report)
    except SourceAdapterError as exc:
        assert "promotion refused" in str(exc)
    else:
        raise AssertionError("review-only source was allowed to promote")


def main() -> int:
    candidate_root = ROOT.parents[1] / "external" / "g2_candidates"

    limongi = adapt_candidate(LIMONGI_ID, root=candidate_root, include_records=True)
    _assert_review_only(limongi)
    assert limongi["verified_acquisition"]["verified_file_count"] == 5
    components = limongi["source_components"]
    assert components["recommended_yields"]["model_count"] == 108
    assert components["wind_yields"]["model_count"] == 48
    assert components["recommended_yields"]["species_count_per_model"] == [333]
    assert components["wind_yields"]["species_count_per_model"] == [333]
    assert len(components["recommended_yields"]["records"]) == 108
    assert len(components["wind_yields"]["records"]) == 48
    assert components["evolutionary_properties"]["row_count"] == 864
    assert components["evolutionary_properties"]["duplicate_model_phase_coordinate_count"] == 10
    assert components["evolutionary_properties"]["all_duplicate_rows_physically_identical"] is True
    assert all(
        value["physical_values_exactly_identical"]
        for value in components["evolutionary_properties"]["duplicate_model_phase_coordinates"]
    )
    assert components["presupernova_properties"]["row_count"] == 96
    assert limongi["source_axes"]["metallicity_mass_fraction_from_source_article"] == {
        "-3": 3.236e-5,
        "-2": 3.236e-4,
        "-1": 3.236e-3,
        "0": 1.345e-2,
    }
    assert "metallicity_feh_to_mass_fraction_requires_approved_solar_mapping" not in limongi["blockers"]
    assert "duplicate_evolutionary_model_phase_coordinates" not in limongi["blockers"]
    relation = limongi["uninterpreted_component_relation_diagnostics"]
    assert relation["overlapping_model_count"] == 48
    assert relation["source_semantics_supported"] is True
    assert relation["runtime_channel_assignment_approved"] is False
    assert "canonical_injected_energy_absent" in limongi["blockers"]
    assert "source_reported_yield_unit_requires_project_approval_because_cds_unit_is_blank" in limongi["blockers"]
    assert limongi["source_semantics_evidence"]["source_author_definitions"]["table8"]
    assert (
        limongi["source_semantics_evidence"]["source_use_terms_evidence"]
        ["candidate_record"]["production_license_status"]
        == "not_approved"
    )

    nugrid = adapt_candidate(NUGRID_ID, root=candidate_root, include_records=True)
    _assert_review_only(nugrid)
    assert nugrid["verified_acquisition"]["verified_file_count"] == 11
    assert nugrid["component_coordinate_sequences_identical"] is True
    for component in nugrid["source_components"].values():
        assert component["block_count"] == 61
        assert component["unique_coordinate_count"] == 60
        assert component["elements_per_block"] == [80]
        assert component["duplicate_coordinates"] == [
            {"initial_mass_msun": 5.0, "metallicity_mass_fraction": 0.01, "multiplicity": 2}
        ]
        assert len(component["records"]) == 61
    assert "duplicate_mass_metallicity_coordinate" in nugrid["blockers"]
    relation = nugrid["uninterpreted_component_relation_diagnostics"]
    assert relation["aligned_block_count"] == 61
    assert relation["source_semantics_supported"] is True
    assert relation["runtime_channel_assignment_approved"] is False
    assert relation["pre_explosion_minus_winds_negative_value_count"] == 0
    assert (
        nugrid["source_semantics_evidence"]["source_use_terms_evidence"]
        ["candidate_record"]["production_license_status"]
        == "not_approved"
    )

    # The generic canonical converter must also refuse an adapter document.
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        source = temporary / "review.json"
        source.write_text(json.dumps(nugrid), encoding="utf-8")
        try:
            convert(source, temporary / "canonical.dat", temporary / "canonical.dat.json")
        except ConversionError as exc:
            assert "source" in str(exc)
        else:
            raise AssertionError("generic canonical converter accepted a review-only adapter document")

    print("G2_SOURCE_ADAPTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
