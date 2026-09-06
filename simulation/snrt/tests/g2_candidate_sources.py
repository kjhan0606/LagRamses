#!/usr/bin/env python3
"""Integration checks for the staged G2 candidate-source audit."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_candidate_sources import audit_candidates  # noqa: E402


def _assert_manifest_mutation_is_fatal() -> None:
    project = ROOT.parents[1]
    source_root = project / "external" / "g2_candidates"
    manifest = source_root / "acquisition_manifest_v1.json"
    with tempfile.TemporaryDirectory(prefix="snrt-g2-manifest-") as directory:
        mutated = Path(directory) / manifest.name
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        payload["candidates"][0]["files"][0]["sha256"] = "0" * 64
        mutated.write_text(json.dumps(payload), encoding="utf-8")
        report = audit_candidates(source_root, mutated)
        assert report["status"] == "candidate_review_blocked_input_integrity"
        assert report["input_integrity_passed"] is False
        assert report["audit_failures"]
        result = subprocess.run(
            [
                sys.executable,
                str(TOOLS / "audit_g2_candidate_sources.py"),
                "--root",
                str(source_root),
                "--manifest",
                str(mutated),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2
        cli_report = json.loads(result.stdout)
        assert cli_report["status"] == "candidate_review_blocked_input_integrity"
        for index, incomplete in enumerate(
            (
                {},
                {"candidates": []},
                {"candidates": [{"candidate_id": "limongi_chieffi_2018_cds", "files": []}]},
            )
        ):
            coverage_manifest = Path(directory) / f"coverage-{index}.json"
            coverage_manifest.write_text(json.dumps(incomplete), encoding="utf-8")
            coverage_report = audit_candidates(source_root, coverage_manifest)
            assert coverage_report["status"] == "candidate_review_blocked_input_integrity"
            assert coverage_report["input_integrity_passed"] is False


def _assert_inline_parser_mutations_are_fatal() -> None:
    project = ROOT.parents[1]
    source_root = project / "external" / "g2_candidates"
    original_manifest = source_root / "acquisition_manifest_v1.json"
    mutations = (
        ("limongi_chieffi_2018_cds", "table8.dat", "limongi"),
        ("nugrid_set1ext", "element_yield_table_MESAonly_fryer12_delay_total.txt", "nugrid"),
    )
    for candidate_dir_name, filename, kind in mutations:
        with tempfile.TemporaryDirectory(prefix=f"snrt-g2-{kind}-parser-") as directory:
            farm = Path(directory)
            for candidate_dir in source_root.iterdir():
                if not candidate_dir.is_dir():
                    continue
                target_dir = farm / candidate_dir.name
                target_dir.mkdir()
                for source_file in candidate_dir.iterdir():
                    (target_dir / source_file.name).symlink_to(source_file)
            manifest_payload = json.loads(original_manifest.read_text(encoding="utf-8"))
            relative = f"{candidate_dir_name}/{filename}"
            target = farm / relative
            target.unlink()
            content = (source_root / relative).read_text(encoding="utf-8")
            if kind == "limongi":
                rows = content.splitlines()
                for row_index, row in enumerate(rows):
                    if row and not row.lstrip().startswith("#"):
                        fields = row.split()
                        fields[-1] = "not-a-number"
                        rows[row_index] = " ".join(fields)
                        break
                content = "\n".join(rows) + "\n"
            else:
                marker = "H Lifetime:"
                line_index = next(index for index, row in enumerate(content.splitlines()) if marker in row)
                rows = content.splitlines()
                rows[line_index] = "H Lifetime: not-a-number"
                content = "\n".join(rows) + "\n"
            target.write_text(content, encoding="utf-8")
            digest = hashlib.sha256(target.read_bytes()).hexdigest()
            for candidate in manifest_payload["candidates"]:
                for entry in candidate.get("files", []):
                    if entry.get("path") == relative:
                        entry["bytes"] = target.stat().st_size
                        entry["sha256"] = digest
            manifest = farm / original_manifest.name
            manifest.write_text(json.dumps(manifest_payload), encoding="utf-8")
            report = audit_candidates(farm, manifest)
            assert report["status"] == "candidate_review_blocked_input_integrity"
            assert report["input_integrity_passed"] is False
            assert any(failure["reason"] == "source_parse_error" for failure in report["audit_failures"])


def main() -> int:
    report = audit_candidates(ROOT.parents[1] / "external" / "g2_candidates")
    assert report["status"] == "candidate_review_only"
    assert report["production_ready"] is False
    assert report["acquisition_manifest"]["status"] == "pass"
    assert report["acquisition_manifest"]["file_count"] == 65

    limongi = report["candidates"]["limongi_chieffi_2018_cds"]
    assert limongi["files"]["recommended_isotopic_yields"]["row_count"] == 3996
    assert limongi["files"]["wind_isotopic_yields"]["row_count"] == 3996
    assert limongi["coverage"]["reported_zams_mass_msun"] == [13.0, 120.0]
    assert "no_age_resolved_cumulative_release_history" in limongi["blockers"]

    nugrid = report["candidates"]["nugrid_set1ext_mesaonly_fryer12_delay"]
    total = nugrid["files"]["total"]
    assert total["block_count"] == 61
    assert total["duplicate_coordinates"] == [[5.0, 0.01]]
    assert total["species_count_per_block"] == [80]
    assert "duplicate_mass_metallicity_coordinate" in nugrid["blockers"]
    assert "no_age_resolved_cumulative_release_history" in nugrid["blockers"]

    huscher = report["candidates"]["huscher2025_agb"]
    assert huscher["source_identity"]["license"] == "cc-by-4.0"
    assert huscher["single_star_grid"]["model_count"] == 120
    assert huscher["population_tables"]["normalization_semantics_pass"] is False
    assert huscher["canonical_rows_emitted"] == 0

    boccioli = report["candidates"]["boccioli_roberti2026_neutrino_ccsn"]
    assert boccioli["source_identity"]["license"] == "cc-by-4.0"
    assert boccioli["grids"]["F23_single"]["model_count"] == 35
    assert boccioli["quality_findings"]["f23_component_mass_closure_pass"] is True
    assert boccioli["quality_findings"]["lc18_readme_consistency_pass"] is False
    assert boccioli["canonical_rows_emitted"] == 0

    doherty = report["candidates"]["doherty2014_sagb"]
    assert doherty["primary_grid"]["model_count"] == 20
    assert doherty["primary_grid"]["tracked_elements_absent"] == ["Ca"]
    assert doherty["mass_closure"]["pass"] is True
    assert doherty["quality_findings"]["source_label_repair_applied"] is False
    assert doherty["canonical_rows_emitted"] == 0

    stockinger = report["candidates"]["stockinger2020_low_mass_ccsn"]
    assert stockinger["model_grid"]["zams_mass_msun"] == [8.8, 9.0, 9.6]
    assert stockinger["yield_mass_closure"]["pass"] is True
    assert stockinger["diagnostic_explosion_energy"]["vsh_quarantined"] is True
    assert stockinger["canonical_rows_emitted"] == 0

    sukhbold = report["candidates"]["sukhbold2016_ccsn"]
    assert sukhbold["z96_grid"]["model_count"] == 13
    assert sukhbold["z96_grid"]["zams_mass_msun"] == [
        9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5,
        10.75, 11.0, 11.25, 11.5, 11.75, 12.0,
    ]
    assert sukhbold["mass_budget_review"]["within_review_bound"] is True
    assert sukhbold["mass_budget_review"]["exact_mass_closure_claimed"] is False
    assert sukhbold["canonical_rows_emitted"] == 0

    transition = report["candidates"]["limongi2024_transition_fates"]
    assert transition["source_identity"]["license"] == "CC BY 4.0"
    assert transition["machine_readable_tp_table"]["data_row_count"] == 963
    assert transition["source_reported_fate_statements"]["potential_ecsn_is_not_a_deterministic_event_assignment"] is True
    assert transition["project_transition_policy"]["unresolved_runtime_edge_interval_msun"] == [8.0, 8.8]
    assert transition["canonical_rows_emitted"] == 0

    roberti = report["candidates"]["roberti2024_ultralowz_ccsn"]
    assert roberti["source_identity"]["license"] == "CC BY 4.0"
    assert roberti["source_grid"]["model_count"] == 34
    assert roberti["source_grid"]["masses_msun"] == [15.0, 25.0]
    assert roberti["source_grid"]["metallicity_mass_fraction"] == [0.0, 3.236e-7, 3.236e-6]
    assert roberti["yield_model_inventory"]["official_mrt_model_count"] == 30
    assert roberti["yield_model_inventory"]["source_only_models_missing_from_official_mrt"] == [
        "015z300", "015z600", "025z450", "025z700",
    ]
    assert roberti["mass_budget_review"]["outlier_models"] == ["025z600"]
    assert roberti["canonical_rows_emitted"] == 0

    heger_woosley = report["candidates"]["heger_woosley2010_popiii"]
    assert heger_woosley["source_grid"]["record_count"] == 660546
    assert heger_woosley["source_grid"]["coordinate_count"] == 5760
    assert heger_woosley["source_grid"]["zams_mass_count"] == 120
    assert heger_woosley["source_grid"]["zams_mass_msun_minimum"] == 10.0
    assert heger_woosley["source_grid"]["zams_mass_msun_maximum"] == 100.0
    assert heger_woosley["physical_semantics"]["canonical_event_energy_selected"] is False
    assert heger_woosley["canonical_rows_emitted"] == 0

    _assert_manifest_mutation_is_fatal()
    _assert_inline_parser_mutations_are_fatal()

    print("G2_CANDIDATE_SOURCE_AUDIT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
