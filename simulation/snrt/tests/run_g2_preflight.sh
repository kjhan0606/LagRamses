#!/usr/bin/env bash
set -euo pipefail

# G2 preflight is intentionally allowed to finish in BLOCKED state.  It proves
# that configuration and conversion controls work while the legacy table and
# one-point fixture remain refused as physical production inputs.

ROOT=/gpfs/kjhan/LRD_JWST
SNRT_ROOT="$ROOT/simulation/snrt"
BUILD_DIR="$ROOT/build/g2_config"
DATA_DIR="$SNRT_ROOT/data"
mkdir -p "$BUILD_DIR" "$DATA_DIR"

mpiifx -O2 -g -traceback -warn all -check all -fpp \
  -module "$BUILD_DIR" \
  -c "$SNRT_ROOT/native/phase0/stellar_enrichment_config.f90" \
  -o "$BUILD_DIR/stellar_enrichment_config.o"
mpiifx -O2 -g -traceback -warn all -check all -fpp \
  -module "$BUILD_DIR" \
  -c "$SNRT_ROOT/native/phase0/g2_configuration_test.f90" \
  -o "$BUILD_DIR/g2_configuration_test.o"
mpiifx -O2 -g -traceback -check all \
  "$BUILD_DIR/g2_configuration_test.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/g2_configuration_test"
(cd "$BUILD_DIR" && ./g2_configuration_test)

"$ROOT/tests/run_stellar_feedback_policy_unit.sh"
"$SNRT_ROOT/tests/run_g2_population_ledger.sh"

"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/stellar_yield_asset.py"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/yield_converter.py"

"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_g2_candidate_sources.py" \
  --json-out "$DATA_DIR/g2_candidate_source_audit.json" \
  > "$BUILD_DIR/g2_candidate_source_audit.stdout"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_g2_source_package_fingerprints.py" \
  --json-out "$DATA_DIR/g2_source_package_fingerprint_audit.json" \
  > "$BUILD_DIR/g2_source_package_fingerprint_audit.stdout"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_source_selection_gate.py"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_g2_source_selection_gate.py" \
  --json-out "$DATA_DIR/g2_source_selection_gate.json" \
  > "$BUILD_DIR/g2_source_selection_gate.stdout"

"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_source_adapters.py"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_source_adapter_closure.py"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_reduced_chemistry_scope.py"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_limongi_decay_projection.py" \
  --json-out "$DATA_DIR/g2_limongi_decay_projection_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_nugrid_channel_projection.py" \
  --json-out "$DATA_DIR/g2_nugrid_channel_projection_review.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_feedback_energetics_sensitivity.py" \
  --json-out "$DATA_DIR/g2_feedback_energetics_sensitivity_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_huscher2025_candidate.py" \
  --json-out "$DATA_DIR/g2_huscher2025_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_boccioli_roberti2026_candidate.py" \
  --json-out "$DATA_DIR/g2_boccioli_roberti2026_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_doherty2014_sagb_candidate.py" \
  --json-out "$DATA_DIR/g2_doherty2014_sagb_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_stockinger2020_candidate.py" \
  --json-out "$DATA_DIR/g2_stockinger2020_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_sukhbold2016_candidate.py" \
  --json-out "$DATA_DIR/g2_sukhbold2016_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_sukhbold_channel_projection.py" \
  --json-out "$DATA_DIR/g2_sukhbold_channel_projection_review.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_limongi2024_transition_fates.py" \
  --json-out "$DATA_DIR/g2_limongi2024_transition_fate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_roberti2024_ultralowz_candidate.py" \
  --json-out "$DATA_DIR/g2_roberti2024_ultralowz_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_heger_woosley2010_popiii_candidate.py" \
  --json-out "$DATA_DIR/g2_heger_woosley2010_popiii_candidate_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_baseline_metallicity_demand.py" \
  --json-out "$DATA_DIR/g2_baseline_metallicity_demand_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_candidate_grid_coverage.py" \
  --json-out "$DATA_DIR/g2_candidate_grid_coverage_audit.json"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tests/g2_limongi_phase_mass_history.py" \
  --json-out "$DATA_DIR/g2_limongi_phase_mass_history_audit.json"

# F-P1 is part of the aggregate G2 feedback preflight.  Its checksum-bound
# package gate intentionally fails if regenerated candidate evidence drifts.
"$SNRT_ROOT/tests/run_fp1_population_fate_contract.sh"

"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/adapt_g2_candidate_sources.py" \
  limongi_chieffi_2018_cds \
  --json-out "$DATA_DIR/g2_limongi_source_adapter_review.json" \
  > "$BUILD_DIR/g2_limongi_source_adapter_review.stdout"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/adapt_g2_candidate_sources.py" \
  nugrid_set1ext_mesaonly_fryer12_delay \
  --json-out "$DATA_DIR/g2_nugrid_source_adapter_review.json" \
  > "$BUILD_DIR/g2_nugrid_source_adapter_review.stdout"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_g2_source_adapter_closure.py" \
  --json-out "$DATA_DIR/g2_source_adapter_closure_audit.json" \
  > "$BUILD_DIR/g2_source_adapter_closure_audit.stdout"
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_g2_reduced_chemistry_scope.py" \
  --json-out "$DATA_DIR/g2_reduced_chemistry_scope_audit.json" \
  > "$BUILD_DIR/g2_reduced_chemistry_scope_audit.stdout"

set +e
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_stellar_yield_asset.py" \
  /home/kjhan/BACKUP/lagRamses/patch/lagRamses/phase0_validation_yields.dat \
  --json-out "$DATA_DIR/g2_phase0_fixture_audit.json" \
  > "$BUILD_DIR/g2_fixture_audit.stdout"
fixture_status=$?
"$SNRT_ROOT/.venv/bin/python" "$SNRT_ROOT/tools/audit_stellar_yield_asset.py" \
  /gpfs/kjhan/Run_JWST/opt_run/yield_table.asc \
  --json-out "$DATA_DIR/g2_legacy_asset_audit.json" \
  > "$BUILD_DIR/g2_legacy_audit.stdout"
legacy_status=$?
set -e

if [[ "$fixture_status" -ne 2 || "$legacy_status" -ne 2 ]]; then
  echo "G2 preflight expected both non-production assets to be rejected" >&2
  exit 1
fi

"$SNRT_ROOT/.venv/bin/python" - "$DATA_DIR/g2_preflight.json" <<'PY'
import json
import sys
from pathlib import Path

output = Path(sys.argv[1])
fixture = json.loads((output.parent / "g2_phase0_fixture_audit.json").read_text())
legacy = json.loads((output.parent / "g2_legacy_asset_audit.json").read_text())
candidate_sources = json.loads((output.parent / "g2_candidate_source_audit.json").read_text())
fingerprints = json.loads((output.parent / "g2_source_package_fingerprint_audit.json").read_text())
selection = json.loads((output.parent / "g2_source_selection_gate.json").read_text())
limongi_adapter = json.loads((output.parent / "g2_limongi_source_adapter_review.json").read_text())
nugrid_adapter = json.loads((output.parent / "g2_nugrid_source_adapter_review.json").read_text())
adapter_closure = json.loads((output.parent / "g2_source_adapter_closure_audit.json").read_text())
chemistry_scope = json.loads((output.parent / "g2_reduced_chemistry_scope_audit.json").read_text())
decay_projection = json.loads((output.parent / "g2_limongi_decay_projection_audit.json").read_text())
nugrid_projection = json.loads((output.parent / "g2_nugrid_channel_projection_review.json").read_text())
energetics = json.loads((output.parent / "g2_feedback_energetics_sensitivity_audit.json").read_text())
huscher = json.loads((output.parent / "g2_huscher2025_candidate_audit.json").read_text())
boccioli = json.loads((output.parent / "g2_boccioli_roberti2026_candidate_audit.json").read_text())
doherty = json.loads((output.parent / "g2_doherty2014_sagb_candidate_audit.json").read_text())
stockinger = json.loads((output.parent / "g2_stockinger2020_candidate_audit.json").read_text())
sukhbold = json.loads((output.parent / "g2_sukhbold2016_candidate_audit.json").read_text())
sukhbold_projection = json.loads((output.parent / "g2_sukhbold_channel_projection_review.json").read_text())
limongi_transition = json.loads((output.parent / "g2_limongi2024_transition_fate_audit.json").read_text())
roberti_ultralowz = json.loads((output.parent / "g2_roberti2024_ultralowz_candidate_audit.json").read_text())
heger_woosley_popiii = json.loads((output.parent / "g2_heger_woosley2010_popiii_candidate_audit.json").read_text())
baseline_metallicity_demand = json.loads((output.parent / "g2_baseline_metallicity_demand_audit.json").read_text())
coverage = json.loads((output.parent / "g2_candidate_grid_coverage_audit.json").read_text())
limongi_mass_history = json.loads((output.parent / "g2_limongi_phase_mass_history_audit.json").read_text())
report = {
    "status": "blocked",
    "gate": "G2",
    "source_selection_matrix": "config/g2_source_selection_matrix_v1.json",
    "physical_contract": "config/g2_physics_contract_v1.json",
    "configuration_test": "G2_CONFIGURATION_TEST_OK",
    "stellar_feedback_policy_test": "stellar feedback policy: PASS",
    "population_ledger_test": "G2_POPULATION_LEDGER_RUN_OK",
    "yield_asset_test": "STELLAR_YIELD_ASSET_TEST_OK",
    "converter_test": "YIELD_CONVERTER_TEST_OK",
    "source_adapter_test": "G2_SOURCE_ADAPTER_TEST_OK",
    "source_adapter_closure_test": "G2_SOURCE_ADAPTER_CLOSURE_TEST_OK",
    "reduced_chemistry_scope_test": "G2_REDUCED_CHEMISTRY_SCOPE_TEST_OK",
    "limongi_decay_projection_test": "G2_LIMONGI_DECAY_PROJECTION_TEST_OK",
    "nugrid_channel_projection_test": "G2_NUGRID_CHANNEL_PROJECTION_TEST_OK",
    "feedback_energetics_sensitivity_test": "G2_FEEDBACK_ENERGETICS_SENSITIVITY_TEST_OK",
    "huscher2025_candidate_test": "G2_HUSCHER2025_CANDIDATE_TEST_OK",
    "boccioli_roberti2026_candidate_test": "G2_BOCCIOLI_ROBERTI2026_CANDIDATE_TEST_OK",
    "doherty2014_sagb_candidate_test": "G2_DOHERTY2014_SAGB_CANDIDATE_TEST_OK",
    "stockinger2020_candidate_test": "G2_STOCKINGER2020_CANDIDATE_TEST_OK",
    "sukhbold2016_candidate_test": "G2_SUKHBOLD2016_CANDIDATE_TEST_OK",
    "sukhbold2016_channel_projection_test": "G2_SUKHBOLD_CHANNEL_PROJECTION_TEST_OK",
    "limongi2024_transition_fate_test": "G2_LIMONGI2024_TRANSITION_FATE_TEST_OK",
    "roberti2024_ultralowz_candidate_test": "G2_ROBERTI2024_ULTRALOWZ_CANDIDATE_TEST_OK",
    "heger_woosley2010_popiii_candidate_test": "G2_HEGER_WOOSLEY2010_POPIII_CANDIDATE_TEST_OK",
    "baseline_metallicity_demand_test": "G2_BASELINE_METALLICITY_DEMAND_TEST_OK",
    "candidate_grid_coverage_test": "G2_CANDIDATE_GRID_COVERAGE_TEST_OK",
    "limongi_phase_mass_history_test": "G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK",
    "fp1_population_fate_contract_test": "FP1_POPULATION_FATE_CONTRACT_OK",
    "fixture_status": fixture["status"],
    "fixture_blocking_reasons": fixture["production_gate"]["blocking_reasons"],
    "legacy_status": legacy["status"],
    "legacy_blocking_reasons": legacy["production_gate"]["blocking_reasons"],
    "candidate_source_status": candidate_sources["status"],
    "candidate_source_input_integrity_passed": candidate_sources["input_integrity_passed"],
    "candidate_source_audit_failures": candidate_sources["audit_failures"],
    "candidate_source_manifest_status": candidate_sources["acquisition_manifest"]["status"],
    "candidate_source_manifest_file_count": candidate_sources["acquisition_manifest"]["file_count"],
    "candidate_source_candidate_count": len(candidate_sources["candidates"]),
    "candidate_source_blocking_reasons": candidate_sources["blockers"],
    "candidate_source_fingerprint_status": fingerprints["status"],
    "candidate_source_fingerprint_input_integrity_passed": fingerprints["input_integrity_passed"],
    "candidate_source_fingerprint_candidate_count": fingerprints["candidate_count"],
    "candidate_source_fingerprint_file_count": fingerprints["file_count"],
    "candidate_source_fingerprint_scheme": fingerprints["scheme"],
    "candidate_source_fingerprint_blocking_reasons": fingerprints["audit_failures"],
    "source_selection_gate_status": selection["status"],
    "source_selection_gate_runtime_activation_allowed": selection["runtime_activation_allowed"],
    "source_selection_review_validation_branch": selection["review_validation_branch"]["candidate_id"],
    "source_selection_review_validation_branch_sha256": selection["review_validation_branch"]["composite_sha256"],
    "source_selection_production_source_id": selection["production_source_id"],
    "source_selection_production_approval_id": selection["production_approval_id"],
    "source_selection_gate_blocking_reasons": selection["blockers"],
    "source_adapter_status": {
        "limongi_chieffi_2018_cds": limongi_adapter["status"],
        "nugrid_set1ext_mesaonly_fryer12_delay": nugrid_adapter["status"],
    },
    "source_adapter_canonical_rows_emitted": {
        "limongi_chieffi_2018_cds": limongi_adapter["canonical_rows_emitted"],
        "nugrid_set1ext_mesaonly_fryer12_delay": nugrid_adapter["canonical_rows_emitted"],
    },
    "source_adapter_blocking_reasons": {
        "limongi_chieffi_2018_cds": limongi_adapter["blockers"],
        "nugrid_set1ext_mesaonly_fryer12_delay": nugrid_adapter["blockers"],
    },
    "source_adapter_closure_status": adapter_closure["status"],
    "source_adapter_closure_gate_blockers": adapter_closure["gate_blockers"],
    "reduced_chemistry_scope_status": chemistry_scope["status"],
    "maximum_observed_omitted_mass_fraction": chemistry_scope["maximum_observed_omitted_mass_fraction"],
    "untracked_ejecta_residual_contract_implemented": chemistry_scope["untracked_ejecta_residual_contract_implemented"],
    "limongi_decay_projection_status": decay_projection["status"],
    "limongi_decay_projection_unresolved_nuclides": decay_projection["source_isotope_coverage"]["unresolved_count"],
    "nugrid_channel_projection_status": nugrid_projection["status"],
    "nugrid_partial_grid_rows": nugrid_projection["row_audit"]["row_count"],
    "nugrid_partial_grid_mass_residual_max_msun": nugrid_projection["row_audit"]["maximum_absolute_source_population_mass_residual_msun"],
    "feedback_energetics_sensitivity_status": energetics["status"],
    "feedback_source_nulls_preserved": energetics["input_projection"]["source_energy_and_momentum_nulls_preserved"],
    "feedback_momentum_semantics": energetics["source_frame_vector_momentum"],
    "feedback_sukhbold_energy_anchor_count": len(energetics["machine_readable_mass_grid_energy_anchors"]["records"]),
    "feedback_sukhbold_exact_launch_momentum_derived": energetics["machine_readable_mass_grid_energy_anchors"]["exact_ejecta_launch_momentum_derived"],
    "feedback_roberti_ultralowz_energy_anchor_count": len(energetics["source_table_ultralowz_energy_anchors"]["records"]),
    "feedback_roberti_ultralowz_exact_launch_momentum_derived": energetics["source_table_ultralowz_energy_anchors"]["exact_ejecta_launch_momentum_derived"],
    "feedback_heger_woosley_popiii_coordinate_count": energetics["source_table_popiii_energy_grid"]["coordinate_count"],
    "feedback_heger_woosley_popiii_exact_launch_momentum_derived": energetics["source_table_popiii_energy_grid"]["exact_ejecta_launch_momentum_derived"],
    "huscher2025_candidate_status": huscher["status"],
    "huscher2025_license_verified": huscher["source_identity"]["license_verified"],
    "huscher2025_single_star_model_count": huscher["single_star_grid"]["model_count"],
    "huscher2025_population_normalization_pass": huscher["population_tables"]["normalization_semantics_pass"],
    "huscher2025_population_integrated_return_range_under_claimed_units": [
        huscher["population_tables"]["integrated_return_minimum"],
        huscher["population_tables"]["integrated_return_maximum"],
    ],
    "huscher2025_blocking_reasons": huscher["blockers"],
    "boccioli_roberti2026_candidate_status": boccioli["status"],
    "boccioli_roberti2026_license_verified": boccioli["source_identity"]["license_verified"],
    "boccioli_roberti2026_model_count": sum(
        value["model_count"] for value in boccioli["grids"].values()
    ),
    "boccioli_roberti2026_f23_component_mass_closure_pass": boccioli["quality_findings"]["f23_component_mass_closure_pass"],
    "boccioli_roberti2026_lc18_readme_consistency_pass": boccioli["quality_findings"]["lc18_readme_consistency_pass"],
    "boccioli_roberti2026_lc18_failed_wind_omission_count": boccioli["quality_findings"]["lc18_failed_models_with_reported_wind_but_zero_wind_table_count"],
    "boccioli_roberti2026_blocking_reasons": boccioli["blockers"],
    "doherty2014_candidate_status": doherty["status"],
    "doherty2014_model_count": doherty["primary_grid"]["model_count"],
    "doherty2014_mass_closure_pass": doherty["mass_closure"]["pass"],
    "doherty2014_tracked_elements_absent": doherty["primary_grid"]["tracked_elements_absent"],
    "doherty2014_source_label_repair_applied": doherty["quality_findings"]["source_label_repair_applied"],
    "doherty2014_blocking_reasons": doherty["blockers"],
    "stockinger2020_candidate_status": stockinger["status"],
    "stockinger2020_model_count": stockinger["model_grid"]["model_count"],
    "stockinger2020_yield_mass_closure_pass": stockinger["yield_mass_closure"]["pass"],
    "stockinger2020_last_finite_diagnostic_energy_erg": {
        model: value["last_finite_diagnostic_energy_erg"]
        for model, value in stockinger["diagnostic_explosion_energy"]["models"].items()
    },
    "stockinger2020_vsh_quarantined": stockinger["diagnostic_explosion_energy"]["vsh_quarantined"],
    "stockinger2020_blocking_reasons": stockinger["blockers"],
    "sukhbold2016_candidate_status": sukhbold["status"],
    "sukhbold2016_model_count": sukhbold["z96_grid"]["model_count"],
    "sukhbold2016_mass_range_msun": [
        min(sukhbold["z96_grid"]["zams_mass_msun"]),
        max(sukhbold["z96_grid"]["zams_mass_msun"]),
    ],
    "sukhbold2016_energy_range_erg": [
        sukhbold["energy_and_fallback"]["final_kinetic_energy_erg_minimum"],
        sukhbold["energy_and_fallback"]["final_kinetic_energy_erg_maximum"],
    ],
    "sukhbold2016_mass_budget_within_review_bound": sukhbold["mass_budget_review"]["within_review_bound"],
    "sukhbold2016_exact_mass_closure_claimed": sukhbold["mass_budget_review"]["exact_mass_closure_claimed"],
    "sukhbold2016_third_party_redistribution_permission_verified": sukhbold["source_identity"]["third_party_redistribution_permission_verified"],
    "sukhbold2016_blocking_reasons": sukhbold["blockers"],
    "sukhbold2016_projection_status": sukhbold_projection["status"],
    "sukhbold2016_projection_record_count": sukhbold_projection["record_count"],
    "sukhbold2016_projection_high_mass_record_count": sukhbold_projection["high_mass_record_count"],
    "sukhbold2016_projection_high_mass_component_counts": sukhbold_projection["high_mass_record_count_by_source_component"],
    "sukhbold2016_projection_source_nulls_preserved": sukhbold_projection["source_nulls_preserved"],
    "sukhbold2016_projection_blocking_reasons": sukhbold_projection["blockers"],
    "limongi2024_transition_fate_status": limongi_transition["status"],
    "limongi2024_transition_fate_model_count": limongi_transition["machine_readable_tp_table"]["model_count"],
    "limongi2024_transition_fate_runtime_edge_msun": limongi_transition["project_transition_policy"]["unresolved_runtime_edge_interval_msun"],
    "limongi2024_transition_fate_interpolation_allowed": limongi_transition["project_transition_policy"]["continuous_fate_interpolation_allowed"],
    "limongi2024_transition_fate_blocking_reasons": limongi_transition["blockers"],
    "roberti2024_ultralowz_candidate_status": roberti_ultralowz["status"],
    "roberti2024_ultralowz_model_count": roberti_ultralowz["source_grid"]["model_count"],
    "roberti2024_ultralowz_metallicity_nodes": roberti_ultralowz["source_grid"]["metallicity_mass_fraction"],
    "roberti2024_ultralowz_official_mrt_model_count": roberti_ultralowz["yield_model_inventory"]["official_mrt_model_count"],
    "roberti2024_ultralowz_source_only_missing_models": roberti_ultralowz["yield_model_inventory"]["source_only_models_missing_from_official_mrt"],
    "roberti2024_ultralowz_mass_budget_outliers": roberti_ultralowz["mass_budget_review"]["outlier_models"],
    "roberti2024_ultralowz_blocking_reasons": roberti_ultralowz["blockers"],
    "heger_woosley2010_popiii_candidate_status": heger_woosley_popiii["status"],
    "heger_woosley2010_popiii_record_count": heger_woosley_popiii["source_grid"]["record_count"],
    "heger_woosley2010_popiii_coordinate_count": heger_woosley_popiii["source_grid"]["coordinate_count"],
    "heger_woosley2010_popiii_mass_hull_msun": [
        heger_woosley_popiii["source_grid"]["zams_mass_msun_minimum"],
        heger_woosley_popiii["source_grid"]["zams_mass_msun_maximum"],
    ],
    "heger_woosley2010_popiii_blocking_reasons": heger_woosley_popiii["blockers"],
    "baseline_metallicity_demand_status": baseline_metallicity_demand["status"],
    "baseline_metallicity_demand_star_count": baseline_metallicity_demand["stellar_population"]["star_count"],
    "baseline_metallicity_demand_observed_range": [
        baseline_metallicity_demand["stellar_population"]["birth_metallicity_mass_fraction"]["minimum"],
        baseline_metallicity_demand["stellar_population"]["birth_metallicity_mass_fraction"]["maximum"],
    ],
    "baseline_metallicity_demand_fraction_below_candidate_domain": baseline_metallicity_demand["candidate_domain_comparison"]["fraction_below_lowest_positive_full_grid_candidate"],
    "baseline_metallicity_demand_blocking_reasons": baseline_metallicity_demand["blockers"],
    "candidate_grid_coverage_status": coverage["status"],
    "candidate_grid_mass_gaps_msun": {
        channel: value["uncovered_runtime_edge_intervals_msun"]
        for channel, value in coverage["mass_coverage"].items()
    },
    "candidate_grid_source_hull_gaps_msun": {
        channel: value["uncovered_runtime_source_hull_intervals_msun"]
        for channel, value in coverage["mass_coverage"].items()
    },
    "candidate_grid_exact_common_metallicity_nodes": coverage["metallicity_coverage"]["exact_common_node_count"],
    "candidate_grid_roberti_baseline_inside_sparse_z_interval": coverage["ultralowz_sparse_candidate"]["baseline_values_inside_zero_to_first_positive_coordinate_interval"],
    "candidate_grid_roberti_production_coverage_approved": coverage["ultralowz_sparse_candidate"]["production_coverage_approved"],
    "candidate_grid_popiii_runtime_mass_gap_msun": coverage["popiii_mass_grid_candidate"]["runtime_channel_3_uncovered_edge_intervals_msun"],
    "candidate_grid_popiii_production_coverage_approved": coverage["popiii_mass_grid_candidate"]["production_coverage_approved"],
    "limongi_phase_mass_history_status": limongi_mass_history["status"],
    "limongi_phase_mass_history_monotonic_violations": limongi_mass_history["mass_history"]["monotonic_mass_violation_count"],
    "limongi_wind_mass_closure_max_abs_msun": limongi_mass_history["terminal_integrated_wind_closure"]["maximum_absolute_residual"]["absolute_residual_msun"],
    "blocking_reason": "No cited, checksummed, approved physical full-grid asset exists for channels 1--3; staged candidates remain review-only.",
}
output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print("G2_PREFLIGHT_BLOCKED")
PY
