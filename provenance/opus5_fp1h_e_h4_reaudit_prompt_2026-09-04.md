# Claude Opus 5 remediation re-audit — F-P1H-E H4

Perform a read-only audit in `/gpfs/kjhan/LRD_JWST`. Do not edit files, run
jobs, build RAMSES, or launch a simulation. Inspect the H4 implementation and
the directly related regression fixtures.

This is a re-audit of H4 in the F-P1H-E fail-closed physical stellar-source
admission bundle for the production/publication-ready lagRamses high-level
hydro project: RT, stellar/AGN feedback, dust, and coupled source terms.
H4 is review evidence only; it must not approve or reconstruct physical
values. The repository must remain at zero physical nodes/canonical rows,
with false promotion/runtime/publication flags and unresolved `[0.8,1.0)` and
`[40,120] M_sun` seams.

The first H4 audit was CONDITIONAL PASS. Its two medium conditions were:

1. Restore the fail-closed source-order check removed by the shared phase
   history refactor.
2. Regenerate the checked-in G2 phase-history report and attest its shared
   implementation hash, with a stored-vs-live freshness test.

The driver also found and repaired two test-fixture interactions caused by the
new H2 source-rights gate: synthetic approved fixtures are now either
review-only when testing the independent yield-auditor projection path, or
use an explicitly scoped synthetic node-auditor seam when testing the
physical-package/converter branches. This must not weaken production code or
permit synthetic sources into physical admission.

Inspect at minimum:

- `simulation/snrt/tools/fp1_limongi_phase_history.py`
- `simulation/snrt/tools/audit_g2_limongi_phase_mass_history.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tools/audit_stellar_yield_asset.py`
- `simulation/snrt/tools/convert_yield_rows_to_canonical.py`
- `simulation/snrt/tests/g2_limongi_phase_mass_history.py`
- `simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tests/stellar_yield_asset.py`
- `simulation/snrt/tests/yield_converter.py`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `simulation/snrt/data/g2_limongi_phase_mass_history_audit.json`
- `simulation/snrt/config/g2_limongi_phase_mass_history_contract_v1.json`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`

Verify specifically:

1. The shared aggregator rejects out-of-contract source phase order,
   duplicate phases, missing PSN, non-positive duration/cumulative age,
   negative cumulative wind mass, and increasing total mass. The LC18
   cross-check still has a controlled non-zero failure path for malformed
   invariants. Its normal unresolved anomaly remains blocked for production,
   publication, conversion, and deposition, and parsed zero values are not
   interpreted as physical zero.
2. Both reports use the same shared implementation and attest its SHA256;
   the checked-in G2 report attests the live audit and shared-code hashes.
   The common aggregation preserves 108 models, 845 unique phase rows, 19
   exact collapsed duplicate rows, and the existing 52/56, 48/4, 53/3, and
   101/7 accounting.
3. Source precision, cross-source residuals, failed-model anomaly, explicit
   review-only semantics, and no-inference policy remain visible. No failed
   wind correction, terminal energy inference, momentum inference, or
   cross-source reconciliation is introduced.
4. The synthetic fixture changes preserve the real H2 fail-closed boundary:
   an approved synthetic source cannot pass source-rights binding, while the
   independent projection and package tests still exercise their intended
   rejection branches. No production path is monkeypatched by the code under
   audit.
5. The driver has independently passed:

   - `python3 simulation/snrt/tests/g2_limongi_phase_mass_history.py`
   - `python3 simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
   - `/gpfs/kjhan/LRD_JWST/simulation/snrt/.venv/bin/python simulation/snrt/tests/stellar_yield_asset.py`
   - `/gpfs/kjhan/LRD_JWST/simulation/snrt/.venv/bin/python simulation/snrt/tests/yield_converter.py`
   - `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`
   - `bash simulation/snrt/tests/run_g2_preflight.sh` with terminal
     `G2_PREFLIGHT_BLOCKED`.

Return severity-ranked findings with file/line evidence and exactly one
verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. Distinguish any remaining
low-severity observations from conditions that prevent H4 closure. Note
read-only limitations. Do not edit files or run jobs.
