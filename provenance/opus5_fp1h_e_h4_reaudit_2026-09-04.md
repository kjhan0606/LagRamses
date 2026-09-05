# Claude Opus 5 H4 remediation re-audit — F-P1H-E

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5 CLI, read-only
Prompt: `provenance/opus5_fp1h_e_h4_reaudit_prompt_2026-09-04.md`

## Verdict

**PASS.** Opus confirmed that both conditions from the first H4
`CONDITIONAL PASS` are closed. H4 remains review evidence only and does not
promote any source package or reconstruct any physical value.

## Closed conditions

1. `simulation/snrt/tools/fp1_limongi_phase_history.py` captures source-order
   phase ranks before canonical sorting and raises on an out-of-contract order.
   Its shared invariant path also rejects duplicate non-collapsed phases,
   missing PSN, non-positive duration/cumulative age, negative cumulative wind
   mass, and increasing total mass. The LC18 consumer propagates a controlled
   error report and exit status 2.
2. The G2 phase-history audit and the F-P1 LC18 cross-check use the same shared
   aggregator and attest its live SHA256. The checked-in G2 report attests
   both its audit-code hash and the shared-code hash, with a stored-vs-live
   test.

## Independent Opus findings

- The common accounting remains 108 models, 845 unique phase rows, and 19
  exact collapsed duplicate rows, with the 52/56, 48/4, 53/3, and 101/7
  checks retained.
- The normal failed-wind anomaly remains blocked for production, publication,
  canonical conversion, and runtime deposition. Parsed zero values are not
  promoted to physical-zero claims.
- Source precision, cross-source residual definitions, unresolved failed-model
  anomalies, review-only rights/publication semantics, and no-inference policy
  remain explicit. No wind correction, terminal energy/momentum inference, or
  cross-source reconciliation was introduced.
- The H2 boundary remains fail-closed: a synthetic approved source cannot pass
  the code-registered source-rights binding. The yield-auditor projection
  fixture is review-only, while converter package/projection branches use a
  test-only, scoped synthetic node-auditor seam; no production module is
  monkeypatched.
- Repository admission state remains zero physical nodes/canonical rows with
  false promotion flags and unresolved `[0.8,1.0)` and `[40,120] M_sun`
  seams.

## Driver verification on GPFS

The driver independently ran and passed:

- `python3 simulation/snrt/tests/g2_limongi_phase_mass_history.py`
- `python3 simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `.venv/bin/python simulation/snrt/tests/stellar_yield_asset.py`
- `.venv/bin/python simulation/snrt/tests/yield_converter.py`
- `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`
- `bash simulation/snrt/tests/run_g2_preflight.sh`, terminating at the
  expected `G2_PREFLIGHT_BLOCKED` state.

The live hashes verified by the driver are:

```text
fp1_limongi_phase_history.py          71ea5453659b15882e5d594c4e47c303d38e7a78e05596857450720a9f4c374c
audit_g2_limongi_phase_mass_history.py 6a2fb242638807f8518c62a838df81f635989348b933da53b07d5b29b122fba2
audit_fp1_lc18_failed_wind_crosscheck.py 7382132a6e6d9b974dff64fc11a8c316c631f311dc244c9b4516624d271ea0f8
```

## Non-blocking follow-up observations

Opus identified these as quality improvements, not H4 closure conditions:

- Expand G2 freshness comparison beyond the two code hashes to include the
  contract/source-adapter hashes and numeric mass-history block.
- Add stored-vs-live freshness checking for the LC18 report as well.
- Carry the existing caveat that source-attested intermediate-burning order
  is unavailable into the G2 report, and distinguish an enforced invariant
  from an observed measurement.
- Narrow the test-only converter seam and restore it with `try/finally`.

## Audit limitations

Opus did not edit files, launch jobs, build RAMSES, or run the Fortran/G2
preflight runners. It did independently run the four focused Python suites;
the driver independently executed the full runners and recorded the expected
fail-closed terminal state above.
