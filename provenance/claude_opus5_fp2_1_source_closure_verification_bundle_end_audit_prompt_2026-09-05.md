# Claude Opus 5 end-of-bundle audit request — F-P2.1 source-closure verification

You are the sole active end-of-bundle auditor under the current project
cadence. Perform a read-only scientific and code-architecture audit of the
implemented F-P2.1 bundle in `/gpfs/kjhan/LRD_JWST`
(`kjhan0606/LagRamses`). Do not edit files, run jobs, build native code, or
launch simulations. Use only file reads, search, and directory listing. Do
not reopen generic old AMR/HDF5/ksection/MPI work outside this bundle.

## Project purpose and decision boundary

The project ultimately targets publication-ready high-level hydro for
LRD/JWST: physically defensible RT, stellar/AGN feedback, and dust coupled to
a future RAMSES production path. F-P2.1 is only a repair and verification
boundary for a single explicit-source spectral/dust candidate path. It must
not be treated as approval of a stellar population, AGN SED,
escape/obscuration model, mixed STAR+AGN source, live RT–hydro feedback,
production run, or publication result.

## Current bundle to inspect

Read these implementation, tests, and records directly:

- `simulation/snrt/snrt_core/sed.py`
- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/snrt_core/dust.py`
- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/build_draine_dust_opacity.py`
- `simulation/snrt/tools/p4_run_transport_pilot.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tools/validate_agn_nine_group_ledger.py`
- `simulation/snrt/tests/source_sed_dust_closure.py`
- `simulation/snrt/tests/agn_nine_group_artifact.py`
- `simulation/snrt/tests/p0_sed_closure.py`
- `simulation/snrt/tests/stellar_photon_ledger.py`
- `simulation/snrt/tests/dust_opacity.py`
- `simulation/snrt/tests/draine_dust_opacity.py`
- `simulation/snrt/tests/agn_photon_ledger.py`
- `simulation/snrt/tests/p4_dust_runner.py`
- `simulation/snrt/tests/p5_dust_runner.py`
- `simulation/snrt/tests/merge_photon_source_ledgers.py`
- `simulation/snrt/P4_DUST_OPACITY.md`
- `simulation/snrt/P4_AGN_RATE_LEDGER.md`
- `simulation/snrt/P4_SOURCE_LEDGER.md`
- `provenance/fp2_1_source_closure_verification_bundle_plan_2026-09-05.md`
- `provenance/fp2_1_source_closure_verification_bundle_implementation_evidence_2026-09-05.md`
- `provenance/claude_opus5_fp2_source_sed_dust_closure_bundle_end_audit_2026-09-05.md`
- `provenance/stellar_dust_asset_baseline.md`
- `provenance/rt_architecture_audit.md`

Also inspect the canonical AGN metadata, static-input metadata, external asset
manifest, transport artifact, and validation JSON referenced by the evidence.
Check the data flow and hashes rather than accepting the evidence document's
claims without verification.

## Required questions

1. Does the explicit AGN path now remain truthful and distinct from the
   parameterized Sazonov reference control, including support, normalization,
   escape-fraction, and intrinsic/escaped photon fields?
2. Is the source-weighted Draine v2 calculation mathematically correct, and
   is the independent offset-grid plus closed-form test meaningful rather than
   duplicating the production algorithm?
3. Are source SED, group-edge, Draine-table, and builder hashes all recorded,
   re-hashed, and propagated through P4/P5 before output creation? Are error
   paths fail closed while the null-identity v1 reference control remains
   intentionally usable?
4. Are the group boundaries, source support, units, empty groups, and source
   identity consistent across SED, H/He closure, dust closure, and runners?
5. Do the canonical regenerated artifacts and validator actually close the
   stated nine-group reference-control gate?
6. Identify any implementation defect in this bundle separately from the
   mixed STAR+AGN admission, astrophysical source selection, dust physics,
   live RHD, and publication decisions explicitly deferred to later gates.

## Verdict format

Return a concise but technically specific report with exactly one overall
verdict: `PASS`, `CONDITIONAL PASS`, or `FAIL`. Order findings by
`BLOCKER`, `MAJOR`, `MINOR`, and `NOTE`, with file and line anchors where
possible. For each finding state whether it blocks F-P2.1 or belongs to a
later science/production gate. State what is genuinely closed, what remains
deferred, and recommend one next coherent bundle.

Do not modify the repository and do not issue a verdict merely because the
plan or evidence says a requirement is met; inspect the implementation and
algorithmic wiring yourself.
