# Claude Opus 5 end-of-bundle audit request — F-P2.2 closure integrity

You are the sole active end-of-bundle auditor under the current project
cadence. Perform a read-only scientific and code-architecture audit of the
implemented F-P2.2 bundle in `/gpfs/kjhan/LRD_JWST`
(`kjhan0606/LagRamses`). Do not edit files, run jobs, build native code, or
launch simulations. Use only file reads, search, and directory listing. Do
not reopen generic old AMR/HDF5/ksection/MPI work outside this bundle.

## Project purpose and decision boundary

The project ultimately targets publication-ready high-level hydro for
LRD/JWST: physically defensible RT, stellar/AGN feedback, and dust coupled to
a future RAMSES production path. F-P2.2 is an engineering integrity and
canonical-coverage bundle for the single explicit-source spectral/dust
boundary. It is not approval of an astrophysical SED, escape/obscuration
model, mixed STAR+AGN source, live RT–hydro feedback, production run, or
publication result.

## Files to inspect

Read the following directly and inspect the data flow rather than trusting the
evidence document:

- `simulation/snrt/snrt_core/provenance.py`
- `simulation/snrt/snrt_core/sed.py`
- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/snrt_core/dust.py`
- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/build_draine_dust_opacity.py`
- `simulation/snrt/tools/p4_attach_pilot_sources.py`
- `simulation/snrt/tools/p4_run_transport_pilot.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tools/validate_agn_nine_group_ledger.py`
- `simulation/snrt/tests/source_sed_dust_closure.py`
- `simulation/snrt/tests/agn_nine_group_artifact.py`
- `simulation/snrt/tests/agn_nine_group_explicit_artifact.py`
- `simulation/snrt/tests/p4_dust_runner.py`
- `simulation/snrt/tests/p5_dust_runner.py`
- `simulation/snrt/P4_DUST_OPACITY.md`
- `simulation/snrt/P4_AGN_RATE_LEDGER.md`
- `simulation/snrt/data/agn_nine_group_external_assets.json`
- `provenance/fp2_2_closure_integrity_bundle_plan_2026-09-05.md`
- `provenance/fp2_2_closure_integrity_bundle_implementation_evidence_2026-09-05.md`
- `provenance/claude_opus5_fp2_1_source_closure_verification_bundle_end_audit_2026-09-05.md`
- `provenance/rt_architecture_audit.md`

Also inspect both canonical pilot and explicit ledger metadata, static-input
metadata, transport artifacts, and validation JSON. Check whether their
recorded hashes and working-tree attestation agree with the current files.

## Required questions

1. Does the code-manifest implementation enforce the complete dependency set
   for the explicit AGN photon closure and source-bound Draine dust closure,
   including exact role/path semantics and current hashes? Is any dependency
   still omitted?
2. Is the canonical payload self-hash unambiguous and fail-closed for modified
   closure arrays, metadata values, or manifest entries? Does the test mutate
   the payload meaningfully?
3. Is the dust status vocabulary fixed and schema-consistent, and are the
   metadata/payload/source-table/builder hashes truthfully propagated into P4
   and P5 artifacts while preserving v1 null-identity controls?
4. Does the explicit canonical nine-group validator actually exercise the
   repaired path rather than merely reusing pilot assumptions? Are the
   synthetic SED and its status clearly non-physical engineering controls?
5. Is the working-tree attestation useful and internally consistent, including
   the fact that the current development tree is dirty? Are canonical asset
   manifests and validator outputs synchronized?
6. Identify defects local to F-P2.2 separately from deferred mixed STAR+AGN
   admission, physical source selection, unimplemented dust physics, live RHD,
   and publication decisions.

## Verdict format

Return one overall verdict: `PASS`, `CONDITIONAL PASS`, or `FAIL`. Order
findings by `BLOCKER`, `MAJOR`, `MINOR`, and `NOTE`, with file/line anchors
where possible. For each finding state whether it blocks F-P2.2 or belongs to
a later science/production gate. State what is genuinely closed, what remains
deferred, and recommend one next coherent bundle.

Do not modify the repository and do not issue a verdict merely because the
plan or evidence says a requirement is met; inspect the algorithms and
artifact wiring yourself.
