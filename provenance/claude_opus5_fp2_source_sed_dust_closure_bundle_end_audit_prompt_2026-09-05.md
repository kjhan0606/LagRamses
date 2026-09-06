# Claude Opus 5 end-of-bundle audit request — F-P2 source SED/dust closure

You are the sole active end-of-bundle auditor under the current project
cadence. Perform a read-only scientific/code architecture audit of the
implemented F-P2 source-SED, photon-ledger, dust-opacity, and mixed-source
closure bundle in `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Do not edit
files, run jobs, build native code, or launch simulations. Use only file reads,
search, and directory listing. Do not ask for or perform a generic re-audit of
old AMR/HDF5/ksection/MPI work that is outside this bundle.

## Project purpose and decision boundary

The project ultimately targets publication-grade high-level hydro for LRD/JWST:
physically defensible RT, stellar/AGN feedback, and dust coupled to a future
RAMSES production path. This bundle is only a prerequisite spectral/provenance
boundary. It must not be mistaken for approval of a stellar population, AGN
SED, obscuration model, dusty LRD interpretation, live RT–hydro feedback, or
publication result. Existing pilot artifacts may remain reference controls.

## Bundle to inspect

Read the implementation, tests, plan, and evidence directly:

- `simulation/snrt/snrt_core/sed.py`
- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/snrt_core/dust.py`
- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/build_draine_dust_opacity.py`
- `simulation/snrt/tools/merge_photon_source_ledgers.py`
- `simulation/snrt/tools/p4_run_transport_pilot.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tests/source_sed_dust_closure.py`
- `simulation/snrt/tests/merge_photon_source_ledgers.py`
- `simulation/snrt/tests/p4_dust_runner.py`
- `simulation/snrt/tests/p5_dust_runner.py`
- `simulation/snrt/P4_DUST_OPACITY.md`
- `simulation/snrt/P4_AGN_RATE_LEDGER.md`
- `simulation/snrt/P4_SOURCE_LEDGER.md`
- `provenance/fp2_source_sed_dust_closure_bundle_plan_2026-09-05.md`
- `provenance/fp2_source_sed_dust_closure_bundle_implementation_evidence_2026-09-05.md`
- `provenance/stellar_dust_asset_baseline.md`
- `provenance/rt_architecture_audit.md`

Also inspect the canonical AGN metadata, static-input metadata, transport
artifact manifest, and validation JSON referenced by the evidence record. Check
whether the source/dust binding attributes are consistent with the runner
logic; do not assume that a passing unit test proves the wiring is correct.

## Required audit questions

1. Is the source SED contract mathematically and dimensionally coherent? Check
   the conversion from energy-fraction density to photon-number density, the
   `L_bol` normalization, support requirements, piecewise-linear integration,
   shared group-boundary ownership, empty groups, and path-free identity.
2. Does the explicit AGN ledger use the same source spectrum for photon groups,
   H/He cross sections, and photoelectron excess energy? Is the built-in pilot
   honestly isolated as a reference control?
3. Does the Draine/WD01 v2 sidecar compute the correct source-photon-weighted
   absorption closure and record enough hashes/units to prevent accidental
   source or group-edge substitution? Identify any gap between recorded and
   independently enforced provenance.
4. Do P4 and P5 fail closed on source/dust identity and edge mismatches, while
   preserving the intended null-identity v1 controls? Are output provenance
   attributes truthful?
5. Does the STAR+AGN merger correctly form aggregate photon-, absorber-, and
   excess-energy-weighted closures and prevent a component-only dust sidecar
   from being used as a mixture closure? Check duplicate IDs and output
   atomicity on rejection.
6. Are the tests and regenerated canonical artifacts sufficient for the stated
   candidate engineering gate? Distinguish implementation defects from
   scientifically deferred decisions.

## Verdict format

Return a concise but technically specific report with:

- `PASS`, `CONDITIONAL PASS`, or `FAIL` as the overall verdict;
- findings ordered by severity (`BLOCKER`, `MAJOR`, `MINOR`, `NOTE`), with file
  and line anchors where possible;
- whether each finding blocks this bundle or belongs to a later science gate;
- an explicit statement of what is genuinely closed and what remains deferred;
- one short recommendation for the next coherent implementation bundle.

Do not modify the repository. Do not issue a verdict merely because the plan
or evidence document says that a requirement is met; verify the algorithms and
data flow from the source.
