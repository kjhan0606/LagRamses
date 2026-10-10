# High-level RT/feedback/dust simulation-readiness closeout

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)

## Decision

The repository now has two explicitly admitted executable smoke profiles:

| Profile | Decision | Evidence |
|---|---|---|
| SNRT/CUDA reference-control | `READY_FOR_BOUNDED_QUALIFICATION` | D4 Slurm job `333272`, two MPI ranks/two A10 GPUs, startup-contract and receiver-rollback PASS |
| CPU legacy comparison | `READY_FOR_BOUNDED_COMPARISON` | two-rank `nstepmax=2` smoke, `feedback_mode='legacy'`, finite-state diagnostics and `Run completed` |
| Physical dusty production | `BLOCKED` | production manifest audit exit 2; live driver is still `ZERO_SCAFFOLD` |

The first two profiles are qualification/comparison runs, not a scientific
production release. The physical profile cannot be enabled by selecting a
different binary or namelist.

## Verification completed

The final bounded checks completed with PASS:

```text
G4_DUST_CLOSURE_PASS tests=3 mapping=explicit thermal=one_pass backend=cpu
SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
SNRT_NATIVE_DUST_MAPPING_RECEIVER_RUN_PASS
SNRT_NATIVE_DUST_TRANSACTION_ALL_OK
FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK
G2_SOURCE_SELECTION_GATE_TEST_OK
PARKED_AGB_SOURCE_CHECKS_RAN
G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK
FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK
```

The reference-control preflight used the recorded SNRT/CUDA binary
`simulation/snrt/build/g5_startup_contract_admission/ramses3d` with SHA-256
`e30a938397781ee02eb060222937007f659bed4a0e5d268c1cee125901e019f0` and the
effective initialized-RAMSES namelist. The CPU comparison used
`bin/ramses_final3d` with SHA-256
`5eaef4bdc321abbe2494a481d229b7170221c0be623d572f722b3f326253bff2`.

## Physical blockers that remain intentional

- The production manifest still lacks approved, checksummed registrations for
  the physical stellar yield table, stellar SED, AGN SED, and dust model.
- `p4_snapshot_hdf5` remains a not-migrated, hash-pending reference payload;
  the existing small P4-derived inputs are validation/control artifacts, not a
  substitute for that physical input contract.
- The thermal atlas has an approved physics provenance record but its license
  status is still pending explicit confirmation.
- The live RAMSES driver has no persistent dust field or nonzero dust optical
  depth, and no live dust thermal, momentum, scattering, or IR receiver.
- The current checkout is dirty because unrelated migration/high-mass work is
  in progress. A physical release requires a clean named commit and a rebuild
  of the selected binary from that commit.
- The 40--120 M☉ source seam and AGB source choice remain review-only; parked
  source readers do not constitute physical approval.

These are physical-release conditions, not reasons to repeat the already
passing native contract smokes. The next implementation work is therefore
physical asset admission and live receiver activation, followed by the G5
hydro/restart/AMR qualification matrix.

## Reproduction entry points

```text
simulation/snrt/runs/simulation_ready_preflight.sh
simulation/snrt/tests/run_g4_dust_closure.sh
simulation/snrt/tests/run_snrt_bundle_gate.sh
```

The latest record-keeping repair was committed and pushed as `7330087`.

After that record update, the integrated gate was rerun from the current
workspace and returned `SNRT_BUNDLE_GATE_PASS` at commit `7330087`.  It
reported production-link PASS in 178.207 s, AGN partition PASS, five native
symbols, dust-ledger receiver PASS, four thermochemistry negative cases, ten
spectral negative cases, eleven transaction rollback/configuration negative
cases with two-rank MPI coverage, A10 CUDA PASS, production-negative PASS,
and diff/conservation checks PASS.  The temporary linked test binary had
SHA-256
`684dec40632e94058f5fd35028aeda20985418eb37a365bd33303f31f7f6ada0` and
was removed by the gate after the run.

## Current G0 audit scope correction — 2026-10-10

The current checkout is `/gpfs/kjhan/chimes-cvode-gpu-20261001` with origin
`kjhan0606/LagRamses`. The external-assets registry is a project-wide
migration/comparison inventory, not a list of assets all selected into one
production run. The G0 auditor previously required production license and
provenance metadata, and rejected blocked/comparison statuses, for every
registry entry. It now applies those admission checks only to IDs in
`required_production_assets`; an absent selected ID reports
`required_asset_not_registered` without a redundant metadata-missing error.
Registry-wide duplicate-ID and absolute-path integrity checks remain.

Verification: `simulation/snrt/tests/production_manifest.py` passes, including
an unselected `not_migrated` legacy comparison asset (must not block), and
selected legacy, missing metadata, missing registry entry and blocked-status
negative cases (must fail closed). `git diff --check` is clean. A read-only
audit of the current manifest still returns `production_gate_pass=false` with
12 concrete blockers: pending thermal-atlas license; missing metadata for the
selected executable, CUDA smoke binary and P4 input; the manifest's recorded
dirty-source status; four unregistered physical stellar/AGN/dust assets; and
the P4 input's missing path/checksum plus `not_migrated` status. These are
not waived by narrowing the audit, and no external registry, physical asset,
license record or source table was changed.

The 2026-09-06 statement that the live driver had no dust state/receivers is
historical and superseded by the integrated relative-motion and effective
population evidence in `dust_remaining_implementation_2026-09-09.md` and
`effective_population_execution_2026-09-11.md`. Those bounded
noncosmological CPU/OpenMP results do not convert this G0 manifest to PASS or
qualify the cosmological production profile.
