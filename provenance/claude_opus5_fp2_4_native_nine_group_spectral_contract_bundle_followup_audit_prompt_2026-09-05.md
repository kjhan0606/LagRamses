# Claude Opus 5 follow-up audit request — F-P2.4 native nine-group remediation

You are the single bundled follow-up auditor for F-P2.4 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.

This is a read-only audit. Do not edit files, run jobs, launch RAMSES, use the
network, or invoke Python merely to reproduce an already recorded native
result. Inspect the actual Fortran/CUDA sources and the recorded native
evidence. The project purpose is production/publication-ready high-level
RAMSES physics for RT, stellar/AGN feedback, and dust; this particular bundle
is only the native SNRT nine-group/source-contract wiring boundary.

The previous bundle-end audit is recorded at:

`provenance/claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_2026-09-05.md`

It returned `CONDITIONAL PASS` with these closure conditions:

1. C1/F3: bound cross sections and absorber-weighted excess energies from
   above, not only from below.
2. C2/F4: execute loader rejection paths and a real version-4 checkpoint
   write/read round trip containing a nonzero state payload.
3. C3: document the quantified emission-mean versus absorber-weighted heating
   residual and the absent native secondary-ionization/recombination channels.

The remediation claims are recorded at:

`provenance/fp2_4_native_nine_group_spectral_contract_bundle_implementation_evidence_2026-09-05.md`

Inspect at least these files:

- `patch/lagRamses/snrt_spectral_contract.f90`
- `patch/lagRamses/snrt_spectral_contract_smoke.f90`
- `patch/lagRamses/snrt_spectral_contract_loader_smoke.f90`
- `patch/lagRamses/snrt_checkpoint_smoke.f90`
- `patch/lagRamses/snrt_state.f90`
- `simulation/snrt/tests/run_snrt_native_spectral_contract.sh`
- `simulation/snrt/config/snrt_group_contract_reference_control_v1.nml`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `simulation/snrt/P4_SOURCE_LEDGER.md`
- `simulation/snrt/P4_TRANSPORT_CONSERVATION_VALIDATION.md`

Verify the following narrowly:

- C1 upper bounds are applied to every supported H I/He I/He II group and do
  not reject the checked-in reference values or admit a plausible decimal
  slip.
- Loader tests exercise actual `snrt_spectral_contract_load` or
  `load_from_environment` failure paths for unset/missing/malformed/version,
  malformed identity, edge mismatch, unknown fraction semantics, candidate
  status, intrinsic-fraction runtime blocking, and reference opt-in.
- The checkpoint smoke calls the real `snrt_state_checkpoint_write` and
  `_read`, uses a nonzero nine-group payload, rejects a candidate identity
  before state mutation, and then restores all payload records under the
  reference contract. Check the serialized field order against the writer and
  reader, including fraction semantics.
- The documentation names and quantifies F1/F2 as later science gates and
  does not claim He, secondary ionization, recombination, physical SED
  approval, live RT+feedback, or production science validation.
- `fraction_semantics` is bound to the checkpoint identity; resolved-domain
  runtime cannot silently accept an intrinsic fraction. A reference control
  cannot run without explicit opt-in.
- No new array-ordering, dimension, unit, source-energy, stale-four-group,
  or CUDA ABI defect was introduced.

Use the recorded GNU Fortran and `mpiifx` runner outputs and the latest
`make -C bin SNRT=1 USE_CUDA=1 ramses` result as evidence. Do not treat the
reference-control contract as a physical approval, and do not require the
later HDF5 restart hooks, He chemistry, secondary ionization, dust, transport
accuracy, or 40--120 M☉ yield seam to be solved in this bundle.

Return exactly one top-level verdict:

`PASS`, `CONDITIONAL PASS`, or `FAIL`.

If conditional or failed, distinguish a blocker for closing F-P2.4 from a
later G3/G4/G5 science or coupling task, and name the smallest concrete next
action. If PASS, state that C1/C2/C3 are closed and list any advisory items
that remain explicitly outside this bundle.
