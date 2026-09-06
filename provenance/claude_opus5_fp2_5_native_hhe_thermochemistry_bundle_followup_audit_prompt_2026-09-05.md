# Claude Opus 5 follow-up audit — F-P2.5 first-audit remediation

You are the single bundled Claude Opus 5 end auditor for F-P2.5 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.

This is a read-only follow-up algorithm, implementation, physics-boundary, and
provenance audit. Do not edit files, launch simulations/jobs, invoke Python or
JAX, use external network tools, or create build artifacts. Inspect the actual
current native Fortran/CUDA source and the recorded native evidence.

The first F-P2.5 audit is recorded in:

- `provenance/claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_end_audit_2026-09-05.md`
- `provenance/fp2_5_native_hhe_thermochemistry_bundle_audit_remediation_2026-09-05.md`

The project goal is production-ready and publication-ready high-level RAMSES
radiative transfer plus stellar/AGN feedback and dust physics. F-P2.5 is only
the native H/He thermochemistry bundle; it is not a physical AGN/stellar SED
approval, live production run, or publication validation.

Read these current files first:

- `provenance/fp2_5_native_hhe_thermochemistry_bundle_plan_2026-09-05.md`
- `provenance/fp2_5_native_hhe_thermochemistry_bundle_implementation_evidence_2026-09-05.md`
- `provenance/claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_end_audit_2026-09-05.md`
- `provenance/fp2_5_native_hhe_thermochemistry_bundle_audit_remediation_2026-09-05.md`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `patch/lagRamses/snrt_cuda_kernels.cu`
- `patch/lagRamses/snrt_cuda_multigroup_interface.f90`
- `patch/lagRamses/snrt_cuda_multigroup_smoke.f90`
- `patch/lagRamses/snrt_transport_step.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_thermochemistry.f90`
- `patch/lagRamses/snrt_thermochemistry_smoke.f90`
- `patch/lagRamses/snrt_state.f90`
- `patch/lagRamses/snrt_checkpoint_smoke.f90`
- `bin/Makefile`
- `simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh`
- `simulation/snrt/tests/run_snrt_native_thermochemistry.sh`
- `simulation/snrt/tests/run_snrt_native_spectral_contract.sh`

Verify specifically whether the first-audit findings are genuinely closed:

1. HIGH-1: Is the FP32 inventory tolerance/guard safe and explicit? Does the
   optional unassigned ledger prevent a tolerance-sized excess from being
   silently assigned or silently lost? Is the live production path protected
   from an inventory overrun after CUDA reduction?
2. HIGH-2: Trace the new species-aware CUDA ABI and memory layout. Confirm that
   `(leaf,group,species)` Fortran storage is read as `(species,group,leaf)` in
   CUDA, that group/species opacity masks are applied, that one shared
   H I/He I/He II inventory is consumed across all groups and all substeps, and
   that the restored photon state and returned group absorption remain
   conservative. Check the driver really uses this ABI and the old scalar
   compatibility entry point is not on the production path.
3. MODERATE-3/4/5: Confirm raw FS2010 channel closure is checked before the
   loaded flag; the hot-temperature He II radiative+dielectronic test is
   meaningful; and saturated-H-II unavailable-target energy is routed to heat.
4. MODERATE-6/7: Confirm authoritative H/He fractions are double precision,
   checkpoint version 6 rejects old/mismatched payloads before mutation, and
   both spectral and FS2010 identities are actually serialized and checked.
5. LOW-10: Confirm the Makefile prerequisites are sufficient for parallel
   production builds. Note any remaining correctness issue in the remediation,
   including the `unassigned` treatment or arbitrary group-order policy.

Separate direct source verification from evidence that is merely recorded in
the evidence document. Do not turn later declared tasks (global implicit
opacity/chemistry fixed point, cooling, dust, physical SED admission, HDF5
restart, convergence, 40--120 M_sun yield seam) into F-P2.5 blockers unless the
current remediation incorrectly claims them complete.

Return exactly one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or
`FAIL`. Then give severity-ranked findings with file/line references, an
explicit disposition for each F-P2.5 acceptance gate, and a clear statement of
whether the bundle may close or must receive another remediation. Do not edit
the plan or source files.
