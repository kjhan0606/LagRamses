# Claude Opus 5 bundled end audit request — F-P2.5 native H/He thermochemistry

You are the single bundled Claude Opus 5 end auditor for F-P2.5 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.

This is a read-only algorithm, implementation, physics-boundary, and
provenance audit. Do not edit files, launch simulations/jobs, invoke Python or
JAX, use external network tools, or generate new build artifacts. Inspect the
actual native Fortran/CUDA source and the recorded native evidence. The
operator pre-approved this coherent bundle and requested continued progress.

The project goal is production-ready and publication-ready high-level RAMSES
radiation transport plus stellar/AGN feedback and dust physics. F-P2.5 is one
native H/He thermochemistry bundle toward that goal; it is not by itself a
physical SED approval, a live production run, or publication validation.

Read these first:

- `provenance/fp2_5_native_hhe_thermochemistry_bundle_plan_2026-09-05.md`
- `provenance/fp2_5_native_hhe_thermochemistry_bundle_implementation_evidence_2026-09-05.md`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `simulation/snrt/config/snrt_secondary_table_contract_v1.nml`
- `simulation/snrt/data/furlanetto_stoever_2010/README.md`
- `simulation/snrt/data/furlanetto_stoever_2010/TABLE_MANIFEST.sha256`
- `patch/lagRamses/snrt_thermochemistry.f90`
- `patch/lagRamses/snrt_thermochemistry_smoke.f90`
- `patch/lagRamses/snrt_thermochemistry_loader_smoke.f90`
- `patch/lagRamses/snrt_nlte_coupling.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_transport_step.f90`
- `patch/lagRamses/snrt_state.f90`
- `patch/lagRamses/snrt_checkpoint_smoke.f90`
- `bin/Makefile`
- `simulation/snrt/tests/run_snrt_native_thermochemistry.sh`
- `simulation/snrt/tests/run_snrt_native_spectral_contract.sh`

Audit the current worktree, not only the historical F-P2.4 snapshot. Treat
unrelated pre-existing dirty-tree changes as out of scope unless this bundle
overwrites or depends on them incorrectly.

Assess at least the following:

1. FS2010 contract and interpolation: Does the native loader genuinely bind
   the 14 vendored tables, dimensions, source identity, manifest identity,
   shared energy grid, finite/non-negative values, bilinear interpolation, and
   low/high energy semantics? Does it normalize the same five channels as the
   reference without hiding an energy defect?
2. Primordial opacity and absorption: Are H I, He I, and He II densities and
   group optical depths in consistent code/physical units? Does the driver pass
   the total absorber budget to the CUDA transport and then partition returned
   absorption using the same start-of-step species inventories, with no H-only
   cap or double counting across groups/substeps?
3. Photoelectron physics: Are primary ionizations, FS2010 secondary H I/He I/
   He II ionizations, unavailable-target energy returned to heat, excitation,
   and the residual ledger physically and numerically closed? Check threshold
   energies and the species-state transitions, including He II -> He III.
4. Recombination: Are H II, He II, and He III case-B rates implemented with
   the stated formulas, including dielectronic He II and exactly
   `alpha_HeIII,B(T)=2 alpha_HII,B(T/4)`? Does the bounded electron-density
   closure remain non-negative and number-consistent? Identify any limitation
   such as split start-of-step opacity versus the later global implicit fixed
   point, but do not turn a declared later gate into a false F-P2.5 blocker.
5. RAMSES coupling: Is temperature extracted from RAMSES total energy with the
   correct kinetic subtraction and unit factors? Is only gas heating added to
   total energy, while ionization potential, excitation/line, recombination,
   background, and metal cooling remain separate? Does a failed chemistry
   partition leave state and energy unchanged? Is the table contract fail-closed
   before source/transport mutation?
6. State/restart: Does checkpoint version 5 serialize and validate all new H/He
   fractions, enforce H I/H II mirror consistency, reject old or mismatched
   payloads before mutation, and preserve nonuniform intensity/state values?
7. Build/evidence: Do the recorded GNU and `mpiifx` native smoke results and
   the `make -C bin SNRT=1 USE_CUDA=1 ramses` link support the stated engineering
   gates? Separate direct evidence from what still needs a live CUDA/RAMSES
   evolution, physical SED admission, HDF5 restart integration, global
   convergence, dust, and the 40--120 M_sun yield seam.
8. Look for algorithmic or wiring defects that Python-only tests would miss:
   Fortran array ordering, optional/allocatable state use, source/transport
   call ordering, numerical overflow/underflow, stale H-only assumptions,
   failure cleanup, checkpoint record order, and Makefile module order.

Return exactly one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or
`FAIL`. Follow it with severity-ranked findings and file/line references,
evidence supporting each finding, and an explicit disposition for every F-P2.5
acceptance gate. Distinguish a concrete blocker for this native bundle from a
later high-level task. Do not grant physical AGN/stellar SED approval, a live
RT+feedback result, or publication science validation.
