# Fable plan audit — F-P2.6 native RT/chemistry transaction and fixed point

You are the primary plan auditor for the next large implementation bundle in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`. This is a read-only
plan audit. Do not edit files, run jobs, build code, invoke Python/JAX, use the
network, or launch RAMSES.

Read:

- `provenance/fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md`
- `provenance/fp2_5_native_hhe_thermochemistry_bundle_plan_2026-09-05.md`
- `provenance/claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_followup_audit_2026-09-05.md`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_transport_step.f90`
- `patch/lagRamses/snrt_thermochemistry.f90`
- `patch/lagRamses/snrt_state.f90`
- `patch/lagRamses/snrt_cuda_kernels.cu`
- `provenance/production_publication_readiness_plan.md`

## Project purpose

The project ultimately targets production- and publication-ready high-level
RAMSES radiative transfer, stellar/AGN feedback, and dust physics. The
proposed F-P2.6 bundle is only a native RT/chemistry correctness boundary. It
must not promote a physical AGN/stellar SED, yield table, dust model, live
hydro run, or publication result.

## Audit questions

1. Is the level-wide transaction boundary physically and computationally
   correct for a neighbor-coupled transport update? Does rollback cover every
   photon, species, and thermal state mutation, including partial CUDA/Fortran
   failures?
2. Is the proposed treatment of `unassigned_absorption_code` conservative and
   fail-closed? Does the plan avoid silently converting an unowned residual
   into gas heat or ionization?
3. Is the proposed fixed-point iteration well-defined? In particular, does
   each trial restart from the same incoming radiation, avoid double consuming
   photons, specify opacity/state relaxation and convergence norms, and commit
   only one converged result?
4. Is the scope feasible without accidentally claiming global MPI/AMR or live
   hydro closure? Are the native tests capable of detecting rollback,
   non-convergence, conservation, and exactly-once thermal updates?
5. Does this bundle materially advance the final RT/feedback/dust goal, and
   are the deferred physical SED, yield, dust, cooling, and G5 tasks clearly
   separated?

## Required output

Return one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or `FAIL`.
Then give severity-ranked findings with concrete plan edits. Distinguish a
plan blocker from a later science/coupling task. If the plan is feasible,
state the minimum conditions for implementation approval. Do not modify the
repository and do not audit unrelated AMR/HDF5/ksection work.
