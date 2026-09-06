# DUST-8 bundle plan: CUDA ledger to FP64 RAMSES trial handoff

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base commit: `490ab27` (DUST-7 fourth-species CUDA boundary)
- Bundle owner: native SNRT RT/feedback/dust integration
- Plan status: implemented; self end audit recorded in the implementation evidence

## Scientific and engineering objective

Complete the next handoff in the high-level RT/dust path: the separate
DUST-7 CUDA H/He/dust photon ledgers must reach the FP64 RAMSES trial path
without reconstructing the H/He partition from an assigned-total value.  A
failed trial must remain non-committing, and a successful trial must retain the
existing H/He chemistry behavior and photon conservation contract.

This bundle is an integration boundary, not a dust-physics approval.  The
project goal remains a production- and publication-ready RAMSES high-level
RT/stellar-AGN-feedback/dust implementation with explicit physical source,
opacity, energy, momentum, and transaction contracts.

## In scope

1. Add a separate prepared-cell transport entry point that invokes
   `snrt_cuda_multigroup_rt_step_species_dust_c` for each existing transport
   substep.  Accumulate raw, returned, assigned-total, H/He-species, and dust
   group ledgers across substeps, together with the depleted H/He inventory.
   Preserve the existing three-species prepared entry point and symbol.
2. Add one small FP64 receiver validator for the two identities
   `raw = assigned_total + returned` and
   `assigned_total = sum(H/He) + dust`, with finite/non-negative/shape checks
   and an FP32-rounding tolerance derived from `epsilon(real32)`.
3. Route the native RAMSES driver through the new ledger-returning entry point
   with an explicitly initialized zero dust optical-depth array.  Feed the
   returned H/He-species ledger directly to native thermochemistry; never feed
   DUST-7 assigned-total output into the old host partition routine.  Retain
   dust/raw/returned ledgers as trial diagnostics and validate them, but do not
   invent a persistent dust thermal, momentum, or abundance field.
4. Add a focused native validator smoke and extend the existing native bundle
   gate to check the new receiver boundary and source-level wiring.  Re-run
   the actual A10 CUDA smoke, GNU/Intel native checks where available, and the
   SNRT production link.  Do not launch a live RAMSES evolution in this bundle.
5. Record source hashes, compiler/device facts, exact tolerance, and the
   intentional zero-dust reference-only status in implementation evidence.

## Explicitly deferred

- Loading a P4 JSON opacity sidecar in Fortran or selecting a cell-level
  dust-to-metal prescription.  Those inputs are not currently present in the
  validated native RAMSES contract, so nonzero live dust would be an
  ungrounded physical activation.
- Persistent dust abundance/temperature state, dust heating into `uold`,
  momentum/radiation-pressure coupling, IR re-emission, scattering, and
  restart/migration of those states.
- Changing the existing transaction snapshot or convergence state solely to
  carry trial diagnostics that are not committed physical state.
- Reworking AMR topology, MPI exchange, source spectra, stellar population
  models, or unrelated legacy feedback code.

## Exit criteria

- The new transport path uses the DUST-7 ABI and preserves the old prepared
  ABI; substep sums and H/He inventory depletion are explicit.
- FP64 receiver validation rejects malformed, negative, non-finite, and
  ledger-inconsistent arrays before chemistry/commit, and accepts the tested
  FP32 rounding envelope.
- The native driver consumes direct H/He CUDA ledgers, reports/validates dust
  trial ledgers, and remains physically zero-dust until a validated opacity
  source-bound gate is available.
- Existing CUDA legacy/zero-dust regression, transaction rollback, native
  bundle link, and diff checks pass.  Any full-link failure caused by unrelated
  pre-existing worktree changes is isolated and reported rather than folded
  into this bundle.
