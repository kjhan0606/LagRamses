# Self end audit: DUST-8 CUDA ledger to FP64 RAMSES trial handoff

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Auditor: Codex implementer, independent final source/evidence review
- Date: 2026-09-06
- Verdict: **CONDITIONAL PASS — bundle scope satisfied**

## Review findings

The DUST-7 fourth-species CUDA ABI is now load-bearing in the native prepared
driver path through a distinct Fortran entry point.  The legacy prepared ABI
remains a separate wrapper, and the production link contains both symbols.
The new path preserves the existing AMR exchange/coarse-flux transaction
boundary and accumulates all DUST-7 ledgers over the same transport substeps.

The driver uses the direct H/He ledger rather than repartitioning the CUDA
assigned total on the host.  It rebuilds the FP32 total optical depth from the
component arrays before the call, avoiding a separately rounded FP64/FP32
component mismatch.  The FP64 receiver validator checks shape, finite state,
non-negativity, and both ledger identities with an explicit `64*epsilon(real32)`
tolerance.  It has GNU and Intel smoke evidence; the A10 CUDA evidence includes
legacy zero-dust bitwise equivalence and a maximum H/He ledger difference of
`3.278255E-07`.

No persistent dust field was added to the transaction snapshot or RAMSES
hydro state.  The named `ZERO_SCAFFOLD` mode is logged and the code explicitly
fails closed if a future dust-only positive absorption reaches the gas
receiver.  Thus the current native status is an integration/reference path,
not a nonzero dust physics activation.

## Scope and residual risks

The consolidated native gate passed production link, receiver smoke,
thermochemistry, spectral, two-rank transaction/MPI, CUDA A10, production
negative, and diff checks.  This is sufficient for DUST-8's declared native
handoff objective.  It does not qualify the deferred P4 opacity sidecar,
cell-level dust-to-metal mapping, dust heat/momentum/IR receivers, persistent
state, restart/migration, or live nonzero-dust RAMSES evolution.  Those remain
explicit follow-on gates and are not silently promoted by this bundle.

## Final disposition

**PASS within DUST-8 scope; conditional for production/publication.** The next
required dust physics bundle must supply and validate a source-bound native
opacity/cell-abundance contract and a persistent dust receiver before changing
`ZERO_SCAFFOLD` or claiming live dust feedback.
