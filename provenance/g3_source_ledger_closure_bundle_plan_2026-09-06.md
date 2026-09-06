# G3 source-ledger closure bundle plan

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Work location: `/gpfs`
Status: operator-approved implementation bundle.

## Boundary

This bundle closes the source-to-photon-ledger merge boundary needed by G3. It
does not select or approve a stellar SED, AGN obscuration model, SNIa DTD,
PISN prescription, physical yield table, or dust mixture. Those decisions stay
candidate-only and fail closed.

## Deliverables

1. Require component ledgers to use identical group edges and internally closed
   CSV/metadata photon totals.
2. Reject globally colliding source IDs before any output is written, and emit a
   deterministic mixed source-kind and source-ID policy.
3. Carry source epoch metadata through the merge. Reject different epochs and
   mark missing epoch metadata as not production-eligible.
4. Preserve the aggregate photon, spectral, and dust-binding provenance. A
   component-only dust sidecar must not be accepted for a mixed source ledger.
5. Run one consolidated CPU/JAX gate covering the merge contract and the
   existing source-SED/dust binding regression.

## Acceptance

The bundle passes only when the normal STAR+AGN merge, aggregate spectral
closure, duplicate-ID rejection, mismatched-edge rejection, mismatched-epoch
rejection, and source-bound dust checks all pass. No runtime source activation
or large RAMSES job is part of this bundle.

## Deferred work

Physical SED selection, escape/obscuration calibration, full SNIa/DTD runtime,
PISN eligibility, approved dust opacity/depletion, live RAMSES coupling, and
the CPU/CUDA build-identity separation remain later gates.
