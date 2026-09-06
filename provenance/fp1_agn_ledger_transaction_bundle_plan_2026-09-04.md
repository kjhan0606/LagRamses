# F-P1.5 AGN ledger transaction bundle plan — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Status: **Claude Opus 5 CONDITIONAL PASS; conditions integrated; implementation
may proceed within the bounded scope below**.

## Why this bundle

The static AGN nine-group source ledger is already closed, but the high-level
feedback path still has an unresolved distinction between a pre-reset sink
state, a deferred feedback accumulator, and a committed SNRT source update.
The current reader rejects duplicate sink IDs but does not canonicalize an
identical restart duplicate keyed by `(nstep_coarse, sink_id)`.  The live SNRT
source loop can also deposit earlier photon groups before a later group fails;
because the accounted-mass marker advances only after the loop, a retry can
re-emit the already deposited groups.

This bundle addresses that bounded engineering gap.  It does not approve a
physical AGN SED, obscuration model, escape fraction, or hydro feedback
prescription, and it does not bypass the blocked G2 physical stellar-yield
gate.

## Scope

### A. Canonical coarse-state ledger boundary

1. Define and implement one canonical reader/audit path for the active
   `agn_coarse_state` JSONL records.
2. Key records by `(nstep_coarse, sink_id)` and deterministically retain one
   byte-equivalent copy of an identical restart duplicate.  Reject a
   same-key duplicate whose physical payload differs.  This bundle does not
   add a run UUID or dump counter; a rewind that produces a different payload
   therefore remains a fail-closed conflict and rewind safety stays open.
3. Preserve both raw `radiative_efficiency` and the
   `effective_radiative_efficiency` actually used for `Lbol`; downstream
   conversion must consume and emit the effective field explicitly.  The
   live SNRT driver currently uses raw `eps_sink`; that F10 convention
   mismatch remains an explicitly open follow-up, not a hidden PASS.
4. Validate the source-owned algebra without inventing missing physics:
   `inflow = min(Bondi, Eddington)`, unit conversion, and
   `Lbol = epsilon_effective * mdot_inflow * c^2`.  Validate that the record
   is explicitly a **pre-reset instantaneous state** and that accumulated
   fields are not silently treated as interval-integrated photon output.
   Use one declared Julian year (`365.25 d`) at every writer/parser boundary;
   accept an idle sink with `effective_radiative_efficiency=0` and zero
   luminosity while keeping raw efficiency in its own valid field.
5. Add a read-only audit tool that checks the compiled production source
   ordering: the coarse-state writer is called once before the AGN feedback
   accumulator reset, and the source record declares its review-only limits.

### B. All-or-nothing SNRT AGN source update

1. Add a transaction-level source routine that computes and validates every
   photon-group increment for one `(nstep_coarse, sink_id)` before mutating
   `snrt_intensity`.
2. If any group, target slot, range, finiteness, or overflow check fails, commit no group
   from that source and leave `accounted_mass` unchanged.
3. Advance `accounted_mass` only after the complete group transaction commits.
4. Key the in-memory accounting by the stable `idsink` identity, not only by
   the mutable sink-array position.  Keep the source loop serial and document
   that `accounted_mass`, the key map, and local intensity mutation are not
   OpenMP-threadprivate; assert the local-leaf/single-MPI-owner contract.
5. Resolve `SNRT_RT_ENABLE` once per process rather than allowing a mid-step
   environment toggle to strand an accounting marker.
6. Keep the current fail-closed runtime controls.  This is a correctness
   hardening change to a runtime-gated path, not authorization to turn
   `SNRT_RT_ENABLE` on in a production RAMSES run.

### C. Evidence

The bundle must include focused tests for:

- exact restart duplicate collapse;
- conflicting duplicate rejection;
- raw/effective efficiency field preservation;
- Bondi/Eddington/Lbol algebra and invalid-value rejection;
- source-order/reset-boundary verification;
- injected mid-group failure proving zero partial state mutation;
- successful multi-group commit proving one accounting advance;
- idle-sink zero-effective-efficiency acceptance;
- one-year convention consistency;
- required-field `null` rejection;
- duplicate-formatting and rewind-conflict rejection;
- stable `idsink` remapping after sink-array reorder;
- failed-transaction `snrt_intensity` checksum identity and same-step
  idempotence;
- an explicit negative regression showing that cross-step deferred
  re-emission remains open and is not falsely claimed closed;
- an OpenMP/MPI ownership/source-order static audit;
- Python syntax, diff whitespace, and the relevant native/production compile
  path where the local toolchain permits it.

The evidence must state which claims are arithmetic/transactional and which
remain unproven physical hydro closure.  No large RAMSES job is needed for
this bundle.

## Explicit non-goals

- no G2 yield-table source selection or 40–120 M☉ promotion;
- no AGN SED/obscuration/escape-fraction science approval;
- no new dust, scattering, IR, or radiation-pressure physics;
- no generic ksection, AMR redistribution, HDF5, checkpoint, or MPI
  hardening unrelated to this source transaction;
- no live production activation and no publication-ready AGN claim;
- no persistent hard-crash journal claim.  A durable journal and the
  cross-coarse-step deferred-accumulator policy remain separate follow-ups;
  the bundle must preserve negative evidence for both.

## Opus 5 plan-audit disposition

Claude Opus 5 returned **CONDITIONAL PASS**. It accepted the scope as a
bounded pre-G2 engineering bundle and confirmed that the transaction is a
real high-level feedback prerequisite. Its conditions were integrated above:
no run identity is added and rewind conflicts remain fail-closed; raw and
effective efficiency fields are separate; idle effective efficiency may be
zero; all paths use `365.25 d`; all inputs are finite; cross-step deferred
re-emission remains open; accounting uses `idsink`; serial/OpenMP/MPI
ownership and one-time environment resolution are explicit. Evidence must be
bound with the same-run freshness pattern rather than a dirty-worktree claim.

## Acceptance gate

The bundle is complete only if the canonical parser/audit and the source
transaction tests pass, the production source ordering is independently
verified, no runtime activation flag is relaxed, and the final evidence is
bound to the current worktree.  Claude Opus 5 is the sole active auditor: it
must first review this plan and then perform the bundle-end read-only audit.
AGY and Grok are not approval authorities for this bundle.

## Implementation checkpoint — 2026-09-04

The bounded implementation is now present in the production source path:

- `snrt_agn_deposit_transaction` prepares every group and commits only after
  all validation/overflow checks pass; the source smoke test covers both a
  successful two-group commit and an injected negative-group failure with
  byte-for-byte state preservation.
- The live source loop uses an `idsink`-keyed accounting map, skips duplicate or
  invalid identities, keeps the mutable source loop serial, and latches
  `SNRT_RT_ENABLE` once per process.  It avoids an unnecessary remap when the
  sink order is unchanged.
- The coarse-state writer emits the explicit pre-reset/instantaneous markers,
  raw/effective efficiency fields, and 365.25-day convention.  The Python
  reader canonicalizes semantic restart duplicates and rejects same-key
  conflicts, malformed values, invalid source algebra, and inconsistent epoch
  metadata.
- `bin/Makefile` now links the complete SNRT module graph, including
  `mpi_mod` and the CUDA ABI interface modules, under the existing
  `SNRT=1`/`USE_CUDA=1` gate.
- The synthetic arithmetic/transaction evidence is recorded in
  `simulation/snrt/data/agn_coarse_ledger_transaction_audit.json`.  It is not
  a physical AGN run and carries no SED, obscuration, hydro-closure, or
  cross-coarse-step journal claim.

The final Opus 5 bundle audit remains the acceptance step.  A real production
ledger must replace the synthetic input before any runtime or publication
claim is made.
