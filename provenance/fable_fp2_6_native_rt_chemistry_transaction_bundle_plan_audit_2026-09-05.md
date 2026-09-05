# Fable plan audit — F-P2.6 native RT/chemistry transaction and fixed point

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Fable
Mode: read-only plan review; no files, jobs, builds, or simulations modified

## Verdict

`CONDITIONAL PASS`

The proposed bundle is the correct next native RT/chemistry correctness
boundary for the project's production/publication high-level RT, stellar/AGN
feedback, and dust goal. The exclusions are clean. Implementation is
authorized after the plan edits recorded in the current plan file.

## Required plan corrections

### B1 — transaction must include coarse/interface flux state

The coarse-level reverse flux correction mutates persistent coarse-leaf
intensity before the prepared transport call completes. A level rollback or a
repeated fixed-point trial would otherwise leave coarse neighbors with flux
from a failed or discarded trial.

Disposition: make the correction transaction-owned, accumulate it in a level
trial buffer, and commit it only with the converged trial. The transport
routine must return the trial field instead of writing persistent coarse
intensity during the trial.

### B2 — commit, rollback, and convergence are MPI collective

Prepared transport already performs collectives and virtual-cell exchange. A
rank-local rollback or rank-local iteration count would consume inconsistent
ghost radiation or deadlock at the next collective. Ranks with no leaves must
also take the same collective decision path.

Disposition: reduce failure and convergence norms globally before every
decision; require the same trial count and branch on every rank. Multi-rank
decision correctness is in scope, while performance/load balance remains G5.

### B3 — native smokes need an explicit-array core and controlled failures

The current RAMSES driver is too stateful for deterministic standalone
rollback testing, and natural CUDA caps prevent the required failure paths.

Disposition: put the transaction/fixed-point core in a new explicit-array
module with a thin driver adapter. Add deterministic smoke-only failure
injection for a named leaf/stage, reject it under approved production status,
and test convergence against a one-cell bisection reference. Test
non-convergence with one iteration or zero tolerance.

### H1 — fixed-point physical definition is incomplete

The plan must specify the closing opacity, frozen quantities, relaxation target,
same-trial commit, and convergence norm. Recommended contract: time-centered
neutral opacity from start and relaxed end states; fixed incoming radiation,
start-of-step inventory, and pre-heating temperature; relaxation only on
opacity-feeding fractions; same-trial intensity/absorption/chemistry/heating
commit; global max of fraction and relative optical-depth changes.

Disposition: the plan now declares these definitions and explicit initial
defaults.

### H2 — rollback continuation is unsafe

The current driver call site has no return path that makes continuing after a
rollback safe. Continuing would freeze chemistry while photons accumulate.

Disposition: route failure through RAMSES clean-stop with failure class,
residual, level, and step. Any continue mode must be an explicit
non-production diagnostic option.

### M1 — inventory tolerance must be scale-aware

The fixed absolute `5e-5` code-density tolerance can be a large fraction of a
low-density inventory and is not a floating-point-sized bound.

Disposition: replace it with a documented multiple of single-precision
epsilon times the larger absorbed/inventory scale, with a double-precision
host floor and a level-summed residual ledger.

### M2 — source-phase ordering must be explicit

The source phase mutates intensity and source-accounting markers before
transport. A pre-source snapshot could erase already-accounted source photons
on rollback.

Disposition: snapshot after source deposit/accounting; keep the source phase as
its own committed transaction.

## Advisory conditions

- Name the H I mirror in the snapshot list.
- Use photon-number closure wording rather than photon-energy closure, because
  the emission-mean versus absorber-weighted heating gap remains open.
- Scope determinism claims to a fixed MPI rank count; reverse virtual-cell sum
  order is not rank-count invariant.
- Restore by copy so a failed call is bitwise equivalent to no call.
- Record iteration cap and its runtime cost multiplier in startup evidence.

## Final assessment

After these edits, F-P2.6 is a feasible, appropriately scoped native bundle
that materially advances the high-level RT/feedback objective. It must not be
used to claim physical SED/yield approval, dust completeness, live production
hydro validation, or publication readiness.
