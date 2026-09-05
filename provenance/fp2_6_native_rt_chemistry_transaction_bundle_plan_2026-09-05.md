# F-P2.6 native RT/chemistry transaction and fixed-point bundle plan

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Repository: `kjhan0606/LagRamses`
Parent: F-P2.5 native H/He thermochemistry bundle

Status: Fable plan audit returned `CONDITIONAL PASS`; its required plan
clarifications are applied below. The operator pre-approved continuation on
2026-09-05. The first bundle-end audit attempt by Claude Opus 5 produced no
verdict, so Fable's fallback audit returned `FAIL`. The identified native
remediation, final hardening, refreshed direct evidence, and production-link
rebuild are complete. The final Fable fallback closure audit returned
`CONDITIONAL PASS` with no blocker; the remaining live-loop evidence and
diagnostic/G5 records are carried forward explicitly.

## Project decision boundary

The project targets production- and publication-ready high-level RAMSES
radiative transfer, stellar/AGN feedback, and dust physics. F-P2.6 addresses a
native RT/chemistry correctness boundary that must be closed before a live
feedback run. It does not select a physical stellar or AGN SED, approve yield
tables, or claim production hydro validation.

## Objective

Make the native `transport -> species absorption -> H/He chemistry -> RAMSES
thermal receiver` path transactional and quantitatively self-consistent:

1. no absorbed photons or thermal/chemical energy disappear when partition or
   chemistry fails;
2. any unsupported/unassigned absorption is a visible, fail-closed condition
   rather than an unowned production loss; and
3. the native path has a bounded, declared local opacity/chemistry fixed-point
   iteration instead of relying only on start-of-step opacity.

## Work packages

### C1 — level transaction and failure rollback

- Define a level transaction around the prepared native RT call, beginning
  after the AGN/stellar source phase has committed its source deposit and
  source-accounting marker. The source deposit remains its own committed
  transaction; a failed RT/chemistry level must not erase it.
- Snapshot the incoming photon state, persistent H/He state (including the
  H I mirror), and RAMSES thermal-energy field after that source phase.
  Commit them only after all leaves pass species partition, photon-number and
  photoelectron-energy closure, chemistry, and receiver checks.
- Treat coarse-parent/interface flux correction as transaction-owned. The
  transport routine must accumulate reverse-communicated coarse corrections
  in a level trial buffer and return the trial field; it must not write those
  corrections into persistent coarse-leaf intensity before the transaction
  commits. This covers both rollback and repeated fixed-point trials.
- On a partition, thermochemistry, non-finite-state, or receiver failure,
  restore the complete level snapshot and return a nonzero error. No failed
  cell may leave a consumed photon packet, updated ion fraction, or partial
  `uold` heating behind.
- Keep per-cell and per-group absorption ledgers available for diagnostics,
  but treat the level rollback as the atomic boundary because transport
  couples neighboring cells.
- Reduce failure status collectively with MPI before any commit or rollback
  decision. Every rank, including ranks with zero leaves, executes the same
  trial count and collective decision path. Production correctness is
  multi-rank collective; load balance and scaling remain G5 work.

### C2 — explicit unassigned-absorption receiver gate

- Promote `unassigned_absorption_code` from a print-only diagnostic to a
  versioned closure field with an explicit numerical-residual tolerance.
- If the partition produces a residual larger than the declared floating-point
  tolerance, fail the level transaction. A tolerance-sized residual must still
  be included in the reported photon/energy closure residual; it may not be
  silently assigned to gas heat or ionization.
- If a cell's chemistry fails after transport, the rollback in C1 must restore
  the removed radiation state. Record the failure class and residual for the
  caller; do not count the failed absorption as delivered feedback.
- Replace the current fixed absolute inventory tolerance with a scale-aware
  bound derived from the measured FP32 reduction error: use a documented
  multiple of single-precision epsilon times the larger of the cell's
  absorbed and available-inventory magnitudes, with a double-precision floor
  for the host ledger. Record the derivation and level-summed residual in the
  commit ledger; do not permit a low-density cell to inherit a fixed
  `5e-5`-scale tolerance.

### C3 — native bounded opacity/chemistry fixed point

- Add a versioned native iteration contract with these initial defaults:
  maximum 32 trials, fraction absolute tolerance `1e-6`, relative optical-
  depth tolerance `1e-5` above an optical-depth floor of `1e-12`, and explicit
  under-relaxation factor `0.5`. The convergence norm is the global maximum
  over leaves of the absolute H II/He II/He III fraction change plus the
  relative group-optical-depth change above the floor. All defaults are
  recorded in native configuration and startup evidence.
- Within the prepared level boundary, restart each trial from the same incoming
  photon state, recompute species opacity and transport absorption, advance
  H/He chemistry, and under-relax only the species fractions that feed opacity.
  Use time-centered neutral fractions formed from the start-of-step and
  relaxed end-of-step states. Hold incoming radiation, start-of-step species
  inventory, and temperature fixed across trials; recombination uses the
  pre-heating temperature. A relaxed trial state is never committed by
  itself.
- The committed intensity, group absorption, species update, and heating must
  all come from the same final trial. A consumed trial photon field must never
  be carried into the next iteration. Only a converged trial may be committed
  by C1; a non-converged trial participates in no state mutation.
- Reduce the convergence norm and failure flag collectively before deciding
  to stop. Every rank runs the same configured number of trial branches, and
  a non-converged result fails closed.
- Keep this first fixed point bounded to the native prepared-level path;
  multi-rank decision correctness is in scope, while distributed performance,
  load balance, and the full hydro timestep limiter remain later G5 work.

The convergence residual is intentionally the residual of successive
under-relaxed opacity predictors (`relaxed_k` versus `relaxed_{k-1}`), with
the declared fraction and optical-depth norms. The raw chemistry map is the
state that is committed from the final trial; a strict raw-map residual is not
silently claimed. The worst-case native cost is therefore up to 32 prepared
transport evaluations per level, each carrying its configured transport
subcycles, in addition to the one-time preparation and commit.

### C4 — native evidence and build boundary

- Implement the transaction/fixed-point core in a new module over explicit
  arrays, with the RAMSES driver as a thin hydro/AMR adapter. Add
  deterministic smoke-only failure injection for a named leaf and stage
  (partition, chemistry, or receiver). The runtime must reject the injection
  control under an approved-production contract and print its disabled status
  at startup.
- Add Fortran/CUDA-native smokes for: successful commit, partition failure
  rollback, chemistry failure rollback, nonzero unassigned rejection,
  converged fixed point, and non-convergent fail-closed behavior. Define the
  non-convergence smoke as maximum-iteration-one or zero tolerance, and use a
  one-cell bisection solution as the convergence reference.
- Check conservation before and after the transaction for photon number,
  species inventory, photoelectron energy, and RAMSES thermal energy. Check
  that only the converged gas-heating receiver changes `uold`.
- Build the changed production module graph with GNU and `mpiifx`, link the
  CUDA binary, run the focused native smokes on `/gpfs`, and record source,
  configuration, binary hashes, iteration cap, and expected cost multiplier.
  On a transaction failure, use the RAMSES clean-stop path with failure class,
  residual, level, and step in the message; continuing is allowed only under
  an explicit non-production diagnostic flag. Do not launch a large RAMSES
  evolution.

## Acceptance gates

- a failed prepared-level call restores the complete snapshot by copy and is
  bitwise equivalent to no RT/chemistry call for photon, H/He, and thermal
  state; source-phase deposits remain committed;
- unassigned absorption cannot pass the native production boundary without a
  declared residual and a successful receiver decision;
- the fixed-point iteration is deterministic at a fixed MPI rank count,
  converges on the reference one-/two-cell controls, rejects non-convergence,
  and never reuses consumed trial photons;
- successful transactions close photon, species, photoelectron, and thermal
  ledgers and update RAMSES energy only once;
- native Fortran/CUDA smokes, GNU/`mpiifx` builds, `git diff --check`, and the
  recorded evidence pass; and
- no physical SED, yield asset, dust model, live hydro run, or publication
  claim is promoted by this bundle.

## Explicit exclusions

This bundle does not resolve the physical AGN/stellar SED decision, the
40--120 M_sun yield seam, SNIa/PISN source approval, dust scattering/IR/grain
evolution, cooling processes outside the native photoheating receiver, HDF5
restart integration, distributed AMR/MPI scaling, or publication convergence.
Those remain separate high-level physics or G5/G6 bundles.

## Audit and approval workflow

Fable is the primary plan auditor for this bundle. Claude Opus 5 is the plan
backup if Fable cannot issue a verdict. Implementation begins only after the
plan audit and operator approval. Claude Opus 5 is the primary end-of-bundle
implementation auditor; Fable is its backup if Opus cannot perform or decide.

## Fable plan-audit disposition — 2026-09-05

Fable returned `CONDITIONAL PASS`. The plan was accepted after adding the
transaction-owned coarse correction, MPI-collective decisions including
zero-leaf ranks, an explicit-array transaction core with smoke-only failure
injection, the fixed-point temperature/inventory/opacity/norm definitions,
source-phase snapshot ordering, scale-aware inventory tolerance, and
clean-stop behavior. The audit report is recorded at
`provenance/fable_fp2_6_native_rt_chemistry_transaction_bundle_plan_audit_2026-09-05.md`.

## Initial bundle-end audit and remediation — 2026-09-05

Claude Opus 5 was attempted as the primary implementation auditor. The
`opus-5` model was rejected by the local CLI catalog and the fallback `opus`
alias timed out after 900 seconds without a verdict; the attempt record is at
`provenance/claude_opus5_fp2_6_native_rt_chemistry_transaction_bundle_end_audit_attempt_2026-09-05.md`.
Under the audit governance, Fable was then used as the fallback. Its initial
verdict was `FAIL`, recorded at
`provenance/fable_fp2_6_native_rt_chemistry_transaction_bundle_end_audit_2026-09-05.md`.

The remediation actions are:

- slot-indexed H/He/H I snapshot and restore, zero-leaf allocated payloads,
  full-capacity coarse trial buffers, collective transport/topology error
  reductions, and collective driver preflight are now implemented;
- scale-aware host inventory tolerance is preserved at the pre-partition cell
  scale, and the CUDA species cap now uses a bounded active-set allocation so
  the capped photon removal equals the actually assigned inventory;
- finite/non-negative trial-field checks, the non-production diagnostic
  failure-injection mode, all three driver failure routes, and global
  unassigned ledgers are now explicit; and
- the direct smoke runner now exercises all named injection stages, mixed
  species inventory caps, GNU/mpiifx/MPI paths, and a static driver-route
  check. Refreshed hashes, cost accounting, and residual semantics are now
  recorded; the final Fable closure audit is the remaining bundle action.

The following remain explicitly record-only or later G5/G6 work: live driver
failure injection without a small initialized RAMSES harness, no-slot coarse
correction accounting, distributed AMR performance/load balance, HDF5 restart
integration, physical SED/yield approval, dust processes, and publication
convergence.

## Final closure disposition — 2026-09-05

Fable's governed fallback closure audit returned `CONDITIONAL PASS` with no
blocking finding. The stale final-binary hash was refreshed in the evidence
record. C1 and C2 are `PASS`; C3 is `PASS` with a record-only live-loop
evidence gap; C4 is `CONDITIONAL` for the same evidence limitation. The
remaining low-severity first-failure/structured-ledger diagnostics, coarse
face/no-slot accounting, rank-uniform pre-collective policy, and distributed
AMR concerns are not required to close this native transaction boundary and
are assigned to later G5/G6 work. This bundle is therefore implementation-
complete conditionally, not a grant of live feedback or publication approval.

Closure report:
[`fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md`](fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md)
