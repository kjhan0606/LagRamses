# Claude Opus 5 plan audit — F-P1.5 AGN ledger transaction bundle

Date: 2026-09-04  
Mode: read-only; no files edited, built, committed, pushed, or jobs launched.  
Repository: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, `main`

## Verdict

**CONDITIONAL PASS.** The proposed scope is bounded, in scope, and suitable
as pre-G2 engineering work. It maps to the existing P1.4/P1.5 priorities and
does not approve the AGN SED, physical hydro coupling, G2 yields, or runtime
activation. Implementation may proceed after the conditions below are folded
into the plan.

## Required conditions

1. A same-key duplicate needs a discriminator for an identical restart versus
   a rewound rerun. Since no run UUID or dump counter is currently present,
   either add one or explicitly keep rewind safety open and reject any
   same-key payload conflict. The latter is the selected bounded disposition.
2. Keep raw and effective radiative efficiency as distinct fields and CSV
   columns. The current reader aliases the effective value to a field named
   `radiative_efficiency`, while the live driver still consumes raw `eps_sink`;
   this F10 convention mismatch must remain open rather than being hidden.
3. Permit `effective_radiative_efficiency=0` for idle sinks and zero luminosity;
   the current open `(0,1)` reader bound is too strict.
4. Harmonize the year convention. The Fortran writer uses 365 days while the
   Python reader uses 365.25 days; otherwise the planned algebra check is not
   round-off level. The bundle selects 365.25 days everywhere.
5. Add finite-value checks. The current Fortran source routine can let NaN
   pass its `<=0` and overflow checks.
6. An all-or-nothing group transaction prevents intra-call partial deposit,
   but does not close the cross-step deferred re-emission path: `accounted_mass`
   resets at a new coarse step while `dMsmbh` can remain nonzero after a
   deferred blast. Keep this as explicit negative evidence.
7. Key the in-memory accounting by stable `idsink`, not mutable sink-array
   index, and test a sink-array reorder.
8. State and test the concurrency contract: the source loop and saved
   accounting arrays are serial/non-threadprivate; local-leaf lookup gives a
   single MPI owner; and `SNRT_RT_ENABLE` must not be re-read mid-step in a
   way that strands accounting state.

## Required tests

Add tests for year consistency, idle sinks, required-field nulls, formatting
duplicates and rewind conflicts, stable sink identity, failed-transaction
state checksum identity, same-step idempotence, and the still-open
cross-step deferred re-emission case. Bind generated evidence with the
same-run freshness pattern; a dirty worktree is not production evidence.

## Evidence boundary

The findings and proposed fixes are arithmetic, wiring, and transaction
claims. They do not establish physical AGN SED, obscuration, escape fraction,
thermal/jet, radiation-pressure, or live hydro closure. Generic AMR/HDF5
hardening is outside this bundle.
