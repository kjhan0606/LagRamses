# F-P1.5-R AGN effective-efficiency convention bundle plan — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Parent: F-P1.5 AGN coarse-ledger transaction bundle  
Status: **Claude Opus 5 plan re-audit APPROVE WITH CHANGES; mandatory
amendments integrated; implementation approved; focused evidence complete;
final bundle-end audit PASS; F-P1.5-R engineering bundle closed**

## Why this is the next candidate

F-P1.2 is closed with a Claude Opus 5 `PASS` as an engineering transaction
bundle.  The remaining high-level AGN source-path mismatch is explicit in
`provenance/opus5_fp1_agn_ledger_transaction_bundle_end_audit_2026-09-04.md`:
the coarse-state writer computes `effective_radiative_efficiency` for `Lbol`,
while the live SNRT driver currently passes clamped raw `eps_sink` to the
photon budget.  In addition, the driver passes retained BH mass (`dMsmbh`)
where the writer's luminosity uses the Eddington-limited supplied inflow
`min(dMBH_coarse,dMEd_coarse)`.  `dMBH_coarse` alone is Bondi supply, not the
writer's inflow.  In a MAD low state the efficiency mismatch can be unbounded,
and even at ordinary efficiency the mass-convention mismatch produces a factor
`1/(1-epsilon)`.  This bundle removes both bookkeeping ambiguities before any
live activation or physical AGN/SED claim.

This is a bounded RT/AGN-feedback wiring task.  It does not select a new AGN
efficiency model, approve an AGN SED, open `SNRT_RT_ENABLE`, or claim hydro
closure.  Physical source and yield gates remain fail-closed.

## Scope

### A. One code-owned efficiency contract

1. First confirm and record the actual initialization branches.  The compiled
   default is `spin_bh=.true.`; `eps_sink` is initialized to zero and is only
   written by `kjhan_growspin` at the end of its first spin update.  With
   `spin_bh=.false.`, the accretion callers deliberately use the model default
   `0.1`; with spin enabled during the uninitialized window, the legacy
   accretion caller reads zero.  The shared helper must expose a distinguishing
   status for these cases rather than silently pretending they are the same
   physical efficiency.
2. Add one pure, RAMSES-independent AGN efficiency helper shared by the coarse
   ledger writer and the SNRT driver.  Its inputs are the raw sink efficiency,
   `spin_bh`, instantaneous Bondi and Eddington rates, `mad_jet`, and
   `X_floor`; its outputs are the resolved base efficiency, effective value,
   inflow rate, and a controlled status.
3. Preserve the existing declared mode transformation exactly:
   `inflow = min(max(Bondi,0), max(Eddington,0))`,
   `edd_ratio = Bondi/Eddington` when the denominator is positive and zero
   otherwise, and for MAD low-state operation
   `epsilon_eff = epsilon_raw * max(edd_ratio,0) / X_floor` when
   `edd_ratio < X_floor`; otherwise `epsilon_eff = epsilon_raw`.
4. Make the asymmetric validity policy explicit at the shared boundary:
   resolved/raw efficiency is strict `(0,1)`, while effective efficiency is
   `[0,1)` because MAD quenching may legitimately produce zero photons.
   Define explicit statuses for `spin_bh=.false.` default, uninitialized
   spin-enabled `eps_sink`, non-finite/non-positive input, and `eps_sink >= 1`.
   `X_floor <= 0` and zero-Eddington handling must be identical in both
   callers.  Remove the driver's hidden `0.99` clamp; no caller-specific
   clamp or fallback is allowed.  An uninitialized spin-enabled efficiency is
   not described as exact legacy parity and must remain a visible review/live
   admission condition.
5. Keep `radiative_efficiency` (raw/resolved base state) and
   `effective_radiative_efficiency` (value used for luminosity and photon
   production) as distinct fields.  Emit a status/mode marker so a review
   record cannot be mistaken for a mode-resolved physical measurement.

### B. Wire the effective value through both consumers

1. Replace the duplicated efficiency arithmetic in
   `patch/lagRamses/sink_particle.kjhan.f90` with the shared helper when
   building `agn_coarse_state_v1.jsonl` and `bolometric_luminosity_erg_s`.
2. In `patch/lagRamses/snrt_ramses_driver.f90`, compute the same effective
   value from the current RAMSES sink/rate state and pass **only the effective
   value** to `snrt_agn_photon_budget`.  Retain raw `eps_sink` only as an
   input to the helper and diagnostic provenance, never as the photon-budget
   coefficient.
3. Declare the first argument of `snrt_agn_photon_budget` as supplied inflow
   mass, not retained BH mass.  Use the per-step increment of
   `min(dMBH_coarse,dMEd_coarse)` keyed by `idsink` for source accounting;
   `dMBH_coarse` is Bondi supply and `dMEd_coarse` is the Eddington limit.
   Use `dMsmbh` only for a one-sided retained-mass consistency check:
   `dMsmbh <= (1-epsilon_eff) * dM_inflow * (1+tol)`, with the existing
   gas-floor clipping and saved-energy suppression recorded as admitted slack.
   The cumulative rate-based supplied ledger is therefore an upper-bound
   source convention when the gas floor clips retained mass; it is not an
   equality claim for actually retained/accreted BH mass.
   Do not require equality and do not disable a source merely because the gas
   floor clips the retained increment.  The singular `epsilon_eff=1` case,
   unexplained excess retained mass, and the uninitialized spin-enabled
   divergence must fail closed or remain explicitly non-promotable; the two
   mass conventions must never be mixed silently.
4. Add the helper object to the unconditional `/gpfs` Makefile `MODOBJ` path,
   with its `amr_parameters.jaehyun.o` dependency as needed.  Add direct
   prerequisites for both consumers and add the missing
   `snrt_ramses_driver.o` dependencies on `amr_commons.o` and `pm_commons.o`
   (the hydro dependency is already transitive through `pm_commons.o`).  The
   source/build report must show one helper implementation, both call sites,
   its provenance hash, and no independent formula.  Include the third
   positional consumer `snrt_agn_source_smoke.f90` in the compile contract.
5. Keep the existing serial `idsink` source loop and all-or-nothing photon-group
   transaction unchanged.  This bundle changes the coefficient provenance and
   supplied-inflow accounting convention, not the runtime activation policy.
6. Record module ownership explicitly: `X_floor`, `mad_jet`, and `spin_bh`
   come from `amr_parameters`; `dMBHoverdt`, `dMEdoverdt`, `dMBH_coarse`,
   `dMEd_coarse`, `dMsmbh`, `eps_sink`, and `idsink` come from `pm_commons`;
   the driver uses explicit `use ..., only:` lists for every one of these
   inputs.

### C. Evidence and regression

Add or update focused evidence for:

- identical raw/effective/inflow results from the writer and driver for
  thermal, MAD-low, MAD-high, zero-inflow, and boundary `X_floor` states;
- retained-versus-supplied mass conversion, including
  `dM_inflow=min(dMBH_coarse,dMEd_coarse)`, the one-sided gas-floor-clipped
  retained-mass bound, photon-budget argument semantics, and a fail-closed
  `epsilon_eff=1`/excess-retained case;
- `spin_bh=.false.` default `0.1`, compiled `spin_bh=.true.` default, and the
  spin-enabled zero/uninitialized `eps_sink` status/divergence;
- invalid/non-finite rates, efficiencies, and floor values failing or falling
  back according to one explicit shared policy, with raw `(0,1)` and
  effective `[0,1)` checked separately;
- exact `Lbol = epsilon_eff * inflow * c^2` and photon-rate scaling;
- a negative static check proving raw `eps_sink` is not passed to the photon
  budget;
- coarse-ledger JSON raw/effective/status field preservation and audit output
  no longer claiming an unresolved writer/driver mismatch;
- source-order/helper identity, helper SHA256 supplied to the audit tool (via
  its explicit helper-path input), and direct Makefile dependencies;
- default non-SNRT build as well as the SNRT/CUDA build graph;
- `snrt_agn_source_smoke.f90` positional API compile and photon-budget mass
  convention;
- unchanged all-or-nothing multi-group transaction, stable sink identity,
  restart-duplicate handling, and failed-transaction state identity;
- deferred-`Esave`/nonzero carry-over case proving no duplicate source emission
  after the mass-accounting switch;
- Python syntax, focused native/helper compilation, production object build,
  SNRT/CUDA dry-run, and `git diff --check` on `/gpfs`.

The tests are arithmetic and source/build-contract evidence.  No large
RAMSES run is needed and no runtime environment flag may be enabled.

## Acceptance gate

The bundle is complete only if:

1. the coarse writer and SNRT driver call the same code-owned helper;
2. both `Lbol` and photon budgets consume the same effective efficiency for a
   given state and the photon-budget mass argument is the limited supplied
   inflow;
3. raw/base and effective fields remain separately auditable, with explicit
   initialization statuses and no silent clamp or caller-specific default;
4. the one-sided retained-mass/gas-floor contract and `Esave` carry-over policy
   are tested without false equality claims;
5. the audit tool verifies the helper path and hash, all existing AGN
   transaction/accounting tests remain green, and both default and SNRT build
   paths pass;
6. Claude Opus 5 completes one read-only bundle-end audit.

The result may close an engineering convention gate only.  It must not set
`SNRT_RT_ENABLE`, approve an AGN SED/obscuration/escape prescription, claim
physical hydro closure, or promote a yield/fate source.

## Explicitly deferred

- AGN SED shape, obscuration, escape fraction, spectral calibration, and
  radiation-pressure/jet coupling physics;
- parity between the new writer/driver convention and the independent legacy
  `accrete_bondi` MAD transform or `AGN_blast` feedback-energy prescription;
  those paths remain outside this bounded helper contract;
- the existing `0.5 * snrt_group_energy_fraction` spectral split remains
  unchanged and unapproved; this bundle does not validate or approve that SED
  assumption;
- live RT-hydro activation, AMR/MPI production coupling, and large runs;
- cross-coarse-step deferred-energy journal and hard-crash exactly-once state;
- stellar physical-source approval, 40–120 M☉ fate promotion, and dust/IR
  closure;
- non-blocking F-P1.2 helper cleanup and `deposit_one_star` live integration
  coverage.

## Governance

This plan is prepared after the F-P1.2 bundle-end `PASS`.  Claude Opus 5 is
the sole active plan/bundle auditor.  GPT-5.6-Sol is called only if Opus
cannot issue a verdict or the operator explicitly requests confirmation.  The
operator approved implementation on 2026-09-04; the implementation and
focused evidence are complete.  The first bundle-end FAIL, conditional-pass
F1 finding, and C1 whole-ledger readability finding were repaired.  The final
Opus re-audit returned `PASS`; only non-blocking follow-ups before live
activation remain.
