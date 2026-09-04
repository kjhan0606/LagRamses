# Claude Opus 5 pre-implementation plan audit — F-P1.5-R

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Auditor: Claude Opus 5 (primary; read-only)  
Prompt: `provenance/opus5_fp1_5_agn_effective_efficiency_convention_plan_audit_prompt_2026-09-04.md`  
Result: completed; no files modified, jobs run, or runtime activated  
Model cost reported: approximately USD 3.45

## Verdict

**APPROVE WITH CHANGES**

Opus confirmed that F-P1.5-R is the correct next bounded high-level task and
that the raw/effective efficiency mismatch is real and can be unbounded in a
MAD low state.  The initial plan fixed only one of two coupled convention
errors and required the following mandatory amendments before implementation.

## Mandatory amendments

1. Reconcile the retained black-hole mass used by the driver with the supplied
   inflow mass used by the coarse ledger.  `accrete_bondi` applies
   `dM_retained = (1-epsilon_eff) dM_inflow`; therefore the photon-budget API
   must declare whether its first argument is supplied inflow or retained BH
   mass and must apply the corresponding `(1-epsilon_eff)` conversion.
2. Make `spin_bh` an explicit input to the shared efficiency helper.  When
   spin is disabled, the legacy caller uses the default `0.1` rather than the
   possibly zero `eps_sink` array value.  The helper must reproduce the same
   branch in the coarse writer and driver.
3. State the asymmetric validity contract: raw efficiency is strict
   `(0,1)`, while effective efficiency is `[0,1)` because a MAD-quenched
   state may have effective zero and therefore zero photons.  Remove the
   driver's hidden `0.99` clamp and make `X_floor<=0`/zero-Eddington handling
   identical.
4. Put the shared helper in the unconditional `MODOBJ` path, not only in the
   SNRT conditional object list; add direct prerequisites for both consumers
   and the missing `snrt_ramses_driver.o` dependencies on `amr_commons` and
   `pm_commons`.  Include a non-SNRT default-build check.
5. Update the stale mismatch assertions in
   `simulation/snrt/tools/audit_agn_coarse_ledger.py:204` and
   `simulation/snrt/tests/agn_ledger_transaction.py:156–158`.
6. Document module ownership for `X_floor`, `mad_jet`, `spin_bh`,
   `dMBHoverdt`, and `dMEdoverdt`, with explicit `use ..., only:` imports in
   the driver.
7. Explicitly retain the existing hard-coded `0.5 * group_energy_fraction`
   spectral split as unchanged and unapproved; this bundle does not approve
   an AGN SED.

## Important feasibility observations

- The writer computes effective efficiency from Bondi/Eddington state, while
  the driver currently clamps raw `eps_sink` and passes it directly to
  `snrt_agn_photon_budget`; the mismatch is material, not rounding.
- `dMsmbh` is retained BH mass (`accrete_bondi:4474`), while the writer's
  instantaneous `inflow_rate` is supplied mass (`sink_particle:2230`).
- `eps_sink` is written by the spin branch (`kjhan_growspin:6126`) and the
  accretion callers use `0.1` when `spin_bh` is false; the coarse writer must
  not silently diverge from that branch.
- The `365.25 d` convention and cgs unit scales already agree; the driver
  operates in seconds and needs no new year conversion.
- The native `sinkprops` conversion has no mode-resolved field and should
  retain its explicit `raw == effective` review-only status.

## Scope accepted by Opus

The bundle remains engineering-only: it may align efficiency/mass
provenance, source/build contracts, and arithmetic evidence, but must not
enable `SNRT_RT_ENABLE`, select an AGN SED/obscuration model, alter the
all-or-nothing photon transaction, or make a physical hydro/dust/yield claim.

Implementation was not authorized by this audit.  The plan must be amended,
then receive the normal pre-implementation approval boundary before source
changes begin.
