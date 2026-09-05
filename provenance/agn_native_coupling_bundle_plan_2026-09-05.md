# AGN native coupling bundle — cell conservation and source ownership

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Status: operator approved entry ("tmddls", Korean-keyboard "승인"); Fable
CONDITIONAL APPROVE; operator explicitly approved both recommended corrections
with "적용하고 진행해"; implementation complete within the declared scope.
Opus end review: CONDITIONAL PASS; bounded B1–B3 repairs: closure PASS.
Evidence: `agn_native_coupling_bundle_evidence_2026-09-05.md`.

Completion qualification: live SNRT with AGN sinks is restricted to serial
fresh starts until source history is durable/migratable. Simultaneous legacy
AGN feedback plus live SNRT is rejected. No full RAMSES link, live evolution
or production/publication qualification is claimed. These restrictions take
precedence over the original work-package availability language below.

## Plan-review disposition

The read-only Fable review completed in about seven minutes; full response:
`fable_agn_native_coupling_plan_audit_2026-09-05.md`. No production code was
changed during this planning turn. The following amendments take precedence
over the initial work packages below.

Both choices below received explicit operator confirmation before implementation:

1. **Correct the resolved-jet normalization**, rather than preserving the
   legacy loss for comparison. Use a volume-weighted kernel so the total
   returned mass equals the loaded mass removed by `average_AGN`. The old
   uniform two-cell example removes mass 1 but returns 0.5. This changes
   effective jet strength relative to old runs; record it as a conservation
   correction, never claim unchanged calibration. The related kinetic budget
   normalization must be checked with the actual allocated mass, not patched
   by an arbitrary efficiency rescaling.
2. **Return unresolved fallback mass at the donor gas velocity**, not the BH
   velocity in the initial draft. With the existing mass-removal operation
   preserving gas velocity, this avoids introducing an unfunded momentum
   change merely to return the same entrained material. The original BH-frame
   resolved jet and its global donor/recoil accounting must not be declared
   closed by this fallback repair alone.

Accepted engineering amendments:

- Treat the midplane as equal opposed lobe contributions: zero directed net
  momentum, retained kinetic energy; membership/radicand handling must match
  in both averaging and deposition passes.
- Check the source time ledger as part of ownership. `AGN_feedback` owns
  the coarse supply reset; the SNRT cursor currently resets each coarse step.
  With legacy feedback disabled, claiming SNRT-only correctness requires
  carrying the cursor across steps (or rejecting that mode until corrected).
- Use rank-uniform namelist flags for conflict checks. A driver-only check
  cannot protect a legacy injection that occurs earlier in the coarse step;
  startup/source-dispatch placement must cover that ordering too. Do not
  blindly treat the auditor's driver-latch placement as sufficient for the
  plan's stronger before-either-source requirement.
- Drop the proposed production-negative-run extension: it cannot reach this
  AGN boundary past the unrelated stellar fate gate. Reuse the existing
  native runner; add the kernel-sum example, fold replay into thermal cases,
  and do not create GPU-parity work for the inert CUDA sink implementation.
- Keep one bundle-end Opus review after implementation and evidence, not
  another plan audit for each bounded amendment.

Simultaneous RT+mechanical physical allocation remains a separate explicit
model decision; these corrections do not approve or activate such a model.

## Goal and limits

The final project goal is production/publication-ready RT, stellar/AGN
feedback and dust in native RAMSES. The operator wants actual hydro changes,
larger bundles and few audits, not proliferating Python instruments. AGB
KL16/CK22 raw-source issues remain parked and do not block this work.

This is the next AGN implementation bundle after the accepted F-P1.5-R
efficiency/source-transaction work. It must improve real cell injection and
prevent silently treating the independently implemented legacy and SNRT
models as an approved combined physical model. It does not authorize a
large run, a new AGN SED, or a new partition of accretion/spin energy.

## Existing boundary and concrete defects

- `sink_particle.kjhan.f90:AGN_blast` injects legacy thermal/kinetic feedback.
  Its energy is based on retained `dMsmbh_AGN` and optional saved `EsaveAGN`.
- `snrt_ramses_driver.f90` emits AGN photons using the supplied-inflow
  convention and the shared effective-efficiency resolver. That convention
  is an engineering approximation, not equality with actually retained mass.
- Thermal cap-crossing branches add the entire capped internal energy to
  the previous total energy rather than adding only available headroom.
  Their saved-energy bookkeeping therefore need not close.
- Single-cell fallback branches use `vol_loc` left over from the AMR loop
  when saving excess energy, instead of the receiving `vol_blast(iAGN)`.
- The jet direction leaves u/v/w undefined when `dzjet == 0`; the radial
  square root also needs a roundoff-safe nonnegative radicand.
- The single-cell jet fallback adds returned density but not its bulk
  momentum or kinetic energy. The approved correction returns it at the
  captured donor gas velocity.

## Work packages (one bundle, no intermediate audit)

### A. Native gas/energy conservation

Use one small pure Fortran cell-injection routine, called by actual
`AGN_blast`, that stages mass, vector momentum, total energy and deferred
energy before mutation. Given the current cell, explicit mass/momentum/
energy increments, cell volume and temperature cap, preserve

`E_old + E_injected = E_new + E_deferred`.

Preserve pre-existing hot gas (the cap limits new internal-energy input,
not a cooling operator), and reject invalid nonfinite/unphysical inputs.
Use the receiving cell volume in all branches. Apply it to distributed
thermal, saved-energy replay, resolved jets and single-cell fallback. Preserve
kernel shape and nominal efficiency parameters, but correct normalization
and explicitly record the resulting change to effective jet strength.

For fallback mass return include `dm*v_donor` and `dm*|v_donor|^2/2`; retain the
existing unresolved-jet decision to thermalize the supplied EAGN budget.
For a resolved jet exactly on its midplane, use zero directed net jet
momentum and retain the local opposed-stream kinetic budget as internal
energy; do not discard that energy or pick an arbitrary lobe. Otherwise
preserve the existing bipolar direction and kernel. Bound the tiny negative
roundoff in `dr_AGN-dzjet**2` before sqrt.

### B. Honest SNRT/legacy source-ownership boundary

Trace the existing dispatch flags and prevent unapproved simultaneous live
AGN photons and legacy AGN feedback before either source mutates state.
Reuse existing configuration/gate machinery; do not add an override that
claims physical approval. Legacy-only behaviour stays available, SNRT with
no legacy AGN feedback stays available subject to its existing runtime gates.
Do not silently suppress thermal feedback, reduce photon luminosity, or
relabel legacy jet energy as independent available spin energy.

This is an interim exclusivity boundary, **not completed simultaneous
RT+jet physical coupling**. The latter needs a declared common time-integrated
mass/energy convention, mechanical versus radiative ownership, SED and
absorbed/escaped energy treatment. Record this next physical decision without
building an unused allocation framework or treating AGB data as its prerequisite.

### C. Proportional evidence and one end audit

Extend the existing `run_fp15_agn_efficiency.sh` native runner and a small
Fortran smoke driver to exercise the production cell routine: cap below/
crossing/already-hot, nonunit fallback volume, moving mass return, opposed
jets/midplane, replay and invalid-input no-mutation. Check conservation of
mass/momentum and gas-plus-deferred energy with explicit synthetic values.
Compile changed production callers separately in a scratch build; do not
overwrite the shared executable. Ownership evidence is compile/static review,
not the unrelated production-negative runner. No new Python validator, schema or broad
checkpoint/AMR/MPI gate. No production or publication claim from a smoke.

Capture one evidence record, then request one Opus 5 bundle-end review
(Fable backup if Opus cannot give a verdict). Fable reviews this plan for
necessity, feasibility, and overinstrumentation, not implementation approval
of an invented physical source. Bounded conditions may be incorporated within
the operator-approved scope; material new physics requires a user decision.

## Physics references and interpretation

- Dubois et al. (2012), dual thermal/jet model:
  https://academic.oup.com/mnras/article/420/3/2662/980280
- Dubois et al. (2021), NewHorizon, Section on MBH feedback:
  https://www.aanda.org/articles/aa/pdf/2021/07/aa39429-20.pdf

The latter describes quasar thermal coupling as representing both unresolved
wind thermalisation and radiation coupling. Thus automatic thermal removal
when RT is enabled is not justified by the word "radiative" alone. The
interim exclusivity rule is an implementation safety decision, not a claim
that these physical processes cannot coexist or that all legacy heat is a
duplicate photon absorption term. No literature coefficient is changed here.
