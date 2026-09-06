# AGN accretion-powered partition reference — approved model, bundle plan

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, baseline `63603fb`.
Authority: operator explicitly approved "별도 비교 모형으로 구현 승인" after
being told the high-state mechanical share is 15%, low-state share 100%,
BH growth is unchanged, luminosity/mechanical normalization change, and MAD
is excluded only from this NEW model. Fable CONDITIONAL APPROVE; driver
amendments below govern implementation, now authorized and starting.
Review: `fable_agn_partition_reference_plan_audit_2026-09-06.md`.

## Final objective and why this package

Production/publication-ready native RT, stellar/AGN feedback and dust remain
the final goal. Two preceding bundles corrected mechanical deposition and
accepted radiative fuel. This bundle actually connects a common accretion
budget to both consumers. It is not another unused ledger or Python gate.
It adds an opt-in comparison model, not a replacement of legacy/MAD and not
approval of a physical AGN SED or cosmological production run.

The separate components of an accretion-powered model have precedent in
Steinborn et al. (2015), https://doi.org/10.1093/mnras/stv072 . This is NOT an
implementation of that paper's fitted efficiency functions. The explicitly
approved 15%/100% split is a controlled reference prescription; in particular
zero low-state radiation and the sharp transition are model limitations,
not observationally calibrated claims.

MAD has an additional BH rotational-energy source. Existing code includes
both a MAD efficiency polynomial and spin evolution, but their joint mass,
angular-momentum and injected-energy closure with RT is not established.
Spin-powered outflows can exceed accreted rest-mass energy (Tchekhovskoy et al.
2011, https://arxiv.org/abs/1108.0412). Do not cap MAD or pretend this
accretion-only partition qualifies it. MAD coupling is the explicit next
physical extension after this comparison baseline, not a discarded capability.

## Approved reference prescription

For EACH actually accepted gas withdrawal, with its contemporaneous epsilon
and code-to-physical mass unit:

```
dm_BH = (1-epsilon_event) * dm_gross
E_release = epsilon_event * dm_gross * mass_unit_event * c_exact^2
chi = instantaneous Bondi rate / Eddington rate
chi >= X_floor: E_heat = 0.15 E_release; E_EM = 0.85 E_release; E_jet = 0
chi <  X_floor: E_jet  = E_release;      E_EM = 0;             E_heat = 0
```

This is a release-energy partition, not 15% of the swallowed rest-mass energy.
Use the existing non-MAD retention efficiency; no BH growth or accretion-rate
recalibration. Relative to fixed-epsilon legacy high mode, mechanical energy
changes from `0.15*epsilon*dm_BH*c_approx^2` to
`0.15*epsilon*dm_gross*c_exact^2`, while EM energy decreases to 85% of release.
Do not silently repurpose `eAGN_T/eAGN_K`: those remain legacy coefficients;
reference shares are named constants of the versioned model. Equality belongs
to the high state. The luminosity discontinuity is deliberately visible in
the reference test, not smoothed without new physical approval.

Nine-group escaped fractions refer to E_EM, once. Omitted radiative energy is
not added to mechanical heat, and resolved photon absorption is not duplicated
by another epsilon_f heating term. Mechanical heat is the chosen non-radiative
release channel, not a second deposition of those same photons.

## One implementation package and actual dispatch

1. Add a process-lifetime, strict opt-in selector
   `SNRT_AGN_MODEL=partition_reference_v1` to the existing small AGN module.
   Empty/legacy leaves current behaviour. Unknown selection is an error, not
   a fallback. Admission requires SNRT build, RT enabled, hydro+sink+Bondi+
   sink_AGN, NDIM=3, NENER=0, mad_jet=false, serial fresh-start, valid X_floor,
   and explicit reference-control opt-in. Check rank agreement before any
   rank-dependent return. Do not allow this profile to masquerade as an
   approved production configuration. Existing legacy+RT exclusion remains
   for every other model; CPU/standalone legacy MAD are untouched.
2. Partition the real accepted-event energy before its existing thread/MPI
   reduction. Keep electromagnetic pending energy; add pending heat, jet and
   jet-associated retained loading mass. Event-time state selects the channel;
   pending heat must not turn into a jet when current chi later changes.
   Use physical erg for all energy reservoirs and code mass for loading.
   Creation starts at zero; real sink mergers sum every pending channel with
   the existing merger helper. No second cumulative cursor, event database or
   restart format. Existing serial/fresh-start restrictions fence these states.
3. Reuse `average_AGN`/`AGN_blast` and their conservative depositors. The
   reference dispatch selects ONE event per sink per coarse step using explicit
   channel identity, source energy and loading-mass arguments. Do NOT fabricate
   Bondi/Eddington rate arrays to trick legacy branch selection. Keep the legacy
   stage as its unchanged recipe and dispatch path. Common geometry, donor
   ownership, two-lobe normalization, cold loading, cap and fallback stay shared.
4. Keep separate pending channels and a physical-erg deferred mechanical
   liability. Replay is thermal and cannot recreate photons or loading mass.
   Replay precedes fresh events; if deferred heat remains, hold fresh channels
   pending (as existing Esave suppresses new accretion). If a fresh event fires,
   consume its requested channel once and transfer its undeposited energy to
   deferred liability, even when the temperature cap accepts none. Do not
   re-load gas on replay. Jet loading uses only retained mass accumulated in
   jet-mode receipts and clears it when that jet event is committed. Other
   channel energy/mass remains pending through mode changes.
5. Retain TAGN/jetfrac timing concepts, now evaluated against the corresponding
   channel's actual energy/loading mass, not a latest-mode reinterpretation of
   all dMsmbh. For a thermal trigger, use energy/mass, gamma and the existing
   temperature unit explicitly. Keep dMBH/dMEd coarse arrays as rate diagnostics;
   neither they nor current eps recompute energy. Existing Esave/dMsmbh legacy
   reset logic must not erase another reference channel. Accretion suppression
   in the reference path consults its deferred liability rather than stale
   legacy Esave. Report this as a distinct temporal prescription.
6. Bind the opt-in model to an explicitly named reference-control spectral
   contract describing fractions relative to the post-partition electromagnetic
   luminosity. Reuse existing group moments, identify their inherited source,
   do not invent a physical SED, and reject other model/contract combinations
   at preflight. Relax source-exclusivity checks only for this admitted profile.
   Actual SNRT source commit/transport rollback ordering stays unchanged.

## Evidence and stopping boundary

Extend the existing native runner and smokes. Exercise actual partition and
reservoir/deposition helpers: high/low/threshold events, varying epsilon,
alternating modes before consumption, merger and failed source retry, cap
defer/replay, and shared-donor overlap. Verify

```
E_release = E_EM + E_heat + E_jet
issued mechanical energy = accepted gas energy + deferred liability
issued EM energy = represented spectral energy + unrepresented energy
```

Track pending channels separately from issued energy; do not count pending
and emitted energy twice. Test that legacy selection retains its old output
and that MAD/unknown/unsupported joint configurations fail admission. Compile
the real changed callers and initialization in scratch, with no shared binary
overwrite. One evidence record and one Opus bundle-end audit cover the package.
Fable reviews necessity, feasibility and overinstrumentation once before edits.

No new Python framework, generic AMR/HDF5 audit, author contact, AGB research,
large job, production activation or new checkpoint support. A short live
multi-source case remains necessary before production/publication claims; the
native tests and compile must not be described as that run. A fully closed
relativistic BH+gas+gravity energy model is not claimed by this subgrid budget.

## Driver disposition of Fable's conditions (governing amendments)

- A1/A4: use a single existing averaging/deposition pair per coarse step,
  inside AGN_feedback. Per sink priority is replay if deferred>0, otherwise an
  eligible jet, otherwise eligible heat. Unselected channels remain pending.
  Heat may wait while jets remain eligible; no fairness/starvation claim is
  made for this reference scheduling rule. Never use or overwrite legacy
  Esave for reference liabilities. Accretion is held while either legacy
  Esave or the reference deferred reservoir is nonzero.
- A2: pass explicit replay/jet identity and actual requested energy/loading
  mass to the existing traversals; relocate the legacy energy recipe intact
  into the caller/helper. Only pass arguments each traversal actually uses:
  average needs mode/loading, blast needs mode/energy. No fabricated rates.
- A3: extend the existing event/thread/MPI pack by three components, not a
  new reduction service. chi uses the actual f_bondi-scaled rates in grow_bondi.
- A5: use the existing reference_control status plus explicit reference opt-in;
  no new spectral contract field/version. Its moments now multiply E_EM and
  therefore retain their electromagnetic bolometric meaning. A documented
  named reference configuration may reuse those moments. Validate reference
  admission at namelist preflight before accretion, not only at the first
  driver call; reject approved_production under this comparison model.
- A6: no cumulative issued counters or expanded source/dump schema. Budget
  closure assertions are native tests on the actual helper outputs.
- A7: heat trigger uses (gamma-1)*E_heat_code/mass_gas*scale_T2 > TAGN;
  jet trigger uses jet-channel retained mass/(BH mass-jet-channel retained
  mass) >= jetfrac, with positive pending jet energy. Include a nonzero-TAGN
  case and an eligible-channel priority/mode-alternation case in the smoke.
- A8: when a jet fires, its loading entitlement is consumed even if the donor
  cap permits less loaded gas; unused entitlement is not a future mass source.
  This matches legacy mass-loading reset semantics. Deferred ENERGY remains.

These amendments implement the approved comparison prescription without
promoting it to a physical MAD/SED model or adding another audit gate.

## Implementation status before bundle-end audit

Implementation is now wired in the production Fortran path. `SNRT_AGN_MODEL`
is latched and admitted only for the serial fresh-start,
`reference_control`/RT configuration; otherwise the existing legacy+RT
exclusion and legacy/MAD path remain in force. Accepted events are partitioned
at event time into the electromagnetic reservoir and three mechanical pending
channels. The channels are reduced through the existing thread/MPI pack,
initialized for new sinks, and merged through the existing sink merger. The
reference `AGN_feedback` branch selects replay → eligible jet → heat once per
sink/coarse step, passes explicit mode/energy/loading inputs to the common
averaging and blast traversals, consumes selected entitlements once, and keeps
deferred mechanical energy in physical erg without writing legacy `Esave`.

The final source was rebuilt with `SNRT=1 USE_CUDA=1`; the consolidated bundle
gate passed after the final loading-mass boundary check. The gate ran the
reference native smoke, production link/symbol checks, thermochemistry and
spectral negative cases, two-rank transaction smoke, CUDA multigroup checks,
production fail-closed checks, and `git diff --check`. The end-audit verdict is
recorded separately. Opus initially returned `CONDITIONAL PASS`; its jet-speed,
replay receiver, rank-consensus, starvation-guard, legacy-parity, admission,
source-hash, and temporary-binary findings were repaired or explicitly fenced
as documented in the implementation evidence. No claim of live evolution,
restart, multi-source MPI qualification, MAD closure, or physical SED is made
by this bundle.
