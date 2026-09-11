# Selective PAH endpoint review — 2026-09-11

CONDITIONAL plan review, not live completion. Driver accepts per-group
admission, exact C/C+ product insertion, single-use captures, particle heat
capacity, material residual and restart identity conditions. For source
formation the driver selected the plan's alternative atomic-donor/local
binding-release convention after the review was commissioned; the review's
precondensed source receipt recommendation is not applied simultaneously.
This explicit change conserves atomic-source energy via +A12*Ninject heat
and -A12*Ninject derived PAH binding, without a new escape ledger. Defaults
and initial existing-cell thermal energy are unchanged.

Only the current869.634-eV fixed-group closure exercises hard atomization
live. This is not broadband stellar/AGN PAH survival. The rejected other
occupied groups remain an explicit scope limitation.

## Original read-only review

Review complete. Verdict and indispensable items follow. The plan file and ExitPlanMode tools are not available in this session, so the review is delivered here.

## Verdict: CONDITIONAL

**Q-GOAL: useful, but as a component qualification, not a simulation-ready channel.** It adds the first carbon-skeleton loss path (currently absent per the header of `patch/lagRamses/dust_pah_live_model.f90`) and the exact C/C+/H return plus chemical-reference bookkeeping that any later yield model needs. The limit is the group binding: CHIMES fixes the nine group edges and means in `patch/lagRamses/snrt_chimes_bridge.c:47-49`. Under the monochromatic closure the only admissible hard event energy is the 500 to 2000 eV group at 869.634 eV. Groups at 17.66, 34.38, 106.64 and 4023.6 eV stay rejected when occupied in PAH-bearing cells (`snrt_ramses_driver.f90:1477-1481`, `dust_pah_mixed.f90:143`). Any stellar EUV or AGN spectrum occupies them, so the integrated run still aborts. The plan says this; it should also say that "208 to 1239.84 eV" collapses to one closure energy, and the worked example should be 869.634 eV, not 250 eV.

**Q-LEAN: no gratuitous machinery.** Pure helper, no new advected scalar, no daughter opacity, no pooling. The 298 K thermal correction is a constant, keep it. Optional: reuse the existing event and cost counters printed at `snrt_ramses_driver.f90:2125-2128` instead of new diagnostics.

## What I verified against source

- Constants match code: B(12)=0, B(0)=48, B(13)=-3.2 eV (`dust_pah_hydrogen.f90:346-349`); IP 7.02 eV in levels (`dust_pah_live_model.f90:147,161`); C, H, C+ enthalpies (`snrt_chimes_atomization.h:10-11`). A12, I_C, the 207.827 eV worst-state threshold and both 250 eV excess values recompute correctly by hand.
- Uvib(298 K): hand sum over the 102 DL01 modes gives about 0.31 eV, consistent with the stated 0.3068 eV to 2 percent. Not recomputed exactly.
- Energy closure: `E_PAH = u + B + q·IP − A12` is consistent with how `pah_level` is built, and the residual `+A12·(N_old − N_next) + I_C·ΔC+` has the same sign structure as the existing H2 correction at `dust_pah_mixed.f90:306-307`.
- Optical support at the actual group energy: both retained tables give positive Qabs of about 4.0e-4 at 1.445e-3 micron (869 eV) in the first radius block, with smooth monotone neighbours. Group 9 is beyond the table and already masked (`dust_pah_live_model.f90:101-105`).
- Plan's premise holds: captures are energy divided by group mean energy (`dust_pah_mixed.f90:241-242`), the old receiver rejects above 13.6 eV (`dust_pah_radiation.f90:67`), and hard primary opacity is masked (`dust_pah_live_model.f90:109-119`). A third gate the plan omits sits at `dust_pah_hydrogen.f90:293`; all three must stay for the residual captures.

## Indispensable changes (conservation, not approximation)

1. **Charge defect unless products are inserted before reconciliation.** `snrt_chimes_reconcile_charged` (`snrt_chimes_bridge.c:634-642`, called from `snrt_chimes_runtime.f90:226-228`) tops up freed carbon as neutral CI and then sets electrons from solid charge plus ion charges. Destroying a cation lowers solid charge by one while C enters neutral, so one gas electron vanishes. The C+ must already be in the state passed as `previous_state` so the free-carbon increment is zero there.
2. **Per-group admission class, not a raised `pah_max_primary_ev`.** Raising the scalar at `dust_pah_live_model.f90:70,116-119` would give groups 5 to 7 real table cross sections and route them to the old receiver. Classify: soft to old receiver, the 869.6 eV group to atomization, all others sigma zero with the driver rejection at line 1481 unchanged. `pah_primary_alpha` then carries the hard opacity into tau automatically.
3. **Single-use photon budget.** Zero the atomization group in `captures(:,q)` before the calls at `dust_pah_mixed.f90:250` and `:254`; the receiver's balance at `dust_pah_hydrogen.f90:155` then closes without edits.
4. **Heat capacity.** Add `(24+h)` particles per event to `cv` in the pattern of `dust_pah_mixed.f90:252,259,268`, or the excess heat gives a wrong temperature.
5. **Residual and receipt.** Add the A12 and I_C residual at `dust_pah_mixed.f90:302-307`; add the `−A12·N_inject` receipt beside `stellar_ramses_runtime.f90:1154-1160` without changing the gas subtraction there.
6. **Mass and restart.** Use `dust_pah_state_mass` and `dust_pah_inventory` (`dust_mass_physics.f90:456-475`) for the phase mass and element gain; bind A12, I_C, and the admission class in `pah_live_identity` (`dust_pah_live_model.f90:182-208`).

**Model approximations, authorized, keep as declared:** unit yield; full local retention of photoelectron and Auger energy as gas heat, which at 870 eV is the largest approximation since both electrons would escape a 24-carbon molecule; monochromatic closure; H-independent hard cross sections; zero fragment kinetic energy. Note that at 869.6 eV every state has Y=1, so the u(j)-dependent threshold is exercised only by the pure-helper tests, not by any run.

This review verifies neither execution nor the exact Uvib value.
