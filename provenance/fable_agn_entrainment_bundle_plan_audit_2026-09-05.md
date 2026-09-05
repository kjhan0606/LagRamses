# Fable — AGN entrainment bundle plan review

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Requested model: `fable`; reported reviewer: `claude-fable-5-1`
Session: `78bb6ce9-4e47-4449-94d6-66c5a3679c93`
Duration: 278.8 seconds
Read-only Read/Grep/Glob; no subagents, permission denials or tool errors.
Verdict: CONDITIONAL APPROVE on the original proposed plan.
Driver disposition and amendments are in
`agn_entrainment_bundle_plan_2026-09-05.md`; they are not a second auditor verdict.

## Verbatim reviewer result

# Audit: `agn_entrainment_bundle_plan_2026-09-05.md`

## Verdict: **CONDITIONAL APPROVE**

Right next bundle, feasible, not overinstrumented. Three plan items need correction before operator sign-off; one of them would otherwise abort the production AGN configuration.

## 1. Necessity — confirmed against code

The three defects the plan targets are real, present at `e427029`, and are exactly the open residuals of the previous bundle (`agn_native_coupling_bundle_evidence:104-124`; Opus closure explicitly declined to re-raise them):

- **Frame mismatch**: `average_AGN:6717-6735` removes `m_load` at donor velocity `(u,v,w)` and captures it in `vloadAGN`; `AGN_blast:7089` passes `vAGN` (BH velocity) to `agn_jet_delta` for the resolved jet. Only the single-cell fallback (`:7142`) uses `vloadAGN`.
- **Lobe imbalance**: one `psy_norm` over the whole cylinder (`:6659`, `:7087-7090`); nothing enforces equal integrated lobe masses.
- **Composition**: only `imetal` withdrawn/returned (`:6725-6729`, `:7179-7187`); `ichem` fields untouched while density changes → abundance mis-mixing.

This is native feedback physics, not tooling. Scope (one helper extension, caller wiring, extend existing smoke + runner, one evidence doc, one end review) matches the operator's cadence.

## 2. Feasibility — checks I performed

- **Depositor contract survives two-lobe input.** `agn_deposit_cell:64` requires `energy_delta ≥ |dp|²/(2dm)`. For a cell receiving both lobes, `dK = (dm⁺|v⁺|²+dm⁻|v⁻|²)/2 ≥ |dm⁺v⁺+dm⁻v⁻|²/(2(dm⁺+dm⁻))` by Cauchy–Schwarz, and stream mixing cannot lower receiver internal energy, so `:65` also passes. No depositor change needed; the plan's "preserve second moment" is the right instruction.
- **Global closure** `Δ(gas energy)+deferred = EAGN` is correct under cold loading: donor loses `½m|v_d|²`, receivers get `½m|v_d|² + f_ek·EAGN` kinetic plus `(1−f_ek)EAGN` thermal (unchanged volume weighting via `p_gas`).
- **MPI path**: only five `nsink` arrays are reduced (`:6771-6783`); `vloadAGN` is never reduced, so remote receivers currently see zero — the plan's item 2 diagnosis is accurate. Extending the packed `MPI_SUM` buffer with 3 velocity components, `W±`, an integer owner count, and `nelt` fractions is mechanical. Owner-uniqueness check is justified: `dr_cell.le.dx_loc/2` (`:6661`) can tie on a cell face.
- **Fallback classification** on globally reduced `W⁺>0 .and. W⁻>0` is rank-uniform; replacing the `vol_gas>0` gate (`:7132`) and `psy_norm>0` (`:6990,:7087`) is coherent. The existing consensus collectives (B2) are reused unchanged.

## 3. Cold-loading convention (§3) — honestly scoped, one omission

The description is accurate: `:6723/:6735` preserve internal-energy *density*; removed mass carries zero thermal energy. It matches legacy stellar loading (`feedback.kjhan3:1836,1850`), so it is project-consistent, not ad hoc. **Missing statement**: because ρ drops while `e_int` density is held, donor specific internal energy (temperature) rises by `ρ/(ρ−m_load/V)`, up to 4/3 at the `0.25ρV` cap (`:6724`), and this heating bypasses the `T2maxAGN` limiter (donor removal never passes through `agn_deposit_cell`). The operator is consenting to that side effect; the plan must say so. Also §3's "do not also return donor thermal energy that was never withdrawn" is correct and prevents a double count.

## 4. Correctness defects in the plan (must fix)

**D1 — Item 4 refusal would abort the production AGN config.** The only production-like AGN namelist (`simulation/snrt/runs/p7_sinkprops_smoke/cosmo.nml:76,78,92,93`) runs `delayed_cooling=.true.` with `sink_AGN`, `spin_bh`, `mad_jet`. "Refuse AGN mass loading in that configuration before mutation" makes the bundle infeasible for exactly the run it serves. The established legacy rule already exists: stellar mass loading withdraws `imetal` and `ichem+ielt-1` only and leaves `idelay` untouched (`feedback.kjhan3:1827-1852`); `feedbk:620` treats `idelay` as a SN-mass reservoir, not a species. Replace refusal with **legacy parity: `idelay/ivirial/ixion/isgs` untouched, documented as a limitation**. Keep a refusal only for `NENER>0` (nonthermal energy inside `row(5)` would corrupt the 5-field `E−KE` reconstruction); production builds have `NENER=0` (`bin/Makefile:25`), so nothing active is blocked.

**D2 — "existing enabled element-density map" is ambiguous.** Two conventions exist: legacy `nelt` (read from the yield table, `init_part.f90:151`, loop `ichem..ichem+nelt-1` in `feedback.kjhan3:1841-1844`) and Phase 0 `runtime_field_map%element_index` (`stellar_ramses_runtime.f90:202-208`, populated only when Phase 0 initializes; `active_element` defaults all-true, `n_stellar_elements=11`). `sink_particle.kjhan.f90` is the legacy path; bind to **legacy `nelt`/`ichem` with the `ichem+nelt-1 ≤ nvar` guard**, and state that the Phase 0 map is not consulted. `nelt=0` (no table) yields an empty loop, so no special case. The `ichem` index is always defined (`read_hydro_params:455`) even when it exceeds `nvar`, so the guard is not optional.

**D3 — §3 donor-heating disclosure** (above).

## 5. Physical choices requiring operator consent (not defects)

- §1 donor-frame injection changes the resolved jet model relative to Dubois-style legacy (jet launched relative to entrained gas, not BH; with nonzero `v_avgptr` the lobes no longer track the BH). Not calibration-neutral, though it removes the `m_load(v_BH−v_d)` recoil error.
- §2 separate lobe normalization with all-or-nothing fallback (no half events).
- §3 retaining cold loading with the uncapped donor heating side effect.

## 6. Recommended cuts (nonblocking)

- Drop "any dust density" from item 4: no dust hydro field exists in `patch/lagRamses` (no `idust/ndust/dust_density` anywhere). State "no dust field; nothing to transport."
- Do not build a two-rank MPI harness. Make pack/allreduce/unpack a pure helper and test pack→sum→unpack serially; the runner (`run_fp15_agn_efficiency.sh`) stays `mpirun`-free.
- "Physical composition bounds": limit to finite, `≥0`, `Z≤1`; do not assert element fractions sum to one (legacy fields carry no such invariant).
- Since `agn_feedback_deposition_smoke.f90:54` is touched anyway, apply the Opus nit (`ierr==agn_deposit_invalid_source`).
- Keep the plan's abundance-change assertion; global species sums alone would pass the old bug.

## 7. Blockers vs. future physics

- **Blockers**: D1 and D2 (plan text corrections before implementation); operator consent on §1–§3. No missing prerequisite in code: depositor, `vloadAGN` capture, consensus stops, and the packed reduction path all exist.
- **Nonblocking future physics** (correctly excluded): enthalpy advection / donor thermal withdrawal, BH recoil, RT+mechanical common energy convention, `NENER>0` support, `idelay` entrainment semantics, live-evolution qualification.

**Verdict: CONDITIONAL APPROVE** — proceed once item 4 is rewritten to legacy parity (refuse only `NENER>0`), item 3 binds to the legacy `nelt/ichem` convention with the NVAR guard, and §3 discloses the donor-heating side effect; then request operator consent on §1–§3.
