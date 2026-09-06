# Opus 5 — AGN native coupling bounded closure review

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Requested/reported review model: `claude-opus-5`
Session: `9597c2b7-2faf-4dc7-b88c-a9d8e9fc52f0`
Duration: 177.1 seconds
Read-only review; no subagents spawned; no permission denials.
Scope: closure of B1–B3 from the initial bundle-end review only.
Verdict: **PASS within that bounded scope**, not production qualification.

## Verbatim reviewer result

## Verdict: PASS (bounded — closes B1–B3 only)

I read the disposition (`agn_native_coupling_bundle_evidence_2026-09-05.md:126-186`) and inspected the changed code myself.

### B1 — closed by quarantine, correctly
`read_params.jaehyun.f90:1705-1712` rejects requested live SNRT with `sink` under `ncpu>1` or `nrestart>0`; the operands (`snrt_requested_max` after MIN/MAX allreduce, `sink`, `ncpu`, `nrestart`) are all rank-uniform, so `nml_ok=.false.` → `clean_stop` at `:1718` is collective-safe. `snrt_ramses_driver.f90:286-289` repeats it on `nsink>0` before the spectral contract load at `:290` and before any state/source mutation; `enabled` is latched and uniform, `nsink` is global, so that `clean_stop` is uniform too. No bypass path exists — no env override, no `if(myid==1)` guard on the condition itself.

The choice not to rebase is right. The auditor's proposed `retained_initialized`-style seed is safe for the *retained* cursor (it is a monotone high-water mark) but not for the inflow cursor, where seeding from `supplied_mass` discards genuine unradiated accretion — an unauthorized silent loss, symmetric to the burst it prevents. Neither branch is derivable from a zero cursor, so refusing the mode is the only honest option.

I checked the residual inside the permitted envelope: sinks created mid-run start with cursor and `dMBH_coarse` both at zero, and `merge_sink` does not transfer `dMBH_coarse` to a survivor (all 14 occurrences in `sink_particle.kjhan.f90` are dumps, grow_bondi accumulation at `:4131-4187`, the AGN_feedback read at `:6275`, and the now-unreachable reset at `:6360`). So serial fresh start is genuinely consistent, not merely undeclared.

### B2 — closed, and the driver's provider correction is right
`clean_stop` is `MPI_FINALIZE` + `stop` (`patch/cuRamses/update_time.f90:351-365`), collective and with no output dump — so a rank-local trip really would have hung, and a global reduction is the necessary fix, not `MPI_ABORT`. All three new paths reduce first: `:6298-6305` (budget/config, before `average_AGN` removes any donor gas), `:6997-7004` (EAGN rejection, after the full local `nAGN` loop, which `cycle`s rather than exiting), `:7151-7159` (cell/deferred, after the level loop, the `PEsaveAGN` gather and the serial fallback loop). Each `MPI_ALLREDUCE` sits at subroutine scope under `#ifndef WITHOUTMPI`, outside every OpenMP region and every `nAGN` conditional; the in-loop accumulation uses `max` and `reduction(max:deposition_error)` (`:7032`). `AGN_feedback`'s early returns (`:6185-6187`) are globally uniform and `AGN_blast` has no early return before the reduces, so no rank can skip a collective.

### B3 — the strict precondition is acceptable as an engineering default
The available options were: abort; skip-and-defer-energy-only, which drops `m_load` and `dm·v_donor` already removed in `average_AGN:6713-6727`; or new defer-whole-event / positivity machinery. Only energy has a carry vehicle (`EsaveAGN`) — mass and momentum have none — so the auditor's option is a silent conservation violation, and the third requires an unapproved prescription. Fail-closed is the correct residual.

The implementation matches the claim: code 2 is set only for nonfinite/nonpositive `row` (`agn_feedback_deposition.f90:32-34`) and materially negative `internal_old` (`:60-63`); everything after `:35` reports code 1; no failing path writes `row` (assignment is at `:77`, after every check). The smoke asserts code 2 plus bitwise non-mutation (`agn_feedback_deposition_smoke.f90:55-58`; that row gives `internal_old = -3.5`). The rejection is also not a roundoff hair-trigger: `tol` at `:59` scales with `max(|E_tot|, |E_trial|, E_kin,input)`, so a supersonic cell's `row(5)-kinetic_old` cancellation noise is covered.

**What this acceptance does not settle.** Accepting the precondition is accepting *fail-closed on invalid receiver state*, nothing more. It carries no judgment on: the value or form of `tol`; whether production-resolution nuclear gas actually produces such cells at a rate that makes aborting untenable (no live evolution was run, so this is untested); or which eventual policy — upstream positivity floor, whole-event deferral with donor restoration, or legacy tolerance parity — should replace it. Any of those is a separate physics/numerics decision with its own authorization.

### Nit (not a condition)
The source side of the new discrimination is untested: `:52-54` asserts `ierr/=0`, not `==agn_deposit_invalid_source`. One token change if the file is touched again.

Production legacy/RT qualification, passive-scalar entrainment and donor/BH-frame recoil remain open exactly as `:104-124` states; I am not re-raising them.
