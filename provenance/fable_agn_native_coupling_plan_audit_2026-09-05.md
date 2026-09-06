# Fable AGN native coupling plan audit

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Requested model: `fable` via `claude -p --model fable`
Reported audit model: `claude-fable-5-1`
Session: `ab33d0ba-8699-451f-a629-dd8bf572e38d`
Execution: read-only Read/Grep/Glob, no Bash, edits, jobs or subagents;
exit 0, no permission denials; duration approximately 414 seconds.
Audited plan: `agn_native_coupling_bundle_plan_2026-09-05.md`

## Auditor response (verbatim)

**Verdict: CONDITIONAL APPROVE.** The bundle is the right next step, the four listed defects are real, and the scope stays inside legacy feedback plus one ownership boundary. It needs one physics decision the plan does not know about, and two ownership repairs before the SNRT-only mode it promises is safe.

I confirmed each claimed defect against the code. In the thermal cap-crossing branches the cell receives the full cap energy plus the deferred excess, so the injected total exceeds the intended increment by the cell's pre-existing internal energy. The single-cell fallbacks save deferred energy with the stale level-loop volume rather than the receiving cell volume. The jet direction is stale when the cell sits exactly on the midplane, and the radial radicand can go slightly negative. The unresolved-jet fallback returns density with no momentum. The proposed cell identity holds in every branch, and the deferred term is provably nonnegative because merging momentum into a moving cell can only lose kinetic energy, so the routine can assert that as an invariant.

## Required physics decisions

1. **Resolved-jet kernel normalization.** The jet weight is normalized to unit sum over cells, but the density increment also divides by the total cylinder volume. Summed over the cylinder, the deposited mass, momentum and kinetic energy are therefore a fraction of what was removed from the host cell:
   deposited / removed ≈ 1 / N_cylinder_cells for uniform cells.
   The host cell loses the full loaded mass in the averaging pass, so this is a net mass and energy loss, not just a shape choice. The CUDA copy carries the same formula. Decide whether this bundle corrects the weight to a volume-weighted normalization or freezes legacy behaviour. I recommend correcting it, because per-cell conservation is meaningless while the kernel sum is off by this factor, but it changes jet strength relative to every prior run, so it is your call.
2. **Fallback mass-return velocity.** Returning the mass at the black hole velocity is consistent with the resolved branch, but the removal took momentum at the cell's own velocity, so the remove/return pair is not closed. Returning at the current cell velocity undoes the removal exactly with no new momentum source. I recommend that. Either is defensible if stated.
3. **Midplane treatment.** Zero net directed momentum with the opposed-stream kinetic budget thermalized is defensible. Implement it as an equal split of the cell's share between the two lobes rather than a new rule, and apply the same membership logic in the averaging pass so the normalization sees the same cells.
4. **Cap headroom fix is a behaviour change.** The plan says it does not change calibration. It reduces injected energy only in cap-crossing cells. Record it as a conservation correction rather than claiming neutrality.

## Implementation repairs

- **Coarse ledger reset ownership.** The supply ledgers are zeroed only inside the legacy feedback routine, while the SNRT driver resets its accounted-inflow cursor at every coarse-step boundary. With legacy AGN feedback off, the driver re-emits the entire cumulative supply every coarse step. The plan advertises that mode as available, so the bundle must either carry the cursor across coarse steps when no legacy reset occurs or fail that mode closed too.
- **Exclusivity gate placement.** Put one fail-closed check in the driver at the enable latch, before any deposit, using the namelist sink flags, which are rank-uniform. No namelist override. The existing enable latch reads a per-rank environment variable, so keep the decision off that value.
- **Standalone module.** The cell routine must depend only on the parameters module and take gamma, temperature scale, cap and volume as arguments, so the existing F-P1.5 runner can compile it exactly as it compiles the efficiency helper. Add the object and the sink-particle prerequisite to the Makefile. Return deferred energy and an error code, and let the caller accumulate into the thread-private saved-energy array.
- **Apply the radicand clamp in both passes**, since the averaging pass uses the same square root to decide membership and weight.
- **GPU path is inert.** The lagRamses sink code never dispatches to the CUDA sink kernels, so the CPU routine is the only live path. Do not touch cuRamses, and do not ask the end audit for GPU parity.

## Gates to cut

- Drop the runtime ownership negative. The production negative runner exits at the stellar fate gate before any AGN step, so the exclusivity gate is unreachable there. Compile plus end-audit reading is enough.
- Fold the replay case into the thermal case. It is the same arithmetic. Fold invalid input into one error-code assertion. Keep cap below/crossing/already-hot, nonunit fallback volume, moving mass return, opposed jets with midplane, and add one kernel-sum check over a synthetic cylinder, which is what exposes decision 1.

Everything else in packages A through C is proportionate. I did not write a plan file since this was a read-only audit.
