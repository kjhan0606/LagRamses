# Opus 5 bundle-end audit — AGN partition reference

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: `claude-opus-5` via read-only CLI session
Scope: final implementation of the approved `partition_reference_v1` bundle; physical prescription, units, event partition, reservoir lifecycle, caller wiring, one-event priority, legacy/MAD isolation, admission, evidence strength, and overclaiming.

## Verdict

**CONDITIONAL PASS.** Opus verified that the approved energy partition is really
wired into the Fortran production path: high state is 85% EM plus 15% heat,
low state is 100% jet, BH retention is unchanged, event-time channel selection
is preserved, deferred replay is thermal, and MAD remains outside the admitted
reference branch. The consolidated native/build gate was considered useful but
not sufficient for a production-path-complete claim.

## Required findings and dispositions

1. **Jet kinematics was unbounded and under-disclosed.** Using retained mass as
   the low-state loading entitlement with the existing Newtonian kinetic
   deposition can produce `v_jet >= c` for high spin efficiencies. The repair
   adds `agn_reference_jet_speed_ok` and a named `0.9c` upper bound. The event is
   rejected before hydro mutation if the bound is violated; energy is not
   silently capped or reassigned. The reference evidence and scope now state
   that this is a Newtonian comparison model, not a relativistic correction.

2. **Deferred replay had a deposit/commit mismatch.** The old unconditional
   `ind_blast<=0` clearing could allow replay deposition while retaining the
   full deferred liability. The repair adds one receiver predicate:
   replay requires positive thermal-bubble `vol_gas`; fresh heat/jet events
   require the unique donor/receiver cell. The same predicate gates both the
   `AGN_blast` call and pending-ledger commit. A smoke assertion covers replay
   with and without a receiver.

3. **Border-list MPI propagation was incomplete.** The reference branch now
   publishes committed pending channels by sink slot using touched-mask plus
   `MPI_MAX` consensus, so ranks outside a sink's local border list do not keep
   stale reference channels. The profile remains serial/fresh-start admitted;
   this repair removes the identified stale-rank write for future qualification,
   but does not claim multi-rank production validation.

4. **Jet no-direction starvation is fail-closed.** A selected reference jet
   without a valid direction now aborts the reference event collectively rather
   than leaving an indefinitely dominant pending jet channel while accretion
   continues. No silent thermal fallback is introduced.

5. **Evidence was strengthened and claims narrowed.** The native smoke now
   includes legacy thermal/non-MAD/MAD recipe parity, expanded admission
   negatives, receiver predicate checks, and the jet-speed guard. The aggregate
   gate records hashes for all changed AGN reference sources. Its production
   link uses a temporary `EXEC` path, so it no longer overwrites the shared
   `bin/ramses_final3d` binary. The caller path is compiled and statically
   checked; it is still not a live RAMSES evolution or a full caller-execution
   test.

6. **Minor cleanup.** `average_AGN` receives only reference mode/loading
   arguments and `AGN_blast` receives only mode/energy arguments. The reference
   heat trigger uses `mass_gas` directly, and partial optional-argument calls
   fail closed. The existing `f_bondi`-scaled chi convention remains explicit
   in the evidence.

## Remaining bounded scope

The reference profile is still comparison-only and is not promoted to a
physical AGN SED, MAD/spin rotational-energy closure, restart-capable model,
live multi-source qualification, or cosmological production configuration.
The aggregate gate is a compile/native contract result, not a live evolution.
