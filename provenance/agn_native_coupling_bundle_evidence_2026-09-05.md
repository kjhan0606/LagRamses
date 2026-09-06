# AGN native coupling — implementation evidence

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Status: implemented; initial Opus CONDITIONAL PASS followed by bounded
safety repairs and Opus closure PASS for B1–B3. No production qualification.
Authority: operator approved the volume-weighted jet correction and donor-gas
fallback velocity with "적용하고 진행해", following the Fable conditional plan review.

## Native changes

1. `agn_feedback_deposition.f90` provides a pure cell depositor and resolved
   jet delta builder. The depositor stages the five conserved fields, limits
   only NEW internal energy, returns excess as integrated deferred energy,
   and rejects invalid input without modifying the row. The cap is not a
   cooling operation on pre-existing hot gas.
2. `sink_particle.kjhan.f90:AGN_blast` calls this exact routine for thermal,
   saved replay, directed jet, zero-axis thermal and single-cell fallback
   deposition. It commits total metal only after the staged row passes.
   Fallback excess uses `vol_blast`, never the final AMR loop's `vol_loc`.
   Fallback writes are serial because several AGN may share a target cell.
3. `average_AGN` sums the Gaussian kernel times receiver volume; the resolved
   density increment is `m_load*weight/sum(weight*volume)`. The jet speed is
   `sqrt(2*f_ekAGN*EAGN/m_load)` so the allocated kinetic budget corresponds
   to the full returned mass. Gaussian shape and nominal feedback efficiency
   parameters are retained, but effective mass/energy injection changes from
   the old nonconserving result. This is NOT calibration-neutral.
4. Both membership passes clamp the radial square-root radicand against tiny
   negative roundoff. On the exact midplane the two opposed lobe contributions
   have zero directed net momentum but retain their kinetic second moment.
   Zero-axis/no-loaded-mass cases thermalize the existing EAGN budget rather
   than losing the nominal kinetic fraction; zero-axis fallback locates the
   donor cell too.
5. `average_AGN` captures donor velocity before mass removal in `vloadAGN`;
   the unresolved fallback restores mass, its original momentum and kinetic
   energy using that velocity, not the BH velocity. This velocity is local
   to the fallback owner (the same rank has its `ind_blast`); the resolved
   jet continues to use its existing BH-frame convention.
6. `read_params.jaehyun.f90` checks the SNRT enable latch only AFTER namelist
   reading and I/O-token release, and before simulation initialization. MPI
   min/max checks require all ranks to agree on the requested mode. Requested
   live SNRT plus `sink && sink_AGN` is rejected, with no override. Both source
   dispatchers have a defensive conflict check. They share the single
   process-lifetime `snrt_agn_rt_requested()` latch in the efficiency module.
7. The SNRT-only driver no longer clears its `idsink`-keyed accounted inflow
   at coarse-step changes. Its cumulative supply has no legacy reset owner
   in that mode. The existing ID-remap path preserves the cursor; new IDs
   begin at zero. The successful all-group transaction remains the only
   source-commit advance. This does NOT implement persistent restart or
   cross-rank cursor migration.

The Makefile links the deposition module and declares direct dependencies.
The existing native AGN runner is extended; no new Python validator, runner,
schema, GPU sink implementation or intermediate external audit was added.

## Reproduction and evidence

Before correction, direct arithmetic from the actual caller gives:

- uniform two-cell jet: mass removed 1, mass returned 0.5;
- cap crossing: old energy 5 + input 8 = 13, but old code gives gas 15 +
  deferred 3 = 18 (unit volume).

Executed from the project root:

```text
bash simulation/snrt/tests/run_fp15_agn_efficiency.sh
python3 simulation/snrt/tests/agn_effective_efficiency.py
git diff --check
```

Results: `SNRT_AGN_EFFICIENCY_OK`, `SNRT_AGN_SOURCE_OK`,
`AGN_NATIVE_CELL_COUPLING_SMOKE_OK`, `SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK`,
`AGN_EFFECTIVE_EFFICIENCY_TEST_OK`; whitespace check clean. Existing Python
regression was run unchanged; it is not new evidence for the full live path.

The native smoke covers below/crossing/already-hot cap, nonunit volume,
gas-plus-deferred closure, donor-velocity fallback, unequal-volume cylinder
normalization and midplane energy, invalid-input nonmutation and idle ratio.
The same production depositor used by AGN_blast is compiled, not an oracle
copy or a duplicated Python formula.

The identical production depositor and smoke were also compiled and run at
`-O3` in `build/agn_native_coupling`; all assertions, including invalid-input
nonmutation, passed. Correction after the end audit: `amr/end.f90` is NOT the
linked `clean_stop` provider. VPATH selects `patch/cuRamses/update_time.f90`,
which calls `MPI_FINALIZE`. The repair therefore globally reduces AGN error
decisions before any call to it; the earlier MPI_ABORT inference was wrong.

Changed production caller compilation was performed under
`build/agn_native_coupling`, using freshly compiled deposition/efficiency
modules and the existing RAMSES module interfaces under `bin/`. Compiler:
`mpiifx -qopenmp -fpp -O0 -g`, with `NDIM=3`, `NVAR=18`, `NVECTOR=500`,
`NPRE=8`, `NENER=0`, `SOLVER=hydro`, `LONGINT`, `QUADHILBERT`,
`OUTPUT_PARTICLE_POTENTIAL`, `USE_FFTW`, `PHASE0_STELLAR_ENRICHMENT`.

- `sink_particle.kjhan.f90`: CPU and `-DSNRT` compile pass.
- `read_params.jaehyun.f90`: CPU and `-DSNRT` compile pass.
- `snrt_ramses_driver.f90`: `-DSNRT` compile pass.

No full RAMSES link or live evolution was performed. The shared executable
was not overwritten. No commit, push or source-model activation was done.

## Claim limits / next physical decision

This bundle fixes native deposition arithmetic and installs an interim
exclusivity boundary; it is NOT completed simultaneous radiation+mechanical
AGN physics, nor production/publication qualification of legacy AGN.

- Resolved-jet donor-to-BH-frame recoil/bulk-energy accounting, imbalance
  between discretely sampled lobes, and non-metal passive-scalar entrainment
  are not closed by the local cell identity or the fallback repair.
- Jet mass-source ownership across full AMR/MPI evolution remains outside
  these focused synthetic tests; GPU sink code is inert in this active patch.
- No common time-integrated radiative/mechanical/spin partition, new AGN SED,
  obscuration or escape prescription has been adopted. An explicit model
  decision is still required to remove the simultaneous-source exclusion.
- The existing supplied-inflow photon convention remains an upper-bound
  engineering convention, not actual retained-mass equality.
- Rank-local cursors are not durable; migration/restart coupling remains
  unqualified and must not be inferred from fixing the coarse-step reset.

These are AGN physics/live-coupling follow-ups, not a reason to reopen parked
AGB raw data, generic HDF5 payload audits or broad AMR infrastructure work.

## Opus conditions and bounded repair disposition

Full initial review: `opus5_agn_native_coupling_end_audit_2026-09-05.md`.
The reviewer confirmed the approved native corrections, caller wiring and
claim limits, then identified three safety conditions. Driver response:

### B1: fail unsupported history ownership closed, do not guess a rebase

Accepted the defect. A zero cursor cannot identify which of a restored or
migrated cumulative ledger has already radiated. Seeding from the current
ledger, however, can silently discard genuine first-interval accretion.
No persistent/migratable cursor is implemented in this bounded repair.
Instead, preflight rejects requested live SNRT with sinks on a restart or
with more than one MPI rank. The driver repeats the check for nonempty sink
arrays before any spectral/state/source mutation. Serial fresh-start SNRT
keeps the coarse-step carry fix; MPI RT without AGN sinks is not disabled.
There is no bypass. Durable history/merger/migration qualification is still
future AGN source work, not a claim of this bundle.

### B2: global error consensus before the actual collective stop

Accepted after confirming the selected `update_time.o` provider. Added
MPI_MAX consensus for all three newly introduced rank-local error paths:
parent AGN budget/config validation, precomputed EAGN rejection (after the
entire local loop, including ranks with zero nearby AGN), and cell/deferred
energy rejection after deposition. No conditional collective sits inside a
rank-local AGN loop. All ranks then enter clean_stop together.

### B3: reject invalid incoming gas explicitly; do not discard its ejecta

Accepted the need to distinguish bad receiver state from bad source input,
but **did not adopt the proposed skip-and-defer-only-energy behaviour**.
Jet mass and momentum may already have been removed from the donor in
`average_AGN`; skipping the receiver while deferring only energy loses that
returned material. Continuing with negative gas internal energy or inventing
a positivity repair is also not an approved AGN prescription.

The explicit policy is strict receiver validation: error code 2 identifies
nonfinite/nonpositive density or materially negative incoming internal energy;
code 1 denotes source/budget/argument rejection. A failed row is unchanged,
and the global stop from B2 prevents continued evolution of the partial event.
The smoke now asserts this negative-internal-energy no-mutation case. This is
a deliberate fail-closed hydro precondition, not legacy numerical tolerance
parity; an approved upstream positivity treatment would be separate work.

### Other review notes

The normalization correction increases both intended returned mass/kinetic
coupling and the scale of the existing resolved-jet passive-scalar mixing
and donor/BH-frame momentum errors. Those are NOT rendered harmless by a
local energy identity; production legacy-jet qualification remains open.
The exact net recoil error is `m_load*v_BH - m_load*v_donor` plus sampled-lobe
imbalance, so its magnitude need not scale uniformly by N in every flow.
The donor still loses the full loaded mass in both old and corrected paths.

Idle zero-Eddington branch selection is now finite zero instead of NaN;
this is an explicit semantic cleanup, not a new radiative efficiency. The
active build requires OpenMP; non-OpenMP support and inert CUDA parity are
not claimed. Unused variables in the rewritten OpenMP loop's private list
were removed. No new audit is requested for each repair: one closure review
covers this entire set and the B3 policy response.

## Closure outcome

The single bounded Opus 5 closure review returned **PASS (B1–B3 only)**:
`opus5_agn_native_coupling_closure_audit_2026-09-05.md`. It independently
accepted unsupported-history quarantine instead of guessed cursor rebasing,
confirmed all three error-consensus collectives, and accepted strict receiver
rejection because skip-and-defer-energy would lose already-removed ejecta.
The numerical tolerance and frequency of invalid receivers in live nuclear
gas remain unqualified; a closure PASS is not a production robustness claim.

After these repairs, the existing O0/check-all native runner, the identical
O3 cell smoke, both CPU and SNRT caller compiles, and the SNRT driver compile
all passed again. No full link or evolution was added. The review's optional
source-error-code assertion refinement is nonblocking and not a new gate.

This approved implementation bundle is complete within its declared scope.
Live SNRT with AGN sinks remains restricted to serial fresh starts, and
simultaneous legacy AGN feedback plus live SNRT remains rejected. The next
physical/live-coupling bundle requires a scoped plan and operator approval;
none was started by this closure.
