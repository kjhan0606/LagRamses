# AGN accretion-powered partition reference — implementation evidence

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Baseline: `63603fb3b4426b5410eee6224a12833df68550d8` (`main`), with unrelated
pre-existing worktree changes preserved. This bundle is not committed or
pushed yet.
Plan: [`agn_partition_reference_bundle_plan_2026-09-06.md`](agn_partition_reference_bundle_plan_2026-09-06.md)
Plan review: [`fable_agn_partition_reference_plan_audit_2026-09-06.md`](fable_agn_partition_reference_plan_audit_2026-09-06.md)
End audit: [`opus5_agn_partition_reference_bundle_end_audit_2026-09-06.md`](opus5_agn_partition_reference_bundle_end_audit_2026-09-06.md)
  — initial `CONDITIONAL PASS`; ranked repairs below were implemented and re-gated.

## Implemented production path

- `snrt_agn_efficiency.f90` adds the strict process-lifetime
  `SNRT_AGN_MODEL=partition_reference_v1` selector and narrow admission
  predicate. Unknown/rank-inconsistent selections fail closed. The reference
  profile requires SNRT+RT, hydro/sink/Bondi/AGN, 3D, `NENER=0`, non-MAD,
  serial fresh start, positive `X_floor`, and explicit `reference_control`
  spectral-contract status.
- `accrete_bondi` calls `agn_partition_release` on the accepted event before
  thread reduction. It records EM erg, heat erg, jet erg, and retained jet
  loading code mass separately. BH retention/growth stays the existing
  `(1-epsilon_event)*gross` result.
- `grow_bondi` extends the existing thread/MPI payload from `5+ndim` to
  `8+ndim`, and accumulates the three mechanical channels into
  `agn_mechanical_pending(1:3,:)`. Its accretion gate also holds while
  `agn_mechanical_pending(4,:)` contains deferred mechanical energy.
- `init_sink`, new-sink creation, and `merge_sink` initialize, sum, reorder,
  and tail-clear all four reference channels. No restart or dump schema is
  added; admission keeps this comparison mode serial/fresh only.
- `AGN_feedback` dispatches the reference mode once per sink/coarse step with
  priority replay → eligible jet → heat. It passes explicit mode, source
  energy, and loading entitlement to the existing `average_AGN` and
  `AGN_blast` traversals, so no fabricated Bondi/Eddington arrays or second
  traversal are used. It never writes legacy `Esave` in the reference branch;
  returned deferred code energy is converted to physical erg and kept in the
  reference deferred channel. A selected event without a valid receiver stays
  pending.
- Legacy and MAD energy selection remain on the original caller path, with
  the existing recipe relocated to the shared `agn_legacy_energy` helper.
  The reference branch is not a MAD or physical-SED implementation.
- The bundle-end repair adds a single replay/fresh receiver predicate, a
  `0.9c` Newtonian reference-jet speed guard that fails before hydro mutation,
  collective no-direction starvation rejection, and touched-sink MPI pending
  state propagation. `average_AGN` and `AGN_blast` now accept only the
  reference arguments each traversal uses; partial optional calls fail closed.

## Native checks

`simulation/snrt/tests/run_fp15_agn_efficiency.sh` passed after the final
boundary-condition repair. The reference smoke verifies:

1. high state: `[E_EM,E_heat,E_jet,m_load] = [0.85,0.15,0,0] E_release`
   for `chi >= X_floor`;
2. low state: `[0,0,1,retained] E_release`/code-mass loading for
   `chi < X_floor`;
3. replay priority and preservation of unselected channels;
4. jet selection from `m_load/(M_BH-m_load)` and rejection when loading is at
   or above the BH mass;
5. one-time jet entitlement consumption, heat trigger above/below `TAGN`,
   and the existing native overlap/deposition closure tests;
6. receiver availability and the `0.9c` jet-speed guard;
7. legacy thermal/non-MAD/MAD recipe parity;
8. narrow reference admission with build/feature/dimension/rank/restart/
   threshold/MAD negatives.

Marker: `AGN_REFERENCE_PARTITION_SMOKE_OK`.

## Consolidated bundle gate

Command: `simulation/snrt/tests/run_snrt_bundle_gate.sh`

Final run result: `SNRT_BUNDLE_GATE_PASS`.

| Stage | Result |
|---|---|
| `production_build` (`SNRT=1 USE_CUDA=1`, temporary `EXEC`) | PASS, 179.354 s |
| `agn_partition_reference` | PASS, 1.089 s |
| native symbol check | PASS, 4 symbols |
| thermochemistry loader | PASS; 4 negative cases |
| spectral contract loader | PASS; 10 cases |
| transaction MPI | PASS; 2 ranks, 11 rollback/config cases |
| CUDA multigroup | PASS; photon and species conservation markers |
| production negative | PASS; fail-closed marker |
| `git diff --check` | PASS |

Production link artifact (temporary and removed on gate exit):

```text
sha256 5dba03f58b784870e2b6351479bac07d191186d0c08113000ce33f1cb9cf90cb
```

Relevant final source hashes:

```text
patch/lagRamses/agn_feedback_deposition.f90 4c725773325bc62e53578c3c78afa1e1d3f310a5e61a0c8c4b99d6f5d439735b
patch/lagRamses/sink_particle.kjhan.f90      630a29c80fa8556dcb13b3067661f9925252a090007ae53209cbd3906e678e03
patch/lagRamses/pm_commons.f90               d7cbee1dc9fe4bcd98f58b125c1ae7e743fecc3a244a449c3c14a5a6481b74f7
patch/lagRamses/init_sink.f90                e972a308765a1c5e032511d7f8a485e0f8a3a03094f71ceef53177e94bbb1b43
patch/lagRamses/read_params.jaehyun.f90      76541ccdd23e57633408e2a620f9b1ceb8136e2f30162dce84732b6f6afc0666
patch/lagRamses/snrt_agn_efficiency.f90      09050e12bab1594ea568be2abe17c118f385be1a857f01a4b68cbce116d43a0a
patch/lagRamses/snrt_ramses_driver.f90       561a26af80a78e863f28900949cf13693f7c12799bf8d065f5ce2fdf9103a8c8
patch/lagRamses/snrt_rt_transaction.f90      c16f8aa476f7de3122f025414585c34f43e12f37e6a313a60ab631551d4ffa14
patch/lagRamses/snrt_transport_step.f90      6a818cf4fa7419f0b87122063b4eb687eb6344107e1c01c12b5f4e4167c2e628
patch/lagRamses/snrt_thermochemistry.f90     b8617d69a08a21e47708674d2965b4aaf5ded3827f36c0d3b95745c005952772
patch/lagRamses/snrt_cuda_kernels.cu         bf95eaffe07ff856f16dfbb76710cd8eee81c19dae40b0c101ec76426c1113fc
```

## Scope and remaining qualification

This is a production-path compile and bounded native contract result, not a
live RAMSES evolution. The reference profile remains comparison-only and is
not approved for restart or production cosmological execution; the MPI ledger
write repair is structural, while multi-rank live qualification remains open.
The Newtonian jet guard is a fail-closed `0.9c` bound, not a relativistic
correction. MAD/spin rotational-energy closure, a physical AGN SED,
smoothing/calibration of the `X_floor` discontinuity, dump/restart
persistence, and a short live single-sink thermal/jet case remain separate
qualification work.
