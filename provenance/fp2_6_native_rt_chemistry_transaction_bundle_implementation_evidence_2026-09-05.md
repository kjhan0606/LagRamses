# F-P2.6 native RT/chemistry transaction and fixed-point bundle — implementation evidence

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent: F-P2.5 native H/He thermochemistry
Work location: `/gpfs`

Status: Fable's final fallback bundle-end closure audit returned
`CONDITIONAL PASS`; there is no blocker. The implementation remediation and
refreshed direct evidence below are complete, including the post-audit fixes
for the validation order, collective preflight clean-stop, non-finite hydro
protection, and the bounded iteration limit. Claude Opus 5 was unavailable
for a decisive verdict. Fable's pre-implementation plan audit returned
`CONDITIONAL PASS`; its required plan clarifications and the initial failure
disposition are recorded in
[`fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md`](fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md).

## Native implementation delivered

- [`patch/lagRamses/snrt_rt_transaction.f90`](../patch/lagRamses/snrt_rt_transaction.f90)
  is the explicit-array transaction core. It snapshots local leaf photon
  slices, H II/He II/He III, the H I mirror, and thermal values; it restores
  by copy on failure and commits only a converged trial. It exposes the
  version-1 configuration, defaults (`32`, `1e-6`, `1e-5`, `1e-12`, `0.5`),
  failure classes, and MPI collective max-failure/min-converged/max-residual
  reduction.
- [`snrt_transport_step.f90`](../patch/lagRamses/snrt_transport_step.f90)
  now has a prepared trial ABI. The incoming local photon field and returned
  local field are explicit arrays, while coarse/fine reverse corrections are
  returned in a full-slot `coarse_flux_trial` buffer. The production path no
  longer mutates persistent `snrt_intensity` during transport. Zero-leaf MPI
  ranks still enter the interface and collective path.
- [`snrt_amr_topology.f90`](../patch/lagRamses/snrt_amr_topology.f90)
  accumulates coarse-parent corrections in the trial buffer rather than in
  persistent radiation state. The correction is therefore replayable for
  every fixed-point trial and is committed only once at the level boundary.
- [`snrt_ramses_driver.f90`](../patch/lagRamses/snrt_ramses_driver.f90)
  now freezes the post-source incoming photon state and start-of-step H/He
  inventories, runs the prepared transport → species partition → native
  thermochemistry chain inside a bounded fixed point, and updates opacity from
  start/end time-centred neutral fractions. Temperature is held at the
  pre-heating value for all trials. H/He state and RAMSES `uold` heating are
  not modified until the collective convergence decision and final commit.
  Failed transport, partition, unassigned absorption, chemistry, receiver,
  convergence, and MPI decisions restore the snapshot and call RAMSES
  `clean_stop`; there is no production continuation after a failed level.
- The fixed point uses the raw chemistry result, its under-relaxed opacity
  predictor, and a same-final-trial photon/absorption/heating set. A trial
  always restarts from the frozen incoming photon field, so consumed photons
  are never carried into the next iteration. Source deposition/accounting is
  intentionally outside this transaction and remains committed if the later
  RT/chemistry transaction fails. The convergence residual is explicitly the
  successive under-relaxed predictor residual; a strict raw-map residual is
  not claimed. The worst-case configured cost is 32 prepared transport
  evaluations per level, each including its transport subcycles.
- [`snrt_thermochemistry.f90`](../patch/lagRamses/snrt_thermochemistry.f90)
  replaces the fixed `5e-5` inventory tolerance with a documented scale-aware
  bound based on a 256-FP32-ULP reduction allowance and an 8-host-FP64-ULP
  floor. The native partition residual is reported and cannot pass the
  production boundary above that bound.
- [`config/snrt_rt_transaction_contract_v1.nml`](../simulation/snrt/config/snrt_rt_transaction_contract_v1.nml)
  records the versioned defaults, frozen fields, time-centred opacity rule,
  collective decision, smoke-only failure injection, and clean-stop policy.
  The failure-injection environment controls are rejected whenever the loaded
  spectral contract is runtime-admissible.
- [`bin/Makefile`](../bin/Makefile) includes the new transaction module and
  prerequisite in the SNRT production module graph.

## Direct native evidence

The new combined runner is
[`tests/run_snrt_native_rt_transaction.sh`](../simulation/snrt/tests/run_snrt_native_rt_transaction.sh).
It performs the explicit-array successful commit, bitwise rollback, coarse
correction commit, converged/non-converged norm checks, all three named
stage/leaf injection selectors with production rejection, and a real MPI
two-rank reduction with a zero-leaf rank. It also statically checks that the
production driver routes partition, chemistry, and receiver failures through
the selector and that diagnostic mode is explicit.

```text
bash simulation/snrt/tests/run_snrt_native_rt_transaction.sh
  SNRT_NATIVE_RT_TRANSACTION_PARTITION_ROLLBACK_PASS
  SNRT_NATIVE_RT_TRANSACTION_CHEMISTRY_ROLLBACK_PASS
  SNRT_NATIVE_RT_TRANSACTION_BISECTION_PASS error=2.3283E-10
  SNRT_NATIVE_RT_TRANSACTION_MAXITER1_NONCONVERGED_PASS
  SNRT_NATIVE_RT_TRANSACTION_MAX_ITER_LIMIT_REJECT_PASS
  SNRT_NATIVE_RT_TRANSACTION_SMOKE_PASS
  SNRT_NATIVE_RT_TRANSACTION_GNU_SMOKE_PASS
  SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS ranks=2
  SNRT_NATIVE_RT_TRANSACTION_DRIVER_FAILURE_ROUTES_PASS
  SNRT_NATIVE_RT_TRANSACTION_DRIVER_HYDRO_PREFLIGHT_PASS
  SNRT_NATIVE_RT_TRANSACTION_SMOKE_RUN_PASS
```

The repeated serial markers are the default plus the deterministic named
partition/chemistry/receiver selector cases. The MPI smoke verifies the
collective max/min/max reduction and zero-leaf snapshot/restore path. The
driver route marker is a source-level reachability check paired with the
selector smoke; no live RAMSES failure-injection evolution was launched. The
only compiler messages are Intel long-global-name warnings generated by MPI
module temporaries; there are no compile or runtime failures.

The existing native thermochemistry, spectral/checkpoint, and CUDA controls
were rerun after this bundle:

```text
bash simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_OK
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK

bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_NATIVE_SPECTRAL_CONTRACT_OK
  SNRT_CHECKPOINT_OK
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK

bash simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
  SNRT_CUDA_MULTIGROUP_OK relative_budget_error=5.365561E-09
  SNRT_CUDA_MULTIGROUP_SPECIES_MIX_OK relative_budget_error=2.483527E-09
  SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK
```

The CUDA control is a real `nvcc`/GPU execution of the nine-group species-aware
ABI. The first case covers unavailable opacity-bearing species and group-9
transport; the second covers mixed H I/He I/He II inventories with
opacity-weighted saturation and redistribution. The scale-aware CUDA residual
guard and host pre-partition tolerance therefore have direct native coverage.

The native production module graph compiles with the same options as the
Makefile:

```text
make -C bin SNRT=1 USE_CUDA=1 \
  snrt_rt_transaction.o snrt_transport_step.o snrt_ramses_driver.o
  exit status 0

make -C bin SNRT=1 USE_CUDA=1 ramses
  exit status 0
  linked bin/ramses_final3d
```

The full link was rebuilt after the remediation and final hardening with
`mpiifx`, CUDA 13.0.2,
OpenMP, and the SNRT/HYDRO module graph. The transaction smoke also compiled
with GNU Fortran using an explicit minimal `amr_parameters` test stub; the
production graph itself was compiled and linked with `mpiifx`.

The linked binary contains the native symbols
`snrt_ramses_advance_level`,
`snrt_transport_absorb_multigroup_prepared_trial`,
`snrt_transaction_reduce_decision`, and
`snrt_thermochemistry_advance_cell`. No RAMSES evolution, feedback run, or
large job was launched.

`git diff --check` is clean for the edited bundle files.

## Captured source and binary identities

```text
c16f8aa476f7de3122f025414585c34f43e12f37e6a313a60ab631551d4ffa14  patch/lagRamses/snrt_rt_transaction.f90
fbe8876a6e49f5adfca649c886dda1c47f5d0d8feb65b327aa0b91b3ce577c20  patch/lagRamses/snrt_rt_transaction_smoke.f90
e7cbb40abeb1dee2df221d61ab0f80cabe503ab75c362f5ba331f70f8bcfe265  patch/lagRamses/snrt_rt_transaction_mpi_smoke.f90
83bc7165d0707ad17780fc81a0c6cd3018ab70efc6f7cb87d97ad2baf20bc9be  patch/lagRamses/snrt_ramses_driver.f90
6a818cf4fa7419f0b87122063b4eb687eb6344107e1c01c12b5f4e4167c2e628  patch/lagRamses/snrt_transport_step.f90
99a52c93aa11674d87b75ffd4d1be4f8ef49b44fe5204f3357be97a5b2b879c3  patch/lagRamses/snrt_amr_topology.f90
b8617d69a08a21e47708674d2965b4aaf5ded3827f36c0d3b95745c005952772  patch/lagRamses/snrt_thermochemistry.f90
544c825ec56c90f1955e392fee48708d0e203a4791272c16ad08a6d1019ad14f  patch/lagRamses/snrt_state.f90
bf95eaffe07ff856f16dfbb76710cd8eee81c19dae40b0c101ec76426c1113fc  patch/lagRamses/snrt_cuda_kernels.cu
86b52d8b9c4bac2bfa85a6d9ea894833a7f7a5d0945ec1c25b43391e1d4621cb  patch/lagRamses/snrt_cuda_multigroup_interface.f90
8957ae2f422a549f991d2a6ee397e756ea12b872eadd629d80b82b2d05cfa0da  patch/lagRamses/snrt_cuda_multigroup_smoke.f90
c04b067e6f99a45e2e35f23d93ca475d9381527ad5f454ef6a4b52d17ffe76cf  bin/Makefile
eb4e5395831fe73090d65cfe16a4a845eecc0a344b07b297b7c81fcbe2542718  simulation/snrt/config/snrt_rt_transaction_contract_v1.nml
1dca124592b9dc3c3d263443b14555ac0dec02531c950004f2667a6e4b22aec7  simulation/snrt/tests/run_snrt_native_rt_transaction.sh
df3c3d15da2e9af68def935dcabe5bd270318082d458c53c848538b911e4e9a8  simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh
b5160fb481bc33b43b7640736f0c8df837b90af11e139131b87f21517c3ba44d  bin/ramses_final3d
```

Configuration and cost record:

```text
contract_version=1
max_iterations=32 (hard upper limit=32)
fraction_absolute_tolerance=1e-6
tau_relative_tolerance=1e-5
tau_floor=1e-12
relaxation=0.5
residual=successive under-relaxed opacity predictors
worst_case_cost=up to 32 prepared transport evaluations per level,
                 each with the configured transport subcycles
```

## Final Fable closure audit — 2026-09-05

The primary Claude Opus 5 audit did not issue a verdict (`opus-5` was rejected
and the `opus` alias timed out at 900 seconds), so Fable performed the governed
fallback closure audit. The full operator-captured report is recorded at
[`fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md`](fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md).
Fable returned `CONDITIONAL PASS` with no blocker.

- The stale binary hash condition was closed by recapturing the final linked
  binary identity above (`b5160fb4...`); the binary exports the four native
  symbols and contains the final-hardening strings.
- The fixed-point loop is correct by source review and the norm/bisection
  controls pass, but a live initialized-RAMSES execution of the complete loop
  remains evidence-level record-only work. The present route check is honestly
  static, not a live failure-injection evolution.
- Low-severity diagnostic/long-term records remain: first-failure detail and a
  structured unassigned ledger field, coarse-face/no-slot accounting and
  non-negativity hardening, rank-uniform pre-collective environment policy,
  and distributed AMR performance. These do not block this native transaction
  boundary and are carried to later G5/G6 work.

Fable marked C1 and C2 `PASS`, C3 `PASS` with the live-loop evidence gap, and
C4 `CONDITIONAL` only because of that evidence gap. It did not approve any
physical AGN/stellar SED or yield table, feedback production run, dust model,
HDF5 restart integration, or publication-scale convergence.

The current linked binary exports the native symbols
`snrt_ramses_driver_mp_snrt_ramses_advance_level_`,
`snrt_transport_step_mp_snrt_transport_absorb_multigroup_prepared_trial_`,
`snrt_rt_transaction_mp_snrt_transaction_reduce_decision_`, and
`snrt_thermochemistry_mp_snrt_thermochemistry_advance_cell_`.

## Scope boundary

This bundle closes the native level transaction and bounded local
opacity/chemistry coupling boundary only. It does not approve a physical AGN
or stellar SED, resolve the 40--120 M_sun yield seam, add DTD/PISN physics,
add momentum/thermal subgrid feedback, implement dust scattering/IR/grain
evolution, connect SNRT state to RAMSES HDF5 restart, or establish distributed
AMR performance and publication-scale convergence. Those remain separate
high-level feedback/dust or long-term G5/G6 bundles.
