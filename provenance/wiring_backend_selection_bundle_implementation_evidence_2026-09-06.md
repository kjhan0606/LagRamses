# Wiring and backend-selection bundle evidence

Date: 2026-09-06
Project root: `/gpfs/kjhan/LRD_JWST`
Base: `85cfa45`
Auditor status: Fable unavailable; see the [Fable attempt record](fable_wiring_backend_selection_bundle_plan_audit_2026-09-06.md)

## Implemented bounded change

`gpu_sink` was previously toggled by an AGN CPU/GPU autotune block in
`patch/lagRamses/amr_step.jaehyun.f90`, but a repository-wide source search
found no consumer of that flag in `AGN_feedback`, `average_AGN`, `AGN_blast`,
or any other AGN routine. The block therefore measured the same path while
claiming to switch backend. It was removed. The namelist symbol remains for
backward compatibility and is now documented in
`patch/lagRamses/amr_parameters.jaehyun.f90` as a non-selecting compatibility
field with default `.false.`.

No automatic selector was added for SNRT or dust because no equivalent
OpenMP implementation exists. No automatic selector was added for stellar or
AGN feedback because no equivalent CUDA implementation reaches the same
source, state, precision, AMR/MPI, and conservation contract.

## Wiring assertions

The following read-only checks were run from the project root:

```text
WIRING_CHECK_BEGIN
SNRT_WIRING_PASS source_to_cuda_to_fp64_to_chemistry
FEEDBACK_WIRING_PASS stellar_cpu_openmp_and_agn_fp64
DUST_WIRING_PASS cuda_ledger_fp64_zero_scaffold_candidate_ir_unwired
BACKEND_DECISION_PASS snrt=cuda_required feedback=cpu_openmp_explicit dust=hybrid_explicit
WIRING_CHECK_PASS
```

The positive SNRT edges are present in
`patch/lagRamses/snrt_ramses_driver.f90`: AGN photon deposition,
`snrt_transport_absorb_multigroup_prepared_dust_trial`,
`snrt_dust_validate_ledgers`, and
`snrt_thermochemistry_advance_cell`. The negative assertion confirms that the
driver no longer calls the old host `snrt_partition_absorption` repartitioner.
The transport adapter reaches
`snrt_cuda_multigroup_rt_step_species_dust` in
`patch/lagRamses/snrt_transport_step.f90`.

The stellar edges are present in `patch/lagRamses/feedback.kjhan3.f90`:
`thermal_feedback` -> `sub_thermal_feedback` -> `phase0_feedback` -> the
stellar runtime, which calls the source/yield/population ledgers and the
transactional RAMSES bridge. The AGN edges are present in
`patch/lagRamses/amr_step.jaehyun.f90` and
`patch/lagRamses/sink_particle.kjhan.f90`:
`AGN_feedback` -> `average_AGN` / `AGN_blast` -> the FP64 deposition
primitives. A repository-wide assertion confirms that `gpu_sink` is not read
by the AGN implementation.

The dust assertions distinguish the live and candidate paths. DUST-8 connects
the CUDA ledger to the FP64 receiver and H/He chemistry while logging
`ZERO_SCAFFOLD`; `snrt_dust_coupling` and `snrt_dust_ir` are not called by the
live RAMSES driver. They remain candidate APIs and are not claimed as live
nonzero dust or IR feedback.

## Consolidated build and native evidence

The existing `/gpfs/kjhan/LRD_JWST/simulation/snrt/tests/run_snrt_bundle_gate.sh`
was rerun after the bounded change. It completed with:

```text
SNRT_BUNDLE_GATE_COMMIT 85cfa457352ee871eb3f9eb51c04bf5e9b21b557
STAGE production_build status=PASS elapsed_s=176.802
STAGE agn_partition_reference status=PASS elapsed_s=1.085
NATIVE_SYMBOLS_CHECK count=5 status=PASS
STAGE dust_ledger_receiver status=PASS elapsed_s=0.194
STAGE thermochemistry status=PASS elapsed_s=0.523
STAGE spectral_contract status=PASS elapsed_s=0.977
STAGE transaction_mpi status=PASS elapsed_s=3.890
MPI_COVERAGE required_ranks=2 marker=SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS
STAGE cuda_multigroup status=PASS elapsed_s=4.452
STAGE production_negative status=PASS elapsed_s=5.770
STAGE diff_check status=PASS elapsed_s=0.069
SNRT_BUNDLE_GATE_PASS
```

The production link used Intel MPI/ifx and CUDA 13.0.2; the CUDA stage ran on
the A10 environment. The gate did not launch a live RAMSES evolution.

## Backend decision matrix

| Subsystem | Implemented path | Equivalent alternative found? | Disposition |
|---|---|---:|---|
| SNRT | CUDA multigroup transport; OpenMP only around host orchestration | No | CUDA required; absence remains fail-closed |
| Stellar feedback | FP64 source/ledger/bridge with OpenMP locks and parallel caller | No CUDA path | CPU/OpenMP explicit |
| AGN feedback | FP64 native/reference routines; legacy `gpu_sink` control had no consumer | No equivalent CUDA feedback path | Explicit model path; dead autotune removed |
| Dust | CUDA primary ledger boundary plus FP64 host receiver; zero scaffold | No equivalent OpenMP transport | Hybrid explicit; no silent switching |

An eventual selector requires differential tests for identical source,
receiver, AMR/MPI ownership, precision, and conservation behavior. Until then,
timing-based switching is not a scientifically valid selection criterion.

## Scope limits

This evidence does not claim nonzero live dust opacity, dust thermal/momentum/
abundance commit, IR transport, HDF5 restart, or production RAMSES evolution.
It also does not claim that the legacy `gpu_sink` name has been removed from
all user namelists; it remains a compatibility input with no new feedback
backend semantics.
