# Fable follow-up bundle-end audit request — F-P2.6 remediation

You are the fallback bundle-end auditor for F-P2.6 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.
The primary Claude Opus 5 audit was attempted twice but produced no verdict:
`opus-5` was rejected by the local model catalog and `opus` timed out after
900 seconds. Fable's first fallback audit returned `FAIL`; the original
report is `provenance/fable_fp2_6_native_rt_chemistry_transaction_bundle_end_audit_2026-09-05.md`.

This is a read-only follow-up audit after remediation. Do not edit files,
launch RAMSES or other jobs, invoke Python/JAX, or use network tools. Inspect
the actual current worktree and distinguish pre-existing dirty-tree work from
the F-P2.6 bundle. The final project goal is production-ready and
publication-ready high-level RAMSES radiation transport with stellar/AGN
feedback and dust physics. F-P2.6 is only the native transport -> absorption
-> H/He thermochemistry -> RAMSES thermal receiver transaction boundary; it
does not approve a physical SED/yield, live feedback evolution, dust physics,
or publication-scale convergence.

Read these first:

- `provenance/fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md`
- `provenance/fp2_6_native_rt_chemistry_transaction_bundle_implementation_evidence_2026-09-05.md`
- `provenance/fable_fp2_6_native_rt_chemistry_transaction_bundle_end_audit_2026-09-05.md`
- `simulation/snrt/config/snrt_rt_transaction_contract_v1.nml`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `patch/lagRamses/snrt_rt_transaction.f90`
- `patch/lagRamses/snrt_rt_transaction_smoke.f90`
- `patch/lagRamses/snrt_rt_transaction_mpi_smoke.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_transport_step.f90`
- `patch/lagRamses/snrt_amr_topology.f90`
- `patch/lagRamses/snrt_thermochemistry.f90`
- `patch/lagRamses/snrt_nlte_coupling.f90`
- `patch/lagRamses/snrt_state.f90`
- `patch/lagRamses/snrt_cuda_kernels.cu`
- `patch/lagRamses/snrt_cuda_multigroup_interface.f90`
- `patch/lagRamses/snrt_cuda_multigroup_smoke.f90`
- `bin/Makefile`
- `simulation/snrt/tests/run_snrt_native_rt_transaction.sh`
- `simulation/snrt/tests/run_snrt_native_thermochemistry.sh`
- `simulation/snrt/tests/run_snrt_native_spectral_contract.sh`
- `simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh`

Re-audit the original findings, not only the prose. In particular verify:

1. Source-phase ordering, complete slot-indexed photon/H/He/H-I/thermal
   snapshot, copy rollback, and transaction-owned coarse corrections.
2. Explicit Fortran/CUDA dimensions and axes; group/species inventory
   consistency; same-level, MPI, coarse-to-fine, fine-to-coarse, and
   zero-leaf collective behavior. Check the species-aware CUDA active-set cap
   for unavailable species, mixed inventories, group order, and subcycle
   accounting.
3. Collective failure/convergence/residual decisions and every prepared-path
   return. Check the driver CUDA/uold/unit preflights and clean-stop on
   non-convergence or coupled failure.
4. Fixed-point frozen inputs, time-centred opacity, under-relaxation,
   residual semantics, maximum 32 iterations, same-final-trial commit, and no
   reuse of consumed photons. Make sure the evidence states that the residual
   is the relaxed predictor residual and records the worst-case cost.
5. Scale-aware pre-partition inventory tolerance, above-tolerance unassigned
   rejection, visible tolerance-sized ledger residual, finite/non-negative
   trial checks, H/He simplex, photoelectron/thermal closure, and one-time
   `uold` update.
6. Native evidence: current thermochemistry, transaction GNU/mpiifx/MPI,
   spectral/checkpoint, CUDA mixed-species, `git diff --check`, full
   `make -C bin SNRT=1 USE_CUDA=1 ramses`, linked symbols, and current hashes.
   Treat the driver-route static check honestly: it is not a live RAMSES
   failure-injection evolution.

Return exactly one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or
`FAIL`. Then list severity-ranked findings with file/line references,
explicitly mark every original B/H/M/L finding as closed, residual, or
blocker, and give a disposition for each F-P2.6 acceptance gate C1--C4 and
each acceptance bullet. Separate record-only/G5/G6 conditions from blockers.
Do not grant physical AGN/stellar SED approval, yield-table approval,
feedback production authorization, dust approval, or publication science
validation.
