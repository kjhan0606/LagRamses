# Claude Opus 5 bundle-end audit request — F-P2.6 native RT/chemistry transaction

You are the primary bundle-end auditor for F-P2.6 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.

This is a read-only audit of algorithms, physical/numerical legitimacy,
native Fortran/CUDA wiring, transaction boundaries, failure behavior, and
evidence. Do not edit files, launch RAMSES or other jobs, invoke Python/JAX,
or use network tools. Inspect the actual current worktree and distinguish
pre-existing dirty-tree work from F-P2.6 changes.

The final project goal is production-ready and publication-ready high-level
RAMSES radiation transport with stellar/AGN feedback and dust physics. F-P2.6
is the native transport → absorption → H/He thermochemistry → RAMSES thermal
receiver boundary. It is not a physical SED/yield approval, a live production
run, or publication-scale convergence claim.

Read these first:

- `provenance/fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md`
- `provenance/fp2_6_native_rt_chemistry_transaction_bundle_implementation_evidence_2026-09-05.md`
- `provenance/fable_fp2_6_native_rt_chemistry_transaction_bundle_plan_audit_2026-09-05.md`
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
- `bin/Makefile`
- `simulation/snrt/tests/run_snrt_native_rt_transaction.sh`
- `simulation/snrt/tests/run_snrt_native_thermochemistry.sh`
- `simulation/snrt/tests/run_snrt_native_spectral_contract.sh`
- `simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh`

Assess at least the following.

1. Transaction boundary: Is the source-phase ordering correct? Does the
   snapshot cover local photon state, H II/He II/He III, the H I mirror, and
   thermal state? Are all trial mutations isolated until commit? Is rollback
   bitwise by copy, and are coarse-parent/interface corrections trial-owned
   rather than persistent?
2. Transport and AMR wiring: Does every fixed-point trial restart from the
   same incoming photon field? Are Fortran dimensions/orderings correct for
   the CUDA ABI? Do same-level, MPI, coarse-to-fine, and fine-to-coarse paths
   preserve the intended state? Do zero-leaf ranks participate in all required
   collectives without dereferencing absent topology arrays?
3. Collective decisions: Are failure, convergence, and residual decisions
   reduced consistently across MPI ranks before commit/rollback? Can a rank
   with a local transport/partition/chemistry/receiver failure deadlock or
   diverge from the collective path? Is non-convergence fail-closed through
   RAMSES clean-stop?
4. Fixed point: Are start-of-step inventory, incoming radiation, and
   pre-heating temperature frozen? Is the opacity feeding state time-centred
   from start and relaxed end fractions? Is under-relaxation, convergence
   norm, optical-depth floor, iteration cap, and final-trial selection
   numerically well defined? Are consumed trial photons ever reused?
5. Species and energy closure: Is group absorption partitioned against the
   same H I/He I/He II inventory used by CUDA, with no cross-group or
   subcycle double count? Does unassigned absorption fail above the declared
   scale-aware tolerance and remain visible at tolerance scale? Does native
   thermochemistry close primary/secondary/photoelectron energy and return
   unavailable secondary targets to heat? Is RAMSES total energy modified only
   once after a successful receiver decision?
6. Numerical robustness: Is the FP32-ULP inventory tolerance derived and
   applied to the correct pre-partition scale? Are nonfinite/negative states,
   density-zero cells, zero-size arrays, overflow, and invalid temperatures
   handled fail-closed? Are H/He simplex and H I mirror invariants preserved?
7. Build and evidence: Do the recorded native smokes, MPI smoke, CUDA smoke,
   GNU/mpiifx module graph, and full CUDA-linked binary support the claims?
   Check the Makefile module order and actual linked symbols. Separate direct
   evidence from claims that still require live RAMSES evolution, HDF5 restart,
   distributed AMR scaling, physical SED/yield approval, dust processes,
   momentum feedback, and publication convergence.
8. Look specifically for defects that Python-only checks would miss: Fortran
   allocatable/assumed-shape behavior, MPI collective ordering, persistent
   writes hidden behind adapters, coarse correction double application,
   wrong species-axis indexing, cleanup paths, and stale H-only assumptions.

Return exactly one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or
`FAIL`. Then provide severity-ranked findings with file/line references,
evidence, and a disposition for every F-P2.6 acceptance gate (C1–C4 and the
listed acceptance bullets). A record-only or later-G5/G6 condition must be
explicitly separated from a blocker for this bundle. Do not grant physical
AGN/stellar SED approval, feedback production authorization, or publication
science validation.
