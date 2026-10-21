# Manual Ada M5 + dust + CHIMES smoke — 2026-10-21

## Scope and result

Ran the bounded `effective_population` live M5/dust/CHIMES profile manually on
lageunha's RTX 5000 Ada. This is a reference-control integration smoke, not a
production-physics qualification or a performance comparison. Two accepted
runs completed four coarse steps: MPI=1/OMP=4 and MPI=2/OMP=4. Both AMR levels
2 and 3 committed M5 order 5; the script's dust, stellar, momentum,
finite-state and nonconvergence checks passed. Neither run produced an
`output_NNNNN` directory.

Accepted evidence is retained at:

- `.m5-effective-amr-lageunha-20261021-r3/result.txt`
- `.m5-effective-amr-lageunha-20261021-r3/run.log`
- `.m5-effective-amr-lageunha-20261021-r3/effective.nml`
- `.m5-effective-amr-lageunha-20261021-r3/inputs.sha256`
- `.m5-effective-amr-lageunha-20261021-mpi2/result.txt`
- `.m5-effective-amr-lageunha-20261021-mpi2/run.log`
- `.m5-effective-amr-lageunha-20261021-mpi2/effective.nml`
- `.m5-effective-amr-lageunha-20261021-mpi2/inputs.sha256`

The effective schedule was `noutput=1`, `aout=2`, `tout=1e30`,
`foutput=fbackup=1000000`, with a four-step bounded run. The run-time checksum
verification passed for the template, IC, effective namelist, executable,
CHIMES library and physical tables.

## CHIMES GPU finding

The first instrumented attempt used two OpenMP threads and completed the
hydrodynamics, but reported `gpu_calls=0`, `cpu_fallback_calls=170020`, and a
zero GPU worker mask. Inspection of `snrt_chimes_rhs_cuda.cu` showed the RT-on
per-cell RHS dispatcher deliberately falls back to CPU when its OpenMP team is
smaller than four. The runner now rejects that configuration instead of
silently calling it a hybrid GPU run.

The accepted four-thread run reported:

```
CHIMES_RHS_COUNTS gpu_calls=20771 cpu_fallback_calls=151126 gpu_batches=5413 device_errors=0 warm_table_uploads=1 peak_inflight_batches=1 streams_used=1
CHIMES_GPU_WORKER_MASK value=15
```

This proves successful CHIMES device RHS work, but not a speedup. Another
lageunha GPU process was active during this first run, so wall time is
deliberately not compared.

The two-rank run used Intel MPI pinning with `I_MPI_PIN_DOMAIN=omp`; it
reported two CHIMES counter records, both with nonzero GPU calls, zero device
errors, `streams_used=1`, and worker mask 15:

```
gpu_calls=16800 cpu_fallback_calls=68100 gpu_batches=4200 device_errors=0 warm_table_uploads=1 peak_inflight_batches=1 streams_used=1
gpu_calls=13200 cpu_fallback_calls=71920 gpu_batches=3300 device_errors=0 warm_table_uploads=1 peak_inflight_batches=1 streams_used=1
```

The integration gate passed with 56 coarse and 64 fine leaf cells and no raw
outputs. RAMSES reported 90.5 seconds total, of which 86.7 seconds (95.8%) was
in cooling. This tiny reference case is not a throughput/scaling comparison;
the one-stream / peak-one-inflight result is not a stream-scaling qualification.

## Runner changes

`simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch` now has a
restricted direct-launch mode for one or two MPI ranks on lageunha, supports
the Ada device only in that mode, captures the CHIMES counters, pins MPI ranks
to separate OpenMP domains, uses portable `grep` checks, and requires the
four-thread minimum for this RT-on hybrid profile. It verifies that each MPI
process emits successful GPU-call and worker-mask records.
`simulation/snrt/tools/chimes_runtime_counts_reporter.c` provides the small
`LD_PRELOAD` destructor that reads existing exported counters; it does not
alter solver state or tolerances.

The separate A40 Slurm smoke request was superseded by the requested manual
lageunha Ada run. No raw output cleanup was needed because both accepted runs
wrote no raw dumps.
