# Manual Ada M5 + dust + CHIMES smoke — 2026-10-21

## Scope and result

Ran the bounded `effective_population` live M5/dust/CHIMES profile manually on
lageunha's RTX 5000 Ada. This is a reference-control integration smoke, not a
production-physics qualification or a performance comparison. The accepted
run completed four coarse steps with one MPI rank and four OpenMP threads.
Both AMR levels 2 and 3 committed M5 order 5; the script's dust, stellar,
momentum, finite-state and nonconvergence checks passed. No `output_NNNNN`
directory was produced.

Accepted evidence is retained at:

- `.m5-effective-amr-lageunha-20261021-r3/result.txt`
- `.m5-effective-amr-lageunha-20261021-r3/run.log`
- `.m5-effective-amr-lageunha-20261021-r3/effective.nml`
- `.m5-effective-amr-lageunha-20261021-r3/inputs.sha256`

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
lageunha GPU process was active, so wall time is deliberately not compared.
The one-stream / peak-one-inflight result is also not a scaling qualification.

## Runner changes

`simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch` now has a
restricted direct-launch mode for the one-rank lageunha smoke, supports the
Ada device only in that mode, captures the CHIMES counters, uses portable
`grep` checks, and requires the four-thread minimum for this RT-on hybrid
profile. `simulation/snrt/tools/chimes_runtime_counts_reporter.c` provides the
small `LD_PRELOAD` destructor that reads existing exported counters; it does
not alter solver state or tolerances.

The separate A40 Slurm smoke request was not used. No raw output cleanup was
needed because the accepted run wrote no raw dumps.
