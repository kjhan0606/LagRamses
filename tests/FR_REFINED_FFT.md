# f(R) FFT corrections on fully covered refined levels

## Scope

This change was built on GitHub `main` commit
`37c05b5a8be1976b6f3a9aa9122b08dddd55c1f2` (2026-09-07).
It extends the existing f(R) spectral Newton correction to refined levels
that cover the entire periodic box. Enable it explicitly:

```fortran
&FR_PARAMS
  ! Retain the intended model parameters and iteration limit.
  fR_fft_refined=.true.
/
```

The default is `.false.`. The existing base-level eligibility and the other
scalar models are unchanged. Refined-level eligibility requires FFTW,
three dimensions, no physical boundary patches, and the global grid count
of a fully covered level. Partial AMR coverage falls back to the existing
solver. All ranks make the same eligibility decision, including empty ranks.
The single-rank transform retains its existing 256-cubed cell limit.

The solver still applies eight spectral Newton corrections, with the existing
fractional step bound, followed by Newton-GS and the original nonlinear
residual test. Neither the nonlinear curvature relation nor `fR_eps` is
relaxed. This is a preconditioner, not a one-pass linear replacement for the
nonlinear f(R) equation.

The active Poisson multigrid implementation provides linear correction MG.
A general nonlinear f(R) FAS solver would additionally need its nonlinear
coarse operator and tau correction. This patch does not implement FAS.

## Reproduction

From this code worktree, build and check eligibility:

```bash
make -C bin -j1 HDF5=1 USE_FFTW=1 EXEC=ramses_fr_refined
python3 tests/test_fr_refined_fft_gate.py
```

The eligibility test compiles the actual Fortran functions with a small
fixture and exercises 14 conditions, including full/partial coverage,
default-off behavior, boundaries, the serial size cap, and overflow rejection.

For the executable comparison, choose a new absolute test directory; the
preparation command refuses to overwrite an existing one:

```bash
python3 tests/run_fr_refined_fft_gate.py prepare --binary bin/ramses_fr_refined3d --root /absolute/new/test-directory --ranks 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tests/run_fr_refined_fft_gate.py run --root /absolute/new/test-directory
```

`prepare` creates a frozen binary, synthetic GRAFIC input, the four effective
namelists, and a manifest with hashes and storage estimates. Review these
before `run`. This is a short solver regression, not production initial
conditions: 16-cubed particles, a sinusoidal displacement, zero velocities,
64 Mpc/h box, starting at a=0.5, and two coarse steps. Refinement reaches a
fully covered 64-cubed force mesh. The fixture requests no reachable full
dump (`noutput=1`, `aout=1.1`, `tout=1d100`, and large periodic intervals).
It writes approximately 110 MiB of first-call force diagnostics across the
four cases instead. Execution is two MPI ranks with one thread each by
default; the runner assigns distinct CPUs on the local host.

`check` reuses existing logs and binary force samples without rerunning:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tests/run_fr_refined_fft_gate.py check --root /absolute/new/test-directory
```

The verifier independently reconstructs the central-gradient force and the
nonlinear n=1 PDE on the periodic lattice. It uses the float32 cosmological
parameters actually read from the GRAFIC header, checks complete unique
cell coverage and physical field sign, and compares on/off fields at the
same epochs. Failed convergence is rejected even if the executable exits
successfully. A single-case check is labeled separately from a full paired
comparison.

## Initial measured result (syntax, 2026-09-07)

Test root:
`/home/kjhan/BACKUP/lagRamses-DE/scalar_fft_checks_20260907_v4`

Authoritative combined report: `CHECK_RESULT.json` under that root.
Each model directory contains `run.nml`, `run.log`, and `force_samples/`.
Build log: `/tmp/paper1-fr-refined-build-20260907.log`.
Executable SHA-256:
`b96bd791a3709d3e9f2f2e7d9ee79706651b5edf8fbb3b5a2e5d35075589a30e`.

Both models used strict convergence, `fR_eps=1e-6`, and `n_iter_fR=6000`.

| Model | Finest-level GS iterations, off / on | Two-step time, off / on | Ratio | Finest-level fifth-force relative RMS difference |
| --- | --- | --- | --- | --- |
| F5 | 1252 / 1 | 6.032 / 0.449 s | 13.4 | 1.885e-6 |
| F6 | 148 / 38 | 1.013 / 0.605 s | 1.67 | 1.955e-6 |

The on cases additionally perform eight FFT corrections per level solve.
Times are the executable's total elapsed simulation time, including those
corrections, not MPI launcher time. Each case was timed once, so these are
initial measurements, not statistically characterized speedups.

All independent nonlinear PDE relative residuals were below 1e-6. The
reconstructed gradient matched the recorded fifth force exactly. The largest
on/off scalar-field relative RMS difference was 6.145e-6 and the largest
density relative maximum difference was 4.571e-12.

These tests validate small full-periodic meshes and the eligibility guard.
The 64-cubed runtime does not exercise the distributed transform branch used
at 256-cubed and above, GPU scalar relaxation, or long cosmological evolution.
It does not establish a force512 production speedup or converged P(k).
Existing grammar production jobs and their frozen executables were not
changed. The distributed follow-up below exercises force512; it does not
replace the need for a representative cosmological production validation.

## Distributed force512 campaign

The 2026-09-07 follow-up is a separate two-node Slurm allocation on grammar,
job `519512`, under
`/gpfs/kjhan/Hydro/DE_nonstd/PaperI_research_extensions_20260905/fr_fft512_gate_20260907_v1`.
The allocation finished after 24 minutes 1 second. The FFT-on cases passed;
both strict FFT-off controls failed to converge, so the overall paired
comparison is correctly marked failed rather than claiming a speed ratio.

- Frozen executable: the same SHA-256 as the small test above; no rebuild.
- Geometry: 128-cubed synthetic particles, full levels 7/8/9, two coarse
  steps from a=0.5 in a 64 Mpc/h box. This tests the distributed force mesh,
  not a 1024-cubed particle production simulation.
- Allocation: two nodes, 16 MPI ranks (eight per node), two OMP threads per
  rank, 320 GiB reserved per node, three-hour Slurm limit. The established
  F5/F6 production nodes are excluded.
- Four cases: F5 on, F6 on, F5 off, F6 off, executed sequentially on the same
  allocation. Each case has a 30-minute cap and the unchanged strict
  `fR_eps=1e-6`, `n_iter_fR=6000` settings. A failed baseline is recorded as
  a failure, not a speed measurement against a converged solution.
- No reachable full outputs: `noutput=1`, `aout=1.1`, `tout=1e100`, and
  `foutput=fbackup=1e9`. Expected force diagnostics total 54.75 GiB.
  GPFS had approximately 170 TiB available at launch.

Each case has `run.nml`, `run.log`, `RUNTIME_RESULT.json`, and, after
verification, `FIELD_RESULT.json`. The combined scheduler/verification log
is `slurm-519512.log`. `PAIR_f5.json`, `PAIR_f6.json`, and `RESULT.json`
are generated only after the respective checks finish. Read those files
for the measured status; a submitted or running job is not a pass.

The runner requires scalar log entries with `distributed=T` at both
levels 8 and 9 for the on cases. It independently reconstructs the
nonlinear PDE and fifth-force gradient from every captured cell. Paired
checks additionally require matching epochs and on/off density and force
agreement. It keeps all raw samples, including partial data from failures.
Scripts and ICs are frozen with hashes in `DISTRIBUTED_MANIFEST.json`.

To reproduce in a new absolute directory, use
`tests/run_fr_distributed_fft_gate.py prepare`, audit the generated effective
namelists and manifest, then submit the copied
`source/run_fr_distributed_fft_gate.slurm` with
`FR_DISTRIBUTED_RELEASE` set to that directory. The reference Slurm script
matches the default 16-rank/two-thread preparation arguments.

### Measured outcome

| Model | FFT-on two-step time | Level-9 FFT-on solve | Independent level-9 PDE residual | FFT-off result at 6000 sweeps |
| --- | --- | --- | --- | --- |
| F5 | 26.864 s | 8 FFT corrections + 1 GS | 1.708e-13 | Failed at level 8, residual 4.829e-3 |
| F6 | 28.418 s | 8 FFT corrections + 1 GS | 6.729e-13 | Failed at level 9, residual 9.753e-6 |

All FFT-on captures passed full lattice coverage, physical field sign,
nonlinear PDE, and fifth-force gradient checks. Both accelerated levels
explicitly reported distributed scalar transforms. The FFT-off times to
failure were 144.396 s for F5 and 1125.040 s for F6, including MPI startup;
they are not successful timings and are not used to compute speed ratios.
No complete force512 on/off field comparison was possible. Neither control
was stopped by the 30-minute wall cap: both stopped for strict nonconvergence.

The aggregate `INCOMPLETE_OR_FAILED_DISTRIBUTED_FORCE512` / Slurm `FAILED`
status reflects the failed controls. Individual FFT-on runtime and
independent PDE statuses are PASS. Existing production jobs `516654_8`
and `519024` remained running on their original nodes after this test.
No snapshots were produced; approximately 30 GiB of test input, samples,
and logs were retained on GPFS.

A local copy of the reports and logs is kept at
`/home/kjhan/BACKUP/lagRamses-DE/fr_fft512_report_20260907_v1`.
See its `REPORT.md` and `RESULT.json` for the scoped conclusions and raw
per-case outcomes. No production speedup or long-evolution claim follows
from this synthetic two-step solver regression.
