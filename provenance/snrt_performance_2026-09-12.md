# SNRT runtime cost and hybrid implementation, 2026-09-12

Scope: existing RT/chemistry/dust implementation performance, not a new physics
gate. Repository `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Operator instruction: stop stream-count fine tuning; retain configurable
defaults. A card's optimum depends on its hardware, MPI sharing and workload.
The stream pool is process-local: streams per rank are not streams per card.

## Implemented changes

- `snrt_hybrid.cpp`: CPU IR batches read the immutable input and halo directly,
  avoiding the seven-neighbor GPU packing buffer and extra output copy. GPU
  stream availability still selects GPU versus CPU; no new autotuner.
- `snrt_ir_cuda.cu`: one block per cell, adjacent threads over contiguous rays,
  instead of one thread serially traversing a cell's entire frequency/angle
  field. FP64 expressions and transactional commit/rollback retained.
- `snrt_dust_live.f90`: OpenMP cell-parallel cosmological material callback,
  cell-major halo packing/unpacking with MPI outside the parallel regions.
  Three rank-local IR phase wall times identify halo, solve and scatter cost.
- `snrt_ramses_driver.f90`: one rank-local cold-chemistry wall timer, included
  in coupling, for subsequent runs. This last timer is NOT in job 541239.

No namelist/default, tolerance, physical model, or Makefile VPATH changes were
made for these performance edits. Other working-tree changes predate this work.

## Native measurements

Evidence directory: `.snrt-performance.jDx9Iz/` (logs, native benchmark source,
before/after binaries). NVIDIA A10 on syntax, one visible card, OMP4,
2048 cells x 136 frequencies x 80 directions, 256-cell batches, two streams.
Nonzero periodic manufactured field. Whole entrypoint timing includes host
packing, allocation and transfers; one warmup, best of two measured calls.
Shared host and short samples: these are local measurements, not universal
speedups or an optimum-stream determination.

| Backend | Before seconds | After seconds | Time reduction |
| --- | ---: | ---: | ---: |
| OpenMP | 0.138062 | 0.109508 | 20.7% |
| CUDA | 0.621927 | 0.479026 | 23.0% |
| Hybrid | 0.187488 | 0.162713 | 13.2% |

After hybrid is still 1.486 times the OpenMP time. GPU offload is not inherently
beneficial for this host-resident IR workload. Current forced-CUDA entrypoint
uses one host worker, so this is not a best-possible asynchronous CUDA pipeline.
All reported mode comparisons have zero output relative difference. Existing
native hybrid busy-pool/failure rollback regression, dust backend parity and
rollback suite, and cosmological dust expansion fixture pass. See
`regression.log`, `backend-regression.log`, `native-regression.log`.

OMP1-to-OMP4 at the same batch size gives 0.313404/0.109508 = 2.86x after
the CPU change. OMP8/batch64 measured 0.074910 seconds but changes two variables
and is NOT a clean thread-scaling measurement. Stream 1/2/3/4 trial logs remain
reference data only; no default change or further tuning follows.

## Integrated comparison

Baseline job 540015, `.galaxy-coupled-pilot.ZGTFWW/`: MPI32 x OMP4 on four
grammar nodes, two 128-cubed cosmological steps; CPU-only build. RAMSES wall
1000.5275 seconds; broad cooling mean 983.730 seconds. Allocation 35.9822 core
hours; actual recorded TotalCPU 28.5575 core hours. These are distinct costs.
Initial grains and photons are zero, with no SF/AGN: this is not a luminous
galaxy benchmark or a GPU end-to-end measurement.

Comparison job 541239, `.snrt-performance.jDx9Iz/full-profile/`, uses the same
NML SHA256 `08ff16e8893ff28d7140885062bae3e9ba7170421944a0bb0f8ced60e4d2283a`
and resources on grammar[031-034]. Binary `ramses_profile_final_3d` SHA256
`26a37d7f84db2317aad9d1e7366408908c8eedadc63a6c1be9c6a8dcec8fedc2`.
Completed 0:0, elapsed allocation 960 seconds, RAMSES wall 947.750337 seconds.
Recorded TotalCPU 1-03:26:02 = 27.4339 core hours; allocation 34.1333 core
hours. Versus baseline, wall falls 5.27%, actual CPU 3.93%, allocation 5.14%.
Per-task MaxRSS 40671576 KiB is approximately unchanged, NOT global RSS.
Cooling mean remains 928.565 seconds (97.9% of the broad timer total).
First step reports coupling 417.377 seconds,
transport 8.045 seconds; IR halo 7.4108, solve 39.6358, scatter 8.2830 seconds.
Transport and IR are INSIDE coupling; do not add coupling to its constituents.
Approximately 354 seconds remain outside these measured constituents, not yet
proven to be exclusively chemistry. The cold-chemistry loop is already OpenMP;
its CVODE dark solver allocates a dense matrix and solver per cell. Solver
reuse or alternative scheduling needs evidence, not a tolerance relaxation.
Different nodes and shared-machine load limit causal whole-run comparisons.

Second step: coupling 478.413, transport 11.611, IR halo 9.7286, solve
40.9833, scatter 8.7110 seconds. Across both steps the three IR timers sum to
114.7525 seconds; transport sums to 19.656 seconds, coupling to 895.790.
About 761.38 seconds of coupling remain outside these subsets. Do not assign
all that residual to CVODE without the chemistry-specific measurement.
Runner reports completed=1, CMB_rank_commits=64, IR_steps=2, no rejection;
final mcons=0 and econs=-1.98e-9 match the baseline reported precision.
Both IR balances and escapes are zero as expected for this initial state.

The subsequent chemistry-timed binary compiled successfully, SHA256
`81fa3e824841817635c5ade95591ac2c851230aa40228a16f94815703ff281ad`.
It was not used for the measurements above. No chemistry algorithm/scheduler
optimization or measured chemistry-only speedup is claimed.

Output policy: noutput=1, aout=1.1, tout=1e100, foutput=fbackup=1000000,
nstepmax=2; no scheduled raw dumps in this test. Retain logs and identities;
post-completion inspection found no `output_*` directories. Cleanup manifest:
zero raw outputs, zero files removed; inputs/logs/binaries retained.

## Tensor Core assessment

The repository ALREADY has a cuBLAS fast-TF32 angular-reduction prototype in
`snrt_cuda_kernels.cu`, used only by the `SNRT_P1_DIAGNOSTIC=1` diagnostic,
not by the physical advance. It is disabled in the integrated pilot.
`tensor_probe.cpp`/`tensor.log` exercise that actual entrypoint with nonconstant
positive inputs: 2048 rows x 80 directions x 16 bins, 0.00079894 seconds,
maximum relative error 5.77739464e-4 against FP64 evaluation of the same FP32
inputs. Requested TF32 compute is verified; hardware Tensor instruction
counters were not measured. This error cannot silently replace strict energy
ledger arithmetic. No Tensor Core physical path was enabled.

NVIDIA's [GA102 architecture whitepaper](https://www.nvidia.com/content/PDF/nvidia-ampere-ga-102-gpu-architecture-whitepaper-v2.pdf)
states GA10x lacks A100-style FP64 Tensor Core acceleration. The
[Ampere tuning guide](https://docs.nvidia.com/cuda/ampere-tuning-guide/)
also emphasizes coalescing and reducing host/device transfers. The current
FP64 IR stencil/exponentials and chemistry integration are not direct TF32
matrix-product substitutions. Existing dense chemistry linear algebra is a
possible future batched solver target, not an implemented GPU speedup.

Next action is to isolate chemistry cost before additional optimization.
The integrated comparison is complete. No stream sweep, new audit gate, or
mixed-precision physics shortcut is required.

## Approved zero-light work (after the operator's RT cost reinspection)

Implemented exact-zero local work avoidance, without bypassing transport,
MPI participation, dark chemistry, or physical input validation:

- Photo solver allocates its large node-ray buffer lazily after a validated
  nonzero ray. Existing zero-light identity receipt and all opacity/phase
  checks remain unchanged.
- Node reconstruction validates every ray before directly publishing zeros
  for the all-zero case; no private node buffers/copies are needed.
- Live grain scattering skips node reconstruction only when BOTH incident
  photon N and energy are exactly zero. It still calls the moving receiver
  to validate phases/directions and produce the transactional work receipt.
- Existing cold-loop wall timing is supplemented by two sums of per-worker
  elapsed time (cold call and grain scatter). These sums are NOT CPU-time
  measurements and are not additive to the encompassing wall timer.

Native 65,536-cell, 80-direction zero-node probe: OMP4 before 12.563927s,
after 0.447718s (~28x for this function only). Post-change invalid-input
rollback and nonzero N/E moment checks pass. `zero-nodes-after.log` and
`zero_nodes_probe.cpp` retain evidence; same-host short benchmark caveats apply.
The live all-zero scatter adapter now avoids even those node-output fills.
The existing nonzero long-interval CHIMES regression passes with nuclear,
charge, photon and energy errors 1.74e-10, 6.88e-14, 5.56e-12, 1.33e-11.
See `zero-light-chem-regression.log`; no tolerances changed.

Full build `ramses_zero_light_3d` SHA256
`e9cf55969ba06c158c043202a64f724ec34bc41657379c331de8a1185cb7dc94`.
Submitted grammar Slurm job **541547**. The new run directory is
`.snrt-performance.jDx9Iz/zero-light-profile/`.
Effective NML is byte-identical to job541239. MPI32 x OMP4, four grammar
nodes, 20-minute limit, two steps, unchanged no-dump policy; 91TiB free at
submission. Job541547 COMPLETED 0:0; RAMSES wall 1081.9171s, allocation
1093s, TotalCPU 28.5994 core hours. This is 14.16% slower in wall and 4.25%
more actual CPU than job541239. Do not claim integrated acceleration.
Same node list does not control co-tenant load. Both IR commits, all 64 rank
CMB commits and mcons=0/econs=-1.98e-9 pass. No raw snapshots were produced.
The native function speedup is not a whole-simulation speedup.

Rank1 two-step nonoverlapping decomposition: cold-loop wall 697.153s,
IR halo/solve/scatter 119.0182s, primary transport 18.817s, remainder of
coupling 192.0298s (coupling total 1027.018s). The cold-call accumulated
worker wall is 2733.235s, not CPU time. Separately timed moving-grain scatter
reports 0.000s: this run does not establish that adapter as a dominant cost.

Structural RT findings retained for the next decision: a 136 x 80 FP64
IR field over 128^3 cells alone is 170GiB (5.3125GiB per rank at MPI32),
before trial/transported/candidate/halo copies. IR halo uses 680 tiles per
substep and remaining serial field scans exist. These are confirmed costs
to investigate, not proof that they explain all unassigned coupling time.
No global empty-IR shortcut has been added: incoming boundaries, sources,
material emission and coarse/fine state must first be accounted for.

## Split attribution follow-up

Grammar job541705 submitted in `.snrt-performance.jDx9Iz/split-profile/`.
Same physical parameters, MPI32 x OMP4; only termination changes to one
step, 12-minute scheduler bound, no raw dumps. Compare to the FIRST step of
previous runs, not their two-step totals. NML SHA256
`74537062aa5ba405b68f49e1d4f07c0dd9501e6c79c4b3d3fde10eb60b927f7d`;
binary `ramses_split_profile_3d` SHA256
`adb9f862a00210e6bbe950d9355d7bda3c69e9c86930b792d72206628b029f20`.
Build and existing CHIMES long-interval regression pass (`split-build.log`,
`split-regression.log`). No physical arithmetic/tolerance changes this follow-up.

Two thread-private elapsed counters separate successful photo and dark
integration calls within cold chemistry; driver sums them over workers.
They are NOT wall-time contributions to add directly to the outer wall timer.
Four sequential driver intervals separate primary loop, pre-IR preparation,
whole live-IR stage and tail/commit. The whole live-IR stage INCLUDES its
previously logged halo/solve/scatter constituents. Job541705 COMPLETED 0:0,
RAMSES wall 440.7924s for ONE step, TotalCPU 12:48:41. Rank1 cold wall
320.992s; summed worker photo 0.519s, dark 1263.376s, whole cold call
1264.771s. Thus 99.89% of the measured cold-call worker time is in dark
integration, not photo chemistry. Primary/preIR/liveIR/tail wall:
332.913/0.001/73.570/5.468s; liveIR includes halo/solve/scatter
7.4469/40.8174/7.8640s. Do not reinterpret these as two-step totals.
32 rank CMB commits and one IR commit pass; no raw dump generated.

## Dark integration optimization

Subsequent [exact-temperature rate reuse](chimes_rate_cache_2026-09-12.md)
job542050: cold52.333s, whole150.466105s, actual CPU03:04:34; same-node
comparison to542047 reduces cold17.9% and CPU11.3%, all checks pass.

Latest: [SUNDIALS build correction and cost attribution](chimes_solver_cost_2026-09-12.md).
The previous numerical library was built without optimization. Precise FP64
Release gives same-node job542047 cold63.739s / whole162.228304s /
TotalCPU03:28:05 versus542031 cold193.346s / whole325.976843s /
TotalCPU08:20:25. All conservation/commit checks pass. Use the documented
release SUNDIALS prefix for future builds; preserve the old baseline library.

Follow-up: [five-reviewer synthesis and P0/P1 bundle](chimes_dark_bundle_2026-09-12.md).
Job542031 completed with native parity and conservation checks passing;
cold wall193.346s, TotalCPU08:20:25, whole wall325.976843s. This is not an
overall speedup claim: node allocation changed and collective time increased.

`snrt_chimes_bridge.c`: for split cold dark modes 1/3 ONLY, the local solver
configuration uses N_spectra=0 when the RT/on-the-spot flags are both one.
Entry requires exactly zero incident photons; these wrappers supply zeros.
The RT flag and on-the-spot choice remain unchanged, preserving case-B rates;
species, thermal/CR/CMB physics and tolerances are unchanged. This eliminates
18 inactive radiation ODE unknowns and photo-rate loops. Global loaded tables
and configuration are not changed; general illuminated mode0 is unaffected.

Initial all-dark-mode candidate failed a hot neutral input's electron
reconciliation (status46, electron -2.8e-48). It is NOT admitted: hot/atomic
mode2 retains the original system. No negative-electron gate was relaxed.
Changing system dimension changes CVODE's error norm/trajectory, so algebraic
zero radiation alone does not establish numerical parity for all inputs.

Native `dark_probe.f90` exercises 100/10000/1000000K and nH=1e-4/.2/100
cm^-3, 32 identical calls per state, neutral H/He/trace C/O, z99 CMB, dt1e10s.
Logs `dark-before.log` and `dark-cold.log`: cold calls are 2.19–3.29x faster;
hot original path is approximately unchanged. Maximum temperature relative
difference 2.14876e-9, species/H absolute difference 1.73005e-17; all nine
cases pass nucleus/charge checks. This matrix does not cover every molecular
mixture, metallicity or long interval. Existing long photo regression also
passes (`dark-regression.log`), but does not replace full cold coupling tests.

Submitted grammar job541870. One-step integrated candidate directory
`.snrt-performance.jDx9Iz/dark-profile/`:
same NML hash74537062aa5ba405b68f49e1d4f07c0dd9501e6c79c4b3d3fde10eb60b927f7d,
MPI32 x OMP4, CPU grammar, 12-minute bound, no scheduled dumps.
Binary `ramses_dark_reduced_3d` SHA256
`9becd02ac0792fd0ee7ae59afed14ac10dd783634be82d0ed2ed9a03e1772aa3`.
Job541870 COMPLETED 0:0: RAMSES wall373.3188s, allocation383s, actual CPU
08:59:34 (8.9928 core hours). Against one-step job541705: wall -15.31%,
actual CPU -29.81%; same physical NML/node list, shared-host caveats remain.
Rank1 cold wall198.752s versus320.992s (-38.08%); summed dark worker
776.150s versus1263.376s. All 32 CMB rank commits and the IR commit pass;
reported mcons=econs=0. No raw output directories exist; cleanup required
zero deletions. Logs/inputs/binaries are retained. No CUDA/stream change.

Following completion, native OMP4 mixed H/HII/H2 validation ran without
another full simulation: HII/H=.01, H2/H=.01, electron/H=.01, residual HI=.97,
same He/C/O totals; 100/10000K x nH1e-4/.2/100, 32 threaded calls plus a
serial checked call per case. All six states pass species/nucleus/charge
checks (`dark-mixed.log`, `dark_mixed_probe.f90`). This is a thread/conservation
check, not an additional before/after state-parity measurement.

Remaining measured liveIR wall73.813s is unchanged versus73.570s, with
halo7.3738/solve40.7998/scatter7.2323s. Primary wall263.674s includes the
198.752s cold wall and7.392s transport, leaving57.530s there. This residual
varies substantially between runs; neither MPI wait nor memory work alone
has been established as its cause. Do not subtract the entire local cold
speedup from the global elapsed time or claim all remaining work is RT.

## IR copy/parallel follow-up

Job542016 submitted, `.snrt-performance.jDx9Iz/ir-profile/`, same one-step
NML hash74537062aa5ba405b68f49e1d4f07c0dd9501e6c79c4b3d3fde10eb60b927f7d.
Binary `ramses_ir_parallel_3d` SHA256
`a54d795a92168302c6845545d929dae3165c29c4d22b04dd8ae99722c28e0c73`.
MPI32 x OMP4 on grammar, 12-minute bound, same no-dump schedule, 91TiB free.

CPU IR absorption now leases before packing and directly reads immutable
input/coefficient arrays when assigned to CPU, writing private transactional
outputs. GPU implementation/lease policy is unchanged. The redundant initial
`transported=energy` copy is omitted when transport dispatch supplies the
whole result. Cell-independent old-energy/interface-flux and candidate-energy
loops now use OpenMP; scalar sums use reductions, so rounding order can
change, but acceptance thresholds do not. Invalid transport trials still
return before public commit. No zero-field transport shortcut is introduced.

Existing CPU backend/rollback and expansion tests pass after rebuilding:
`ir-parallel-verified.log`, `ir-expansion-verified.log`. A10 hybrid mixed/busy
pool/rollback test passes (`ir-hybrid-regression.log`). A single timer around
the existing primary decision collective distinguishes rank wait from local
work; it is a subset of primary wall, not an additive new stage.
Job542016 COMPLETED 0:0 and evaluated in the same working turn. RAMSES
wall312.5697s versus373.3188s; actual CPU08:23:17 versus08:59:34, a 6.72%
CPU reduction. All32 CMB rank commits and one IR commit pass; reported
mcons=econs=0. No raw snapshots exist; zero files deleted.

Rank1 liveIR65.828s versus73.813s (-10.82%), with halo7.2470,
solve35.6311 versus40.7998 (-12.67%), scatter9.1364s. Cold wall197.777s
and summed dark worker776.141s are essentially unchanged. Primary wall
213.491s includes cold197.777s, transport7.530s, decision collective7.674s;
remaining local primary work is0.510s. The earlier57.530s primary residual
does not recur. The old run did not measure decision wait separately:
do NOT assert that all its residual was MPI wait, or attribute the entire
60.749s global reduction to IR edits. Rank timing/load variability remains.
The scoped IR implementation/test bundle is complete; no active job remains.

## M5 128-cubed preparation follow-up (2026-09-24)

Job 400507 (`.galaxy128-m5-d03-gpu4-hybrid-mpi16-o2/`) timed out after four
hours on H200 x4, MPI16 x OMP2. Rank 1 completed 135 material tiles of 256
cells during its first RT subcycle. Direct aggregation of those 135 log lines
gives 82.114 s/tile for serial per-cell preparation (range 69.347–91.122 s),
18.456 s/tile for the hybrid GPU/IR stage (range 15.349–21.222 s), and
0.095 s/tile for commit. This is an observed rank-local split, not a
whole-run speedup or a CHIMES-only attribution. The effective first-step
parameters imply approximately ten RT subcycles from the unchanged 1/12
transport CFL. Consequently the pending eight-hour MPI32 x OMP1 job 402942
is a scaling probe, not a credible full-step completion gate.

The M5 tile now OpenMP-parallelizes only its independent per-cell preparation;
the native IR/material tile transaction and commit remain sequential. Cell
failures, including negative codes, are checked after the parallel region.
CHIMES split times use its existing thread-private counters, summed over cells;
these worker elapsed sums must not be subtracted from tile wall time. Existing
OMP4/serial CHIMES parity evidence above supports native-cell concurrency,
but this new M5 tile path still requires an integrated OMP2+ comparison.
Intel Fortran built `ramses_m5_preppar3d` in
`.galaxy128-m5-d03fix-serialtile-build/` (SHA256
`0353a870a27842ceb64f6678b610e34ab284684a3ff65d762f3fc93103e20e7d`).
Job 402942 was initially submitted with the earlier binary. While it was
still `PENDING` on 2026-09-24, its run-directory binary was replaced with
the above SHA256 build after preserving the old copy as
`ramses_m5_3d_pre_d03fix_cc4e30c2`. The new build contains the corrected
cosmological D03 material-temperature branch and the per-cell preparation
OpenMP loop; the submitted MPI32 x OMP1 layout does not exercise preparation
parallelism. The existing Slurm script reads the run-directory binary and
`launch.sha256` at execution time; the updated checksum and five pinned
SUNDIALS libraries pass preflight. Its sourced run environment now purges
inherited modules and explicitly loads the pinned Intel MPI/compiler and
CUDA stack; a non-interactive `set -euo pipefail` preflight passes. Job
402942 is still a scaling probe, not
the integrated OMP2 comparison. No full-step result or speedup is claimed.

Resource check at 2026-09-24 21:18 KST: H200 node syn104 has 64 CPUs and
eight GPUs; five GPUs were allocated. A separate active reservation holds
16 CPUs and two GPUs for another user through 2026-10-19, so a future
single-node MPI32 x OMP2 (64-CPU) proposal is not schedulable outside that
reservation. The pending four-GPU job 402942 cannot start with only three
physically unallocated GPUs; its scheduler reason remains `Priority`.

The first small integrated M5 preparation OMP1/OMP2 comparison, job 403046,
failed before evolution in its OMP1 leg (`FAILED 1:0`, syn101, 9 s). Its
effective batch environment exported `SNRT_RT_LEVEL=2`, and the runtime
explicitly rejected that level restriction with `SNRT live IR AMR requires
all levels: unset SNRT_RT_LEVEL`, followed by `Too many errors in the
namelist`. The executable's `srun` step itself exited 0; no M5 parity or
physics result follows. The partial run and log are preserved. Corrected job
403205 uses a new `.m5-prep-omp-compare-20260924-r2/` run directory, the
same SHA256 binary and byte-identical OMP1/OMP2 effective namelists, but
unsets `SNRT_RT_LEVEL` after sourcing the shared environment. It remains a
small integrated parity check, not a 128-cubed performance result.

Job 403205 then advanced through its first OMP1 M5 commit and main step with
`mcons=econs=0`, but failed on the second M5 material tile with status 280;
its OMP2 leg was not entered. The tile IR stage had succeeded. In the
`material_cell` packet-restore path, the Fortran `intent(out)` status was not
initialized before skipping the standalone IR call and testing `status/=0`.
The first tile happened to see zero, while the second read an undefined value.
The source now sets `status=0` after validating the packet. This explains the
observed rejection without relaxing any material admission criterion; a
rebuilt, rerun binary is still required to verify the correction. Job 403213
completed the rebuild (`0:0`) in the existing pinned 187-variable M5 build
profile, producing SHA256
`0e09eee5586a2a0044883c03b74539c78422e40421e2cac451d0b08cacc8613f`.
Job 403214 submits the same two effective namelists in a fresh
`.m5-prep-omp-compare-20260924-r3/` directory, with the corrected binary
and all-level IR environment. Both sequential legs completed (`0:0`):
OMP1 in 18m38s and OMP2 in 17m57s, with two M5 commits and `Run completed`
in each log. Both final HDF5 files were 101,018,160 bytes, and whole-file
`h5diff -q` returned zero (exact equality of all compared HDF5 objects).
The photon/energy/state ledger lines also match to printed precision, and
both runs report main-step-2 `mcons=-3.00E-07`, `econs=-3.27E-04`.
On this small integrated case, the stiff second-tile preparation fell from
101.240s (OMP1) to 44.998s (OMP2), while the enclosing M5 commit fell from
352.622s to 301.316s. These are measured timings, not a 128-cubed or
multi-GPU throughput claim. Effective namelists and logs are retained; the
two evaluated raw `output_00001` trees are recorded separately for cleanup.
The same tile's measured CHIMES photo worker sum was 99.892s (OMP1) and
87.770s (OMP2), versus dark-chemistry sums of 0.705s and 0.622s. Thus the
small illuminated case points to the spectral photo/CVODE path, not the dark
solver or IR stage, for the next bounded cost investigation. Worker sums and
tile wall times are different quantities; these numbers cannot be projected
to the 128-cubed cosmological run without its own split counters. In
particular, pending H200 MPI32 x OMP1 job 402942 measures rank scaling but
does not exercise the new per-cell OpenMP preparation loop.
A safeguarded-Newton replacement for the 128-node spectral maximum-entropy
reconstruction was tried locally and discarded: across 3,240 reconstructions
on one login-node core, the existing bisection took 0.153s and Newton 0.148s.
Its distributions matched to L1 <= 5.6e-12, but this ~3% kernel-only
difference is too small to justify changing the production closure. No such
solver change or trial test file remains in the worktree.
While H200 job 402942 was still `PENDING`, its run-directory binary was
atomically replaced by the same status-fix build and `launch.sha256` updated;
`sha256sum -c` and pinned CHIMES/SUNDIALS linkage pass. The preceding binary
is preserved as `ramses_m5_3d_pre_statusfix_0353a870`. This does not turn the
pending H200 scaling probe into an OMP parity or full-step result.
The Slurm request is explicitly `--partition=h200` and
`--gres=gpu:H200:4`; `squeue` reports `gres/gpu:H200:4` and `PENDING
(Priority)` as of 2026-09-25 00:29 KST.

Prestart input/output recheck on 2026-09-25: the effective namelist is
`/gpfs/kjhan/LRD_JWST/.galaxy128-m5-d03-h200x4-hybrid-mpi32-o1/run.nml`
(SHA256 `141eb841a742c320a5e9284ffebeee50d176efcf009bfd993e6308c6e9c95be0`).
It requests one intentionally unreachable full dump at `aout=1.1`, `tout=1e100`,
`foutput=fbackup=1000000`, and `nstepmax=1`: zero routine dumps are expected
from this short scaling probe. An unexpected 128-cubed base dump could still
be of order 200--300 GiB using the earlier M5 estimate, with AMR growth not
bounded by that figure; GPFS showed 253 TiB free (shared, not reserved).
The pending batch already checks `launch.sha256` at execution. That manifest
now covers the executable, effective namelist, both sourced environments,
stellar/SNIa/yield inputs, three radiation/dust contracts, CHIMES main/atomic/
molecular banks and all nine group HDF5 files; all 26 entries pass. The
pre-enriched IC is the previously checked 22-field manifest, not rehashed on
the login node. The group and SED contracts explicitly say
`reference_control`, so this run is a wiring/scaling probe, not an approved
science-production SED or calibration result.
The production-input audit found no `approved_production` stellar/group
namelist among the active effective-population inputs or
`simulation/snrt/config/` controls. The active mixed PARSEC v5 SED has
`status='reference_control'` and an explicitly zeroed unmeasured late tail;
the nine-group contract has `contract_status='reference_control'`, and the
active DL01 dust contract has reference opacity/thermal/IR statuses. The
independent BPASS alternative is likewise a reference comparison, not a
matching effective-SSP replacement. Native
`snrt_spectral_contract:runtime_status_allowed` requires the explicit
`SNRT_ALLOW_REFERENCE_CONTROL=1` opt-in for these controls. A later science
run therefore needs a physically justified, mutually consistent SED/group/
dust contract and actual approval/identity evidence; relabeling these files
would not satisfy that requirement. This is a physical-input blocker for
production claims, separate from job 402942's code/scale measurements.

## 2026-09-25 H200 large-step feasibility result

Job 402942 became runnable on syn104 and was explicitly confirmed as the
operator-authorized `--gres=gpu:H200:4`, MPI32xOMP1 128^3 M5 run. Effective
NML `/gpfs/kjhan/LRD_JWST/.galaxy128-m5-d03-h200x4-hybrid-mpi32-o1/run.nml`
has SHA256 `141eb841a742c320a5e9284ffebeee50d176efcf009bfd993e6308c6e9c95be0`;
binary SHA256 is `0e09eee5586a2a0044883c03b74539c78422e40421e2cac451d0b08cacc8613f`.
The launch manifest passed in scheduler stdout. It uses an optimized SUNDIALS
runtime and the existing CHIMES rate cache, but predates the 2026-09-25
halo/regrid MPI tiling fix. This is a reference-control wiring/scaling probe,
not an approved stellar/group/dust science-production input set.

The first coarse step chose `dt=1.782e-1` code units under
`aexp_step_limit=.01`, unlike the old tiny-step preflights at `1e-6`.
Level 7 has 262144 grids globally, 8192 per rank (65536 base cells per rank,
about 256 material tiles of 256). Rank 1 completed only tiles 1--12 in
1h42m; each completed tile consistently reported about 203--205 s material
preparation, 45--48 s material stage, and 197--201 s CHIMES dark worker
time. Even an optimistic 250 s/tile without startup, communication, or
variation is ~17.8 h for 256 tiles, versus the 8 h allocation. The dark
number is worker wall within preparation, not additional CPU time or an
exclusive diagnosis of CVODE arithmetic. M5 commit, main step 1, and a
checkpoint were absent. This run cannot establish 128^3 throughput or
physics completion. Do not compare its large-dt chemistry rate directly
with prior tiny-step fixed-level tests.

To avoid consuming the remaining six-plus hours of four H200 cards without
a reachable first-step result, the owned job was canceled at 1h42m40s:
Slurm `CANCELLED by 10396`. No `output_*` or `jobcontrol.txt` was created;
thus there is no checkpoint to restart and no raw simulation dump to delete.
Effective inputs, launch manifest, binary, scheduler output, and
`live-402942.log` are preserved. A future large-step attempt needs a
measured chemistry/material cost reduction, a binary with the new regrid
exchange, and a wall-time estimate based on its actual timestep. The
physical-input approval blocker above remains independent.

### Substep correction and bounded next diagnosis

The ~17.8 h projection above is only for **one material substep**, not for
the first coarse step. `snrt_moment_ramses.f90` sets
`nsub=ceiling(12*chat*dt/(mesh%dx*sl))` and calls `material_tile` once per
substep. At `a=0.01`, `H0=67.4 km/s/Mpc`, a 12.5 cMpc/128 base cell, logged
`dtnew=0.1782`, and `chat=0.01c`, `units.f90` gives `scale_t=4.578e13 s`,
`dt=8.158e12 s`, physical cell width `3.013e21 cm`, and a CFL quotient
`9.74`. Thus the first coarse step has **at least ten** material substeps
at level 7, with `subdt≈8.158e11 s`. At unchanged measured 250 s/tile,
ten passes over 256 tiles/rank would take ~178 h even before startup,
transport, AMR and later-level work. This is an extrapolation, not a measured
completed-step runtime; level/refinement and chemistry evolution can change
the cost. It strengthens the cancellation decision, not a production timing
claim. The existing native 1e10 s fixture is not a matched cell comparison.

Fable read-only plan advice returned Q-GOAL/Q-LEAN and recommended a single
representative-cell solver-cost experiment before a Jacobian implementation.
The driver adopts the diagnostic ordering, **not** its unmeasured conclusion
that a new Jacobian is the next fix. The advice compared a 1e10 s neutral
fixture with the live large-step cell and mistakenly treated the logged
`dtnew` as the CHIMES interval; the actual CHIMES interval is `subdt` above.
Likewise `prep_s` includes `chimes_dark_s`: at OMP1, ~199 s of ~203 s prep
is already attributed, while the separate ~45--48 s material stage is the
next measured cost. GPU dust/material batches ran in that stage; the GPU was
not wholly idle. The bridge's three-attempt tighter-tolerance retry exists,
but its live frequency is unmeasured. An exact-interval single-cell profile
should count network/CVODE calls, steps, Jacobian RHS, LU setups and failures
before any physics, tolerance or solver-code change. Do not launch another
whole-box H200 run merely to gather those counters.

The first bounded native profile is now complete, but it is **not** a live-cell
reproduction. Jobs 403605 and 403606 ran on an A100/H100 allocation with one
CPU worker, using the existing diagnostic `chimes_cost_probe.so` and an
optional target interval in `chimes_dark_mixed_test.f90`. Both completed 0:0
with native element/charge admission. Job 403606 used `subdt=8.158275e11 s`,
`nH=0.2 cm^-3`, `T=165.4 K`, no dust/radiation, and the actual IC's eleven
element mass fractions converted to number ratios (but a manufactured neutral
species state). It measured 0.01552 s for the receiver call and 0.01166 s
inside one `chimes_network` call; CVODE accepted 197 steps, 217 ordinary and
632 numerical-Jacobian RHS evaluations, four Jacobians, 31 LU setups, one
error-test failure, zero nonlinear-convergence failures. There was one
network invocation, so no bridge retry in this cell. Job 403605's less exact
C/O-only 100 K variant took 0.01055 s. The first submission, 403602, failed
before computation because `LD_PRELOAD` was inherited by `srun` itself;
403605 moved preload to the target process and passed. Logs and binary hash
remain in `.snrt-performance.jDx9Iz/target-subdt-20260925/`.
Job 403613 rebuilt after the final input/guard edit, repeated the target
profile with unchanged solver counters, and passed the ordinary
`DARK_MIXED_THREAD_PARITY_CONSERVATION_PASS` regression; it completed 0:0.
No RAMSES output was scheduled or generated by these native probes.

These 0.01--0.02 s diagnostics are ~50x cheaper than the live average
~0.78 s/cell; matching `subdt` and initial elemental mix is insufficient to
explain the live hot path. In particular, four Jacobians and zero nonlinear
failures in the manufactured neutral state do **not** justify implementing
Fable's proposed analytic Jacobian yet. The live pre-photo species state,
temperature/density distribution, and bridge retry frequency remain
unmeasured. Preserve existing solver tolerances and conservation admission;
the next profile must sample an actual live input before selecting a solver
optimization or an effective-chemistry approximation.

Source inspection narrows one hypothesis: `snrt_chimes_photo.cpp:photo_step`
copies the species and ray state unchanged when validated input radiation has
zero total energy. The z=99 run has no initial luminous particles/AGN in its
IC and initializes source arrays to zero; the CMB is a separate thermal bath,
not a nine-group stellar photon payload. Subject to the actual stored M5
payload being zero as expected, the photo pre-stage cannot explain a ~50x
dark-chemistry difference by itself. Do not reinterpret the native neutral
probe as the actual live cell: its gas temperature, density and advected
species were manufactured, and a bridge retry can still change cost.

The first-step call order also passes through RAMSES `cooling_fine` before
SNRT, but this is **not** evidence of duplicate gas cooling. The actual
Makefile object is `cooling_fine.kjhan.o`, resolved from
`patch/cuRamses/cooling_fine.kjhan.f90` (not `hydro/cooling_fine.f90`). Its
`coolfine1` branch sets `delta_T2=0` when dust cooling is
`chimes_neq_v1` or `snrt_hhe_cie_metals`; the ordinary `solve_cooling`
branch is skipped. The input namelist's `cooling=.true.` is required by
the dust admission gate, but the thermal sink belongs to the native SNRT
receiver in this mode. Any post-hydro temperature difference from the
manufactured native cell remains possible; a second legacy cooling sink
is not the explanation for the measured CHIMES cost.

### H200 CPU accounting anomaly (unresolved)

Slurm `sacct -j 402942` reports 32 allocated CPUs but only `01:48:43`
`TotalCPU` for the 1h42m46s `srun` step (`AveCPU=00:03:22` across 32 tasks).
If the rank-local chemistry really occupied all 32 allocated cores during
that interval, aggregate CPU time would be far larger. This is a **question**,
not proof of single-core execution: Slurm accounting might omit Intel MPI
child CPU time, or only one rank might dominate while others wait. The
native single-cell/live-cell ~50x discrepancy makes this worth distinguishing
before changing CHIMES chemistry or its Jacobian.
The earlier H200 job 400507 (MPI16xOMP2) likewise reports `TotalCPU=04:04:37`
for a 4h `srun` step despite measured rank-local progress through 135
material tiles. Repetition across two distinct MPI layouts weakens any claim
that job 402942 alone was accidentally bound to one CPU; it also raises the
possibility that H200 step CPU accounting here reports only a subset of rank
processes. Neither explanation is established without affinity/per-rank CPU
measurements.

A three-CPU-second-per-rank affinity/accounting probe first requested four
H200 GPUs as job 403622. It remained pending under `QOSMaxGRESPerUser`:
`sacctmgr` reports the `normal` QoS per-user GPU cap as eight, and six were
already in use. This is a CPU-binding probe, so the owned unstarted job
403622 was canceled and replaced with job 403645 requesting **one** H200,
the same MPI32xOMP1 CPU shape and CPU binding. The reduced GPU count means
it will not reproduce four-GPU contention, but the CPU-affinity/accounting
question remains testable. At last observation 403645 was pending for
priority; no result or hardware conclusion exists yet. Its bounded script
and output target are `.snrt-performance.jDx9Iz/h200-affinity-20260925/`.
Neither probe reads nor modifies the other running jobs.

### 2026-09-25 bounded H200 diagnosis result

The one-GPU replacement 403645 completed 0:0 in five seconds. All 32
`srun --cpu-bind=cores` tasks saw the same **32-core allowed set** rather
than an erroneous one-core mask; each consumed 3.000 CPU seconds in about
3.01 wall seconds. Slurm recorded `01:37.653` step TotalCPU, consistent
with 32 simultaneous busy workers plus launch overhead. The common allowed
set is not exclusive per-rank affinity, but it has enough cores for these
32 tasks and they demonstrably ran in parallel.

The Intel-MPI-linked Fortran probe 403652 then reproduced the original
`srun --mpi=pmi2 --cpu-bind=cores` launch path on H200. It completed 0:0
in seven seconds. `MPI_Reduce` measured 96.0007 summed rank CPU seconds
for 32 ranks, maximum rank wall 3.4051 s; Slurm recorded `02:04.465`
step TotalCPU including MPI startup/finalization. This rules out a general
H200 or Intel-MPI accounting failure as the explanation for the
`TotalCPU≈wall` records of jobs 400507/402942. It does **not** identify
which application phase left most ranks idle.

Finally job 403665 repeated the existing representative *manufactured
neutral* CHIMES cell 256 times on each of 32 H200 ranks, with exactly the
same `8.158275e11 s` substep interval, 165.4 K, nH=0.2 and zero radiation.
It completed 0:0 in eight seconds; all 32 rank records passed element/charge
admission, each reporting approximately 3.86--3.89 s for 256 solves
(~0.0151 s/cell). Slurm step TotalCPU was `02:12.307`. Thus CHIMES itself
does not exhibit a 50x slowdown merely from 32 concurrent ranks, H200 CPU
placement or reading its pinned tables. The *actual integrated* cell path
still measured ~0.78 s/cell inside `chimes_last_split_wall(2)`; this
discrepancy is now specific to live state/path, bridge retries, or
application-level waiting rather than generic MPI/CHIMES concurrency.

All three diagnostic jobs used a single H200 only to satisfy the cluster
GPU-job requirement; no CUDA computation or RAMSES dump was requested.
Logs/scripts are under `.snrt-performance.jDx9Iz/h200-affinity-20260925/`
and `.snrt-performance.jDx9Iz/target-subdt-20260925/`. The old four-GPU
402942 run should **not** be relaunched unchanged; the next evidence must
capture one real live cell's receiver controls/state or process CPU around
that receiver before selecting a chemistry solver change.

### 2026-09-25 actual-cell capture queued

The 256-repeat H200 fixture above re-solves the *same manufactured initial
state* each time; it is a concurrency/CPU-accounting control, not an evolved
cell trajectory or a reproduction of the full molecular receiver. Source
inspection confirms the live `chimes_last_split_wall(2)` encloses the dark
transition and, when a thermal root is crossed, atomization plus a second
atomic dark solve. The native fixture calls only one direct dark transition.
No chemistry-tolerance or conservation rule has been changed on this basis.

An opt-in `SNRT_CHIMES_CAPTURE_ONE=1` diagnostic now records exactly the
first rank-1 live molecular receiver input (nine controls, eleven normalized
elemental abundances, 157 normalized species abundances, and incident photon
number/energy sums) before the chemistry call. The capture is OpenMP-guarded,
flushed, and inactive by default; its Fortran module compiled successfully.
Job 403698 started on syn104 with a single H200, MPI32xOMP1,
320 GiB allocation. The run uses the existing 128^3 one-step reference input
in a new directory and stops as soon as the first capture appears, with no
scheduled dump (`noutput=1, aout=1.1, tout=1e100, foutput=fbackup=10^6`; current
`aexp=0.01`, `nstepmax=1`). The effective namelist is
`.galaxy128-m5-d03-h200x4-hybrid-mpi32-o1/run.nml`; the job script/log target
is `.snrt-performance.jDx9Iz/live-cell-20260925/`. At submission the GPFS
scratch mount had ~253 TB free. No capture or runtime result exists yet.

The existing native dark test now accepts `SNRT_CHIMES_REPLAY_FILE` to parse
these four capture records and run the *same* species/controls in the direct
dark solver when incident radiation is exactly zero and the cell is below
the molecular temperature ceiling. It refuses other inputs rather than
claiming they are a valid direct-dark reproduction. Its unchanged default
parity/conservation path passed as job 403704 (COMPLETED 0:0 in two seconds);
the test-binary SHA256 was
`ee8543eacc3f335d556a6a57964ca9c2293b267f08f05389b8f844fbde934072`.
No actual-cell replay result exists until job 403698 reaches its first cell.

While 403698 remained live, source inspection found one avoidable setup cost:
the active `adaptive_loop.jaehyun` built the RAMSES CIE table at startup and
the active `cooling_fine.kjhan` rebuilt it each coarse step even when
`dust_cooling='chimes_neq_v1'`. In that mode `coolfine1` deliberately sets
the RAMSES CIE `delta_T2` to zero, while the native CHIMES receiver owns the
thermal update. The two table-build call sites now skip only this admitted
CHIMES mode; initial gas temperature setup and the `cooling_fine` energy/floor
pass are unchanged. Both changed Fortran objects compiled with the active
HDF5/SNRT/dust-dynamics/CHIMES/CUDA flags. This is **not yet an integrated
runtime validation** or a measured speedup; job 403698 runs its already
linked earlier binary, and a later bounded coupled run must verify the
optimized path before making such a claim.

### 2026-09-25 real-cell replay and first-tile follow-up

Job 403698 completed 0:0 in 41m20s after its script observed exactly one
rank-1 live molecular receiver input and stopped the `srun` step on purpose.
Its step is therefore `CANCELLED`, not a simulation completion claim. The
captured cell 18447 had nH=0.1653613771 cm^-3, Tgas=158.8543648 K,
Tdust=1 K, subdt=8.307003378659e11 s, zero grain ratio, and exactly zero
incident photon number/energy. Its 157 species were neutral (with the actual
eleven-element inventory); no output dump was made. The captured log SHA256
is `3672825ee02619ab4823a832f575cf322a5bbf7e6b26e4064d21d11420b638c3`.

Job 403735 replayed **that exact input** in the direct molecular dark solver
on A100/H100. It completed 0:0: plain and CVODE-profiled receiver wall were
0.01485 and 0.01471 s respectively. The profiled call had one CHIMES network
integration, CVODE wall 0.01440 s, 187 accepted steps, four Jacobians, one
error-test failure and no convergence failure. The accepted temperature was
158.8543648 K. This directly rules out the first live cell's neutral state,
temperature, dust temperature and exact substep as sufficient explanations
for the old tile-average ~0.78 s/cell. It does **not** prove that all cells
are cheap: rare costly cells or application-level contention remain possible.

The existing per-cell dark wall array now emits one first-tile summary:
maximum dark-cell wall and counts over 0.1 s and 1 s. The amended Fortran
object compiled. Job 403742 started with one H200, MPI32xOMP1 and the same
one-step reference namelist, but its first-tile timing is **invalid as a
parallel benchmark**: a read-only `/proc` inspection on syn104 found all 32
`ramses_final3d` processes confined to CPU 12, despite 32 CPUs being
allocated (`cpuset.cpus.effective` listed 32 CPUs; `cpu.max` was unlimited
and `nr_throttled=0`). Slurm `sstat` showed about 1.5 minutes of average CPU
per rank after 45.5 minutes elapsed; the H200 was sampled at 0% utilization.
The old Slurm/OpenMP binding combination is thus consistent with a one-core
placement error, rather than evidence for slow CHIMES chemistry or a GPU
bottleneck. No first-tile timing appeared, and own job 403742 was canceled
to avoid further wasted allocation. No scheduled dump was made; the log
did not show the old `Computing new cooling table` message, but absence
alone is not full runtime validation.

The bounded first-tile script now overrides OpenMP binding and launches
`srun --cpu-bind=none`, preserving MPI32xOMP1 and the identical physics
input. Job 403780 completed 0:0 in 1m42s (srun step was intentionally
canceled after the diagnostic appeared). Five live RAMSES processes sampled
after startup each had access to the entire 32-CPU allocation, rather than
CPU 12 alone; Slurm recorded 35m35s of total CPU in 1m42s wall. Rank 1's
first 256-cell material tile took 3.931 s preparation, 4.786 s stage,
0.005 s commit. Summed CHIMES dark wall was 3.725 s; the slowest dark cell
was 0.017739 s, with zero cells over 0.1 s or 1 s. The next three rank-1
tiles had preparation times 3.949, 3.948 and 3.943 s, with dark sums
3.746, 3.744 and 3.740 s. The log contained no `Computing new cooling
table` message, no unexpected raw dump was produced, and the bounded script
printed `FIRST_TILE_DARK_DISTRIBUTION_PASS 403780`. The log SHA256 is
`0b72bb676aec153f3834468798488c17f0a0902da981b514a287d1d0ecb273a7`.
This validates the corrected parallel affinity and first-tile chemistry
cost for this one reference input; it is not a completed coarse step, a
strict isolated CIE-table speedup measurement, or production-readiness proof.

### 2026-09-25 bounded full 128^3 coupled step

The corrected affinity now permits a meaningful full-step test. Job 403794
was submitted on syn104 (one H200, MPI32xOMP1, 320 GiB, eight-hour cap) with
a private executable copy SHA256
`906f4bca3135ae4f4ba729854879a9701d8ed1ddd9c84d5c9cc860e19c229152`
and effective namelist SHA256
`141eb841a742c320a5e9284ffebeee50d176efcf009bfd993e6308c6e9c95be0`.
Run directory: `/gpfs/kjhan/LRD_JWST/.galaxy128-m5-d03-h2001-m32-affinity-20260925/`.
The job uses the proven `--cpu-bind=none`/`OMP_PROC_BIND=FALSE` combination,
and disables per-tile performance logging. The source/IR contract remains
`NONPRODUCTION` reference control. This is a one-coarse-step coupled
engineering gate, not a z=6 science trajectory. Its dummy future `aout=1.1`,
`tout=1e100`, and million-step output/backup intervals request zero raw
dumps; `/gpfs` had 252 TB free before launch. Acceptance requires batch
completion, `Run completed`, no critical error/MG/NaN marker, and no
`output_*` directory. Even a pass would not approve the physical reference
sources or longer evolving/restart behavior, and this run itself is not an
illuminated-source test. Existing bounded source-active M5/CHIMES/D03 and
restart evidence is retained in
[`snrt_m5_full_connection_2026-09-15.md`](snrt_m5_full_connection_2026-09-15.md);
it need not be rerun merely to compensate for this deliberately source-dark
cost gate.
The end evaluation will also read the existing `SNRT_MN_COMMIT`,
`SNRT_MN_LEDGER`, `SNRT_MN_ENERGY`, `SNRT_MN_STATE`, `NaN_CHK`, and stage-timing
lines: require finite/nonnegative physical inventories and energies, xHII in
[0,1], zero NaN counters, and no native accounting rejection. The native M5
precommit number/actual-energy residual limit is 1e-9 of its global inventory
scale; no post hoc tolerance change or extra diagnostic framework is needed.
If an unexpected full dump occurs, the earlier M5 estimate is roughly
200--300 GiB for the 128^3 base state, before additional AMR growth; the
run-directory postcondition treats even one such dump as a failed gate.
The effective namelist explicitly disables `agn`, `smbh`, `sink_AGN`, CR,
and delayed cooling. It enables channel-resolved stellar enrichment and
star formation, but a single first step at `aexp=0.01` does not itself prove
that a star formed or any stellar photons/feedback were injected. The
acceptance scope is therefore gas/dust/RT/CHIMES execution and cost, not
AGN or source-active feedback validation.
At 1m46s elapsed the job was RUNNING on syn104. Its 32-rank `sstat`
average CPU was 1m42s (consistent with restored parallel execution);
average RSS was about 7.3 GiB/rank. The integrated log had entered SNRT
RT/dust-scattering work but had not yet reached `Run completed`.
At 39m35s elapsed, all ranks still accumulated CPU during the first fine
step: average CPU 38m51s, minimum CPU 38m26s, maximum CPU 39m02s, with
MaxRSS 7.73 GiB. This rules out the prior one-core affinity failure, but
CPU time alone cannot distinguish useful chemistry/transport work from
active MPI or GPU wait polling. Exact late-tile cost remains unmeasured.
No new step marker or error had appeared at that point.
At ~53 min, one brief read-only GDB attach to an owned rank on syn104
(PID 654000) sampled `snrt_moment_transport::mn_project` at the weighted
harmonic `matmul`, called by CPU fallback of `mn_dispatch_project` from
`snrt_dust_live_moment_tile` and `snrt_ramses_advance_level`. GDB detached
normally; the simulation remained running. This **one call-stack sample**
confirms useful material/M5 projection work at that instant, not the
fraction of total wall spent there or a controlled CPU-vs-GPU comparison.
The original four-H200 reference run script was also corrected for future
reuse: it now overrides OpenMP placement and uses `srun --cpu-bind=none`.
Its environment was not changed; `launch.sha256` now additionally pins the
modified script (27 entries). The script passed `bash -n`, and its recorded
SHA256 was checked directly; the full asset manifest is rechecked inside
any future batch, not rehashed on the login node. This does not alter job
403794, which uses its separate frozen binary and script.
At 1h15m, a second brief read-only sample of the same owned rank (PID
654000) found it in `snrt_dust_ir_advance`, called from
`snrt_dust_live_stage` and `snrt_dust_live_moment_tile`. GDB detached and
the job remained RUNNING, with average CPU 1h14m36s and maximum per-rank
RSS 7.73 GiB. The optimized frame did not expose a tile index, so these
two call-stack samples establish changing active operators, not a measured
tile progression or ETA. No new step marker or error had appeared.
The effective first fine-step `dt=1.782e-1`, `aexp=0.01`, 8.425 h^-1 Mpc
box, 128-cell width and `reduced_c=0.01` imply
`chat*dt/(mesh%dx*scale_l)≈0.812` from the actual cosmological `units`
routine. The live M5 code chooses `nsub=ceil(12*ratio)=10`, so this one
hydrodynamic step repeats transport and the CHIMES/dust material tile sweep
ten times. Extrapolating the measured first 256-cell tile across 256 tiles
per rank and ten substeps yields about 6.2 hours of tile work before other
transport, communication or load imbalance. That is a rough lower-scope
cost indication, NOT a measured full-step ETA or grounds for changing
physics/timestep settings mid-run; the eight-hour wall cap remains unchanged.
Code-path check while job 403794 remains active: the M5 tile passes
`incoming_radiation` to `snrt_dust_live_stage`, selecting its
`local_material` path. That path sets the internal IR substep count to one
and skips its separate halo/transport sweep; the tenfold repeat above is
the outer conservative M5 subcycling, not an additional IR subcycle
multiplier. The per-cell preparation and several IR loops have OpenMP
pragmas, but this run requests MPI32 x OMP1, so any speedup from those
pragmas is not measured here. The first-tile 3.9 s preparation and 4.8 s
stage are still the relevant bounded cost observations. Do not infer a
multi-thread or multi-GPU speedup from this static inspection.
The effective `run.nml` does not set `n_cuda_streams`; the RAMSES default is
one stream **per MPI rank**, passed to `snrt_backend_initialize` during
parameter reading. With MPI32 and one allocated H200, this must not be
described as one or three streams per card. `SNRT_BACKEND=hybrid` may send
supported 384-direction M5 closure/projection and 32-direction face tiles
to an available shared stream, otherwise to CPU. The rank sample in
`mn_project` demonstrates one CPU fallback, not a measured fallback rate
or GPU utilization. No stream tuning is inferred or applied to the live run.
At 2026-09-25 09:56:05 KST, a single `nvidia-smi` query launched as an
overlapping read-only step inside job 403794's own allocation reported its
H200 NVL at 89% GPU utilization and 17,903 MiB used. This confirms GPU
activity at that instant, not a time-average, M5-only attribution, hybrid
CPU/GPU dispatch fraction, or a speedup. The simulation step was not changed.
Code inspection of the live configuration clarifies the hybrid admission
semantics. `cuda_stream_pool.cu` owns one stream pool per MPI process and
`cuda_acquire_stream()` tests only that process's local `busy` flags;
`snrt_hybrid_try_acquire_c()` also checks device free memory divided by the
number of local ranks sharing the GPU. Neither check tests whether other MPI
ranks currently occupy the H200's compute units. Thus this is rank-local
stream/memory admission with CPU fallback, not the requested card-wide
"GPU idle, otherwise CPU" arbitration. The present one-H200/MPI32 result
cannot establish dispatch efficiency or card-wide scheduling. Treat this as
an open runtime-policy limitation, distinct from the conservation/physics
gate; assess it with the completed step's cost evidence before changing the
live allocation or adding cross-rank coordination.
During the live first coarse step, rank 1 reported
`SNRT_MN_COMMIT level/order/substeps=7/5/10 wall=12396.753 s`, a zero photon
and actual-energy residual in this source-dark control, max xHII
`8.1730786833908e-45`, and positive gas energy `1.5551740961900e-08`
code units. The reported gas-heat ledger term is negative
`-1.9295081588235e47`, consistent in sign with cooling but not by itself
a full gas-energy audit. The rank-1 dispatch counter was 89,128,960
closure cells on CPU and zero on GPU; face counters were zero. A separate
single-time H200 utilization sample therefore cannot be credited to GPU M5
closure. The run was still active after this marker, so this is not a pass.
An `output_00000` directory appeared at star-formation time despite the
dummy future snapshot schedule. Its 32 listed entries were only
`stars_00000.out00001` through `.out00032`, all observed zero bytes; no
AMR/hydro/particle snapshot file had appeared at that point. The source
`star_formation.kjhan.f90` opens birth-record files when
`sf_birth_properties` is enabled even if no stars form. The submitted runner's
postcondition rejects *any* `output_*` directory, so it will report a false
"Unexpected raw dump" failure if this remains the only output. Do not
delete the birth files while the run is active. At terminal state, classify
the actual directory contents and physics markers separately from that
overbroad shell postcondition; amend the next launcher to distinguish birth
records from full snapshots. No current-job script or namelist was changed.

### Terminal evaluation of job 403794

Slurm ended `FAILED 1:0` after 03:30:51, but its `403794.0` simulation step
ended `COMPLETED 0:0` after 03:30:50. The application logged `Run completed`
at 12605.041 s, one main step with printed `mcons=econs=0`, one M5-L7 commit,
zero reported photon/actual-energy residual, finite state values and zero
printed `NaN_CHK` counters. The batch exit 1 came from the shell's
`Unexpected raw dump` test. Exact inspection found 32 zero-byte
`stars_00000.out00001`--`.out00032` birth-record files in `output_00000`
and no AMR/hydro/particle snapshot there. This is an output-classification
error, not a failed simulation process or a 200--300 GiB raw snapshot.

The run is **not** a coupled AMR/production pass. At the next level-8 entry
all 32 ranks emitted `diag_check_eint` negative-internal-energy reports,
with local minima roughly `-1.866e-10` to `-4.696e-10` code units and tens
of thousands of affected cells per rank. This is an actual invalid gas
state despite `Run completed` and the L7 radiation ledger. The effective
`run.nml` omitted `interpol_var` and `interpol_type` in `REFINE_PARAMS`;
the compiled hydro defaults are `interpol_var=0`, `interpol_type=1`.
Earlier producer isolation in
[`m5_l8_interpolation_boundary_evidence_2026-09-16.md`](m5_l8_interpolation_boundary_evidence_2026-09-16.md)
found this `interpol_var=0` path negative at L8 entry and the explicit
`interpol_var=1, interpol_type=1` path clean at that checkpoint. The new
run therefore repeated a known inadmissible effective input; it does not
show that the documented producer-side repair failed. The next full-coupled
attempt must pin both interpolation values in its *effective* namelist,
check the L8 entry and complete the L8 M5/material transaction before any
128^3 science or larger-resolution readiness claim. Preserve job 403794's
original input/script/log identities as failed engineering evidence; do not
reinterpret the source-dark, reference-control run as calibrated physics.

### Corrected bounded L8 submission

Job 404118 was submitted 2026-09-25 10:16 KST and allocated immediately on
syn102 (`a100`, `gpu:A100:1`, MPI30 x OMP1, 512 GiB, 48 h cap). This is a
bounded hydro/M5/dust/CHIMES L8-entry/transaction check, not an illuminated
science trajectory or an AGN-feedback test. The new, non-overwriting run
directory is
`/gpfs/kjhan/LRD_JWST/.galaxy128-m5-l8-a1001-m30-20260925.4Ia0pF/`.
Its effective `run.nml` SHA256 is
`576c6929e1417f0cb24141edb0ebd00f0318f47336542427084088b7e6320a8a`;
the frozen executable remains
`906f4bca3135ae4f4ba729854879a9701d8ed1ddd9c84d5c9cc860e19c229152`.
The copied executable contains `sm_80` and `sm_90` cubins, including the
A100-capable target. This run changes only the bounded coarse-step count to
2, limits AMR to levels 7--8, and explicitly sets
`interpol_var=1, interpol_type=1`; the source-dark reference physics contract
and other namelist choices are retained. `SNRT_HYDRO_ENTRY_STOP=2` requests
normal bounded stop after a successful level-8 M5 transaction.

The effective output schedule has `noutput=1`, unreachable `aout=1.1` and
`tout=1e100`, `foutput=fbackup=1000000`: zero full dumps expected. A surprise
full 128-cubed dump would be of order 200--300 GiB before AMR growth, and
shared GPFS had 252 TB free at preflight. The new runner accepts only
`output_00000/stars_00000.outNNNNN` birth records, rejecting other output
artifacts; unlike 403794, it does not equate a birth-record directory with
a full snapshot. The batch pins binary, effective namelist and sourced
environment hashes, and checks the bounded L8 completion marker and fatal
diagnostics. One H200/MPI32 alternative had a scheduler test-start on
2026-10-01 even at an 8 h cap; A100/MPI30 fit syn102's available CPU and
started immediately. The A100 run is for physics/throughput qualification,
not a same-hardware performance comparison to H200 job 403794. Final
numerical and raw-output evaluation is pending its terminal state.

### M5 GPU dispatch root-cause check

The frozen executable used by 403794 and 404118 contains the CUDA device
symbols `snrt_mn_cuda_closure_c` and `snrt_hybrid_try_acquire_c`, but this
alone does not wire them into the M5 Fortran caller. `nm -u` on the current
`bin/snrt_moment_dispatch.o` shows no reference to the CUDA acquire or M5
closure/project/flux functions. Disassembly of the frozen executable's
`snrt_moment_dispatch_mp_mn_dispatch_closure_` OpenMP worker contains CPU
`mn_reconstruct`/`mn_project` calls and no CUDA acquire/closure call;
projection and flux entry points likewise contain no CUDA operator calls.
The source wraps those calls in `#ifdef HYDRO_CUDA`, whereas `bin/Makefile`
adds `-DHYDRO_CUDA` only under `USE_CUDA=1`. This proves that the linked
M5 dispatcher in this binary is CPU-only even though CUDA objects are also
linked and dust/IR batches use the GPU. The precise build-history mistake
(e.g. reuse of a CPU-built object across flag changes) is not established
from this binary alone. It is **not** evidence that card-wide occupancy
arbitration sent 89 million closure cells to CPU. Do not alter job 404118's
frozen binary; it remains a valid CPU-dispatch M5 physical boundary test,
not a GPU-M5 throughput test. Before a later GPU performance claim, rebuild
the dispatcher in an isolated CUDA-enabled profile, verify its undefined
CUDA references before linking, and exercise a bounded GPU/CPU parity case.
The first part of that correction was checked without changing the live
job: an isolated Intel `mpiifx -qopenmp -fpp -O0 -DHYDRO_CUDA` compilation of
the current dispatcher source succeeded in
`.mn-dispatch-cuda-check.4HcC6k/`. Its object has undefined references to
`snrt_hybrid_try_acquire_c` and all three M5 CUDA operators, as required.
No full executable was rebuilt or benchmarked by this isolated compile.
To prevent a repeat of this mixed-profile link, `bin/Makefile` now checks
the M5 dispatcher object immediately before a `SNRT=1 USE_CUDA=1` RAMSES
link. It requires the `snrt_mn_cuda_closure_c` undefined reference and
fails with an explicit rebuild instruction if the object is CPU-only. No
VPATH order or compiler flag was changed. `make -n` shows this check before
the link; the exact predicate rejects the current stale object and accepts
the isolated CUDA-built object. This guard alone does not establish GPU
runtime behavior; the subsequent complete build and bounded parity result
are below.
The targeted search found `snrt_moment_dispatch.f90` as the only current
`patch/lagRamses/snrt*.f90`/`dust*.f90` source whose Fortran CUDA calls are
conditional on `HYDRO_CUDA`; `snrt_runtime_backend.o` already references
the hybrid/CUDA runtime symbols. The guard is therefore scoped to the
observed mixed-build gap rather than a new general instrumentation layer.

### Exact-vacuum CPU closure shortcut (future binaries only)

The source-dark H200 step called rank-1 CPU closure for 89,128,960 cells.
Inspection showed `mn_reconstruct` already returns zero angular intensity
and zero cache for an exactly zero moment vector, without changing its
warm-start hint; `mn_project` of that zero intensity is zero. The CPU arm of
`mn_dispatch_closure` now performs that same assignment directly for an
all-zero tile or an all-zero cell within a mixed tile. Nonzero, non-finite
and unrealizable moments still take the existing solver/validation path;
forced-CUDA prelaunch failure behavior is untouched. The counter continues
to count CPU-dispatched cells, including exact-vacuum shortcuts, not matrix
projections actually executed.

The existing direct Fortran policy fixture was extended to check both an
all-zero tile and a mixed zero/isotropic tile, including projected moments,
angular output, cache and unchanged warm hint. It prints
`MN_DISPATCH_POLICY_PASS` in both CPU-only and `-DHYDRO_CUDA` builds with
`-fcheck=all`. The current source also compiles with Intel `mpiifx` and
retains the three M5 CUDA operator references under `HYDRO_CUDA`.
Artifacts are in `.mn-vacuum-policy-test.BAzEv8/` and
`.mn-dispatch-cuda-check.4HcC6k/`. No live job or numerical tolerance was
changed; job 404118 uses the earlier frozen binary and cannot measure this
shortcut. A measured speedup, if any, awaits a separately rebuilt bounded
run; do not extrapolate from the skipped algebra alone.

### Clean CUDA M5 build and bounded source-active CPU/GPU comparison

An isolated full `SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1`
`NVAR=187` build in `.m5-cuda-cleanbuild.jHpfd5/` completed as Slurm job
404126 (`COMPLETED 0:0`). The executable SHA256 is
`ef26fad594922f9972148b7d713388d95cfe64f298b73e88b2e1a0d07d20e61a`.
The freshly compiled dispatcher object has unresolved references to the CUDA
acquire and M5 closure/project/flux operators, all resolved in the linked
binary; `ldd` reports no missing libraries. The `patch/lagRamses` SNRT and
`patch/cuRamses` hydro VPATH selection was preserved. Build success alone is
not evidence of device execution.

The first paired job, 404128, failed its CPU step because its reused effective
input was from a `NVAR=199` dust-dynamics build; the clean binary was
`NVAR=187`. It was not a CUDA numerical failure. The failed CPU log is
retained under `.m5-cuda-parity.FnciQO/cpu/` and was not overwritten.

Corrected job 404130 used the same `NVAR=187` binary and the same effective
namelist (SHA256
`9e4c89f69918c96879666d01c65753c4626e8368d6e8242dbdf2cca5289f5009`)
for consecutive `SNRT_BACKEND=openmp` and `cuda` runs on one A100. The
effective input, `.m5-cuda-parity.FnciQO/effective.nml`, derives from the
existing source-active CHIMES/D03/AGN control, selecting `tensor_mn` order 5,
two steps and no periodic full dump. `noutput=1`, unreachable `aout=2.0`,
`tout=1e30`, and `foutput=fbackup=1000000` yielded zero `output_*` directories.
Both steps and the shell completed `0:0`; logs are under `cpu187/` and
`cuda187/` in that parity directory. Step 2 had nonzero source, absorption,
gas/dust exchange and M5 commit in both runs. The CUDA log identifies the
A100 and reports 24,064 closure cells and 368,640 face cells dispatched to
GPU at step 2 (CPU arm zero); the OpenMP run reports the same counts on CPU.

At source-active step 2, source photons match exactly at printed precision.
The relative CPU/GPU differences are `3.55e-11` for absorbed photons,
`3.15e-11` for absorbed radiation energy, `3.60e-13` for maximum xHII, and
`3.90e-10` for gas energy. Photon residual divided by injected photons is
`7.15e-16` on CPU and `2.36e-15` on GPU; energy residual divided by absorbed
energy is `8.45e-15` and `2.14e-15`, respectively. The active input is a
`reference_control` comparison, not calibrated source/IR/dust approval.
This is a two-step small-grid GPU-path and numerical-parity result, not a
GPU speedup claim, AMR L8 result, or full production-readiness certificate.
The CUDA run's in-simulation elapsed time (73.23 s) exceeded the CPU run
(68.32 s), and both spent about five minutes in startup, so this case is
especially unsuitable for a throughput claim.

### NVAR=199 dust-dynamics CUDA build and material-tile flag correction

The 128-cubed L8 job 404118 uses a `NVAR=199`/`DUST_DYNAMICS` executable, so
the `NVAR=187` CUDA pair above cannot be substituted into that trajectory.
An isolated `SNRT=1 DUST_LIVE=1 DUST_DYNAMICS=1 CHIMES=1 HDF5=1 USE_CUDA=1`
build was made in `.m5-cuda-dustdyn-build.wuQP7D/`. The first attempt, job
404135, encountered an Intel `icpx` 2025.3 optimizer crash (exit 139) while
compiling `snrt_hybrid.cpp`; the unchanged resumed build, job 404136,
completed `0:0`. The initial binary SHA256 was
`7d1743e35df5b3382e5afbd1dc9e2bf57c39816445b0991b09e6ef6920979f4b`.
Build flags show `-DNVAR=199 -DDUST_DYNAMICS -DHYDRO_CUDA`; the M5 dispatcher
object references all three CUDA operators and the linked binary resolves
them without missing shared libraries. This verifies the build profile but
not the physics behavior.

The existing source-active `NVAR=199` M5/D03/sublimation comparison input
(`.rt-m5-full.oDbped/cold-mpi-final/effective.nml`, SHA256
`747f30709c15238651a4f0e4e5eac1857f747a889f75cbfaaf07f58023ca1484`)
then exposed a real first-tile regression before CUDA dispatch. Job 404137
failed on one MPI rank; job 404142 reproduced the same failure on both ranks
with the old successful two-rank geometry. A failure-only diagnostic build,
job 404143, localized job 404144's `code=1` to the tile's supposedly
uniform feature-flag check, not a CHIMES calculation or GPU kernel. The
check combined `.neqv.` comparisons with `.or.` without parentheses. Because
Fortran evaluates `.or.` before `.neqv.`, it did not mean an OR of pairwise
flag mismatches and could reject uniform valid tiles. Each pairwise `.neqv.`
is now parenthesized; temporary failure diagnostics were removed. The frozen
404118 executable was not altered.

The fixed `NVAR=199`/CUDA executable in `.m5-nvar199-diagnostic.CA2eyv/`
(job 404145 `COMPLETED 0:0`) has SHA256
`b72ed7ae2a272e585e0a78b603dde8379c4ced78ff49b0d7002e7dd475f5f8d4`.
Paired two-rank OpenMP and forced-CUDA runs of the unchanged effective input
on one H100 completed as job 404147 (`COMPLETED 0:0`), with both M5 commits
and nonzero source/absorption at step 2; neither wrote a full output
directory. At step 2 the CUDA rank-1 counters report 21,632 closure cells
and 290,304 face cells on GPU, zero on CPU, with reciprocal CPU counts in
the OpenMP control. Source photons, final photon inventory, max xHII and
absorbed energy agree at printed precision. Absorbed photons differ by
`1.33e-14` relative and gas energy by `3.26e-13` relative. Photon residual
relative to source is `8.96e-17` on CPU and `2.56e-17` on CUDA; energy
residual relative to absorbed energy is `9.99e-15` and `4.40e-15`.

This closes the small source-active `NVAR=199` GPU-routing/CPU-parity defect;
it does not establish a speedup, an L8 GPU run, or approval of the
`reference_control` spectrum, group and dust contracts. Failed job logs and
pre-fix binaries remain preserved for reproducibility.

### Scientific-input approval boundary (operator confirmation, 2026-09-25)

The operator confirmed that no scientifically approved replacement files or
versions currently exist for the `reference_control` stellar SED, nine-group
contract, or dust opacity/IR inputs. Keep runs using these inputs classified
as technical/comparison controls even if their integration checks pass. This
is an explicit scientific-input blocker for a production-physics claim, not a
reason to relabel the controls as approved or to stop bounded code validation.
