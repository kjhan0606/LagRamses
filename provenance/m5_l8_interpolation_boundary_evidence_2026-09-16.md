# M5 level-8 interpolation producer evidence — 2026-09-16

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Bundle: `M5-L8-HYDRO-BOUNDARY`
Purpose: isolate the first negative hydro internal-energy state before the
fine-level Godunov/M5 transaction.

## Controlled comparison

Both jobs used a fresh copy of the same 12.5 cMpc cosmological IC and the
same bounded executable (`ramses_m5_boundary3d`, SHA256
`f82b957717aea84d9c5f1c4de536b31d6136a8e32c5ece137c5cdd60d296e973`).
The focused `SNRT_HYDRO_ENTRY_STOP=1` mode stops at the level-entry probe;
the temporary value `2` is reserved for a later M5-enabled run and stops only
after a successful level-8 M5 transaction. Neither mode clips or repairs a
state or changes the receiver semantics.
The test binary was built with the repository Makefile and current VPATH
order, using the coupled NVAR=199/NVECTOR=32 CPU configuration.

The two bounded runs deliberately disabled RT/dust/CHIMES so that the
hydrodynamic producer could be isolated. They therefore do not reproduce the
full 670735 physics contract. In both test namelists `interpol_type` was
omitted and consequently used the compiled hydro default `1` (limited
MinMod); it was not the old generator value `0` (straight injection).

| Job | `REFINE_PARAMS/interpol_var` | Purpose | Result | Elapsed |
|---|---:|---|---|---:|
| 362179 | 0 | baseline producer isolation | all 8 ranks reported negative `eint` at level 8 entry | 00:24:00 |
| 362215 | 1 | RAMSES internal-energy interpolation candidate | no level-entry negative diagnostic was emitted; clean bounded stop | 00:22:52 |

The baseline level-entry minima were approximately
`-3.683e-10` to `-4.694e-10` code units, with 214107--337985 negative
cells per rank. This places the producer upstream of fine-level Godunov and
the M5 material receiver for the tested trajectory. The `interpol_var=1`
run reached the same level-entry checkpoint without a reported negative
internal-energy cell.

## Interpretation and limits

This is evidence for the existing RAMSES internal-energy interpolation path
with a limited MinMod slope as the smallest producer-side repair candidate.
It is not evidence that the full M5 transaction, gas/dust energy ledger, or
z=6 science trajectory has passed: both bounded jobs intentionally stop before
M5, so no ledger drift was measured. The earlier production job 670735
remains a failed result;
its material receiver rejection is retained as a correct fail-closed guard.

No clipping, M5-side floor, receiver relaxation, refinement change, or new
profiler is justified by this comparison. The cancelled job 362821 is not a
science result: it used `interpol_var=1` with the full M5 contract but was
stopped before reaching the level-8 probe because level-7 M5 was still
running.

## Full-coupled bounded attempt

Job `363017` was rebuilt from the current source and run with the full
670735-like RT/dust/CHIMES contract, explicitly
`interpol_var=1, interpol_type=1`, and `SNRT_HYDRO_ENTRY_STOP=2`. It was
intentionally bounded to one node, 8 MPI ranks × 4 OpenMP threads, and 80 GB.
The environment selected `SNRT_BACKEND=openmp` and
`SNRT_DUST_BACKEND=openmp`; this was therefore an OpenMP correctness/entry
test, not a GPU performance test.

The run reached the RT/M5 initialization and emitted 16 dust/IR hybrid batch
reports, but produced no level-8 entry, `godunov+sync`, M5 transaction, error,
or completion marker through `00:27:38`. The sampled per-step maximum RSS
stabilized at approximately 26.1 GB; no OOM or NaN was reported. It was
cancelled at that reproducibility point because the earlier full-coupled job
`362821` showed the same plateau and was cancelled after approximately 27
minutes. This is an **inconclusive execution/performance stop**, not a
physical pass or a physics failure. The OpenMP full M5 path must be made
observable/tractable (or run on its intended GPU backend) before claiming the
M5-enabled boundary gate.

Run artifacts: `../.m5-l8-coupled-interp1-run/live-363017.log`, SHA256
`308d9b6f457a7c603686d18c8da336ccbb450ef2463e81c2a7628ea2cee03439`;
Slurm state `CANCELLED`, elapsed `00:28:51`, with no completion marker.
This attempt does not authorize a 256^3 or 512^3 launch.

The corrected CUDA-enabled retry was job `363143`. It requested and received
one A10 (`gres/gpu=1`); the rebuilt executable linked CUDA 13.0.2,
cuFFT/cuBLAS, and the SNRT device objects (executable SHA256
`aeed5827f41930fe33216641aea9a9f162e0af5fcacffc8dd2b99a5c39f58a62`). The
runtime reported all eight local MPI ranks mapped to GPU 0 and the intended
hybrid counts were `CPU=0 GPU=1` for the dust/IR batches. A direct utilization
sample reached 100% during kernel execution.

This corrected run still did not reach the level-8 probe, `godunov+sync`, M5
transaction, error, or completion marker by `00:30:28` of useful execution
(`00:31:30` Slurm elapsed including cancellation). Sampled MaxRSS stabilized
near 26.3 GB and no OOM/NaN was reported. The run was cancelled after the
same post-batch plateau persisted; its status is **inconclusive execution
performance**, not a physical pass/fail. The last log contains the CUDA
hybrid dust/IR reports with `CPU=0 GPU=1` and no transaction verdict; SHA256
`482c62fcdc44459e4e187392d9cc5a6a6ed5ae1e5364624d505625371ba1c067`.

The two attempts separate the failure modes: `363017` was a CPU-only binary
and environment path, while `363143` proved CUDA allocation and kernels but
still stalled before the coupled transaction gate. The next work item is a
targeted stage-boundary/performance diagnosis of the full material/CHIMES
operator (with no physics or refinement changes), followed by one more
bounded gate run. Broad-resolution launches remain prohibited.

## Transport-only separation result

To locate the plateau without adding a broad profiler, the temporary
`SNRT_HYDRO_ENTRY_STOP=3` probe was moved to the exact boundary immediately
after the number, energy, and IR M5 transport loops and immediately before
the per-cell material/CHIMES loop. CUDA-enabled job `363257` then completed
the probe successfully in `00:13:02` Slurm elapsed (`00:12:08` application
step, `17:34.180` total CPU seconds, Slurm MaxRSS `25.04G`). Its log contains
the marker:

`SNRT bounded M5 transport stage completed; stopping before material/CHIMES`

and no level-entry negative, `godunov+sync`, transaction error, or completion
error. The transport-only result is therefore a conditional M5 transport
pass for this boundary trajectory. It is not yet a full coupled RT/dust/
CHIMES pass because material and ledger publication were intentionally not
executed. The earlier 30-minute plateau is localized to the subsequent
material/CHIMES portion (or its coupling), not to M5 angular transport.

The temporary stop value `3` and its source hook must be removed after the
next full-coupled diagnosis; it is not a production control. Log artifact:
`../.m5-l8-coupled-interp1-run/live-363257.log`, SHA256
`88a6d90122aa21783f04a6c56a397e2d495e8b655a73c62d64404c7e5ab0b327`.

## Required next action

The next effective cosmological hydro/M5 namelist must explicitly contain
`REFINE_PARAMS: interpol_var=1, interpol_type=1`; the relative-dust comparison path remains
explicitly pinned to `interpol_var=0`. Run-directory generation already
defaults the cosmological wizard and `hr5_production` preset to the tested
`interpol_var=1, interpol_type=1` pair, while
the relative-motion generator deliberately overrides it to `0`.

Before resubmitting the 128^3 science trajectory, retain this evidence and
resolve the OpenMP full-coupled plateau, then perform one bounded M5-enabled
level-8 transition with the explicit setting and
`SNRT_HYDRO_ENTRY_STOP=2`. Record finite states, absence of a rejected
transaction, energy-ledger drift, elapsed/TotalCPU/MaxRSS and output policy.
No 256^3 or 512^3 launch is implied by this evidence.

## OpenMP material-initialization fix and H100 check

The first full-coupled retry with the new per-cell OpenMP material loop was
job `373163` on one A10. It failed at `snrt_dust_ir_init` with
`forrtl: severe (151): allocatable array is already allocated`, while the
stack showed concurrent entries from `material_cell`. This was a shared lazy
initialization race introduced by parallelizing the loop, not an OOM or GPU
capacity failure. The job used 25.1--25.6 GiB sampled RSS and was stopped by
the application error at `00:12:01`.

The source was changed so `snrt_dust_live_prepare` is called once, serially,
after the M_N payload contract is accepted and before the OpenMP material
loop. The per-cell guard remains for serial/restart callers; it is no longer
the first initializer in the M_N worker region. The CUDA+CHIMES rebuild
completed successfully with executable SHA256
`c6e37acba755db120a3db3cc75c8d340554ebc9ae0407851b263f536b30c009f`.

Job `373224` reran the same full-coupled contract on one H100 NVL with 8 MPI
ranks × 4 OpenMP threads. It reached M_N transport and the material/dust
batch loop without reproducing the allocatable-array race, material rejection,
or NaN/error diagnostic. H100 utilization was sampled at 95--100% with
4.3--4.4 GiB device memory; host RSS stayed near 24--25 GiB. However, the
same 16 material batch reports remained the last visible progress after
`00:31:36`, so the job was cancelled as a performance diagnostic rather than
left to consume the full six-hour allocation. Its accounting state is
`CANCELLED`, and its log SHA256 is
`dbd706213ffefa08ca600cabedfdcc1d3e3250c7f9f985c765f4dc1a8f543b2f`.

This is a conditional software-concurrency result: the initialization race
is not reproduced, but the full material/CHIMES stage is still not a bounded
production pass. The abundant A100/H100/H200 capacity should therefore be
used for a targeted tile-level material benchmark and one controlled
multi-GPU transport comparison, not for a larger science-resolution launch.

## Tile material/dust dispatcher result

The first tile attempt (`373244`) stopped at smoke-test code `407` before the
material section. This was a fixture mismatch, not a material or GPU failure:
the sourced production contract is `chimes_transition_d03`, whose runtime
path intentionally supplies zero direct D03 grain absorption so CHIMES owns
the absorption. The fixture nevertheless compared against the pure-D03
absorption table. Production code was not changed for this mismatch.

For the actual dispatcher measurement, the smoke executable was run in its
bounded material-only path with the supported `SNRT_SPECTRAL_MODEL=fixed`
control (the test still uses the same `nc=1031` heterogeneous material
payload and conservation/rollback checks). H200 job `373259` completed all
three independent backend cases with the CUDA+CHIMES build
(`bin/snrt_dust_backend_smoke`, SHA256
`a635ddc87d655f4c7c2ae6eae5122449bd1fdab295fc470ff840996ec6f13316`):

| backend | elapsed application time | material batch routing | result |
|---|---:|---|---|
| OpenMP | 2.937 s | CPU=5, GPU=0 | PASS |
| CUDA | 2.611 s | CPU=0, GPU=5 | PASS |
| auto/hybrid | 3.104 s | CPU=4, GPU=1 | PASS |

The numbers include the fixed smoke preamble and are not a kernel-only
benchmark; they are a contract/routing comparison. All cases reported
`DUST_BACKEND_FORTRAN_PARITY` at `9.138381e-12` and `9.095450e-12`, and
`DUST_MATERIAL_TILE_BACKEND_CONTRACT_PASS`. The auto case proves that the
5-batch, 1031-cell tile can mix CPU and GPU leases while preserving the
material/IR checks. The small auto slowdown is expected at this tile size and
does not justify stream-count tuning yet.

Artifacts are under `../.material-tile-h200-20260916/`; Slurm job `373259`
completed in `00:00:09` and used no RAMSES output. The next bounded experiment
is a four-GPU H100 transport/level-8 comparison using a fresh run directory;
it remains separate from the production trajectory.

## Multi-GPU H100 boundary check

The approved bounded comparison was submitted as job `373262` from the fresh
directory `../.m5-l8-multigpu-h100-20260916/`. Its contract was the same
CUDA+CHIMES M5 coupled boundary as the preceding H100 check, with
`SNRT_HYDRO_ENTRY_STOP=2`, 8 MPI ranks, 4 OpenMP threads per rank, and no
requested RAMSES output. The executable, namelist, and environment hashes
were respectively

`c6e37acba755db120a3db3cc75c8d340554ebc9ae0407851b263f536b30c009f`,
`d38a2c629ee2a088eeddc12c9e7bfca2b8a093cbbe179b0c993cf63385f99b68`, and
`a83c920750ec582fb64c79e54fe419f03b37fa4937c561d6270ef129a9360f9e`.

The scheduler reported `gres/gpu=4`, but the application processes received
`CUDA_VISIBLE_DEVICES=0,2,3` and the CUDA pool consequently enumerated three
devices. The fourth physical H100 was unavailable to this step because the
node also carried another user's one-GPU job (`362141_1`); the node's GRES
inventory also reports five H100s while `nvidia-smi` exposes four. This is a
cluster allocation/inventory issue, not a SNRT device-selection failure. The
application mapped its eight local ranks over the three visible devices as
0,1,2,0,1,2,0,1 (local CUDA indices, corresponding to physical 0,2,3).

The run completed initialization, passed Morton-hash and pre-RT NaN checks,
entered M5 transport and the material/CHIMES operator, and showed active GPU
work on the three visible devices. It did not reach the level-8 transaction
completion marker by the predeclared 30-minute diagnostic boundary. No
`ERROR`, material rejection, NaN, multigrid nonconvergence, or MPI abort was
observed; the last visible records were the 16 per-rank IR-absorption and
dust-scattering batch reports. It was therefore cancelled at `00:30:27` as
an **inconclusive execution/performance result**, not as a physics failure.
Slurm accounting recorded `TotalCPU=01:31:01`, MaxRSS approximately 25.16 GiB,
and no `output_*` directory. The final log SHA256 is
`0e2f36f8b3d17f72b2d8b905344ab83a86b176f2bea8b5ef14a67b6378399690`.

This check confirms that the OpenMP material-initialization race fix remains
effective under concurrent MPI/OpenMP/CUDA execution, but it does not provide
a four-GPU scaling number and does not promote the full coupled 128^3 path to
production-ready. The transport-only pass and the 1031-cell tile pass remain
the valid bounded evidence. The next implementation target is the material/
CHIMES coupling throughput path: preserve per-cell CHIMES correctness while
reducing launch/dispatch overhead and exposing a genuinely batched work unit;
only after that should a clean, uncontended multi-GPU boundary rerun be used
for scaling. No 256^3 or 512^3 launch is justified by this result.

## H200 coupled boundary check

To separate GPU generation from the coupled-material bottleneck, job `373277`
was launched from the fresh directory `../.m5-l8-multigpu-h200-20260916/`
with the same input, CUDA+CHIMES executable, `SNRT_HYDRO_ENTRY_STOP=2`, and
no requested RAMSES output. The resource request was 4 H200 GPUs, 4 MPI
ranks, and 4 OpenMP threads per rank; the 4-rank layout was chosen to avoid
MPI-rank oversubscription on the requested GPUs. The run directory hashes
matched the preceding controlled inputs:

`c6e37acba755db120a3db3cc75c8d340554ebc9ae0407851b263f536b30c009f`,
`d38a2c629ee2a088eeddc12c9e7bfca2b8a093cbbe179b0c993cf63385f99b68`, and
`a83c920750ec582fb64c79e54fe419f03b37fa4937c561d6270ef129a9360f9e`.

Although Slurm recorded `gres/gpu=4`, the CUDA pool enumerated only two H200
devices (PCIe `00d8:00` and `00b4:00`) and mapped ranks `0,1,2,3` as
`0,1,0,1`. The node simultaneously carried three other one-GPU jobs, and
its live GRES accounting was 7/8. This is the same scheduler/device-cgroup
discrepancy seen on H100; the run is not a four-GPU scaling measurement.

The run passed initialization, Morton-hash, and all pre-RT NaN checks, then
entered M5 transport and the material/CHIMES stage. GPU and CPU activity were
observed, and no error, material rejection, NaN, multigrid nonconvergence, or
MPI abort was emitted. The last visible records were the eight per-rank
IR-absorption/dust-scattering batch reports. No level-8 transaction completion
marker appeared by the predeclared boundary; the job was cancelled at
`00:29:37` as an **inconclusive execution/performance result**, not a physics
failure. Accounting recorded `TotalCPU=01:36:06`, MaxRSS `48.32G`, and no
`output_*` directory. The final log SHA256 is
`30bcbbc307f4aae659ae2ed9cc7e390a3490255be339ee38b1e24ef405583d2f`.

Together with H100 job `373262`, this shows that changing H100 to H200 and
reducing the MPI rank count does not remove the coupled-stage plateau. The
valid conclusions remain: M5 transport-only is conditionally passed; the
1031-cell material tile is passed for OpenMP/CUDA/auto routing; the full
per-cell CHIMES/material transaction is not yet a bounded production pass.
The next implementation should therefore target a batched, transaction-safe
material coupling path, while treating Slurm GRES/CUDA visibility as a
separate cluster configuration issue. No broad-resolution launch is justified.

## Local material-only transport fast path

The first scoped throughput change removes spatial work that is impossible in
the M5 material receiver. `snrt_dust_live_stage` supplies the already
transported incoming angular field through `incoming_radiation`; it sets
`material_only`, rejects nonzero neighbors/ghosts/transport callbacks, and
therefore performs no second spatial transport. Previously
`snrt_dust_ir_advance` still allocated the six-face `remote`/`blocked` payload,
validated spatial ghost metadata, checked reciprocal neighbors, and replayed
the face loop for every local cell. The new branch skips only those operations
and computes the unchanged per-cell optical response (`transmit`, `loss`, and
`response`) required by the material fixed point. Absorption/emission
iterations, thermal/material callbacks, gas exchange, closure thresholds and
transactional publication are unchanged. The non-local path is unchanged.

The modified source is
`433aa42adc73053c2e94c6b730b2a5505fae29dd4ba7748a1c52dc25c116ca17`.
The initial Intel ifx CPU smoke binary was
`285ca852f7cb683ebaac618728e60086c715e81026a42d775b7d3bb7518b7a9b`;
the subsequent CUDA-enabled H200 smoke rebuild is
`da75984e40c4c064a8d92d8a015150b06117bfdc0d837def6c8a4302f8a1f749`.
Both rebuilds used the existing Makefile/VPATH and unchanged FP64 physics.

The material-only 1031-cell smoke passed with OpenMP and auto dispatch,
including `DUST_MATERIAL_TILE_BACKEND_CONTRACT_PASS`, gas/dust conservation,
rollback, and both material-energy parity checks. The complete non-local
OpenMP smoke also passed the existing IR transport, scattering, sublimation,
stiff gas-coupling and rollback assertions. The local logs are retained under
`.material-tile-h200-20260916/`.

This is a bounded kernel-path improvement, not a full coupled speedup claim:
the standalone smoke was rebuilt and parity-tested, but the complete
M5+CHIMES level-8 run has not yet been resubmitted. The next implementation
remains a genuinely batched transaction-safe material coupling unit; this
change establishes the cheaper local operator that such a batch can call.

The CUDA-enabled H200 check (Slurm job `373298`) reran the same 1031-cell
material-only contract independently for OpenMP, CUDA and auto/hybrid routing.
All three returned zero and emitted
`DUST_MATERIAL_TILE_BACKEND_CONTRACT_PASS`; the elapsed application times were
3.029 s, 2.823 s and 2.792 s, respectively. CUDA enumerated one H200 and the
auto run mixed CPU=4/GPU=1 material batches. The exact material-energy parity
receipts remained `9.138381e-12` and `9.095450e-12`. No RAMSES output was
created.

## Batched M_N material transaction adapter

`snrt_dust_live_moment_tile` is now the batch counterpart to the existing
one-cell adapter. For a tile of `nc` cells it reconstructs all incoming M_N
moments, invokes one material-only `snrt_dust_live_stage`, and performs the
output projection/reconstruction checks before publishing any energy,
material, temperature, gas, PAH, sublimation, or phase state. Spatial
transport remains owned by the M_N operator; the adapter does not reintroduce
the legacy six-face stencil. Mutable optional fields are copied into a local
write set, preserving all-or-nothing publication on rejection.

The production module and native regrid fixture compiled with the CUDA+CHIMES+
HDF5 build contract. The bounded native run used the v3 control contract
required by this fixture (the first v4-environment attempt stopped at the
fixture's intentional contract-version gate). The successful run emitted:

`SNRT_M5_NATIVE_REGRID_PASS signed_moments=1 no_SN_allocation=1 refine_restrict=1 reject_atomic=1`

`SNRT_M5_LIVE_DUST_CELL_PASS material_only=1 no_halo=1 rollback=1 balance=-8.81621E-39`

`SNRT_M5_LIVE_DUST_TILE_PASS batch=2 one_stage=1 parity=1 rollback=1`

The tile test matched the prior one-cell result for both cells and verified
that a negative capacity in one member leaves both tile states and both
projections unchanged. The successful log SHA256 is
`21265709556c6bc84cbe6a2de2f1dbb5540db3764851ca5a7dcaeb8a0a66bfc4`, with an
empty stderr log (`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`).
The source hashes are `174eafb80b76515cf6b1828f40c57312b5ef941b44c4be5ed57f36acbc473b96`
for `snrt_dust_live.f90` and
`b803aa98e1ea4dd5e3f6f565680beb562e809bed9a022bb2583e70c953d9935a` for the
native fixture; the linked test executable is
`1fe05da2fc46cd00e5c855a8a5f59be053d1cdb8fb5b5af4710c7aa4f99e8249`.

This establishes a reusable batch transaction boundary and a native parity/
rollback gate. It does not yet replace the RAMSES `material_cell` call site:
the upstream CHIMES per-cell preparation/post-publication code still needs a
separate, bounded tile collector before production runs can claim end-to-end
batched throughput.

## L8 coupled rerun and unresolved-item disposition (2026-09-17)

The remaining code-side failures in the first boundary reproduction were
repaired and tested in a fresh run directory.  The repairs were deliberately
limited to the demonstrated path:

1. the RAMSES material transaction width was increased from 32 to 256 cells,
   reducing the 262,144-leaf reproduction from roughly 8,192 to 1,024
   material dispatches;
2. the restore path no longer reads the unallocated `packet%gas_momentum`
   component (the earlier width-256 run's genuine SIGSEGV at line 529);
3. a CUDA stream that is still leased is never reclaimed merely because
   `cudaStreamQuery` succeeds; an exhausted pool falls back to CPU in auto
   mode, preserving stream-local scratch ownership;
4. the M5 closure dispatch can replay the GPU-exported dual cache to produce
   the identical angular state without a second CPU Newton solve; and
5. the post-material angular projection is dispatched through the same
   OpenMP/CUDA lease policy, with the existing CUDA `live_project` kernel.

The CUDA+CHIMES+HDF5 build completed with bounds checking and Intel MPI.  The
final executable was `bin/ramses_m5_width256f3d`, SHA256
`94f75470abe9a038eead55f0740f7bdfdad43a4326820002d04e1f124f4ecc81`.
The source build included the new `snrt_moment_dispatch` and
`snrt_moment_live_cuda` objects through the existing Makefile/VPATH order.

The final bounded run was Slurm job `373873` in
`.m5-l8-width256f-rerun-20260917/`, using the exact coupled L8 input contract,
`n_cuda_streams=3` in the validation namelist only, `SNRT_BACKEND=auto`,
`SNRT_DUST_BACKEND=auto`, eight MPI ranks, four OpenMP threads per rank, and
one A40 allocation.  It passed initialization, the 46,425-grid Morton checks,
hydro/Poisson, pre-RT NaN checks, and the M5/material transaction entry.  The
log recorded both CPU and GPU material batches (`batch_cells=256`) and no
SIGSEGV, OOM, non-finite state, rejected transaction, MPI abort, or MG
nonconvergence.  MaxRSS was `28,803,236 KiB` and TotalCPU was `02:32:17`.

The job reached its declared 45-minute wall limit (`TIMEOUT`, `00:45:02`)
before the bounded completion marker could be emitted.  This is an
**inconclusive cost-bound result**, not a correctness failure: the full
per-cell CHIMES/material transaction is still too expensive to certify within
this short boundary allocation.  No `output_*` directory was created.  The
earlier diagnostic jobs `373789`, `373832`, and `373855` remain historical
cancelled/inconclusive runs; they are not relabeled as passes.

Consequently the active unresolved list is now split cleanly:

- **Closed engineering defects:** tile over-fragmentation, the unallocated
  restore read, CUDA stream ownership race, and duplicated CPU M5 closure/
  projection work.
- **Still open performance/science qualification:** a full L8 coupled
  transaction completion and the previously planned 128^3/256^3/512^3
  science campaign.  These require a separately budgeted production
  allocation and are not silently converted into a pass by the bounded run.
- **Out of this M5 boundary:** the already recorded physical-source gates
  (including the 40--120 M_sun and external AGB-source admission issues) remain
  fail-closed/parked according to their own provenance records; they are not
  changed by this runtime repair.

No new scientific parameter or acceptance gate was introduced.

## Binary-identity correction and latest resource check (2026-09-17)

The SHA256 shown above belongs to the executable used by the formal `373873`
45-minute run.  It must not be called the latest `width256f3d` executable.
The rebuilt current executable is `bin/ramses_m5_width256f3d`, SHA256
`2dc3985abb9f9511043d8a3bc1510aa7dc835c6b4f8697c60b76e7341905e458`.

The first 4-GPU retry (`373941`) was stopped after detecting that its copied
executable had the older `94f75470...` identity.  It is excluded from evidence.
Two fresh checks used the current SHA and the same coupled L8 input: `373942`
used four OpenMP threads per rank, and `373943` used eight OpenMP threads per
rank to consume the allocated eight CPUs.  Both reached initialization,
Morton checks, time integration, M5/material hybrid calls, and showed no
SIGSEGV, OOM, non-finite state, rejected transaction, or MPI abort.  They were
operator-stopped before completion to avoid repeating the already established
long CHIMES cost plateau; neither is a completion claim.  No `output_*`
directory was produced.

The focused hydro-restriction wiring test and the current CUDA projection
symbol check pass.  Therefore the recheck found no new code-side defect.  The
only unresolved item in this boundary remains the separately budgeted full
L8 completion/performance qualification, not a failed correctness gate.
