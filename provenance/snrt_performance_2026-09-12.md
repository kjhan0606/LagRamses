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
