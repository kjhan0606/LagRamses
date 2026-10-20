# M5 radiation-momentum coupling — 2026-10-11

## Scope and implementation

The M5 material interaction now carries the receiver-resolved signed angular
first moment of absorbed UV energy from the same ray/node absorption used by
the scalar photon and energy ledgers. H/He contributions use the accepted
atom-budget cap; the dust receiver uses its uncapped accepted absorption. The
RAMSES update combines UV and IR first-moment changes, divides by physical
`c` (not reduced `c`), updates mixture momentum, and debits midpoint kinetic
work from the gas thermal reservoir. In the optional relative-dust path, UV
impulse is partitioned by phase opacity and phase midpoint work is removed
from dust heating; the IR receiver stages phase momentum/work transactionally.

This bundle qualifies the configured **static/co-moving dust + M5 gas/AGN
profile**. The sink/AGN dust profile explicitly rejects
`dust_relative_motion=.true.`; that separate state-layout combination is not
silently claimed here.

## Evidence

- Native receiver smoke: `SNRT_MOMENT_BAND_MOMENTUM_OK`, using two unequal,
  non-collinear rays and checking all four receiver first moments.
- Native contract command:
  `SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1 bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh`
  passed through `SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK`.
- Amended integrated H100 job `419549` completed with Slurm exit code 0 in
  41m27s; the strengthened runner emitted
  `M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS`. No `output_NNNNN` raw-output
  directory was created. The run used M5 order 5, live dust mass, CPU CHIMES,
  and a nonzero AGN luminosity (`3.0561937668561012e40 erg/s` at coarse steps
  1 and 2). Both source-active ledger rows had positive source, absorption,
  and dust heating; photon-ledger closure residuals were
  `1.6847628294300e-18` and `-4.9402324234015e-18`. The corresponding
  energy-ledger residuals were `7.5234439923036e-17` and
  `2.3066940033601e-16`.
  The actual stored-state `Delta P = J_rad` assertion was exercised for
  16,896 successful material-cell calls (tile spans summed over levels and
  substeps, not unique cells) without a closure failure. The acceptance
  script itself reported this count. Effective namelist SHA256:
  `6c8e0a0be400110576835fdc218c0124f8394f7670b93629590d56db4c5772fa`.
- The preliminary H100 job `419495` also passed the earlier acceptance
  version in 41m43s; it is retained as repeatability evidence, not used in
  place of the amended-run verification.
- Backup audit found three early finite-value returns after a successful dust
  receiver could retain `status=0`. They now emit a specific diagnostic and
  return `status=1`. The acceptance runner additionally requires positive
  source, absorption, and dust heating, rejects those diagnostics, and derives
  successful momentum-assertion call coverage from committed tile spans. The
  reviewer confirmed status propagation to the enclosing M_N level transaction
  and corrected its initial ledger-column misread. This amended code/runner is
  rebuilt and exercised by job `419549` below.
- Amended isolated build job `419546` passed; binary SHA256
  `d1144a3f680b485df11d5bb292621984ffe5d2141627155662fa1f8fae1166fd`.
  Integrated H100 rerun: job `419549`, Slurm `COMPLETED 0:0`, 41m27s. Result:
  `.m5-source-dust-chimes-419549/result.txt`.

## Audit disposition and limits

Claude Opus 5 could not audit in this environment: its CLI rejected model name
`opus5` and reported no login. The authorized GPT-5.6-sol backup audit first
returned conditional pending additional source context, then revised to
**Pass with verification caveat** after inspecting the complete UV/IR phase
work and runtime-admission wiring. No algorithmic, wiring, or dimensional
blocker was found for the scoped co-moving profile.

The assertion's units and tile-span coverage were independently reviewed;
the revised finite-error status path is structurally fail-closed through the
M_N level transaction. The amended integrated rerun passed its positive
source/absorption, tile-call-coverage, finite/error-marker and no-raw-output
checks, closing this bounded bundle for the admitted co-moving profile.
Relative-dust phase-state checks are present in code, but the sink/AGN profile
rejects that layout, so it is not qualified here. This reference-control
integration is not broad production qualification, parameter calibration, or
convergence evidence.

## Selected effective-population connection — in progress, 2026-10-12

Job `419903` reuses `m5_source_active_dust_chimes_reference.sbatch` with
`M5_SOURCE_PROFILE=effective_population`. Its initial Slurm state was
`PENDING (Resources)`, not a completed qualification. Requested resources:
one A100/H100 GPU, MPI2/OMP2, 64 GiB, one-hour limit. Existing cosmological
Grackle calibration job `419850` is untouched; it is RT/dust-off science,
not evidence for this high-level connection.

This is the already selected effective SSP2/Kroupa mixed-v5 and empirical
Maoz/N100 Ia profile documented in
`effective_population_execution_2026-09-11.md`, not a new physical-model
choice. The existing source wrapper and runtime Ia approval are retained.
The wrapper explicitly declares truncated mixed evolution/late-SED omission
and `reference_control`; this execution does not promote it to a fully
approved physical source release or change the global default Chabrier IMF.

The evaluated prior short input is reused from
`/gpfs/kjhan/LRD_JWST/.effective-coupled.Q9Kw29/short/run.nml`
(SHA256 `1cc6534578fcaffa84c02d0dbc037973d32691aa059efc738bc9249dfa2cd358`).
Only `rt=.false.` is replaced by the canonical M5/order5 selector, and
periodic `foutput=2` is replaced by `1000000` for this no-dump execution.
The old environment's deprecated RT-enable/level restriction is not inherited.
All-level transport is selected; matched v5 SED, yields/history and the current
`fp2_snia_effective_population_runtime_v1.nml` are bound explicitly.
Actual existing BH seed is reused, not a synthetic photon injection. Model
settings remain: 64 owned leaves, bounds level2--4, four coarse steps,
courant3e-6, stellar birth, C/silicate coadvected two-size grains with
growth/sputtering/coagulation/shattering and condensation, CPU CHIMES157,
no CR/MHD/Fe/PAH/drift/shock carriers. GPU hydro remains explicitly disabled.

The pinned M5 executable is the `419546` binary above. It has NVAR187,
NENER0; the old effective-profile binary had NVAR199 only to reserve twelve
optional relative-dust momentum fields. The shared reader places CHIMES157
after the same element/dust window; those additional fields are inactive for
this co-moving profile. No chemical field is dropped or remapped here.

Effective input destination:
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-effective-dust-chimes-419903/effective.nml`.
The runner records actual input/binary/library/data hashes. Output policy is
one future dummy clock (`aout=2`, `tout=1e30`), foutput/fbackup1000000,
four-step evolution: expected raw dump count/size/total are zero. GPFS free
space at submission was approximately254TiB. Unexpected raw dumps or runtime
errors are retained for inspection and cause failure; no active checkpoint or
science output is deleted.

Acceptance retains live M5 source/absorption/dust heating, stored-state
momentum assertions, finite states and no-error/no-dump checks, and adds actual
stellar births, population ID2, matched common-IMF source binding and admitted
Ia accounting. As in the previous ~0.2Myr integration, a due Ia event is not
expected; nonzero event/debit evidence is the prior native/live population
test, not fabricated in this young-source run. Aggregate radiation source
positivity alone does not distinguish stellar photons from AGN photons;
that attribution must be checked before claiming the missing stellar-to-M5
connection verified. No pass is claimed at submission.

### Live launcher correction — 00:52 KST

419903 started on A100 node syn101. Its first interval was unexpectedly
slow; read-only process snapshots showed MPI progress/halo exchange on one
rank and M5 angular slope/transport on the other, rather than a CHIMES
integration stall. A second snapshot showed movement into a later transport
stage, so silence was not treated as a dead job.

The actual thread affinity identified a launcher error: both rank main
threads were pinned to CPU4 and both OpenMP workers to CPU5, although Slurm
allocated CPU4--5 and CPU10--11. `--cpu-bind=none` plus
`OMP_PROC_BIND=close,OMP_PLACES=cores,I_MPI_PIN=0` let both ranks choose the
same initial OpenMP places. Allocated CPUs10--11 were unused by those workers.
The reusable launcher now requests `srun --cpu-bind=cores`, giving each rank
its own CPU domain before OpenMP constructs core places. Bash syntax and
diff checks pass. No physical input, tolerance or transport formula changes.

Without restarting or cancelling the active calculation, the driver moved
the validated rank2 main PID4140587 from CPU4 to CPU10 and its OpenMP worker
TID4141077 from CPU5 to CPU11 inside the same allocation. `taskset` confirmed
both changes at00:52:26. Rank1 remains on CPU4--5. This live correction must
be reported with the eventual runtime; the original job's spooled script
still contains the old unbound launch. Later affinity/forward progress must
be rechecked because runtime binding can be reapplied on team transitions.
This is not yet a measured speedup or physics PASS. Read-only debugger helper
steps exited; both simulation processes were detached and left running.

### Selected effective-population result — bounded PASS, 00:53 KST

419903 terminated `COMPLETED 0:0` in12m56s; the actual simulation step419903.0
terminated `COMPLETED 0:0` in12m50s, MaxRSS825796KiB. `Run completed`, coarse
step4 and the runner's profile-specific acceptance marker are present.
The auxiliary debugger probe419903.2 failed127 because `/usr/bin/gdb` is not
installed on syn101; this is not a simulation failure. Later Intel debugger
helpers completed and detached. A final affinity observation could not start
because the simulation had already completed; it was not restarted.

Effective namelist SHA256:
`82a664e4f80a44a7f31442a26425e767a77ee8a259c31ee29cd219d85aa68919`.
Diff against the prior evaluated profile confirms exactly the two declared
selector/output changes. The matched SED, history, Ia contract, binary,
CHIMES library and spectral/dust data hashes are retained in
`.m5-effective-dust-chimes-419903/inputs.sha256`.

The admitted stellar population is ID2 with `common_imf_mixed_evolution`,
actual stellar births64 per step (256 by step4), and empirical Ia accounting.
M5 order5 commits four intervals; the initial interval is dark and the later
three have positive source, absorption and dust heating. Actual AGN
luminosity is positive, maximum2.5728170830984776e34erg/s. The stored-state
radiation momentum assertion covers128 successful material-cell calls in
the source-active tile-span counter, not128 unique cells; initial dark
commits have their separate diagnostic coverage. No runtime error,
nonconvergence, rejected transaction or nonzero NaN counter was found.
The runner confirms raw_outputs0: no raw dump was produced or removed.

**Stellar-to-M5 source attribution is independently discriminating.**
The native M5 ledger integrates photons in code volume and density units.
Multiplying the summed source by `units_length**3 * scale_nH`, with
scale_nH=`0.76 * 1.66e-24 / 1.6605390e-24`, gives
4.471853742514118e67photons. The existing AGN JSON records cover steps0--4.
Even converting *all* logged Bondi-supply rest mass to photons at the lowest
group edge0.01eV gives at most7.823335301941491e61AGN photons. This purposely
ignores accretion caps, radiative inefficiency, partition losses and higher
photon energies; it is an upper bound, not an AGN emission estimate. The
initial existing BH has zero pending fuel, with no new/merged sinks or
diagnostic injection. The only source producers in the inspected M5 path
are actual AGN receipts and stellar intervals. Total injected photons exceed
that conservative AGN bound by571604.51; hence the stellar contribution is
at least4.471845919178816e67photons and cannot be an AGN-only false pass.
No extra diagnostic framework, perturbed physics or comparison rerun was
needed for this attribution.

Source-active photon closure residuals are8.7458799751319e-25,
2.5539211509994e-22 and-3.4553773722912e-23; energy closure residuals are
5.7276178772706e-22,-6.2110520348292e-22 and6.4030572955502e-21 in their
native ledger units. They satisfy the existing scale-relative checks.
Gas/dust/IR states remain finite and nonnegative; xHII grows from9.988e-5
in the dark interval to3.582e-4. Final printed global mcons=-8.95e-7 and
econs=5.11e-4 are **not** claims of exact global conservation in a cooling,
accreting, star-forming run. All four printed mcons/econs rows match the
previous selected S_N-profile run at its printed precision. The explicit
source/receiver ledger and stored-momentum checks are the narrower
conservation evidence; no global tolerance was relaxed.

Timing is not a clean scaling benchmark: unbound launches shared CPU cores,
and live affinity was changed during the trajectory. M5 commits took9.372s
(dark),310.581s,300.373s and31.999s; only the final interval was after the
affinity correction. Closure dispatch used both CPU and CUDA; this tiny run's
face work stayed on CPU under auto policy. The marked last-interval reduction
is consistent with correcting starvation but is not an independently
controlled speedup measurement. Future launches use disjoint Slurm CPU
domains from startup.

This closes the bounded missing **selected effective stellar source + M5 +
AGN + co-moving dust + CPU CHIMES** execution connection. Its scope remains
the declared noncosmological, young-source effective model. The separate
nonzero-Ia event/debit test remains the Ia evidence. This no-dump run is not
a new M5 restart proof, arbitrary high-Z/late-age source coverage, physical
release approval, cosmological dust qualification or galaxy calibration.

### Selected M5 restart continuation — submitted, not yet evaluated

Job419934 uses the same pinned binary/source/material profile on A100 syn101,
MPI2/OMP2,64GiB,one-hour limit, with corrected disjoint Slurm core binding.
`M5_VERIFY_RESTART=1` extends the same existing runner: continuous steps0--4
write checkpoints at2 and4, then a second process restores the actual step2
checkpoint and advances to4. The checkpoint is linked, not duplicated.
No scientific input/physics parameter changes; only nrestart/output controls
differ. Output policy: noutput1,aout2,tout1e30,foutput2,fbackup1000000.
Estimated64MiB/dump, three unique dumps<=192MiB; GPFS free space254TiB.
Raw data are retained until terminal success and driver evaluation.

Actual inputs below `.m5-effective-restart-419934/`:
`continuous/effective.nml` SHA256
`f3c7e4bb8f51f60895e3136be12911aa88678aa80417de9795916c1d19f859a4`;
`resumed/effective.nml` SHA256
`b367b29be16a95cf79e27501795b054634cf1f40ac8ce18a5123abc297f48e18`.
The existing `p0_hdf5_stellar_restart.py` gains an optional
`--compare-m5-coupled LEFT RIGHT --source-checkpoint MID` path that bypasses
its synthetic/reader campaign and compares the actual integrated files.
It requires a populated M5/order5/format1000+ checkpoint, actual ordinary
stellar mass return, all hydro/chemistry/dust/radiation/IR/gravity/particle/sink
datasets and source/contract attributes. Signed higher SH moments are allowed,
but initial photon monopoles must be positive. The prior effective-profile
comparison tolerances remain rtol1e-12/atol1e-25; integer fields and identity
attributes match exactly, clocks retain rtol1e-12/atol0. This does not qualify
changed rank count, remapping/refinement, all scientific source ages, or new
physical source approvals.

The comparator unit check admits matching signed-moment data and rejects
changed source identity or changed moment state. The existing default P0
synthetic/legacy-fail-closed mode still passes. These validate the checker,
not actual restart success. The live job was RUNNING at01:09 KST; no result
or pass is inferred yet. The current launcher is an existing local untracked
file, not part of the recorded HEAD. No commit/push performed.

#### Affinity follow-up and restart resubmission — 01:14 KST

Direct thread inspection of419934 disproved the assumption that generic
`srun --cpu-bind=cores` alone preserved distinct MPI/OpenMP domains here:
both rank main threads used CPU26 and both workers CPU27, although the
allocation also included CPU60--61. The inspected stacks/progress show a
live computation, not a dead job; however this is an incorrect launcher
placement and wastes the assigned resources. The earlier one-line launch
fix and its documentation are superseded by the following verified recipe.

The driver reused the exact CPU-slice/explicit-OMP-places binding already
used by Grackle replay419714 and calibration419850, now retained as
`simulation/snrt/runs/slurm_rank_bind.sh`. Under unbound srun, the wrapper
reads the actual allocated CPU set, takes the local-rank disjoint slice,
sets explicit `OMP_PLACES`, unsets inherited KMP affinity, and execs taskset.
The M5 runner copies and hashes this wrapper with its inputs. Local wrapper
tests produce distinct rank0/rank1 places and shell/diff checks pass; actual
compute-thread affinity remains to be observed in the replacement job.

419934 was deliberately cancelled after6m52s, while its first step2 HDF5
output had been started. Its inputs/logs/partial-or-complete raw checkpoint
are preserved pending evaluation/cleanup, not used as restart qualification
or relabelled a physics failure. Replacement job419937 has identical physics,
output schedule, binary and MPI2/OMP2/resources, but the corrected explicit
rank wrapper from startup. The previous job cannot be called the final
restart result; no comparison PASS is claimed yet.

#### Actual restart failure and migration repair — 01:28 KST

419937's continuous and resumed MPI simulation steps both completed0:0
(1m40s and1m14s); the overall batch failed1:0 in3m00s at the comparison.
Both runs reached coarse step4 with finite states and the native restoration
marker. This is NOT a solver crash or OOM. Explicit rank places in the logs
are CPU26--27 and60--61. Direct live worker affinity was not sampled before
termination, so the wrapper logs alone are not a thread-affinity proof.

The first reported `gravity/level_2/f_1` discrepancy was a checker row-order
error: the distributed Hilbert reader reorders three octs. After aligning
exact level/xg coordinates, all eight gravity datasets are bitwise equal.
The checker now requires a bijection of exact grid positions and matching
son/owner topology before comparing the corresponding per-cell payloads;
its predeclared numerical tolerances remain unchanged. Signed moments remain
allowed; row permutation is not mistaken for a physical force difference.

There is also a real coupled-state failure, not explained by that ordering:
83 coarse and84 fine hydro datasets exceed the fixed tolerance, both primary/
IR checkpoint level payloads differ, and the first restarted photon ledger
implies4.208173823562504e-8 fewer initial CODE photons than the continuous
trajectory. The stellar source totals themselves agree at printed precision.
Inspection confirms that restart forces a3/5 to4/4 oct-ownership remap;
`load_balance.kjhan.f90` transfers uold/phi/f but had no transfer of SNRT's
separate primary/IR store. Its defragmentation hook only renumbers local slots.
That missing ownership-transfer connection must be repaired, not waived by
loosening a tolerance or classifying the restart as successful.

`snrt_regrid_migrate` now uses the expanded OLD emission/reception maps and
cpu_map2 to copy only octs whose ownership is changing, before the linked lists
and old domain maps are retired. It carries the existing typed primary,
chemistry and IR payload without an FP32 round trip; receives are validated
collectively before publication, and errors abort through the existing regrid
guard. Payload scratch is peer-batched, not ngridmax*full-angular storage;
the target budget is64MiB with at most8 octs per peer/batch (the minimum one
oct per active peer sets the floor). Non-RT runs do not enter this branch.
No domain-reader rewrite, new physics, source data change, or expanded audit
campaign is introduced.

Isolated incremental build419944 completed0:0 in49s from the preserved419546
build tree; NVAR187 and all original physics options/VPATH remain unchanged.
New executable SHA256:
`20300e7b70173d9a36bea375eaa78407c50ede16d2d2a51cf87062a18ffe166d`.
Same-profile MPI2/OMP2 restart rerun419946 is RUNNING on A100 syn101 at01:28KST;
actual input branches are below `.m5-effective-restart-419946/`.
It uses the same three-dump schedule/<=192MiB estimate; failed419937 raw data
are retained until the diagnosis and corrected-run evaluation finish.
Default P0 synthetic/legacy-fail-closed mode still passes. Neither build
success nor submission is coupled-restart qualification.

#### Migration result and numerical residual isolation — 01:40 KST

419946 continuous/resumed steps completed0:0 in1m42s/1m12s; the batch
failed1:0 in2m58s on the unchanged strict numerical comparison, not a solver
failure. Actual live `/proc/.../task/.../status` inspection during the resumed
step confirms distinct MPI main/OMP worker placement26/27 and60/61. Both
CUDA helper sets also remain within their corresponding two-CPU domains.
The diagnostic helper419946.2 completed0:0 and left no process running.

The missing-photon inventory discrepancy is repaired: the resumed total
photon and absorption receipts now match the continuous trajectory at the
printed13-digit precision, including9.1494096591907e-7 photons after the first
resumed interval and1.7475173505024e-6 after the second. All eight gravity
datasets remain bitwise equal at corresponding coordinates. The remaining
fine-cell hydrogen-fraction difference is2.154720087250983e-13 absolute,
maximum pointwise relative6.1542e-10;59 hydro fields still exceed the strict
comparison tolerance. This is a smaller, distinct residual, NOT continued
loss of one oct's radiation. The strict restart gate remains FAIL.

Particle comparison also needs physical-event alignment: checkpoint stars'
IDs are preserved;48 subsequently born stars have different sequential labels
because grid/particle traversal was reordered. An exact bijection of birth
time and position matches all256 stars; initial/current mass, metallicity,
release time and progress are bitwise equal in that alignment. Label-set
inventory is unchanged. The comparator preserves all pre-existing/nonstellar
IDs, permits only that declared bijection of new labels, and records its
count separately from bitwise datasets. It compares all numerical particle
fields at the same unchanged tolerances; new-star velocity differences up to
9.93e-22 still FAIL rather than being waived. Sink clouds retain exact
ID/level/position matching.

Read-only internal diagnosis found no omitted persistent closure history:
M5/IR packing is FP64, warm closure data are reset each advance; CHIMES rates/
CVODE state are freshly initialized per cell and the thermal cache is reset.
No new cache/physics patch is justified by those findings. An otherwise
identical OpenMP-only diagnostic419969 also completed both simulations0:0
but fails the same strict hydro140 comparison; CUDA alone does not explain
the residual. This control does NOT replace the default hybrid qualification.

One concrete remaining perturbation is source accumulation through reordered
particle linked lists. M5 now sorts each grid's stellar source indices by
actual source-input tuple `(tpp,xp[1:3],mp0,zp)` before ordinary FP64 addition.
The bounded heap sort is NlogN, stores only one grid's particle indices, does
not reorder the particle tree, and preserves all source rates/intervals.
Equal-key sources have identical contributions; new star IDs are not a
physical tie-breaker. This is numerical reproducibility, not a new model.
Migration communication also now checks MPI completion collectively before
validating a received buffer; the clean Makefile dependency includes the new
efficiency-module use without altering VPATH. New isolated build419974 and
default-auto same-profile restart419975 have been submitted; no PASS is yet
claimed for them. Backend choices are optional runner diagnostic environment
flags, not a namelist/schema/default change. All failed-case raw files remain
needed for the ongoing evaluation and are not deleted yet.

#### Source-order hypothesis rejected; bounded disposition — 01:48 KST

Build419974 completed0:0 in30s, executable SHA256
`674ecf7ddca040fbad7e133b2892a96c747a17532d50ee0fbff14558ddc217a6`.
Default-auto test419975 again completed both simulations0:0 (1m41s/1m12s)
but its batch failed1:0 in2m58s at the strict hydro140 comparison.
The canonical hydro failure counts and fine HII difference are unchanged
from419946, including2.154720087250983e-13. Deterministic source ordering
does not reduce this residual in the actual case. That experimental sorting
patch was therefore reverted; it is not accumulated as an unproven fix.
The real migration repair and MPI-completion check remain. A unit test of
coordinate/particle-event alignment passes and explicitly rejects changing
a pre-existing stellar ID; default P0 schema/legacy checks remain PASS.

The most useful current evidence is419946's actual same-profile continuous,
midpoint and resumed outputs. These remain retained for the unresolved
strict-equivalence diagnosis.419934 (cancelled incomplete dump),419937
(pre-repair loss),419969 (CPU-only control), and419975 (ineffective source-sort
control) have finished their bounded evaluations. Their exact HDF5/cooling
files are scheduled for removal in the one-off419981 closeout job, retaining
input/build identities, logs, text metadata and a per-file checksum/size/
deletion manifest. No science archive, calibration output,419946 checkpoint,
external source data or whole run directory is included.

Measured HDF5 sizes are101272496 bytes at step2 and101286832 bytes at step4
(~96.6MiB each), rather than the preliminary64MiB estimate. The reusable
runner's prospective three-dump estimate is corrected to<=312MiB; its dummy
future clock and step2/4 schedule are unchanged. This is an estimate correction,
not a numerical-acceptance change.

419981 also rebuilds the retained current code with the current ROOT Makefile
in an isolated copy, after removing the ineffective source-sort patch.
VPATH/compiler/physics options are unchanged; Grackle trace remains inactive.
No new integrated PASS, current-source release approval or commit/push is
claimed from this build. The remaining M5 restart issue is the strict small
numerical residual, not missing primary/IR ownership transfer, wrong CPU
placement, lost source mass or failure to evolve.

419981 terminated `COMPLETED0:0` in22s. Current retained-source build SHA256
is`bd31e7e25467575707a9f76eeea50557065b79a85879c5eeca792b56d39646cc`;
its matching object/build record is in`.m5-ramses-momentum-build-419981/`.
The cleanup completed for all20 validated files,911641472bytes (~0.849GiB).
Per-file original absolute paths, hashes, sizes and removed=true are retained
in`evaluated-raw-cleanup.json` beneath that build root. These raw data are
deleted (regeneration requires rerunning the recorded inputs), not archived
or recoverable from the metadata.419946's three raw files remain intact for
the unresolved comparison;419850 remains live and untouched on H100 syn08
at01:52KST. The preserved source/input roles and strict FAIL are unchanged.

#### First-divergence localization — 02:14 KST, 2026-10-12

The existing 64-cell MPI2/OMP2 profile is unchanged. A temporary environment-
enabled trace records only level2, coarse-step counter2 (the first interval
after the midpoint checkpoint), with all187 gas fields, three ion fractions,
9x36 primary number/energy moments and136x36 IR moments. Default execution
does not enable it; no new dump clock, physics parameter or acceptance gate
was added. Coordinates use the actual blocked-AMR `IGRID_OF/ICHILD_OF`
macros. The first diagnostic419993 accidentally used the legacy arithmetic
cell-index formula in its logger; its coordinate records are invalid and are
not used as position-matched evidence. That logger error did not alter the
simulated state. Corrected trace build420001 completed0:0 in26s, SHA256
`b4b84784d2909590e3f14fa8c6e807a92109221674e091f646018c0ce4ff817a`.
Both simulations in420002 completed0:0; its unchanged strict comparison still
failed on hydro/level_1/uold_140 (batch3m05s). Actual first-interval trace
contains256 unique records,64 cells at each of four stages, aligned by exact
oct coordinates and child identity:

- Input: clocks, all gas/chemical fields, fractions, primary and IR state
  match exactly at corresponding cells.
- After source injection: source number/energy and every other field still
  match exactly. This excludes an initial ownership/state loss or reordered
  source addition as the first divergence in this case.
- After transport: primary number and energy first differ by up to
  `6.617444900424222e-24` absolute;163/197 individual primary coefficients
  exceed the unchanged1e-12/1e-25 comparison. IR differs by up to
  `5.877830486518305e-39`, with no strict IR mismatch. Gas/fractions remain
  exactly equal because transport has not yet changed them.
- After matter coupling: max ion-fraction difference is
  `1.0797732066109011e-13`;28 individual fractions exceed the same criterion.
  Gas field29 differs by up to`8.149955142111825e-17`. This establishes
  transport as the first divergent operator and matter evolution as the
  subsequent amplification, not a missing saved CVODE history.

A bounded transport correction is under test: integrate each face in fixed
quadrature order; retain six physical face channels through reverse MPI SUM;
combine these channels in -x,+x,-y,+y,-z,+z order at the recipient. The existing
shared-face owner still computes the flux once with equal/opposite increments.
Fine/coarse contributions still SUM into their corresponding coarse face.
The128-column donor workspace is reused; no new mesh-capacity-sized array or
persistent angular state is introduced. Per-face moment integrals do add
active-face scratch, and reverse halo traffic is six channels instead of one;
these costs are explicit, not a claimed free performance improvement.
Build420007 completed0:0 in24s, SHA256
`404b2fdd85c82adbe293a18497bd62588c9a829164babe2bd55dab8dbb3a5b9c`.
Same-profile default-auto restart420008 is submitted with the temporary trace;
no PASS is claimed before evaluating it. It does not alter physical inputs,
source approvals, the numerical tolerance or the S_N path. This static-mesh
case cannot by itself prove decomposition-independent coarse/fine rounding.

420008 completed0:0 in3m06s; continuous/resumed steps completed0:0 in1m45s/
1m16s. The unchanged exact-coordinate/event comparator passes449 datasets:
427 bitwise,21 floating differences inside rtol1e-12/atol1e-25, plus the
separately counted new-star ID relabelling.39 contract attributes match
exactly, as do time/aexp/dtnew/dtold/nstep/nstep_coarse. All actual stellar
mass/progress and chemical abundance fields remain exact after alignment.
The largest reported primary/IR packed-state difference is1.410593220986745e-37;
particle velocity differences are at most5.048709793414476e-29.

The temporary stage trace also verifies the causal change: all primary and IR
moments now match exactly immediately after transport, and all three ion
fractions match exactly after matter coupling. Only5 gas momentum entries
differ by up to1.232595164407831e-32, and IR after material differs by up to
4.70197740328915e-38; none exceeds the fixed criterion. This is actual
same-profile default-auto evidence, not replacement by the CPU-only control.
The original and corrected traced jobs both take about3m05–06s, but this tiny
case does not measure large-MPI communication cost or strong scaling.

The temporary trace code and its AMR-index import have been removed from the
retained source; only the real conservative face-summation change remains.
Trace-free current build420011 completed0:0 in22s, executable SHA256
`1a4a3240482d3d827e85add2d68c3c439c9418cf1a0af1b63027d54b9a15fd65`.
Same-profile trace-free restart420012 is running; its result must be inspected
before closing the current-source numerical issue. No physical asset role,
release approval, clean commit or push is manufactured from these results.

#### Trace-free current-source closure — 02:25 KST, 2026-10-12

420012 completed0:0 in3m05s; its continuous and resumed simulations completed
0:0 in1m45s and1m16s. Current binary is the trace-free420011 executable above.
The final comparator again passes449 datasets/427 bitwise/39 exact contract
attributes with the original numerical tolerances and exact clocks. Both
branches contain256 real stars and identical total ordinary stellar return
`2.1016243221234917e-12` code mass. Native error/nonfinite/nonconvergence guards
remain clear. Photon accounting residual is at most`1.6238579302762e-22`;
energy residual is at most`3.7192011356687e-21` in the printed native units.
Each illuminated interval logs rank1 closureCPU/GPU=15232/1920 and
faceCPU/GPU=215040/0: the test really includes both closure arms in default
auto mode, not an accidentally CPU-only replacement. No forced-CUDA face
performance or large-MPI scaling claim follows from these counts.

The strict numerical restart issue for this declared same-MPI/static64-cell
selected effective model is CLOSED. Core clock/cosmology/CPU-tree audits,
new source sorting, modified physics and loosened tolerances were not needed.
The original S_N path is untouched. Generic coarse/fine conservation uses
the same shared-face and receiver-SUM rules; all-rank/coarse-interface bitwise
invariance and the six-channel communication cost remain outside the claim.

After inspecting the final comparison, source-return/nonfinite/ledger evidence,
cleanup420015 is submitted for42 exact evaluated HDF5/cooling/temporary-trace
files from419946,419993,420002,420008 and420012. Inputs, logs, binary/build
identities, compact comparison JSON, raw hashes/sizes and deletion manifest
remain. No active checkpoint from these five completed cases is still needed.
The long calibration419850 and its preserved step380/science outputs are
excluded; it is still RUNNING on H100 syn08 at02:25KST, step400/a0.08613
(z~10.61). Physical reference-control/source-use roles and active-checkout
release approval are not resolved by this successful numerical closeout.

Cleanup420015 completed0:0 in2s. All42 validated files are removed,
`1749234116`bytes (~1.63GiB); their recorded paths, sizes, checksums and
removed=true are in420011's`evaluated-raw-cleanup-420015.json`. Exact absence
of all42 files and retention of420008/420012 comparison JSON are verified.
These raw outputs are deleted, not archived; reproducing them requires the
retained input/build recipes. The preserved calibration checkpoint and live
run inputs remain intact. No simulation process or scientific archive was
cancelled, moved or deleted. No commit/push was performed in this closeout.

#### Mixed-level consumer regression — 2026-10-12

Job420024 completed0:0 in2m57s on A100 syn101, MPI2/OMP2, using the same
trace-free420011 binary and selected stellar/AGN/dust/CHIMES contracts.
The existing runner's optional `M5_AMR_INTERFACE_CHECK=1` uses existing
refinement keys: eight central level2 cells refine to level3, leaving56
coarse and64 fine leaves. Both levels retain eight octs throughout the
reported trajectory; the driver inspected that topology and successful M5
commits on levels2 and3, not only a process exit status.

The four-step run includes real star formation (417 final stars), positive
radiation injection/absorption and dust heating, dark and illuminated
material transactions, and the native cell momentum assertions. Runtime
error/nonfinite/nonconvergence guards are clear. The maximum printed photon
accounting residual is9.1356580416678e-23 in native units. This is bounded
mixed-level execution evidence for the changed transport consumer, not an
independent analytic coarse/fine flux-error measurement, mixed-AMR restart
equivalence, cosmological qualification or all-rank invariance.

Effective input SHA256:
`dbd00a880b477d0bdd4b905c68dae69582d81774cae317d0f100ecd5c1999dc9`.
Input-inventory SHA256:
`df20743aee4cc8fe1385df964e67609cabb1ef178f932f79fe7520e23c13909d`.
Result SHA256:
`74288455ae230d133b0089066c8c84733eebd6a07f370e4530ad8c163dba0241`.
These are retained beneath `.m5-effective-amr-420024/`; the future dummy
output clock and million-step output/backup intervals produced no raw dumps,
so no deletion is needed. Inputs, logs, AGN/sink history and result remain.

The actual-root `manifests/native_m5_release_candidate_2026-10-12.json` now
binds this evidence together with the static restart result. It is explicitly
not an approved launch-admission manifest. Selected reference-control/source
roles, exact release provenance/license records and committed source identity
remain unresolved; old P4/thermal-atlas/JAX inputs unselected by this native
profile are not additional prerequisites. The Grackle cosmological science
run remains separate and untouched.

Selected-data provenance inspection also finds an existing explicit pending
license record for the KL16 research mirror inside the unchanged mixed-v5
manifest (SHA2568117b5000e0d1e3e2f509f2d33a2ffa22d028aff8c5e09981af0165c0c422a74).
It is recorded as pending, not replaced by a license inferred from public
availability. Local CHIMES source includes LGPL headers and COPYING.LESSER;
that software notice is not evidence of the external atomic/dust/yield data
terms. No redistribution, upstream contact or raw-input modification occurs.

The operator has been asked about candidate use of the exact stellar SED,
AGN partition and dust approximations with their limits retained. This is
not a repeat request for the already approved galaxy-effective-model scope;
nor would such approval confer an upstream license or universal physics
qualification. No response or approval ID is manufactured while waiting.

#### Selected upstream notices checked — 2026-10-12

Read-only checks against the official Draine distribution give SHA256
`b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c`,
exactly matching the selected dust contract's original optical source.
The [official model description](https://www.astro.princeton.edu/~draine/dust/dustmix.html)
and [original table](https://www.astro.princeton.edu/~draine/dust/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all)
establish its provenance; neither inspected document provides an explicit
data-license statement. That is unresolved metadata, not a declaration that
scientific calculation is prohibited.

The PARSEC v2 VMS database still fails verified TLS access. No insecure
fallback or archive download was made. The retained author readme matches
its recorded SHA256`d4dcbc8db36e6d71780a16d9eb21a260d9ec5db3fbc23e8bcab67be247525e31`.
Its MIT notice belongs to embedded Microsoft formatting CSS, not stellar
tracks. The [Costa et al. paper](https://arxiv.org/html/2501.12917v2) has a
CC BY4.0 article notice and identifies the public database; this inspection
does not establish that notice's scope over the separate track archives.

The GitHub root-tree API for the pinned COLIBRE mirror revision
acecc180347d9505780a7985ec2c4a8a11f5c61f returns an untruncated listing with
no root license file; its pinned README provides usage instructions but no
license statement. This agrees with the existing KL16 manifest's pending
original-terms status rather than resolving it by inference.

The [official CHIMES download page](https://richings.bitbucket.io/chimes/download.html)
distinguishes source and reaction-data repositories. Its
[acknowledgements page](https://richings.bitbucket.io/chimes/acknowledgements.html)
requests the Richings2014a/b citations. The currently inspected data root and
README have no explicit license statement; the local source's LGPL-3.0-or-
later notice is not substituted for those data terms. This is a current
metadata check, not a newly claimed data revision; exact runtime bytes remain
bound by the accepted compute-node inventories.

These bounded findings are recorded in the existing release candidate, with
no new licensing framework, scientific source replacement or admission waiver.
No author was contacted, no mail sent, no raw source modified and no new
simulation submitted. Candidate-use approval remains pending as asked above.

#### Operator candidate-use approval — 2026-10-12

The operator subsequently replies **승인** to the exact pending request for
conditional candidate use of the selected mixed stellar SED, AGN partition
and D03/DL01 dust approximations. Tracking ID
`LRD-EFFECTIVE-SOURCE-USE-2026-10-12` records that actual decision; it is not
an external certification. The candidate now records
`operator_approved_effective_model_pending_release`, preserving the source
mass/age/spectrum/material limitations and unchanged raw input identities.

The historical contracts keep their original reference-control status and
approval fields so their accepted evidence remains reproducible. The operator
approval is an explicit candidate-use overlay, not a rewrite of the authors'
data or a claim of universal source accuracy. It does not resolve upstream
data terms, provide a committed clean source identity, or satisfy the full
physical-production preflight by itself. The former candidate-use approval
blocker is closed; the remaining release conditions stay explicit.
