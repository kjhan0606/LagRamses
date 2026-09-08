# Native runtime controls (implementation, not physical approval)

The selected single-rank/OpenMP comparison implementation is now closed.
For its fixed executable, exact inputs, reproduction commands and limits, start
with the [closeout handover](../../provenance/rt_feedback_dust_comparison_closeout_2026-09-07.md).
The wider capabilities described below are not all qualified by that one profile.
`mkrun.py` now offers this reference comparison in the CLI/GUI **Run mode** menu;
it prepares the namelist, input copies and runtime environment without launching
a job. See the [wizard guide](../../patch/cuRamses/aux/README_mkrun_gui.md).

Work/build in `/gpfs/kjhan/LRD_JWST`; repository identity is
`kjhan0606/LagRamses`. Keep isolated build directories at the same depth as
`bin` so the Makefile VPATH resolves the intended lagRamses patch first.

The live-dust build is `SNRT=1 DUST_LIVE=1 HDF5=1 USE_FFTW=0`, with NVAR=30.
Select `USE_CUDA=1` for CUDA/hybrid support, or `USE_CUDA=0` (also the omitted
option) for a toolkit-free CPU/OpenMP build. Always use separate build
directories: make does not detect changes in preprocessing flags in existing
objects. Use `-j1` for an initial clean build; do not rely on the legacy core
module ordering being a complete parallel-build dependency graph.

The CPU build uses the same host RT/dust kernels, feedback physics and HDF5
identities. `auto` sends all batches to CPU workers; `openmp` remains the
explicit host path. Forced `cuda` fails at startup. Unavailable CUDA diagnostic
entry points return an error without modifying output arrays; they do not
emulate GPU results. No CUDA headers, nvcc, driver or runtime libraries are
needed by this build. MPI, Intel Fortran/C++/OpenMP and HDF5 dependencies
remain. This is execution portability, not added or newly approved physics.

## Execution placement

| Operator | Available implementation / selection |
| --- | --- |
| Primary SNRT transport, H/He inventory cap and primary dust absorption | OpenMP or CUDA; automatic selection available |
| Stellar/AGN source evolution and mechanical deposition | Existing native CPU/OpenMP implementation; no equivalent CUDA implementation claimed |
| Dust material implicit energy/temperature solve and group emissivity | Shared FP64 OpenMP/CUDA kernel; automatic selection available |
| IR transport and local absorption/source response | Shared FP64 OpenMP/CUDA kernels; stream-availability hybrid |
| IR outer iteration, conservation sums, MPI exchange and coarse/fine correction | Existing native host implementation |

`SNRT_BACKEND=auto` (default) now means **stream-availability hybrid**, not
whole-call cell-threshold selection. `openmp` and `cuda` remain explicit
whole-call overrides. Legacy diagnostic GPU entry points remain GPU-only.

`SNRT_DUST_BACKEND=auto|openmp|cuda` separately places dust material emission
and IR transport/local absorption; when unset it inherits `SNRT_BACKEND`.
`SNRT_HYBRID_BATCH_CELLS` sets the work-batch size (default 256, admitted
range 1--1048576), not a GPU eligibility threshold. The former
`SNRT_GPU_MIN_CELLS`/`SNRT_DUST_GPU_MIN_CELLS` controls are no longer used.

For each batch an OpenMP worker tries the common cuRamses CUDA stream pool
once: a free slot sends that batch to its stream; a busy pool immediately
executes that batch on the same CPU worker. Other workers continue processing
their own batches while the GPU worker waits for its stream, not the device.
The existing `n_cuda_streams` namelist field controls pool size. Slot leases
are per MPI process/device context, NOT a GPU-wide utilization test or an
inter-process lock. Different ranks may share a physical GPU.
GPU admission also checks batch memory plus 16 MiB headroom, a 20% reserve,
and rank sharing. Forced primary/material CUDA retains whole-call admission.
Forced IR CUDA processes batches on one worker, requires a lease for every
batch and rejects missing resources rather than silently using the CPU.
Stream-ordered allocation/copies/kernels/free replace device-wide
synchronization in these operators. Initialization remains collective at
startup, never deferred until only dust-owning ranks enter.
Both processors use the same 80-iteration FP64 material solve. The original
Fortran implementation remains the independent test reference. Trial outputs
publish only after success; CUDA errors are not silently replayed on the CPU.
This adds no dust physics, namelist field, or persistent device checkpoint.
SNRT batches gather the six-neighbor snapshot from immutable old state; all
groups/directions and sequential atom-inventory consumption of a cell stay
together. Owned outputs are staged until every batch succeeds; ghost workspace
is not published. Dust likewise stages material and IR results. IR gathers
old local/remote neighbor snapshots, honors coarse-owned blocked faces and
retains the original host boundary ledgers and global conservation sums. Busy-slot
fallback happens BEFORE GPU execution; a device/physics error rejects the
transaction rather than replaying a partially executed batch on the CPU.
Timing-dependent CPU/GPU assignment can produce the measured FP32 rounding
differences, so hybrid restart is not promised bitwise deterministic.

Hydro selection remains independent: initialization, hybrid entry, local
CUDA flux dispatch, CFL, restriction and hydro synchronization honor
`gpu_hydro`; the force-gradient stream path honors `gpu_poisson`.
Initializing the shared pool for SNRT/dust must not enable these other sectors.

Device assignment follows the existing cuRamses local-rank modulo visible-device
mapping. Node-local UUID exchange identifies ranks sharing a GPU, including
rank-specific visibility masks. The forced primary memory estimate includes the wrapper's
arrays, 64 MiB headroom, a 20% reserve, and division between sharing ranks.
It is a preflight estimate, not a reservation against unrelated processes.
An allocation or kernel failure still triggers the enclosing RT transaction's
failure path; no automatic CPU replay occurs after partial GPU execution.

Host photon state is authoritative at every call; backend switches need no
unsaved persistent device state. OpenMP honors `OMP_NUM_THREADS`; when absent,
the primary host operator uses the scheduler per-task CPU budget if available,
otherwise divides available CPUs among node-local MPI ranks. This does not
change the global hydro/feedback OpenMP team settings.

## Stellar photons and restart

### Separate binary-history source work (not runtime activation)

`tools/build_cosmic_binary_histories.py` runs the **Fortran BSE engine** through
an isolated `cosmic-popsynth==4.2.0` installation. It generates full linked
histories for an explicit grid of 1--256 ZAMS binaries, rather than trying to
invert BPASS's marginal component weights or repair its nonfinite mixed ages.
It requires the unmodified `examples/Params.ini` from COSMIC commit
`f9b90f451bca014e9e0adb3a410bee8752e30e53` (SHA256
`68b2c0a7a90cd179e4f98c27428935472713e161ea0a19e5c8bf7eef8d3d04e9`).
The complete effective initial conditions, source parameters, shared IDs,
ordered phase records and surviving companions are retained. Same-age rows
must not be averaged or deduplicated.

Example (use an isolated Python environment containing that version):

```sh
python simulation/snrt/tools/build_cosmic_binary_histories.py \
  --initial-binaries simulation/snrt/tests/fixtures/cosmic_binary_grid.csv \
  --params /absolute/path/to/COSMIC-v4.2.0-Params.ini \
  --output-dir /absolute/path/to/new-binary-history-output --seed 20260907
```

This is an **unweighted comparison grid**, not a sampled SSP, a calibrated
SNIa rate, or a new simulation namelist option. BPP CO-WD disappearance is
reported separately from merger-induced massless components; the preceding
phase row does not supply instantaneous explosion mass. No N100 yields or
BPASS spectra are attached. The BSE helium-accretion destruction prescription
also differs from BPASS's Chandrasekhar channel; a current software release
does not make that prescription modern or publication-validated SNIa physics.
Channel-resolved event ejecta, population weights and an atmosphere/SED
calculation consistent with this evolution remain required for runtime use.

### Existing native photon inputs

`SNRT_STELLAR_SED` optionally supplies a native `&snrt_stellar_sed` table.
Photon rates are per **initial** solar mass, age is proper time in Myr,
and metallicity is a mass fraction. The table declares IMF identity and mass
limits, single/binary population identity and fraction, common transport
identity, and `linear_age_linear_Z` interpolation. No age/Z extrapolation is
performed. Integration is piecewise analytic in age, including crossed knots.

The production call uses actual STAR particles' `mp0`, `tpp`, `zp` and positions;
it does not turn SN mechanical energy into a spectrum. Each recursive AMR step
carries its own starting proper time into the RT source transaction. Source
failure restores the pre-injection photon state. Existing HDF5 particle fields
preserve birth mass/time and feedback progress; HDF5 additionally binds the
enabled switch and the exact stellar photon table. Tables must agree across
MPI ranks and cannot silently change on restart.
Stellar-enabled checkpoints use format 3 (primary) or 4 (primary plus IR),
so pre-stellar executables reject them instead of silently dropping future
stellar emission. Formats 1/2 remain unchanged when the stellar source is off.

`config/snrt_stellar_sed_reference_control_v1.nml` is **synthetic**, not BPASS,
and is admitted only against the matching reference transport closure. It is
not a Chabrier physical SED simply because the IMF selector is Chabrier.
The native stellar adapter admits a physical table only when its
`approved_production` status, `approval_id` and `fraction_semantics` match the
admitted common transport contract. No such joint physical table/closure is
supplied by the reference fixture; scientific input selection is still required.

AGN HDF5 checkpoints preserve radiation, thermal, jet, mass-loading and deferred
receipts. MPI source injection has one spatial owner; an owner-count reduction
rejects duplicates and publishes the commit marker to every replica. Consumption
occurs only after the global RT/chemistry/dust commit. Mechanical receipt
publication continues to use the existing collective sink-slot implementation.

New HDF5 sink checkpoints mark `sink_stat_global_schema=1`. Cloud drift
statistics are rank-local, unlike replicated BH positions: the writer now
sums them across MPI ranks, and restore contributes that sum on one rank to
the next position-update reduction. Previously only rank 0's partial sum was
saved and then restored everywhere, biasing the restarted BH/AGN position.
Old multi-rank sink dumps without this marker are rejected as incomplete;
old single-rank dumps remain readable. No old output has been rewritten.

Primary photons and same-level IR use private, at-most-16-component halo
packets on the established RAMSES emission/reception grid maps. This reduces
message count without changing the coarse/fine flux calculation or generic
hydro communication/autotuning. Temporary full-cell storage is bounded to
16 FP64 columns, in addition to the persistent radiation and packet buffers.

## Minimal native checks

The existing Makefile now has `snrt_backend_smoke` and
`snrt_stellar_source_smoke` targets. The former compares CUDA/OpenMP using the
shared cell cap, checks photon/atom ledgers, and checks atomic bad-input rejection.
The latter checks analytic photon production, split-step equivalence, pre-birth
zero emission, and age/Z bounds using the explicit reference table.
These are implementation regressions, not publication-validation gates.

`snrt_dust_backend_smoke` compares v3/v4 material **thermal** evolution against the
original native solver, including zero-density cells and error rollback.
CPU/GPU worst relative difference was 6.19e-16. The separately authorized
parallel comparison in `.parallel-runtime.luzQV6` exercises actual stars,
BH accretion, primary RT and live dust with one rank/OpenMP, one rank/CUDA,
and two ranks/two GPUs using auto placement. It does not qualify multi-node
scaling, general AMR geometry, or long integrations. See the
[parallel extension record](../../provenance/real_source_integration_progress_2026-09-07.md#operator-authorized-mpigpuopenmp-extension)
for executable identity, restart results and limits. No overall speedup is
claimed from the short, partly concurrent runs.
With the pre-hybrid `ramses_parallel_dispatch3d`, two-rank CUDA-to-CUDA restart matches all RT and dust
arrays exactly (91/94 hydro/RT arrays exact; remaining momentum differences
at most 4.14e-25 absolute). CUDA-to-OpenMP restart has a maximum
array-normalized difference of 1.61e-7; dust energy differs by 9.62e-8.

The replacement `ramses_hybrid3d` in `.hybrid-runtime.Vb2XNr` exercises mixed
CPU/GPU batches on two MPI ranks. Each rank's first primary/material call
processed three batches on CPU and one on GPU (64 cells/batch, four threads,
one stream). Hybrid-to-OpenMP restart has maximum array-normalized difference
2.07e-7 (RT), 6.55e-8 (dust energy); dust/IR balance remains below 1e-9.
`snrt_hybrid_smoke` checks mixed execution, an externally held stream forcing
all batches onto CPU, absent GPU, six-neighbor batch/ghost crossings, and
whole-call rollback on late-batch errors. This is not a speedup benchmark.
See the [hybrid implementation record](../../provenance/real_source_integration_progress_2026-09-07.md#stream-availability-hybrid-replacement).

The IR extension `.ir-hybrid.dYEXir/ramses_ir3d` additionally exercises IR
transport and local absorption on both devices in the same two-rank profile:
each rank's first call has CPU=3/GPU=1 for all four operators. Hidden-GPU
OpenMP restart has 94 finite hydro/SNRT arrays (73 exact), with maximum
array-normalized differences 2.30e-7 for RT and 5.95e-8 for dust energy.
Maximum IR balance is 6.46e-10. Independent Fortran material/IR parity is
within 1.88e-15 for thin/finite optical depth, cross-batch neighbors, a remote
ghost and a blocked coarse-owned face. This qualifies neither multi-node
scaling nor a performance advantage. See the
[IR extension record](../../provenance/real_source_integration_progress_2026-09-07.md#ir-stream-hybrid-extension).

For this NVECTOR=500/NVAR=30 executable, use `OMP_STACKSIZE=512M` (and ensure
any Intel `KMP_STACKSIZE` agrees). The compiled `godfine1`+`unsplit` stack
frames alone consume 336 MiB: an inherited 128 MiB worker stack crashed when
hydro ran on a worker rather than the main thread. The mkrun comparison
environment supplies 512M explicitly. Budget memory for the chosen team size.

The live IR solve uses a relative energy tolerance of `1e-9`. Separately
advected dust mass/energy can arrive slightly below the bath temperature.
Only a relative bath-energy deficit within that same solve tolerance is
admitted (with a 64-ULP numerical minimum). The original material energy is
retained in the closure calculation, so restoring the bath temperature uses
the numerical energy-error budget; it is not unrecorded background heating.
Larger deficits still reject the entire coupled transaction. Native tests
cover charged corrections and rejection without partial state publication.

Physical production remains blocked where the selected yield/fate package,
stellar/AGN SED or dust material data have not been scientifically admitted.
No missing physical approval is synthesized by this runtime implementation.

## Reusable integrated reference input

`config/snrt_stellar_agn_dust_reference_smoke.nml` reproduces the actual-star,
BH-accretion and live-dust wiring control. Copy it to `run.nml` in a NEW
directory, together with `config/snrt_agn_accretion_smoke.ic_sink` as `ic_sink`.
Use the NVAR=30 build above. The environment is:

```sh
OMP_NUM_THREADS=2
SNRT_RT_ENABLE=1
SNRT_BACKEND=auto
SNRT_AGN_MODEL=partition_reference_v1
SNRT_REDUCED_C=0.01
SNRT_RT_LEVEL=3
SNRT_ALLOW_REFERENCE_CONTROL=1
SNRT_GROUP_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
SNRT_SECONDARY_TABLE_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_secondary_table_contract_v1.nml
SNRT_STELLAR_SED=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_stellar_sed_reference_control_v1.nml
SNRT_DUST_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/dust_native_reference_control_v3.nml
```

Export these values for the executable and unset any `SNRT_DRIVER_TEST_SEED_SOURCE`
or RT failure-injection controls. This profile intentionally uses CPU hydro and
Poisson, with only primary SNRT placed automatically. Before launching, check
the actual copied namelist and free space: this short evolution control has
`nstepmax=3`, `noutput=1`, `aout=2`, `tout=1e30`, `foutput=1`, `fbackup=1000000`;
three periodic dumps of about 4.6 MiB each are expected. No scheduled time is
reached. To reproduce the demonstrated serial backend-switch restart, copy
checkpoint 2 to another NEW directory and change only `nrestart=2`, optionally
hiding GPUs with `CUDA_VISIBLE_DEVICES=''` and retaining `SNRT_BACKEND=auto`.
This integrated stellar profile was exercised with one rank; the separate
two-rank reference used AGN without star formation. These are not the same test.

## Draine optics in the native live-dust input

The existing thermal builder can also emit the native v3 namelist. From the
project root, choose NEW output paths and run:

```sh
JAX_PLATFORMS=cpu simulation/snrt/.venv/bin/python simulation/snrt/tools/build_draine_dust_thermal.py \
  --source external/draine_wd01_rv31/kext_albedo_WD_MW_3.1_60_D03.all \
  --output /path/to/new/dust_thermal.json \
  --native-output /path/to/new/dust_native.nml \
  --native-source-ledger simulation/snrt/data/p4_pilot_agn_photon_ledger.json \
  --native-background-k 10 --reference-heat-capacity-per-h 1e-24
```

Set `SNRT_DUST_CONTRACT` to that native file and retain
`SNRT_ALLOW_REFERENCE_CONTROL=1`. The raw opacity table supplies absorption
and dust mass per H. Primary opacity is sampled at each ledger group's mean
energy (an explicit grey approximation, not a full SED-weighted average).
Native IR uses 136 quadrature nodes across the supplied raw energy domain and
65 temperatures including the 10 K bath with these defaults; the separate JSON
sidecar retains its older primary-group thermal representation.

The constant heat capacity above is a **test value**, not inferred from the
opacity data. The exporter therefore never emits production approval. It
checks ledger group edges, ordering, representative energies and edge checksum;
raw-source and ledger checksums are retained. Scattering and stochastic grain
heating are not supplied by this conversion.

For a particle-free cooling control based on `config/dust_live_ir_smoke.nml`,
leave `sf_birth_properties=.false.`: its particle-only element list is not
initialized in a no-PIC run. The demonstrated Draine run used CPU hydro,
Courant factor .005, dust mass per gas mass .006024096385542169 and initial
dust energy per gas mass 9.01356598145909e-16 in that fixture's code units
(20 K for the explicit test capacity). Do not reuse these energy units with
another material capacity or units definition.

## SNIa population identity

The physical SNIa handoff is checked against the actual configured IMF,
population model and binary fraction before runtime initialization completes.
`config/fp2_snia_runtime_contract_v1.nml` is the approved Kroupa/binary/.5
baseline; default Chabrier/single-star particles do not match it. The handoff's
IMF conversion factor does not authorize an undeclared change of target IMF.
This check does not replace the existing AGB white-dwarf reservoir requirements.

## High-mass models: native source history and HDF5 restart

`STELLAR_ENRICHMENT_PARAMS` now accepts `high_mass_preset`:
`source_consistent` (default), `wind_only_collapse`, or `mixed_remnant`.
The latter requires `high_mass_remnant_adjust_max_fraction > 0` and <=1,
defined relative to INITIAL stellar mass. Other presets require this value
to be zero. Missing or invalid choices do not inherit a previous override.
There is no conservation-disable option. Legacy mode rejects an override.

The native connection is implemented:
namelist -> source-node history/endpoint resolution -> common IMF source cells
-> cumulative differences -> existing RAMSES mechanical deposition and stellar
mass/progress commit -> HDF5 identity and restart checks. It is not a Python
runtime. `source_consistent` requires identical declared source identities and
mass closure; `wind_only_collapse` preserves winds and suppresses terminal
ejecta/energy/momentum; `mixed_remnant` preserves ejecta and adjusts the remnant
only within the explicit limit. The remnant is a baryonic residual, not a
gravitational black-hole mass. Raw source rows remain unchanged.

Activation requires `fate_policy='user_selected_model_v1'` together with
`high_mass_history_path`. This separate user-responsibility route does NOT
create a scientific approval: leave `fate_map_sha256` and `fate_approval_id`
empty. The original approved-package path and its restrictions remain intact.
The route requires channel-resolved, per-star cumulative yields, wind and SNII,
PIC and HDF5 output/restart. It now supports two explicit populations:
single-star/binary fraction zero/SNIa off, or an effective binary SSP with
positive binary fraction, SNIa on and AGB enabled as the WD supplier (AGB upper
mass <=8 Msun). History and SNIa handoff must match the configured IMF,
population and binary fraction. No second binary-fraction factor is applied
to an empirical DTD whose normalization already includes it. This is not a
binary stellar-evolution calculation; actual yield inputs must justify that
effective population approximation. PISN remains disabled.

For the approved N100/Maoz input, explicitly select Kroupa (`imf_id=1`),
`population_model='binary_ssp'`, `binary_fraction=.5`, `use_snia=.true.` and
`use_agb=.true.` and set `PHASE0_SNIA_RUNTIME_CONTRACT` to
`config/fp2_snia_runtime_contract_v1.nml`. Default Chabrier/single-star behavior
is unchanged and does NOT silently inherit that approval. SNIa uses the real
per-event input and DTD, not canonical per-star channel-4 rows. Existing
population accounting debits only the AGB WD reservoir, retains NS/BH remnants
and rejects a WD shortfall. The supported SNIa metallicity factor remains the
explicit approved constant-unity model, not an invented Z-dependent DTD.

`PHASE0_YIELD_TABLE` supplies the canonical per-star cumulative source table.
The history file is `/stellar_high_mass_history/` version 1 (Fortran namelist),
with up to 512 nodes. See `tests/fixtures/phase0/high_mass_history.nml` for the
complete schema: declared model/wind/terminal identities and coordinates,
IMF/population/binary/support identity, mass, Z, lifetime in years and raw
terminal outcome (0 failed / 1 exploded). These labels describe user-supplied
inputs; they are not independent proof of source authenticity. Every exact-Z
branch must contain ascending nodes from 40 through 120 Msun, with matching
wind/SNII rows and an endpoint at each declared lifetime.

The explicit timing approximation is `wind_linear_terminal_step`: linearly
interpolate cumulative wind rows, emit terminal ejecta at the node lifetime,
and freeze both afterwards. Pre-terminal SN leakage is rejected. Mass uses
nearest-source-node cells, lower-node ties, and budget scaling by M/M_node
(constant mass fractions inside each cell), not interpolation of discrete
explosion outcomes. IMF integration splits at 40 Msun and source-cell edges
and uses analytic IMF mass weights shared by the channels. The history namelist
now accepts `metallicity_policy='exact_nodes'` (backward-compatible default) or
`'linear_Z_cumulative_mixture'`. The latter evaluates BOTH neighboring source-Z
histories at the requested age, then combines all cumulative budgets with fixed
linear birth-Z weights. This is an SSP expectation mixture, not interpolation
of individual explosion fates or terminal lifetimes. Each branch retains its
own event time, resolved remnant and wind budget; fixed nonnegative weights
preserve mass closure, monotonic gross release and split-step equivalence.
At least two Z branches are required. Z=0 is supported when supplied; no
logarithmic-Z assumption, physical extrapolation or edge clamp is introduced.
Only relative 32-epsilon roundoff may match an endpoint. Rotation/engine
alternatives still require a separate declared source branch.

The selected source-cell convention also scales low-mass AGB payloads by
M/M_node so WD return uses the same per-initial-mass normalization. Constant
age payloads remain exactly constant instead of acquiring an artificial
negative energy increment from interpolation roundoff. Legacy table behavior
outside this selected route is unchanged. Energy columns are non-bulk source energy;
the existing bridge separately accounts for directed kinetic energy.

HDF5 particles save the enabled switch and `stellar_source_model_values`:
actual canonical numerical rows, source labels, resolved remnant/correction,
timing/mass policy version and population/channel settings. Every MPI rank
must agree. Restart refuses activation/deactivation, altered values or altered
preset before restoring/evolving particles. Moving identical input files is
allowed; a path or user-typed hash alone is not treated as an identity check.
The identity additionally binds the Z policy, revised low-mass interpolation
semantics, and every consumed SNIa population/event/coupling value, including
DTD normalization and energies. Existing exact-Z/high-mass-only/SNIa-off
identity layout is retained; changed low-mass semantics require a fresh run.

Use `mkrun.py --mode gui` or the CLI full parameter editor's feedback section
to inspect the choices. Both share `ramses_nml_generator.py`; hydro generation
now emits the complete required stellar population fields in ONE namelist
group instead of the old feedback-mode-only stub. Generated inputs still
need real source assets and runtime admission; `CHANGE_ME` is not runnable.

The reusable `config/high_mass_feedback_reference_smoke.nml` is a SYNTHETIC
four-step, level-3 CPU hydro/stellar-feedback wiring control, not a physical
yield package. Copy to `run.nml` in a NEW directory and use the NVAR=30 build
above with `PHASE0_STELLAR_ENRICHMENT=1`. Unset `SNRT_AGN_MODEL`,
`SNRT_STELLAR_SED`, `SNRT_DRIVER_TEST_SEED_SOURCE`; set `SNRT_RT_ENABLE=0`,
`OMP_NUM_THREADS=2`, and
`PHASE0_YIELD_TABLE=/gpfs/kjhan/LRD_JWST/simulation/snrt/tests/fixtures/phase0/high_mass_history_yields.dat`.
The profile contains its absolute history path. Its tiny legacy reader table
is only for existing particle initialization, not the mechanical source.
Output policy is `noutput=1,aout=2,tout=1e30,foutput=1,fbackup=1000000`;
four measured dumps total about 1.62 MB. Review output/storage before launching.
The native fixture `high_mass_history_test.f90` checks all three presets,
event timing, split-step equivalence, ledger closure and invalid inputs.
Live single-rank/MPI2 and restart evidence is recorded in
`provenance/real_source_integration_progress_2026-09-07.md`.

`config/high_mass_snia_z_reference_smoke.nml` exercises the combined path:
real approved N100/Maoz SNIa plus SYNTHETIC wind/SNII/AGB on Z=0,.02 branches,
with initial gas Z=.01. Use `PHASE0_YIELD_TABLE` pointing to
`tests/fixtures/phase0/high_mass_snia_z_yields.dat`, and the SNIa environment
above; keep the RT environment disabled as in the previous control. The profile
uses 1 Gyr/code-time so the four-step run crosses the 40 Myr SNIa delay.
Its intentionally synthetic early WD formation is not physical AGB timing.
Native tests and live single-rank/MPI2/restart checks passed; DTD and Z-policy
changes reject on restart. The SNIa namelist's two source-code-style `&`
continuations were removed for standard GNU/Intel input compatibility; all
physical numbers remain unchanged.

### Actual LC18 wind with explicitly selected approximations

`tools/build_lc18_native_wind.py` now connects checksum-verified LC18 Set R
table8 wind isotope ejecta and table7 phase-summed lifetimes to the native
`user_selected_model_v1` / `wind_only_collapse` route. Python only prepares
small input files offline; Fortran evaluates and deposits feedback in RAMSES.
Select one rotation (0, 150 or 300 km/s); each branch contains 40, 60, 80 and
120 Msun at source Z=0.00003236, 0.0003236, 0.003236 and 0.01345. There is no
extrapolation to primordial or super-solar metallicity. IMF defaults to
Chabrier; `--imf-id` also accepts the existing 0, 1 and 4 alternatives.

Example local comparison input (1000 km/s is a test choice, NOT a calibrated
LC18 wind velocity or a hidden default):

```sh
python3 simulation/snrt/tools/build_lc18_native_wind.py \
  --output-dir /gpfs/kjhan/LRD_JWST/NEW_LC18_INPUT \
  --rotation 0 --wind-speed-km-s 1000 \
  --timing uniform_until_terminal --composition as_tabulated_mean \
  --energy isotropic_thermalized
```

The required options explicitly select uniform release until the actual
lifetime, constant time-averaged ejecta composition, and wind energy
E=0.5*M*v^2 thermalized under an unresolved isotropic-wind approximation.
Net directed wind momentum is zero; the bridge separately adds stellar bulk
ejecta kinetic energy. Isotopes are assigned to their tabulated parent
elements without decay. Total mass is the sum of all isotope ejecta, including
untracked elements; rounded evolutionary masses are not forced to match it.
There is no terminal explosion and the remaining mass becomes the remnant.
The exporter supports only a single-star population; it does not silently
enable the separately approved effective-binary SNIa model.

The new optional history field `net_yield_policy` defaults to `supplied`.
This exporter selects `unavailable_diagnostic_zero`: all net columns must be
zero placeholders, explicitly marked unavailable rather than interpreted as
zero physical net production. Actual gross ejecta drive feedback. Native
startup logs this distinction, and restart identity binds it along with all
source values, approximations and wind-speed-dependent energy.

Copy `config/lc18_wind_local_smoke.nml` into a NEW run directory and replace
its `CHANGE_ME_LC18_INPUT/history.nml` path; point `PHASE0_YIELD_TABLE` to the
generated `yields.dat`. Retain the disabled RT environment described above.
The profile uses 100 Myr/code-time, four fixed level-3 steps and the same output
policy `noutput=1,aout=2,tout=1e30,foutput=1,fbackup=1000000`.
Review effective namelist, output count and storage before launching.
Measured four-dump storage is 2,189,868 bytes. Native optimized Intel and GNU
bounds-checked source tests, live mass/energy checks, same-input restart and
changed-wind-speed restart rejection passed; details are in the progress log.

**Scope limit:** the actual LC18 wind connection is implemented under these
explicit comparison approximations. It is not a phase-resolved wind spectrum,
composition/velocity model, a physical terminal-explosion model, or a complete
AGB/WD source package. The previous combined SNIa control still has synthetic
AGB inputs. No all-channel production/publication approval follows from this
smoke run. Raw sources are unchanged, historical review-only contracts are
not promoted, and locally generated numerical files are not added for source
redistribution without resolving the applicable permission.

### KL16 selected gross normalization

The operator-selected KL16 source reader now provides normalized gross ejecta
in each active row's `selected_ejecta` payload. It uses the sum of ALL 78 raw
element masses as the denominator, scales to the selected expelled mass, and
leaves the selected remnant unchanged. Reduce to eleven tracked elements only
after this step; use `untracked_ejecta_msun` for the remaining metal budget.
The payload records policy `all_listed_elements_to_selected_expelled_mass_v1`,
its normalization factor and original sum. Raw data/fingerprints, initial
compositions, raw net-yield diagnostics and the two excluded nodes are retained
unchanged. Normalized net yields remain unavailable, not copied from raw net
diagnostics. This resolves the gross mass-budget mismatch at input preparation;
it does not by itself provide the missing lifetime/release/energy inputs or
activate a complete physical AGB/SNIa run. No runtime namelist default changes.

### Actual KL16 lifetime, envelope/WD and wind-energy connection

The subsequent operator-approved implementation connects these missing inputs
through `tools/build_kl16_lc18_native.py`. `data/kl16_stellar_lifetimes.csv`
transcribes 71 durations from [Karakas 2014 Table 1](https://arxiv.org/html/1408.5936v1#S3.T1)
and [Karakas & Lugaro 2016 Table 1](https://arxiv.org/html/1604.02178v1#S2.T1).
All 62 active yield nodes match exactly in initial mass, Z, initial helium and
overshoot. M_mix is a post-processing choice, not a different stellar lifetime.
The source's total stellar duration to the AGB endpoint is used, not an
empirical lifetime fit or a post-AGB/WD cooling time.

The combined exporter selects 58 ordinary CO-core AGB nodes over the common
1--6 Msun support at Z=.007,.014,.03. Hybrid CO(Ne), ONe and commented 8 Msun
yield models do not supply N100 WDs. The history option
`agb_release_policy='terminal_step'` releases the normalized envelope and
creates the WD together at each actual node's endpoint. Before that endpoint,
both are zero. This is an explicit terminal-envelope approximation, not a
phase-resolved AGB wind. The native evaluator uses nearest source-mass
fractions with fixed linear-in-Z mixtures; sparse, independently timed source
nodes do not need manufactured rectangular age grids. Other channel grid
checks and the default `cumulative_linear` policy are unchanged.

Example isolated comparison input (both speeds are explicit test choices,
NOT velocities measured or calibrated by KL16/LC18):

```sh
python3 simulation/snrt/tools/build_kl16_lc18_native.py \
  --output-dir /gpfs/kjhan/LRD_JWST/NEW_KL16_LC18_INPUT \
  --rotation 0 --massive-wind-speed-km-s 1000 --agb-wind-speed-km-s 15 \
  --agb-release terminal_envelope --agb-energy isotropic_thermalized \
  --population single --imf-id 1
```

AGB energy is `0.5*M_return*v_wind**2`, thermalized as an unresolved isotropic
wind; directed wind momentum is zero. Stellar bulk momentum/kinetic energy
remain the native deposition bridge's responsibility. Net yields remain
explicitly unavailable; normalized gross ejecta and untracked metals drive
the return. IMF defaults to Chabrier; this comparison explicitly uses Kroupa.
The builder refuses existing output directories and binds all source hashes,
model coordinates, approximation choices and energy values to restart identity.

Use `config/kl16_lc18_local_smoke.nml` in a NEW run directory, replace the
history placeholder and set `PHASE0_YIELD_TABLE` to the generated table.
Unset `PHASE0_SNIA_RUNTIME_CONTRACT` for this SNIa-off control. Disable live RT
as above. Review the effective namelist/output budget before launch. The
four-step native CPU hydro run and same-input restart pass, with conserved
gas+star mass and positive internal energy; details and binary/input hashes
are in the existing progress record. No Python is called by RAMSES.

**SNIa-on qualification is NOT complete.** `--population snia_baseline
--imf-id 1` prepares the existing effective-binary population without changing
its approved DTD. This DTD starts at 40 Myr, whereas the selected WD supplier
first returns WDs at 68.89 Myr. The native combined run rejects the causal WD
deficit (exit 82), rather than borrowing future WDs, clipping the rate or
silently delaying the DTD. Checks immediately before WD events also prevent a
large timestep from hiding an earlier deficit. Resolving this requires an
explicitly selected, consistent SNIa population/WD-supply model; the current
AGB+LC18 positive run is not an all-channel production approval.

### Optional Fishlock low-Z AGB extension

The default KL16/LC18 export is unchanged. To extend the common active-channel
metallicity hull from `.007--.01345` to `.001--.01345`, add
`--low-z-agb fishlock2014_raiteri96` to the combined builder above. This selects
15 Fishlock et al. (2014) gross-yield nodes at Z=.001, initial mass 1--6 Msun;
the 7-Msun ONe model is excluded. The pinned `yield_z001.txt` is taken from
the existing COLIBRE source snapshot. Unlike the KL16 rows, these nodes use
the sum of **all** tabulated gross element masses without renormalization;
the remnant is the complementary initial mass. Tracked elements are a subset,
not the total expelled mass. Net-yield diagnostics remain unavailable.

This is an explicit **cross-model comparison**, not source-matched Monash
evolution: [Fishlock Table 1](https://arxiv.org/html/1410.7457v1) does not give
total stellar lifetimes. Low-Z event ages use the published Raiteri et al.
(1996) Padova fit as given in
[Valiante et al. (2009), equations 3--6](https://arxiv.org/html/0905.1691v1).
The fit is used only for these Z=.001, 1--6 Msun nodes. KL16 source lifetimes
are unchanged. The terminal-envelope/WD approximation and native linear-Z
mixture retain each source node's own event age; no yield extrapolation or
interpolated individual stellar track is claimed.

Select the exported history with `high_mass_history_path` and the table with
`PHASE0_YIELD_TABLE`, as above. No new RAMSES namelist field is required.
For the existing SNIa comparison also use `--population snia_baseline
--imf-id 1` and its explicitly selected effective-SSP runtime contract.
This does not solve the strict WD-supply/40-Myr DTD mismatch, establish a
microscopic binary population, or extend the LC18 upper Z boundary.

The existing native source fixture accepts an optional fourth argument
`fishlock2014_raiteri96` for the extended input. It checks all 73 AGB nodes,
time splitting at Z=.004, full effective-SSP DTD mass closure at six Z values,
and rejection outside the common hull. A four-step native RT/feedback/dust
run and restart at Z=.004 are recorded in the population bundle progress.
These short runs test wiring, not a galaxy calibration or late AGB evolution.

### Effective SSP SNIa accounting (implemented; bounded native verification)

The operator-approved phenomenological SSP route retains the empirical DTD
and N100 event source but debits SNIa ejecta from the remaining stellar-particle
mass AFTER the interval's generic AGB/SNII/wind returns. This is a new accounting
choice, not proof that KL16 supplies early N100 WDs. The strict WD-supply route
and its causal rejection remain unchanged and remain the default.

Select the separate example contract explicitly:

```sh
export PHASE0_SNIA_RUNTIME_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml
```

Its `&snia_population_realization` group adds
`mass_accounting='effective_ssp'` and
`accounting_approval_id='SNIA-EFFECTIVE-SSP-2026-09-07'`. Omitting these fields
in the original contract retains `strict_wd` with an empty accounting approval.
Unknown modes, missing effective approval, and a strict mode carrying the
effective approval reject. The three original DTD/event/thermal group approval
and source bindings remain unchanged; the new approval covers accounting only.
The existing Kroupa/binary population applicability is NOT broadened to default
Chabrier or another IMF by this selector.

Native startup reports the mode and its unresolved-population limitation.
Effective mode does not call the WD reservoir setter or alter the generic
living/remnant partition. The existing N100 event-budget arithmetic is reused
with total SSP capacity as its mass limit; its historical `wd_reservoir_debit`
field is the equal N100 ejecta mass, not evidence of an actual WD inventory.
Persisted mass and generic cumulative returns reconstruct prior Ia return,
which is also checked against the cumulative DTD. All channel deltas and the
particle mass/progress are staged before commit; overdraw rejects the interval
rather than clipping its rate. The zero-terminal-remnant invariant still applies.
Restart binds the new accounting identity and rejects switching modes in either
direction. The historical strict identity layout is preserved.

These are external SNIa contract fields, not additions to the main lagRamses
namelist; no `mkrun.py`/`ramses_nml_generator.py` field or default changed.
Actual KL16+LC18 CPU hydro with effective SNIa and same-input HDF5 restart pass
within the selected comparison. Common source Z support is .007--.01345,
not the broader AGB-only grid; wind speeds remain explicit comparison choices.
Evidence is in the effective-SSP section of the
[source integration record](../../provenance/real_source_integration_progress_2026-09-07.md).

Microscopic binary evolution is a separate medium-term study, not a prerequisite
for using this effective route. These bounded checks do not qualify all source
physics or a simultaneous RT/AGN/dust production run. See the operator decision in the
[completion plan](../../provenance/production_completion_bundles_2026-09-07.md).

### Effective SNIa / AGN RT / real-opacity dust integrated control

`config/kl16_lc18_snia_agn_dust_smoke.nml` combines the real KL16/LC18 history
and effective SNIa with accepted BH accretion, reference AGN transport and live
Draine-opacity dust. Copy it into a NEW run directory, replace the history
placeholder, and copy `config/kl16_lc18_snia_agn_dust_smoke.ic_sink` as `ic_sink`.
Use the integrated reference environment above, with these changes:

- Unset `SNRT_STELLAR_SED`: the existing Chabrier/single-star control does NOT
  match this Kroupa/binary population. Stellar mechanical energy is not a SED.
- Set `PHASE0_YIELD_TABLE` to the actual `--population snia_baseline --imf-id 1`
  KL16/LC18 export and `PHASE0_SNIA_RUNTIME_CONTRACT` to the effective contract.
- Set `SNRT_DUST_CONTRACT` to the native Draine optical/thermal export (the
  recorded local example is `.physical-inputs.A0XtMt/dust_native.nml`). Its
  constant TEST heat capacity is still a reference approximation.
- Use `SNRT_BACKEND=openmp` for this demonstrated profile; all GPU hydro flags
  are false. The AGN spectrum remains the existing reference group contract.

Four fixed level-3 steps use noutput=1, aout=2, tout=1e30, foutput=2,
fbackup=1000000. Two periodic dumps are about 55 MB each. Output numbering is
sequential: checkpoint **1 is coarse step 2**. To resume to step 4, copy
`output_00001` to another NEW directory, set `nrestart=1`, and keep source inputs
and all environment contracts identical. Canonical BH coordinate sync now uses
the freshly rebuilt sink map rather than stale fine-level particle lists.

This measured control reaches 102.39049 Myr, with oldest processed stellar age
51.36384 Myr: wind/SNIa, AGN primary absorption and dust IR are exercised; AGB
input is loaded but its first terminal release is not yet reached. The standalone
AGB evidence above remains separate. The seed is an explicit numerical control,
not a realistic galaxy initial condition. Mass accounting must include the AGN
radiative rest-mass loss and count sinks once, not sum their cloud tracers.
See the combined-control section of the existing integration record for results,
the retained failed runs and remaining physical-SED/heat-capacity limits.

### Tabulated dust material energy (native v4)

The native IR material solve now accepts a monotonic internal-energy table
`U(T)` instead of the v2/v3 constant heat capacity. In v4 the conserved volume
energy is `nH * relative_dust * U(T)`. Both the initial temperature inversion
and implicit radiative cooling use **linear U in log(T)**, with no physical
extrapolation. This is not `C(T)*T`. Primary absorption is accounted once in the
coupled IR/material solve. The full U array is bound to the HDF5 checkpoint;
changing it rejects restart. Versions 2/3 keep their original representation.

Use the existing `tools/build_draine_dust_thermal.py` with
`--material-energy-table material.json`, replacing (not combining with)
`--reference-heat-capacity-per-h`. The material JSON requires:

```json
{
  "schema": "snrt_dust_material_energy_v1",
  "source_id": "YOUR_VERIFIED_MATERIAL_SOURCE",
  "source_url": "YOUR_PRIMARY_SOURCE_URL",
  "composition": "YOUR_EXPLICIT_GRAIN_COMPOSITION",
  "energy_zero": "U(0)=0; no zero-point energy",
  "temperature_k": [10, 20, 100],
  "internal_energy_erg_g": [100, 800, 100000]
}
```

Numbers above are **synthetic format examples**, not physical data. Temperature
and specific energy must be positive, finite and strictly increasing, covering
the requested emission temperature interval and bath. The exporter preserves
material knots, converts erg/g to erg/reference-H using the optical table's
reference dust mass, and writes `material_sha256`,
`internal_energy_per_h_erg_input`, and zero constant capacity in the v4
namelist. Material/Planck knot union must fit the 256-temperature bound.

The exporter still writes **reference-only** status and empty approval.
Supplying a thermodynamic table does not establish that its composition and
grain-size approximation match the optical mixture. No physical table or new
population assumption is selected automatically. The existing coupled
KL16/LC18/SNIa/AGN fixture remains v3 until such input is explicitly selected.
The v4 wiring fixture is `config/dust_native_reference_control_v4.nml`.

### DL01 bulk material connected to the coupled control

`data/dust_dl01_bulk_030_material_v1.json` now provides a literature-based
graphite/silicate internal-energy model, not the synthetic cubic fixture.
`tools/build_dl01_dust_material.py` integrates the normalized Debye vibrational
mode distributions from [Draine & Li (2001)](https://arxiv.org/abs/astro-ph/0011318).
Graphite has one 2D mode at 863 K and two at 2504 K per atom; the silicate
prescription has two 2D modes at 500 K and one 3D mode at 1500 K. Energy excludes
zero-point energy. The implementation recovers three kB per atom for the
high-temperature heat capacity and uses MgFeSiO4 as the explicit representative
silicate mass per atom. These are bulk limits, not finite-grain PAH modes.

The supplied **comparison** selects graphite mass fraction 0.30 and silicate
0.70. This choice is not inferred from WD01 total opacity and is not an exact
reconstruction of its size-dependent PAH mixture. All material shares one
temperature. Stochastic heating, sublimation, and material-dependent separate
temperatures remain outside this closure. No production approval is assigned.

To select another mass fraction or temperature grid, generate new files:

```bash
simulation/snrt/.venv/bin/python simulation/snrt/tools/build_dl01_dust_material.py \
  --graphite-mass-fraction 0.30 --output NEW/material.json
simulation/snrt/.venv/bin/python simulation/snrt/tools/build_draine_dust_thermal.py \
  --source external/draine_wd01_rv31/kext_albedo_WD_MW_3.1_60_D03.all \
  --native-source-ledger simulation/snrt/data/p4_pilot_agn_photon_ledger.json \
  --material-energy-table NEW/material.json \
  --output NEW/thermal.json --native-output NEW/dust.nml
```

The tracked export is `config/dust_dl01_bulk_030_reference_v4.nml`.
For the demonstrated coupled run use
`config/kl16_lc18_snia_agn_dl01_dust_smoke.nml` with the environment of the
KL16/LC18/SNIa/AGN control above, except point `SNRT_DUST_CONTRACT` at this v4
file. Its primitive dust-energy IC has already been converted to U(20 K).
Do not reuse the old constant-capacity IC energy or an old v3 checkpoint.
Copy the shared `.ic_sink` fixture (exactly 12 numeric fields) and resolve the
history placeholder as before. The reference opt-in remains necessary.
The existing profile with a constant test capacity is retained unchanged.

Stellar radiation is disabled unless `SNRT_STELLAR_SED` is explicitly supplied.
The user-approved independent BPASS comparison below is now available; the
existing exact-population v1 control must not be relabeled to bypass matching.

### BPASS independent radiation population (reference comparison, v2)

The existing feedback population remains Kroupa 0.08--120 Msun, binary SSP
fraction 0.5 in the coupled profile, including the effective-SSP SNIa model.
BPASS instead retains its own `bin-imf135_300` identity: slopes -1.30/-2.35,
break 0.5 Msun, mass range 0.1--300 Msun and BPASS's built-in binary recipe.
It is not Kroupa and its binary recipe is not represented by the feedback's
scalar fraction. In the native v2 table, IDs/fraction -1 mark this distinction.
No feedback mass, momentum, thermal energy, IMF or DTD parameter is changed.

For a **fresh**, deliberately selected comparison run with the coupled profile
above, add this environment setting (retain the reference-control opt-in):

```bash
export SNRT_STELLAR_SED=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_stellar_sed_bpass_independent_v2.nml
```

The explicit `population_binding='independent_radiation_reference'` is accepted
only with `status='reference_control'`, never `approved_production`. Version 1
retains exact feedback population matching and its original restart identity.
Version 2 binds the rates, source SHA, population, escape fraction and policies
to MPI consensus and HDF5 restart. It cannot be enabled on an old source-free
checkpoint, or changed on restart. Startup logs name the independent model.
No new main RAMSES namelist field or generator option is introduced.

`tools/build_bpass_native_sed.py` converts the pinned local Galacticus BPASS
HDF5 into a ~247 kB native namelist (52 ages, 13 Z, 9 groups). The 531 MB HDF5
is an **offline input**, not loaded by lagRamses. Its Lsun/Hz spectra already
have the upstream 1e6-Msun burst normalization removed; the exporter integrates
photon rates per initial Msun and applies the requested escape fraction once.
The supplied comparison uses escape fraction 1.0. BPASS's log-solar coordinate
is converted with its own Zsun=0.02; this does not change feedback solar Z.

Explicit limitations, recorded in `data/snrt_stellar_sed_bpass_independent_v2.json`:

- Birth to 1 Myr holds the first 1-Myr spectrum. This is a named comparison
  approximation, not an additional BPASS measurement. No other age/Z
  extrapolation is allowed (source range Z=1e-5--0.04, age up to 100000 Myr).
- The missing wavelength tails are zero. Integration clips to the measured
  1--100000 Angstrom domain before inserting boundaries, so it does not invent
  a triangular IR tail from the 0.01-eV group edge.
- Common reference AGN/stellar group mean energies and cross sections remain
  the transport closure. The new SED supplies photon counts, not a new
  spectrum-specific energy/cross-section calculation.
- This deliberately mixed population is not a common-population consistency
  claim or a production/publication approval of the combined physical model.

Native smoke tests can be run with the existing Makefile
`snrt_stellar_source_smoke` target and
`python3 simulation/snrt/tests/bpass_native_sed.py --native-smoke PATH_TO_BINARY`.
The latter checks the actual parser/integrator, legacy matching, spectral
partition/no-tail leakage, age splitting/bounds and nine invalid inputs.
Measured live/restart evidence is recorded in
`provenance/real_source_integration_progress_2026-09-07.md` (repository root).
# Explicit ordinary CCSN source extension (2026-09-07)

The existing offline builder now accepts `--massive-source lc18_set_r` and
requires `--snii-energy-erg` for that choice. It supplies LC18 13--120-Msun
source nodes: table9 winds plus table8-minus-table9 terminal ejecta at
13/15/20/25; Set R wind-only nodes at 30 and above. Default `wind_only` output
is unchanged. Python prepares input files; RAMSES consumes them in Fortran.

```sh
python3 simulation/snrt/tools/build_kl16_lc18_native.py \
  --output-dir NEW_INPUT_DIRECTORY --rotation 0 \
  --massive-wind-speed-km-s 1000 --agb-wind-speed-km-s 15 \
  --agb-release terminal_envelope --agb-energy isotropic_thermalized \
  --population snia_baseline --imf-id 1 --low-z-agb fishlock2014_raiteri96 \
  --massive-source lc18_set_r --snii-energy-erg 1e51
```

These speeds/energy are explicit comparison parameters, not LC18 measurements.
Use the generated history v2, `PHASE0_YIELD_TABLE`, and channel lower masses
`13,1,13,3,140` (upper `120,6,120,8,260`). Existing high-mass presets apply
only at >=40. The nearest-source-cell model uses mass-fraction scaling,
including energy; it is not discrete SN sampling. No 8--13-Msun yield/fate is
invented. No extrapolation outside the common active Z hull .001--.01345.

`mkrun.py` terminal/GUI Run mode `comparison_ccsn` packages the locally tested
input and new executable from `.ccsn-source.lobKc9/`, without launching it.
It leaves the original comparison modes and radiation/AGN/dust defaults alone.
The selected BPASS SED is still an independent radiation population, not a
new same-population proof. A four-step smoke run does not guarantee SN ages.
Full scope/evidence: `provenance/population_sed_z_snia_bundle_plan_2026-09-07.md`.

## Physical-input extensions (2026-09-08; explicit comparisons)

Sections #1 physical-source/native coupling, #2 additional RT/AGN/dust physics,
and #3 backend/operations expansion are preapproved in full (operator update,
2026-09-08). Earlier approval holds below/in historical records are superseded;
this authorization is not a declaration that missing physical inputs exist.
The following are optional source selections, not changes to the defaults or
claims that the whole physical-input programme is complete.

The combined builder additionally accepts:

- `--massive-wind-timing phase_mass_loss_or_uniform`: actual LC18 Table5
  cumulative mass-loss timing, normalized to the unchanged integrated wind
  budget. Six of 36 rotation-zero nodes have zero printed phase loss but
  nonzero isotope winds; these explicitly keep uniform timing. Phase-specific
  composition and speed are **not** supplied by this option.
- `--agb-source-scope kl16_envelopes_to7`: KL16 1--7 Msun, including the
  Z=.007, 7-Msun hybrid CO(Ne) envelope. Its remnant remains in the total
  ledger but is excluded from the strict CO-WD SNIa supply. With
  `--low-z-agb fishlock2014_raiteri96`, the Fishlock 7-Msun ONe envelope is
  also included; its remnant is likewise excluded from strict CO-WD supply.
- `--agb-net normalized_initial_MZY`: signed AGB net = normalized gross
  ejecta minus normalized source initial elemental fraction times returned
  mass. Initial composition is matched by initial M/Z/Y; duplicate records
  must have identical abundances. Both normalizations are explicit source
  model choices. Other yield-table channels remain net-unavailable; their
  placeholder zeros do not mean zero nucleosynthesis. Gas deposition remains
  gross. For Fishlock, use its own same-model `X0(i)` column instead of a
  KL16 initial-composition match. Fishlock gross ejecta are not rescaled.

Example: append these three options to the previous command and use
`--low-z-agb none`. Set channel upper masses to `120,7,120,8,260`.
The common active Z hull is then .007--.01345, not a new extrapolation.
Effective SSP SNIa remains separate from a microscopic progenitor model.
The unmodified 40-Myr DTD still fails a strict source-CO-WD causality check:
the earliest selected CO-WD appears at 48.44 Myr.

Histories use sparse per-star age knots; no cross-product of every other
star's phase ages is required. Prepared native wind/SN and terminal AGB
maps validate their own complete source histories. Source coordinate
comparison is exact, to preserve real hour-separated late burning phases;
physical budget tolerances are unchanged.

`mkrun.py` terminal/GUI `comparison_ccsn` now offers `baseline`, `phase`,
`agb7`, `agb7_net`, `agb7_lowz_net`, and `agb7_pulses`. It copies the selected local inputs and selects the
corresponding tested binary; no simulation is launched by generation. The
new inputs/binaries are under `.physical-extension.7rcxv4/`. A fresh clone
does not include these local scientific assets or executables.

The `agb7_lowz_net` selection combines 77 AGB nodes (61 KL16 + 16 Fishlock),
LC18 phase winds/ordinary CCSN and the existing effective SSP SNIa. It extends
the common active Z hull to .001--.01345 for AGB1--7 and signed AGB net.
Generate it with the combined-builder options above but select
`--low-z-agb fishlock2014_raiteri96`. The Fishlock7 endpoint returns
5.855482487283754 Msun, leaving 1.144517512716246 Msun, consistent with the
1.145-Msun core printed in [Fishlock Table1](https://arxiv.org/pdf/1410.7457).
Its Raiteri96 age is 53.248752 Myr, an explicit cross-model approximation.
Do not interpret this as an extension to 8 Msun or as an ECSN prescription.
Source initial fractions sum to .9999998846994237 before normalization.
The earlier statement that Fishlock initial composition was unavailable was
incorrect: the already staged yield file includes the required `X0(i)` column.
Native GNU/Intel endpoint, mixed-Z CO-inventory, signed-net and full effective
SSP closure tests pass; existing source inputs/defaults remain byte-identical.

### Fishlock thermal-pulse wind timing

Select `--agb-release fishlock_tp_mass_loss` instead of `terminal_envelope`,
together with `--low-z-agb fishlock2014_raiteri96`. This uses the pinned
[Fishlock CDS TP structure table](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/797/44/table2.dat)
(`external/g2_candidates/fishlock2014_pulses/table2.dat`, local scientific asset).
The 772 printed rows contain two identical duplicates; the reader collapses
only those duplicates and retains 770 pulses across 16 stellar models.

The cumulative AGB wind follows initial mass minus tabulated TP total mass,
at the relative times supplied by the preceding-interpulse periods. The last
TP **left limit** is aligned to the selected Raiteri96 terminal age. Pre-TP
loss is explicitly uniform; unresolved final-envelope loss is a terminal
jump. This is not a measurement of Monash total lifetimes or new post-table
pulses. Composition remains the integrated mean and speed the selected
constant; KL16 nodes retain their terminal-envelope prescription.

Native `wind_history_terminal_remnant` separates cumulative wind from WD
formation. The history's `agb_wind_jump_*` arrays specify each source-node
terminal jump, and are bound to restart identity. No remnant is created by
the earlier wind. Gross/net/energy totals and non-CO exclusions are unchanged.
At 7 Msun/Z=.001 the last TP left limit returns 4.8883 Msun with no remnant;
the terminal event brings the total to 5.855482487283754 Msun and creates the
1.144517512716246-Msun ONe remnant (not strict SNIa CO-WD supply).

`mkrun.py` selection `agb7_pulses` packages the local `agb7-pulses` inputs and
`ramses_physical_pulses3d` binary. Its standard four-step comparison is still
Z=.01; the focused `live-pulses`/`restart-pulses` checks instead use Z=.004
to include the Fishlock/KL16 mixture, with matching initial H/metal fractions
and unchanged dust-to-metal ratio/temperature. This short live check does
not reach terminal AGB ages; the existing native source fixture tests those.

An independent solar-only 9--13-Msun comparison is also available:

```sh
python3 simulation/snrt/tools/build_sukhbold_native_lowmass.py \
  --output-dir NEW_INPUT_DIRECTORY --solar-coordinate .02 \
  --wind-speed-km-s 100 --lifetime-model raiteri96_padova \
  --remnant-policy stable_segment_residual --max-remnant-difference-msun .1
```

It uses the local Sukhbold Z9.6/W18 stable wind and terminal segments and
PHOTB final kinetic energies. The .02 solar label, Raiteri96 lifetimes,
uniform fixed-speed winds and residual remnant are explicit approximations;
they are not recovered KEPLER lifetimes or a non-solar yield grid. Selected
radioactive sidecars are not added to stable yields (overlap/incomplete
inventory). Source residual remnants differ from PHOTB values by at most
.094343 Msun here. History v3 declares [9,13]; old v1/v2 contracts remain
unchanged. Native consumer/SSP tests cover this input; it is not yet packaged
as a full RAMSES comparison by mkrun. Use single population, Chabrier id=2,
binary fraction=0, wind/SN mass range [9,13], and disable absent channels.
Do not merge it into LC18's non-solar branches or extrapolate it to 8 Msun.

Source descriptions and redistribution restrictions:
[Sukhbold archive](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/index.html).
Actual builds, endpoint/SSP and live/restart evidence, along with remaining
physical-input dependencies, are recorded in the project plan cited above.

### Optional native primary dust scattering (2026-09-08)

`mkrun.py` comparison modes offer `Primary dust scattering=isotropic_elastic`.
The default `none` leaves the absorption-only model and old sidecar identity
unchanged. The opt-in uses
`config/dust_dl01_bulk_030_scattering_reference_v4.nml`, regenerated from the
same Draine raw optics and DL01 material with the existing thermal exporter
and `--native-scattering isotropic_elastic`. Native v4/reference status and
`SNRT_ALLOW_REFERENCE_CONTROL=1` remain required. This does not grant
production/publication approval to the approximation.

For each owned cell and primary photon group, direction-integrated bins obey
`q'_d = exp(-tau_sca)*q_d + (1-exp(-tau_sca))*sum(q)*w_d/sum(w)`, with
`tau_sca = nH * relative_dust * sigma_sca * reduced_c*c*dt`.
The coefficient is `C_ext*albedo`, sampled at the same representative photon
energy as the absorption/energy ledger. The native CPU and CUDA kernels share
this cell operator. `SNRT_BACKEND=auto` assigns independent cell batches to
available shared CUDA streams or OpenMP workers; `openmp`/`cuda` force a leg.
This is primary RT, not the separate material/IR backend selection.

The update is first-order split after transport and absorption. Every
nonlinear trial starts from the same incoming state; the scattering result
is committed once with the existing RT transaction. All batches must succeed
before publishing. Each group conserves photons and energy to FP32 rounding;
scattering contributes no absorption or thermal source. Subsequent transport
uses the redistributed angular state. Actual coefficients and the model
marker are bound into the existing HDF5 dust identity. Changing the model
or coefficients on restart is rejected, including switching on/off.

Limits: isotropy is a chosen comparison approximation, **not** the measured
Draine phase function (the source also supplies anisotropy moments). No
HG/anisotropic scattering, grain/gas recoil or radiation pressure, frequency
redistribution, IR scattering, stochastic heating or sublimation is added.
The split spatial scheme is not asymptotic-preserving in the unresolved
optically-thick diffusion limit; accurate transport needs spatial/time/angular
resolution studies. See [Draine source definitions](https://www.astro.princeton.edu/~draine/dust/dustmix.html).

Local tested builds: `.snrt-cpu.OKoz9T/ramses_scatter_cpu3d` and
`.physical-extension.7rcxv4/ramses_scatter3d`. Existing native hybrid and
Fortran tests cover analytic isotropization, photon conservation, half-step
composition, zero opacity, busy streams, forced CPU/CUDA, late-error rollback
and ABI/layout. The bounded live/restart evidence is in
`.dust-scatter.CWEvh3/`; no external audit or full-domain MPI/GPU qualification
is implied by the two-rank local-operator test.

### Native gas/dust/IR thermal coupling (2026-09-08)

Comparison mode now offers `Dust gas thermal exchange=hydrogen_accommodation`.
Default `none` preserves old inputs. The selected v4 sidecars are
`dust_dl01_bulk_030_exchange_reference_v4.nml` (absorption-only primary) or
`dust_dl01_bulk_030_scattering_exchange_reference_v4.nml` (with isotropic
primary scattering). They declare geometric collision area per reference H
`3.495e-22 cm2` and thermal accommodation `0.5`. The former is the explicit
effective-sphere comparison `3*(1.398e-26 g/H)/(4*3 g/cm3*0.1 micron)`,
**not** a grain-size distribution inferred from Draine optical opacity.
Other choices require explicit `--collision-area-per-h` / `--accommodation`
with the existing exporter `--native-gas-exchange hydrogen_accommodation`.

The hydrogen-equivalent geometric law follows equations20--21 of
[McKinnon et al. (2021)](https://academic.oup.com/mnras/article/502/1/1344/6067372):
`K = (nH*relative_dust)*area_per_H * nH * sqrt(8*kB*Tgas/(pi*mp)) * 2*kB*alpha`.
Gas loses `K*(Tgas-Tdust)` to grains per physical volume/time. This is an
explicit approximation, not an electron/ion Coulomb or molecular collision
network; gas composition/Cv comes from the existing accepted H/He state.
Accommodation and effective collision area are model inputs, not newly
measured properties or production approval of the adopted dust mixture.

The implementation here retains finite tabulated grain U(T); it does not
assume instantaneous dust equilibrium. Each existing IR substep solves gas,
material and radiative emission **jointly**, with gas Cv/composition fixed
during that substep. Thermal speed follows the implicit final gas temperature:
`K(Tg)=K0*sqrt(Tg/Tg0)`. For gas energy Eg, eliminating the gas equation gives
`Q = r*(Eg-Cg*Td)`, `r=K(Tg)*dt/(Cg+K(Tg)*dt)`, `Tg=(Eg-Q)/Cg`.
The positive root is found using a scaled cubic in sqrt(Tg) and bracketed
Newton iterations; failed roots reject the entire material trial. The
previous fixed-K mode remains available in the low-level comparison ABI,
but is no longer the selected live solver.
The same monotone material solver then solves
`Ed(Td)-Ed_old + dt*P_net(Td) = dt*H_abs + Q`.
IR reabsorption remains in the existing fixed-point solve; Q is subtracted
from the gas and included in the existing total-energy residual. Gas,
material, IR photons and their temperatures publish only on success. The
ordinary no-exchange ABI layout is unchanged; material modes2/3 add gas data
and transfer to the existing CPU/CUDA/hybrid material kernel. Mode2 is the
frozen-speed comparison and mode3 is the variable-speed live solve.

Do not apply a gas/dust kick after a whole IR step: the first trial of that
split failed in the real control because a stiff grain reservoir could not
cool radiatively during the kick. Its log and old binary are preserved, but
mkrun selects only the replacement joint solver. No temperature clamp, source
range widening, or reduction of the collision coefficient was used to pass.
The existing bath/material temperature domain remains enforced (10--300 K
for this sidecar); cold-gas solutions below the bath, sublimation, charging,
grain drift/size evolution and radiation-pressure force injection are not
supplied by this model. Coupling to primary chemistry remains operator split,
so this is not proof of time-resolution convergence in arbitrary environments.

The code runs natively; `SNRT_DUST_BACKEND=auto` uses the existing nonblocking
stream/OpenMP material and IR dispatch. Collision model, coefficients and
joint-solver marker are bound into HDF5 restart identity. Never switch them
on/off or edit their values during restart. No new main RAMSES namelist key.
Frozen-speed executables remain `.snrt-cpu.OKoz9T/ramses_exchange_coupled_cpu3d`
and `.physical-extension.7rcxv4/ramses_exchange_coupled3d`. The new mkrun
selection is `ramses_exchange_nonlinear_cpu3d` / `ramses_exchange_nonlinear3d`
in the same respective build directories. HDF5 exchange algorithm marker3
rejects marker2 checkpoints instead of silently changing the integrator.
No-exchange restart identities and sidecar physical coefficients are unchanged.
Evidence: `.dust-exchange.WKZl1j/` (old) / `.dust-nonlinear.TtLIOc/` (new)
and the existing population/source bundle record.

### Optional trapped CR comparison (2026-09-08)

In `PHYSICS_PARAMS`, `cr_enabled=.true.` selects the NENER=1 advective
gamma=4/3 CR fluid. `cr_sn_fraction` and `cr_snia_fraction` partition the
existing SNII and coupled SNIa energy budgets; both default to zero and
must lie in [0,1]. `cr_transport='advective'` is the only admitted closure.
Optional `cr_sf_support=.true.` adds CR effective pressure support to
`sf_virial=.true.` models 1,2,4; the thermal SF threshold excludes CR energy.
Diffusion, streaming, losses and CR ionization chemistry are not included.

Use `python3 mkrun.py --mode cli` (or `--mode gui`), select the comparison
run mode, then explicitly enable CR and/or dust mass. These profiles require
the documented local physical inputs and CPU binaries; a fresh Git clone
alone does not include the source datasets or compiled executables. Setup
checks their presence, writes a new run directory and never launches a job.
CR requires periodic noncosmological CPU hydro, HDF5, channel feedback,
no sinks/AGN/delayed cooling, and hllc/hll/llf. The general cosmological
wizard does not enable this reference. HDF5 binds the CR on/off state,
fractions and SF coupling, so these cannot change across restart.
For build commands, executable identities and measured scope, see
[the CR/dust evidence](../../provenance/cosmic_ray_dust_evolution_status_2026-09-08.md).

### Optional bulk dust mass evolution (2026-09-08)

`PHYSICS_PARAMS` now accepts the following explicit reference model (off
by default). `mkrun.py` comparison setup and the shared GUI namelist
generator expose the same parameters; setup never launches a calculation.

```fortran
  dust_mass_enabled=.true.
  dust_growth=.true.
  dust_sputtering=.true.
  dust_condensation=0d0,.2d0,.15d0 ! wind, AGB, SNII; no SNIa condensation
  dust_grain_radius_cm=1d-5
  dust_grain_density=3d0
  dust_sticking=.3d0
  dust_growth_max_temperature=300d0
  dust_metal_atom_mass=24d0
  dust_injection_temperature=20d0
```

These are comparison coefficients, not calibrated efficiencies or measured
element-specific depletion fractions. Dust is a subset of the advected
total-metal and mixture-density fields, not an additional mass source.
Stellar injection condenses the selected fraction of each channel's
returned mass minus H/He; its tabulated material energy at injection
temperature is taken from the available source energy after CR partition.
Existing source locking, MPI exchange and progress accounting are reused.

At each leaf level, before the SNRT step, a frozen-rate analytic update
solves `dD/dt = (A-B)D - A*D*D/Z`. Here D and Z are dust and total-metal
densities, `A=3*S*rho_Z*v/(4*rho_grain*a)` below the growth temperature cap,
`v=sqrt(8*k*T/(pi*metal_atom_mass*mp))`. Thermal sputtering uses
`B=3*3.2e-18*(rho/mp)/a / (1+(2e6/T)^2.5)` in cgs
([McKinnon et al. 2018, section 3.4](https://arxiv.org/html/1805.04521v1)).
The update preserves a zero seed and `0<=D<=Z`; total metal mass is not
changed. Material energy follows mass at fixed grain temperature for this
split update, with an equal opposite gas-energy transfer. CR and kinetic
energies are excluded from the available thermal reservoir. Latent heat
and element/size-resolved grain evolution are not modeled.

Requires SNRT/DUST_LIVE/HDF5, metal hydro and channel-resolved feedback,
active RT and a v4 material table admitting the injection temperature.
Current scope is periodic/noncosmological, no sinks/AGN and no external
metal cooling: gas-phase element depletion is not yet wired into that
cooling model. Primary H/He chemistry and the existing gas/dust/IR thermal
solver remain active. Do not describe this bulk fixed-size/composition
closure as shattering/coagulation, unresolved SN-shock destruction, dust
drift or radiation pressure. Rates and chemistry are operator split;
the reference trials below are not arbitrary-resolution convergence.

The 13-value `header/dust_mass_values` HDF5 attribute binds the entire mass
model and its coefficients. Restart requires exact identity, including
on/off state; absent attributes remain valid only with evolution off.
The tested combined CR/dust executable is
`.cosmic-ray.kyySgK/ramses_dust_mass3d` (NENER=1, CPU hydro). Existing CR
and no-mass-evolution executable files were preserved. Dust mass evolution
is a native Fortran level/source update, not a Python postprocessing step.
RT/material/IR retain their existing stream/OpenMP dispatch; no new GPU
dispatch or cosmological CR admission is claimed for the mass update.

Evidence: `.dust-mass.w1Jo9A/` and
`provenance/cosmic_ray_dust_evolution_status_2026-09-08.md`.

### Optional composition budgets and depleted scalar cooling (2026-09-08)

The ordered medium-term extension starts with an explicit comparison:

```fortran
  dust_mass_enabled=.true.
  dust_mass_model='carbon_olivine_v1' ! default remains bulk_v1
  dust_cooling='depleted_scalar'    ! default none; opt-in approximation
  cooling=.true.
  cooling_method='original'
  haardt_madau=.false.
  J21=0d0
  self_shielding=.false.
  gpu_hydro=.false.
```

Carbon and MgFeSiO4 olivine are separate passive mass densities appended
after dust material energy. All eleven stellar element fields must be active
and continue to represent gas+dust. The NENER1/virial comparison uses fields
22/23 for C/olivine; `NVAR=30` has room. Initialize aggregate dust and its
constituents consistently; the generated composition comparison starts all
three at zero, not by guessing the composition of an old bulk seed.

Condensation is computed on each physical source release segment before
age/Z/IMF mixing. AGB ejecta first reserve the smaller C/12 versus O/16 CO
inventory; only excess carbon/oxygen is eligible for dust. CO is not an
extra dust carrier or a time-dependent molecular abundance. The existing
wind/AGB/SNII efficiencies 0/.2/.15 now apply to eligible carbon and
stoichiometrically limited olivine, rather than all returned metals. Fe
inside olivine is not counted again as metallic Fe. Aggregate metal return
still uses total returned mass minus H/He, not the sum of tracked metals.

Growth/sputtering use separate element-limited capacities with the existing
fixed-radius, shared-density/rate approximation. Gas-phase elements are
total elements minus grain stoichiometry; they are not separately advected
or added to mixture density. Aggregate dust face flux and prolongation are
the sum of constituent fluxes/prolongations. SF removes the appended fields
through its existing all-passive astration path. Source energy, reverse MPI
exchange and progress transactions remain shared with feedback; the local
mass step uses one collective per level, not one per species.

`depleted_scalar` passes `(rho_Z-rho_dust)/rho/0.02` to the actual original
cooling solver. It retains that solver's equilibrium H/He/Compton thermal
closure alongside SNRT photoionization/recombination, not a single unified
non-equilibrium cooling closure. Both J21 and the legacy residual UV floor
are zero in this admitted isolated comparison. Zero rates in logarithmic
cooling tables use the existing metal-zero convention of 1e-100; no physical
UV field is reintroduced to avoid log(0).

This is **not individual-element metal cooling**: selective depletion cannot
be represented exactly by a solar-mixture curve. Physical element-rate data
and a qualified receiver remain part of bundle 1, not a completed feature.
Composition-dependent opacity/material properties, evolving sizes, CO/H2
chemistry, and their model-selection tests belong to the subsequent bundles.
Bulk remains selectable and default; its original no-cooling setup is unchanged.

The HDF5 `dust_composition_values` attribute additionally binds the model,
cooling choice and species index. A changed model is rejected even when
NVAR is unchanged; old scalar-only checkpoints are not silently converted.
The tested CPU executable is
`.cosmic-ray.kyySgK/ramses_dust_composition_zero_uv3d`; `mkrun.py` and the
shared GUI generator expose these choices. Existing stream/OpenMP RT/dust
dispatch remains unchanged; GPU hydro is not admitted for these new carriers.
Evidence and limits: `provenance/dust_composition_medium_term_plan_2026-09-08.md`.

### Individual-element CIE comparison (2026-09-08 continuation)

With `dust_mass_model='carbon_olivine_v1'`, select `dust_cooling='wss09_cie'`
to replace the original cooling update with an embedded, author-published
[Wiersma/Schaye/Smith CIE table](https://local.strw.leidenuniv.nl/WSS08/).
Keep `cooling=.true., cooling_method='original'` for this explicitly selected
dispatch, and the same no-HM/J21/self-shielding admission above. The mass
closure and `none`/`depleted_scalar` comparison options are retained.
`mkrun.py`/GUI select `.cosmic-ray.kyySgK/ramses_dust_cie3d` for this option.

The native routine subtracts grain stoichiometry from actual element fields
and uses gas H/He and nine individual gas-phase metals. It does not add the
old scalar-Z contribution on top. H/He net rates and their CIE electron
density interpolate in the published He/H grid. Metal weights are number
abundances relative to the table's solar values, WSS09 equation (3). The
ASCII table omits the solar electron reference for equation (4), so that
correction is **not** reconstructed from an unrelated H/He-only column.
Untracked metal mass is not relabelled Fe or renormalized into tracked metals.

The thermal coordinate is T/mu, with CIE H/He particles/electrons and one
particle per gas-phase metal nucleus; metal-supplied electrons are neglected
under the stated trace-metal approximation. Physical gas H, not a hardcoded
primordial H fraction, sets the density conversion in this new update. The
mass/temperature constants match the active RAMSES `units()` convention.
A piecewise-linear signed net rate is integrated analytically in T/mu at
fixed density/composition. Equilibrium zeros cannot be crossed by a large
cooling timestep. This replaces the original thermal rate, not SNRT's
photoheating; optically thin radiative losses are not deposited into dust IR.

Published support is T=100--9.5907e8 K and nHe/nH=0.0786528--0.106898;
outside support, including evolution leaving it, the comparison rejects
the step instead of extrapolating, clamping or switching secretly to scalar Z.
It is low-density CIE, **not LTE** and **not a local-radiation/NEQ metal
network**. Metal-electron effects, dense/molecular gas, super-solar abundance
accuracy and consistency with time-dependent SNRT ion fractions are not
qualified by this implementation. In particular it is a baseline for bundle
3, not the final thermochemistry for an AGN radiation field or cooling shocks.

Source `simulation/snrt/data/wss09_cie_z_collis.txt` has 352 rows/25 columns;
the numerical import is in `patch/lagRamses/dust_wss09_cie_data.inc`. The
source's signed element contributions are retained. Two small HDF5 attributes
`dust_cie_hhe_values` and `dust_cie_metal_values` bind the full compiled table
and conversion constants, sharing the existing restart-identity collective.
Legacy `cooling_*.out` remains the original table diagnostic, not a dump of
the active WSS09 curve; the CIE checkpoint attributes contain that table.

### Two-size carbon/olivine comparison (2026-09-08)

The native `dust_mass_model='carbon_olivine_2size_v1'` option adds four
transported grain masses, not just diagnostic size fractions. It requires
the same all-element, channel-feedback, CPU-hydro, noncosmological HDF5/SNRT
v4 contract as composition mode. `mkrun.py` and the shared namelist/GUI
generator expose it together. Bulk/none remain the defaults.

```fortran
  dust_mass_enabled=.true.
  dust_mass_model='carbon_olivine_2size_v1'
  dust_size_radius_cm=5d-7,1d-5
  dust_size_density=2.2d0,3.3d0
  dust_small_injection_fraction=0d0,0d0
  dust_coagulation=.true.
  dust_shattering=.true.
  dust_sn_shocks=.false.
```

Size/material definitions and collision reference coefficients are motivated
by [Dubois et al. 2024](https://arxiv.org/html/2402.18515v2), sections 3.2,
3.5 and 3.6. This is **not a full reproduction** of that model. Existing
wind/AGB/SNII condensation coefficients and no-Ia-condensation policy remain.
Each composition's injected mass is split before the existing source commit;
default injection is entirely in large grains. Total metals include dust.

With the NENER1/virial comparison layout:

| Native field | Meaning |
| --- | --- |
| 20, 21 | Total dust mass density, mixed-reference material energy |
| 22, 23 | Carbon and MgFeSiO4 total dust densities |
| 24, 25 | Carbon small and large mass densities |
| 26, 27 | Silicate small and large mass densities |
| 28 | Optional transient coupled SN energy density |
| 29, 30 | Optional transient freshly injected C/silicate dust densities |

Fields 22/23 are sums of their two size bins; field 20 is their sum. Face
fluxes and AMR prolongation derive those dependent sums from bins, rather
than independently limiting each aggregate. SF uses the existing passive
mass-removal path. Initial masses must already obey the sums/element bounds;
there is no implicit conversion of scalar seeds or old checkpoints.

Growth shares each composition's finite element reservoir between both bins.
Radii and material densities scale geometric accretion; the effective
accreting atom mass remains the configured 24-mp reference. Thermal erosion
retains the shared Tsai-Mathews fit scaled by 1/radius, not the Hu species
fits. Symmetric analytic reservoir subflows and exact quadratic-donor
collision transfers preserve nonnegative masses and element limits. The
growth/erosion rate resolution is 0.1 per substep, capped at 4096 substeps;
an unsupported stiff interval is rejected, not silently frozen or accepted.
No automatic hydro-timestep retry is claimed.

Coagulation uses actual local nH>=1000 cm-3 and T<10000 K, cloud fraction
one, small-grain dispersion 0.1 km/s. Shattering acts below that density.
There is no hidden unresolved-cloud density replacement, turbulent-PDF
growth multiplier, grain charge enhancement, ice mantle, or ISM CO network.
Thermal growth temperature excludes CR/kinetic energy; adding CR pressure
does not secretly heat grains or supply a turbulent Mach number.

`dust_sn_shocks=.true.` additionally enables an **energy-equivalent ambient
SN comparison**, not a measured SN event-count model. It uses coupled
SNII+thermal-Ia energy after CR partition, canonical 1e51 erg events and
6800 Msun swept mass, with size/material-dependent destruction. Wind/AGN
energy is excluded. The three transient source rows use the existing cell
lock, progress commit and reverse-MPI exchange. On the owning cell, after
hydro/cooling and before growth/RT/SF, expose only pre-existing dust and
protect all same-step fresh ejecta. Clear all three rows in the same staged
mass update; no extra MPI collective. They must be zero in uniform ICs and
at completed checkpoints. Dust energy lost/gained in the combined mass
update is transferred oppositely to gas thermal energy.

This option is off by default: overlap with resolved shock sputtering,
event clustering, physical-event counting and spatial/timestep dependence
are not calibrated. Fresh-ejecta survival depends on the declared timestep
source split; it is not a separately resolved reverse-shock calculation.
The 16-value `dust_size_values` checkpoint attribute binds model parameters
and field indices using the existing restart check. Changed size, injection,
collision or shock settings cannot silently resume the same checkpoint.

With `dust_optics_model='fixed_mix'`, size/composition-dependent optics are
not used. The optional composition-dependent internal energy and collision
area below replace only the material/geometry part of that fixed reference.
The D03 option described at the end of this document now connects local
optics as well. WSS09 CIE remains a comparison, not radiation-dependent
NEQ cooling. The 8/16-bin accuracy/cost reference remains a model-selection
task; four evolving masses do not establish multibin equivalence.

### Local composition material and grain collision area (2026-09-09)

`dust_material_model='dl01_composition_v1'` is an explicit two-size comparison
option; default `fixed_mix` is unchanged. It uses local carbon/silicate mass
fractions to construct U(T), retaining the bulk graphite and olivine-like
Debye-mode assumptions of `build_dl01_dust_material.py`. Both sizes of a
composition use its same bulk specific energy. All four bins share one
temperature: no PAH/finite-size vibrational modes, stochastic heating,
separate grain temperatures, latent heat or sublimation.

The compiled 162-node 5--300 K curves are generated by the existing material
builder via `--fortran-composition-output`, stored in
`patch/lagRamses/dust_dl01_composition_data.inc`, and resampled linearly in
log T onto the active v4 temperature grid. Out-of-range grids are rejected.
No molecular/atomic NLTE rate is supplied by these solid heat-capacity data.

Native wiring uses the same material definition at each receiver:

- stellar source energy = sum of the two source dust masses times their
  U(T_injection), charged to the existing source energy;
- mass evolution decodes the old common T using old composition, evaluates
  new-composition U at that T, and exchanges the energy difference with gas;
- IR's joint gas/dust solve receives each cell's U(T), and uses it for
  material bounds and all temperature/energy conversions;
- hydrogen-equivalent gas accommodation uses actual geometric area per
  volume, sum_bin(3*rho_bin/(4*solid_density*radius)), with the existing
  accommodation coefficient. This is not an optical cross section.

The per-cell U input is supported by the shared OpenMP/CUDA material cell
operator and existing free-stream-or-CPU hybrid scheduler. Transport/emission
opacity and the reference-H normalization remain unchanged; do not interpret
this as composition-dependent optical radiative transfer. Native hydro for
this comparison is still CPU-only. CUDA material kernel support does not
authorize GPU hydro or demonstrate full-GPU hydro integration.

`mkrun.py`/GUI expose the option and select
`.cosmic-ray.kyySgK/ramses_dust_composition_material3d`. Starting with zero
dust avoids inventing a composition/energy conversion for old scalar ICs.
The HDF5 `dust_material_composition` attribute binds model version and all
compiled U/T values using the existing restart validation. Enabling,
disabling or changing this material closure on restart is rejected rather
than reinterpreting saved Ed. The old fixed-material checkpoint path remains.

Optical source limitation found while preparing the next connection:
the actual [Draine Gra_81 / smoothed-silicate tables](https://www.astro.princeton.edu/~draine/dust/dust.diel.html)
contain 81 radii and 241 wavelengths spanning 0.001--1000 microns. The
high-energy end is about 1.24 keV, below the current 2--10 keV group; their
long-wavelength limit also excludes part of the existing IR grid. They were
inspected, not substituted for the existing WD01/D03 whole-mixture table.
Do not zero missing groups, repeat endpoint values or silently extrapolate.
Full-band source-backed size optics and the RT/CIE/NEQ connection remain open.

### Prepared D03 composition/size optics (2026-09-09)

The narrow Gra_81/suvSil_81 data are not the only available source. Actual
D03 dielectric files extend to 18--19 keV before their anomalous endpoint
rows, despite the old web-page 2 keV description. They cover all nine
source groups and 136 current IR nodes. The offline builder
`tools/build_d03_grain_optics.py` computes sphere Mie coefficients from five
cached originals in `data/draine_d03/`, retaining Qabs, Qsca and g. The full
byte identities and approximations are in `data/dust_d03_optics_generation_v1.json`.

Native `dust_composition_optics.f90` supplies four optical bases and local
mass-weighted absorption/scattering/g, with an explicit physical unit
conversion and spectral/size/density checks. This D03 comparison requires
0.01/0.1 micron radii and graphite/silicate densities 2.2/3.8 g/cm3; it must
not silently replace the current 0.005/0.1 micron, 2.2/3.3 evolution model.
See [D03 material assumptions](https://arxiv.org/html/astro-ph/0308251).

Four pure components and one mixture passed the existing native IR
initializer/advance with identical absorption/emission opacity and matching
DL01 U(T), including CPU, hybrid and forced CUDA. Energy conservation and
stationary-background tests passed; reference/backend differences <=6.09e-16.
The existing offline optical test entry `tests/draine_dust_opacity.py --d03`
also passes Mie amplitudes, Rayleigh limits, full spectral coverage and angular
weighting. Generated include and manifest rebuild identically from originals.

**Historical preparation status, superseded by the live option below.**
At this point the binary linked the new module, but the live driver retained
fixed-mixture optics. No unconnected namelist/mkrun option was advertised.
The subsequent implementation needed to switch primary
absorption/scattering and cell-dependent IR absorption/emission together,
bind the optical identity on restart, and preserve the existing GPU/CPU lease
policy. Hard-X g is nearly one: retain the angular information rather than
misrepresenting the full Qsca as accurate isotropic transport. No stochastic
heating, separate grain temperatures or RT-dependent metal NEQ is claimed.

### Live D03 size/composition optics (2026-09-09 closeout)

The implemented explicit comparison selection is:

```fortran
! Within the existing dust-enabled channel-feedback comparison:
dust_mass_model='carbon_olivine_2size_v1'
dust_material_model='dl01_composition_v1'
dust_optics_model='d03_transport_v1'
dust_size_radius_cm=1d-6,1d-5
dust_size_density=2.2d0,3.8d0
```

Default `dust_optics_model='fixed_mix'` is unchanged. D03 requires the v4
material contract with scattering enabled and exact spectral binding; use
`dust_dl01_bulk_030_scattering_reference_v4.nml`, or
`dust_dl01_bulk_030_scattering_exchange_reference_v4.nml` for gas exchange.
The contract retains its reference-H normalization, bath and spectral grid;
native D03 replaces its fixed optical coefficients at all receivers.
`mkrun.py` (including GUI) sets these together and selects the tested CPU
NENER=1/NVAR=30 `.cosmic-ray.kyySgK/ramses_dust_d03_live3d`. Existing comparison
guards, including noncosmo/no sinks/CPU hydro, remain in force.

Local four-bin mass fractions now feed primary absorption and scattering,
IR absorption/emission and IR angular scattering. Absorption is distinct
from scattering and equals the Kirchhoff emission opacity. Scattering uses
the labelled delta-isotropic transport approximation Qsca*(1-g) for both
primary and IR, preserving photon number/weighted IR energy and the first
angular moment, **not** a resolved Mie phase function or full angular X-ray
halo. IR scattering is FP64 and uses quadrature-weighted bins, not an
unweighted sum of intensities. It is first-order operator split within the
existing transport substeps, with no new MPI synchronization layer.

Four shared optical/emission bases and four weights/cell extend the existing
joint material/IR receiver. No per-cell T-by-frequency emissivity cube or
simulation-time Python is needed. Shared OpenMP/CUDA kernels and the existing
free CUDA stream -> GPU / busy stream -> CPU policy handle these fields.
Live CPU MPI/restart and forced/automatic CUDA kernel tests passed; this does
not establish multi-GPU hydro or cosmological production qualification.

Restart attribute `dust_optics_d03` binds 1900 doubles, version 2 (primary
and IR transport scattering), including all source data and physical sizes.
Switching optics on/off or changing it on restart is rejected. Earlier
fixed-mixture checkpoints still use their original unchanged identity.
The existing local DL01 U(T) identity remains separate and unchanged.

Evidence and limitations are recorded in
[dust implementation closeout](../../provenance/dust_composition_medium_term_plan_2026-09-08.md#bundle-3-live-optical-connection--driver-closeout-2026-09-09):
1031-cell heterogeneous native CPU/GPU parity/energy/rollback tests;
MPI2 x OMP2 continuous/restart equality (90 hydro and 4 RT datasets);
positive dust/element/thermal states; changed-model restart rejection; and
mkrun/GUI consistency. Fixed-vs-D03 sensitivity is measured, not calibrated.

This closes the local optical wiring, not full many-bin dynamics validation
or radiation-dependent metal NEQ/H2/CO chemistry. The first native bounded
8/16/32-bin size-shift reference is now implemented in `dust_mass_physics`;
its mass/area/number comparison and explicit limiter errors are recorded in
the provenance document above. It is not a selectable live multibin model.
WSS09 CIE,
common grain T, 20 K dielectric and uncalibrated size/SN prescriptions remain
explicit approximations. PAH stochastic heating, sublimation, separate Fe,
dust drift/AGN force and cosmological deployment are not supplied by this
option. No new default activation or blanket publication-readiness claim.
