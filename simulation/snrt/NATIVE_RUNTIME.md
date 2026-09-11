# Native runtime controls (implementation, not physical approval)

## Explicit effective SSP population (2026-09-11)

`population_model='effective_ssp'`, `binary_fraction=0`, `use_snia=.true.`
selects ordinary mixed-v5 histories plus a full-initial-mass empirical DTD.
This is population ID2, not a resolved binary population or a zero Ia rate.
The SED wrapper must bind the actual IMF and population ID2, while the v5
ordinary history retains its original single-star input identity. It requires
all four ordinary channels, their2--600Msun source windows, full.08--600 IMF
normalization and `user_selected_model_v1`/`source_consistent` settings.
The DTD sidecar must use `mass_accounting='effective_ssp'` and a baked-in
binary-incidence rate. Strict-WD and frozen-contact modes are not admitted
for ID2. Existing ID0/ID1 behavior and default Chabrier selection are unchanged.

`config/fp2_snia_effective_population_runtime_v1.nml` supplies the existing
Kroupa Maoz/N100 reference under ID2; it is not a Chabrier rate conversion.
Explicitly select matching material/SED IMF identities. The old zero-binary
single-star setup does not become Ia-enabled merely by setting use_snia.
[Native and actual MPI2/OMP2 evidence](../../provenance/effective_population_execution_2026-09-11.md).

## Explicit approximation comparisons (2026-09-11)

Defaults are unchanged; comparison admission is not a physical-production
approval. Select a matched source package as a unit, not independently mixed
yield/history/SED files.

- `dust_pah_model='pah_h2_catalytic_v1'`: existing H0--13 neutral/cation
  states plus H13+H -> H12+H2, shared finite atomic-H donor and explicit
  local thermal accommodation. Same3584 PAH carriers; CHIMES/fixed groups,
  no Fe or drift. [Bounded native/live/restart PASS](../../provenance/pah_catalytic_comparison_implementation_2026-09-11.md).
- `parsec_mixed_lowmass_truncated_v1`: version5 material/history/SED
  package, sources2--600Msun with the full.08--600 IMF denominator. Different
  evolutionary/terminal grids and omitted post-track light are explicit;
  no SNIa or binary population. [Actual-source and MPI2/OMP2 restart PASS](../../provenance/parsec_mixed_lowmass_evidence_2026-09-11.md).
- Converter `phase_escape_f22_bistability_v1`: wind-velocity comparison,
  not a new mass-loss/LBV history. [Native binding and source evidence](../../provenance/parsec_bistability_comparison_2026-09-11.md).
- `radioactive_model='lc18_al26_fe60_transparent_v1'` with
  `radioactive_companion_path`: matched prompt-projected LC18 only,
  `RADIOACTIVE=1`, two additional tail fields. Al26/Fe60 are subsets of
  existing gas metals, not extra mass. Transparent MeV photons, no decay
  heating; noncosmo CPU/OpenMP MUSCL, no dust/CHIMES/MHD/SGS/GPU hydro.
  [Native and MPI2/OMP2 live/restart evidence](../../provenance/stellar_radioactive_implementation_2026-09-11.md).
- Frozen-He-contact SNIa replaces empirical event counts only under
  explicitly effective SSP accounting and exact birth Z=.01. It is not a
  disjoint microscopic binary population. [Contract and native/MPI restart evidence](../../provenance/snia_frozen_he_hybrid_plan_2026-09-11.md).

- `dust_iron_model='fe_uv_cycle_v1'`: metallic-Fe UV stationary-charge
  comparison through13.6eV, CHIMES electron/ion cycles and secondary-electron
  partition; trace-charge/energy and relaxation-domain checks remain strict.
- `dust_iron_model='fe_thermal_limit_v1'`: full thermal-retention upper
  comparison through10000eV, NOT a photoelectric/Auger cascade model. Both
  Fe choices retain the300K grain-temperature domain and exclude PAH/drift.
  [Physics, native and MPI2/OMP2 live evidence](../../provenance/fe_photon_comparisons_plan_2026-09-11.md).
- `dust_pah_model='pah_atomization_limit_v1'`: catalytic H2 plus single-photon
  complete C24 atomization into H/C/C+, with unit energy-allowed yield and
  local excess-energy retention. No daughter-grain network or higher charges.
  Only explicitly admitted fixed monochromatic groups are supported; the
  current hard group is869.634149eV, not its entire broadband interval.
  Atomic source condensation releases the same derived binding energy used
  by destruction; selecting the model does not heat pre-existing material.
  [Physics, native and actual comparative MPI/restart evidence](../../provenance/pah_single_photon_atomization_plan_2026-09-11.md).

Fe photon comparisons pass bounded three-step MPI2/OMP2 live evaluation;
this is not a Fe restart or full hard-photon cascade qualification. PAH
atomization passes the weak/strong mono-group8 live comparison and7591-dataset
restart match, with measurable carbon return distinct from stellar uptake.
Neither result is a universal broadband-source or calibrated-galaxy approval.
Frontends currently pass57 tests (one display-dependent Tk skip); native
configuration passes48 assertions (`.frontend-final.NGzCu0`).
See the [execution record](../../provenance/remaining_physics_bundle_execution_2026-09-11.md)
for current scope, source limitations and completed raw-output cleanup.

## PARSEC low-Z printed-precision source (2026-09-10)

Explicit feedback conversion with `--metallicity-grid precision_eleven
--wind-km-s 1000` supplies495 nonrotating14--600Msun nodes at
Z=1e-11,.0001,.002,.004,.006,.008,.01,.014,.017,.02,.03. Rebuild the SED
from that package and select all three existing yield/history/SED paths.
The native binder requires the new `parsec2025_w17_hw02_precision_rate_v1`
identity on both sides. No new namelist key or default change.

This model preserves gross elemental yields, reconstructs baryonic masses
only within actual printed intervals, and uses positive tabulated RATE for
wind timing with an explicit fixed-speed comparison. Phase-speed extreme-Z
extrapolation is not admitted. Omitted author Z=1e-6/.001 branches are
replaced only by the declared neighboring-Z cumulative mixture, not a claim
to their original fates. Z=0 and <14Msun remain unsupported; this is not a
full population. [Details and evaluation status](../../provenance/parsec_low_z_precision_implementation_2026-09-10.md).

SNRT Intel builds now require `-no-ftz` on the Fortran MAIN (Makefile
enforces this). Otherwise subnormal FP32 photon increments can disappear
while paired FP64 energy remains, corrupting weak hard-band means. This is
gradual underflow within the existing representation, not an unlimited
dynamic-range claim or a relaxed spectral admissibility test.

## PARSEC five-metallicity common source (2026-09-10)

Add `--metallicity-grid metal_rich_five` to the PARSEC feedback converter
to select actual Z=.008/.014/.017/.02/.03 (225 nonrotating14--600Msun nodes).
Default `solar_pair` and its old two-Z physical payloads are unchanged.
Rebuild the matched SED from that feedback package; its converter derives
the exact same Z list. Native material/radiation share the linear-Z bracket
and relative exact-node tolerance and require actual source binding before
publication. Select the existing yield/history/SED paths together; no new
RAMSES namelist key or full-SSP claim. Other Z with unresolved source budgets
are not silently admitted. [Source, review, native tests and MPI2 exact-restart evidence](../../provenance/parsec_metallicity_extension_implementation_2026-09-10.md).

## LC18 prompt elemental decay (2026-09-10)

The LC18 converter accepts explicit
`--decay prompt_t12_le_100yr_baryonic_v1`; the combined KL16/LC18 converter
accepts `--massive-decay` with the same value. Both retain the unchanged
`as_tabulated_no_decay` default. This completes short nuclear chains at
stellar release and stops at half-life>100yr; Ni56/Co56 reach tracked Fe,
but Al26/Fe60 remain in their parent elements. It is not a finite-horizon
or live ISM decay calculation, and adds no decay heating or photons.
Total return, remnant, ages, energies and AGB sources stay unchanged.
Use the ordinary `PHASE0_YIELD_TABLE` / `high_mass_history_path` connection;
no new NML controls or online Python. Physical source identity binds the
selected projection. [Prescription, review disposition and native evidence](../../provenance/lc18_prompt_decay_implementation_2026-09-10.md).

## PARSEC phase-dependent wind energy (2026-09-10)

`tools/build_parsec_pair_feedback.py --wind-model phase_escape_f22_v1`
generates an opt-in cumulative phase-energy table for the existing PARSEC
v4 population. Do not supply `--wind-km-s` for that named model; default
`fixed` still requires an explicit speed and preserves the old source bytes.
Select the output through existing `PHASE0_YIELD_TABLE` and
`high_mass_history_path`, and rebuild its matched `SNRT_STELLAR_SED` using
`tools/build_parsec_native_sed.py`. Actual feedback and SED model IDs must
agree; restart checks all consumed energy rows and photon moments.

Wind energy thermalizes through the existing native source bridge, with no
second radial kick or SN-class CR partition. Cool10km/s, hot H-poor d1.6
and other hot d2.6 branches are explicitly simplified escape-speed closures,
not a resolved atmosphere/LBV/bi-stability model. Mass/element endpoints and
explosion energy stay unchanged. Both PARSEC converters now use the source's
Lsun3.846e33, correcting the initial SED converter's nominal-unit assumption.
No RAMSES namelist key, carrier or default changed.
[Prescription, invocation, source-unit correction and native/live evidence](../../provenance/parsec_phase_wind_implementation_2026-09-10.md).

## Matched PARSEC high-mass stellar radiation (2026-09-10)

Existing `SNRT_STELLAR_SED` also accepts explicit external source v4 from
`tools/build_parsec_native_sed.py`. This is **HIGH-MASS ONLY**, nonrotating
14--600Msun atZ=.008/.014 by default (five-Z extension above), with the same native IMF cells/full .08--600
denominator and terminal ages as the PARSEC five-fate v4 feedback source.
It is not a full SSP, default stellar source, or original NLTE spectrum.
The ordinary low-mass/AGB radiation contribution is not supplied.

The converter uses actual five cumulative-above-edge Q constraints with
a within-band Planck prior and sub-LW bolometric residual. Missing late
photon times use actual full-track L/Teff Planck spectra, explicitly NOT
measured Q. Native code reads positive interval Q/E, integrates the same
piecewise-linear cumulative law without subtracting large endpoints, and
binds all consumed moments/weights in MPI and HDF5 restart identity. No
per-cell source history or online Python is used. No post-terminal photons.

Requires active channel-resolved PARSEC v4 wind/CCSN/pair feedback with
matching source windows and a supported individual-star IMF. Startup loads
and checks the actual prepared feedback table before admitting radiation.
Requires an existing band-spectrum mode and explicit reference-control
opt-in. It cannot be silently mixed with an independent BPASS source.
Old v1/v2/v3 sources and all defaults are unchanged. See the
[implementation, scientific limits and measured tests](../../provenance/parsec_common_radiation_implementation_2026-09-10.md).

## Optional H/He band spectrum (2026-09-10)

`SNRT_SPECTRAL_MODEL=fixed` (also unset) retains fixed source-weighted group
coefficients. `hhe_maxent64_v1` enables native intragroup maximum-entropy
N/E reconstruction, Verner H/He spectral absorption, survivor hardening and
actual accepted photoelectron energy into chemistry. Requires `SNRT=1`,
`DUST_LIVE=0`, `CHIMES=0`; use the separate D03 option below for live grains.
`SNRT_BACKEND=auto|openmp` uses OpenMP; forced CUDA is unsupported. The selector
must match across MPI ranks and restart; no implicit fixed/spectral migration.
64 quadrature nodes are re-derived, not stored in each cell. Source N/E
injection is consistently FP32-quantized with a bounded total energy rounding
check. This does not reconstruct a unique SED; v1/v2 stellar sources use
nominal energies while v3/v4 supply actual source Q/E. This original mode's
secondary ionization retains a mean-energy closure; use the FS2010 option
for node-resolved secondary deposition. CPU cost is appreciably higher, so
the default is unchanged.

### D03 grain spectral comparison

`SNRT_SPECTRAL_MODEL=hhe_d03_maxent128_fs2010_v1` selects128-node positive
N/E reconstruction, Verner H/He and node-resolved FS2010, with actual D03
graphite/silicate absorption and elastic delta-isotropic scattering. It
requires `SNRT=1 DUST_LIVE=1 CHIMES=0`, `dust_optics_model='d03_transport_v1'`,
the DL01 two-size material and a v4 scattering-enabled dust contract. Radii
are .01/.1 micron and densities2.2/3.8 g/cm3, as for the existing D03 mode.
Fe, PAH, sublimation and relative dust motion reject. CPU/OpenMP only;
automatic backend selection stays on CPU, forced CUDA rejects.

Four actual grain mass columns compete with H/He at every spectral node.
Rejected atom-limited captures retain their own E/N; unlike the fixed
operator they are not transferred through a second dust sink. Captured
grain energy goes directly to the transactional material/IR solve, not a
nominal group mean. Elastic angular relaxation of the reconstructed survivor
spectrum follows absorption; this split does not resolve multiple internal
scatter/absorb histories or apply a grain force. No grain photoelectric
escape is added. Old fixed/64-node behavior and namelist defaults remain.

The original D03 representative-energy/IR tables are unchanged. The separate
`dust_d03_band_data.inc` and `data/draine_d03/band128_manifest.json` bind the
source grid, actual opacity and transport-scattering coefficients. Native
checkpoint version10 and the disjoint HDF5 format include the nodal content
SHA256. Nodes are temporary quadrature, not additional cell carriers. Source
preparation uses `tools/build_d03_grain_optics.py --band-nodes 128`; runtime
uses no Python. Fixed20K dielectric data and1/3--2/3 graphite remain physical
approximations;128 nodes do not resolve every near-edge feature.
[Physics, finite-inventory convention, numerical bounds and MPI2 restart evidence](../../provenance/snrt_band_implementation_2026-09-10.md).

### Static Fe spectral comparison

`SNRT_SPECTRAL_MODEL=hhe_d03_fe_maxent128_fs2010_v1` adds two actual Fe
columns to the preceding four-bin operator. Requires the explicit
`dust_iron_model='fe_electric_compare_v1'`, CPU/OpenMP, CHIMES=0 and a fresh
`SNRT=1 DUST_LIVE=1 DUST_IRON=1` build. The NENER1/virial hydro profile is
NVAR32: Fe31:32 follows the reserved dust window. With CHIMES=1 the old
Fe188:189/NVAR189 profile is unchanged, but spectral modes remain rejected.

Without CHIMES, only static seeds are admitted: `dust_cooling='none'`,
zero `dust_condensation` and `dust_fe_condensation`, growth/sputtering/
coagulation/shattering/SN shocks/Fe kinetics all false, Fe sticking zero,
no relative motion. Both namelist generators implement these prerequisites.
The generated setup selects no incompatible hard stellar SED by default.

This is NOT general stellar/AGN Fe heating. The existing 300 K and4 eV
limits stay: nodal absorption above4 eV in an Fe-bearing cell rejects the
entire trial. A broad1--5.6 eV packet may violate this even with mean2.4 eV;
no spectral tail is removed. The completed live test uses an explicitly
synthetic sub-eV source, not a physical SED claimed to pass these limits.
Electric/eddy mu=1 optics omit spin absorption and photoelectron escape;
six components share the existing material temperature. No new CUDA path.

`tools/build_fe_grain_optics.py --band-nodes 128` creates the separate
`dust_fe_band_data.inc`; `data/dust_fe_band128_manifest.json` binds sources.
Both Fe/D03 node hashes bind native version11 and HDF format40, with no
additional per-cell spectral carriers. Old modes/identities are unchanged.
[Bounded implementation, Fe flux repair and exact MPI2 restart evidence](../../provenance/snrt_fe_spectrum_implementation_2026-09-10.md).

## Angular resolution (2026-09-10)

Build `SNRT=1 SNRT_ANGULAR_LEVEL=0|1|2` for respectively 80, 320, or 720
discrete-ordinates directions (8x10, 16x20, 24x30 Gauss-Legendre/product
rules). Default0 preserves the old nodes, weights and ordering exactly.
Use a separate clean build directory when changing this flag: Make does not
track compiler-option changes. This is not a namelist or runtime switch;
existing mkrun profiles still use their existing binaries.

Directional storage and work scale by roughly 1:4:9. Source deposition,
transport, scattering, IR and checkpoint dimensions share the selected count.
Restart rejects a different angular resolution; there is no implicit remap.
The [bounded implementation evidence](../../provenance/snrt_angular_implementation_2026-09-10.md)
includes 320-ray MPI2 hydro+feedback+dust/restart and 320/720-ray native CUDA
transport. It does not establish angular convergence for every source/mesh,
full GPU evolution, or new intragroup spectral physics.

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
A native collision-aware reference now also implements same-material
coagulation and fragmentation, with explicit noninteracting mass outside
the radius grid. The [group-6 comparison](../../provenance/dust_collision_comparison_2026-09-10.md)
finds two-size/multibin differences but retains the live two-size default;
the collision reference does not add NML options or live hydro carriers.
WSS09 CIE,
common grain T, 20 K dielectric and uncalibrated size/SN prescriptions remain
explicit approximations. PAH stochastic heating, sublimation, separate Fe,
dust drift/AGN force and cosmological deployment are not supplied by this
option. No new default activation or blanket publication-readiness claim.

### Coupled atomic cooling with the dust comparison (2026-09-09)

`dust_cooling='snrt_hhe_cie_metals'` is an opt-in native alternative to
`wss09_cie`, not a rename of that model. It requires composition dust and
the existing active SNRT/v4/HDF5/noncosmo comparison. mkrun and its GUI select
`.cosmic-ray.kyySgK/ramses_dust_atomic3d` (SNRT/DUST_LIVE, NENER=1, NVAR=30,
CPU hydro). No simulation-time Python or new hydro fields are introduced.

H/He collisional ionization, case-B recombination and thermal losses use
the actual SNRT ion fractions, with gas-phase element H/He densities. The
same number inventories feed absorption and the optical-depth predictor;
the same ion-dependent heat capacity feeds temperature, dust evolution and
gas/dust exchange. The photo receiver defers recombination to the atomic
receiver; cooling_fine does not apply a duplicate CIE thermal sink. Atomic
evolution also runs in cells without absorbed photons or with dust-only
absorption. Kinetic and CR energy are excluded from its thermal budget.
Photoheating remains its own nonnegative ledger; escaped collisional/
recombination cooling is not claimed to remain in the stored RT field.

Fits follow [Rosdahl et al. 2013, Appendix E](https://academic.oup.com/mnras/article/436/3/2188/1247446).
A positive implicit H/He substep with frozen temperature/electron density
is controlled by step doubling (1e-3 local error, 1e-5 fraction scale floor,
10% relative state-change bound). It charges collisional ionization with
the same rate and threshold as the reactions. Final publication is atomic
on success; the existing level transaction handles failure. This is
operator splitting, not a fully implicit joint radiation/metal network.

Only metal columns of [WSS09](https://arxiv.org/abs/0807.3748) enter this
closure, weighted by depleted individual abundances (eq. 3). They retain
CIE ion populations; there is no invented eq. 4 solar-electron correction.
Metal electrons are not added to the H/He heat capacity. The atomic domain
is 1--1e9 K; nonzero metals also require the tabulated 100--9.5907e8 K range.
No extrapolation, H2/CO, molecular cooling or full metal NEQ is supplied.
The old CIE option and its He-abundance domain remain unchanged.

Restart stores a distinct closure code, the actual WSS09 tables and the
18-value `dust_atomic_cooling` identity. Changing closures on restart is
rejected. Native one-zone and coupled MPI2 x OMP2 restart results are in
the dust provenance document. Neither short-test convergence nor this
H/He correction qualifies cosmological production or galaxy calibration.

### Hot bulk material extension (2026-09-09, not sublimation)

`dl01_composition_v1` preserves the existing 5--300 K table exactly. For
300--3000 K it now evaluates the same normalized DL01 bulk vibrational mode
integrals in native Fortran (16-point quadrature). It does not extrapolate
the last tabulated slope. This is the harmonic solid's internal energy,
not a claim of phase stability, latent heat, sublimation, or finite-size PAH
physics. The actual runtime domain remains the intersection with the
selected material/IR contract; the default generated run stays at 300 K.

A matching hot contract can be built with the existing material and thermal
builders using `--temperature-max-k 3000`. Their combined knots must remain
within the native 256-node bound, and the optical spectral grid must still
match. All primary/IR/emission receivers and the existing OpenMP/CUDA cell
layout continue to use the same per-cell curve. No runtime Python is added.
Cold material checkpoint identity is unchanged; hot contracts use version 2
with curves sampled through 3000 K. Material/IR contract changes on restart
remain forbidden. `mkrun.py` documentation reflects this distinction, but
does not silently select hot contracts or old binaries for them.

The existing native dust test passes independent 256-node integral values,
300 K continuity, mixture and domain checks. A 1000 K D03 material/IR trial
conserves energy and agrees between Fortran and OpenMP to 1.03e-16 relative.
The hot-contract MPI2/OMP2 live/restart pair agrees bitwise on all 565
hydro/SNRT datasets and its material identity. Optical functions remain
frozen at 20 K; this is not temperature-dependent optics, a new GPU test or
sublimation validation. See the remaining-physics provenance for run paths.
Source: [DL01 mode model](https://arxiv.org/html/astro-ph/0011318).

### Native radiation-dependent chemistry comparison (2026-09-09)

`dust_cooling='chimes_neq_v1'` is a separate opt-in 157-species CHIMES
receiver, including H2/CO, molecular cooling and all ionization stages of
eleven elements. It requires two-size DL01 composition dust, the existing
noncosmological periodic CPU-hydro profile, and gamma=5/3. This does not
promote the comparison to a calibrated production model.

Build a pinned, double-precision CHIMES library with
`data/chimes_native_receiver.patch` (adapter ABI5), SUNDIALS 5.8.0 and HDF5.
ABI4 libraries remain supported for the existing neutral-solid live receiver;
the new charge-aware cell interface requires ABI5. Build the new library in
a separate directory: do not overwrite a library needed by retained binaries.
From a new build directory immediately below the project root:

```sh
make -f ../bin/Makefile -j1 SNRT=1 DUST_LIVE=1 NENER=1 HDF5=1 \
  CHIMES=1 USE_FFTW=0 EXEC=ramses_chimes \
  CHIMES_DIR=/absolute/patched/chimes \
  SUNDIALS_DIR=/absolute/sundials-install
```

This profile reserves NVAR=187 and defaults to NVECTOR=32. Retaining the
legacy 500-grid batch with 157 extra species exceeds a 512 MiB OpenMP
worker stack in the current MUSCL implementation. The output binary is
`ramses_chimes3d`; set `SNRT_CHIMES_BINARY` to its absolute path for mkrun.
Set `SNRT_CHIMES_MAIN_DATA` to the pinned main table and
`SNRT_CHIMES_GROUP_DIR` to the directory containing all nine group files.
Use `tools/build_chimes_group_tables.py` and patched official chimes-tools;
do not substitute broad-spectrum UV background tables for these bins.
Source versions, hashes, build artifacts and execution evidence are in
`provenance/dust_remaining_implementation_2026-09-09.md` at the project root.

The actual SNRT photon densities and exact group edges/means feed chemistry.
The old H/He receiver and generic CIE cooling are bypassed, avoiding duplicate
absorption and thermal evolution. Species are transported nuclear-number
carriers, not extra mass. The chemical/dust face fluxes share a common
normalization; total elemental fluxes follow their nuclei. Grain growth
reserves molecular C/O and other bound elements. Feedback/destruction adds
neutral atoms while preserving existing molecules. HDF5 stores all species
and binds the tables and closure identity; changing them on restart fails.

Limits remain explicit: mean-matched within-bin spectra, local shielding,
thermal-width cell columns for molecular line escape, translational EOS,
no additional CR ionization inferred from CR energy, and no transported
fluorescent/Auger-cascade emission. Photoelectric gas heating is disabled
until its energy can be removed from the existing dust absorption ledger.
The existing pinned FS2010 table now partitions primary photoelectron
energy inside the CHIMES chemical RHS. The corresponding nonthermal share
is subtracted from primary heating, not added on top of it. Atomic base
and Auger-producing photoionization channels use the same admitted shell
energies for both budgets; secondary electrons consume no extra primary
photons. Atomic target shortages continuously limit the ionization and HI
excitation fractions, returning the unavailable share to heat. This is a
declared extension of a primordial atomic-gas table, not a molecular or
arbitrary-composition electron-degradation model. Excitation photons escape
in this comparison. The chemical restart identity is now version 4. It
includes up to two retries of nucleus/charge-failing cells from their
original input, with solver tolerances tightened by 100 at each retry;
the conservation acceptance tolerances are unchanged. Version-2 chemistry
checkpoints require the preserved ABI3 binary/library; version-3 checkpoints
require their preserved pre-retry binary. The neutral-solid receiver accepts
ABI4/ABI5; its existing restart identity is unchanged. ABI5 adds a private
per-cell gas+solid charge constraint for the new local PAH receiver, not an
automatic change of the live `pah_neutral_absolute_v1` model. The latter still
rejects its unsupported photon/charge domain. Native component evidence and
live wiring evidence are in
[`medium_physics_implementation_2026-09-10.md`](../../provenance/medium_physics_implementation_2026-09-10.md).
The absorbed-photon diagnostic alone is not a radiation/material energy
closure test. Metallic Fe grains, stochastic PAH emission, sublimation and
relative dust dynamics are not implemented by this chemistry selection.

### CHIMES atomic spectral building block (2026-09-10)

`snrt_chimes_spectrum.cpp` supplies an immutable native 128-node bank for
311 photoionization/Auger reactions and 682 partial shells, using the same
pinned Verner95/96 and KM93 inputs as the grey model. Fortran bindings are
in `snrt_chimes`; the CHIMES Makefile profile links the implementation.
Each direction's photon number/energy reconstructs separately. The returned
moments are `sum N*sigma`, `sum N*sigma*E`, and
`sum N*sigma*(E-shell binding)` per reaction and band, before multiplication
by reduced light speed and target density. They are coefficients, NOT
accepted finite-inventory captures or net gas heat. A shared loaded handle
can be read by parallel cell calls; load/free remain outside those calls.

Build the separate input with `tools/build_chimes_group_tables.py
--chimes-tools <patched-tools> --main-data <pinned-main.hdf5>
--band-nodes 128 --output <new-directory>`. Native loading requires 128
nodes; 64/256 builder choices are data-comparison outputs, not live modes.
Setting `SNRT_CHIMES_BAND_TABLE=<new-directory>/atomic_shells.h5` adds its
checks to the existing `snrt_thermochemistry_smoke`; the table path alone is
NOT a RAMSES runtime selector. The explicit hot-only selector is below.

The moment-provider API alone is not an evolving finite-inventory receiver
and its primary-electron moment must not be treated as net heat. See
[implementation and native evidence](../../provenance/snrt_chimes_spectrum_implementation_2026-09-10.md).

### Conservative atomic photo operator and hot-cell split (2026-09-11)

`chimes_band_photo_step` now evolves CHIMES157 atomic photoionization/Auger,
shell/node FS2010 secondaries, and finite directional photon N/E in one
native CVODE solve. Each occupied energy node has one survival fraction
shared by its rays; rays are reconstructed independently at entry. The
operator returns heat, three secondary-ionization energies, excitation,
unresolved binding/cascade energy, absorbed energy and primary photon count.
It does not itself include collisional chemistry, recombination, molecular
dissociation, dust or a gas-temperature update. Failures publish nothing.

`chimes_cell_band_hot_atomic` connects that result to actual nonradiative
CHIMES with zero photons, updating heat using the changed particle count.
It is a first-order split with an explicitly selected hot-only live adapter
(below). Its ninth ledger entry is signed dark-step gas thermal change, NOT
transported cooling radiation. It requires dust-free atomic gas with zero
solid charge and T>1e5K (the CHIMES hot-network criterion). Molecular inputs
and a photo/dark endpoint crossing that threshold reject the staged step;
H2+ is never silently erased to make the split work.

The existing CHIMES build links this operator and SUNDIALS SPGMR; no external
CHIMES library replacement, changed defaults or new HDF5 input is needed.
The same `SNRT_CHIMES_BAND_TABLE` extends the existing native smoke. Tests
cover photon exhaustion, H-/H2+, metal/Auger/secondary budgets, hot native
chemistry/cooling, dt convergence and rollback. See
[implementation, review disposition and evidence](../../provenance/snrt_chimes_photo_implementation_2026-09-11.md).
Cold molecular spectra, competing grain absorption and full escaping-energy
closure remain outside this mode; it does not admit a general CHIMES
spectral simulation.

### Hot atomic spectral live wiring (2026-09-11)

Set `SNRT_SPECTRAL_MODEL=chimes_hot_atomic_maxent128_fs2010_v1` and
`SNRT_CHIMES_BAND_TABLE` to the pinned 128-node atomic shell bank
(`998970ed5bb4cfeda01913a72cc3fc62af8ce2fdb55dc924f2704e3abb5391de`).
Keep the existing `dust_cooling='chimes_neq_v1'`, `cooling=.true.`,
CHIMES main/group table environment and CHIMES=1/DUST_LIVE=1 build.
The dust layout is infrastructure only: all dust carriers must be zero,
two-size/DL01 material layout selected, condensation fractions zero,
growth/sputtering/coagulation/shattering/SN shocks off, and no sublimation,
Fe, PAH or relative-motion model. Gas must remain atomic and above 1e5 K
at operator boundaries. Invalid cells reject the collective staged update.
Neither cold input nor threshold crossing is silently converted to a grey
model. The default remains the existing grey receiver.

Paired N/E transport has zero gas absorption; this receiver owns all atomic
captures. It receives per-direction physical photon number and actual
energy, performs the photo + dark split once, and commits species, gas
energy and outgoing N/E together through the existing RT transaction.
FP64 energy corrections are rebased to outgoing FP32 photon counts.
The absorbed-energy diagnostic uses actual shell absorption, not the grey
reference group mean. Excitation/unresolved cascade budgets are not added
to gas heat or invented as an advected chemical reservoir; escaping
cooling/cascade radiation remains unresolved as in the native split.

Restart preserves the same directional state width but binds new model
versions (native 12, HDF5 base version +40), the atomic bank SHA256 and
chemical identity version5. Changing model or bank rejects the restart;
fixed/old spectral formats retain their previous interpretation.
No new namelist field is introduced: the generic generator/GUI explains
the explicit environment selector. `mkrun.py` cold/dusty comparison bundles
explicitly export `SNRT_SPECTRAL_MODEL=fixed`, preventing shell inheritance,
unless the deliberate cold comparison opt-in below is supplied.

### Cold molecular / competing D03 spectral live wiring (2026-09-11)

`SNRT_SPECTRAL_MODEL=chimes_cold_d03_maxent128_fs2010_v1` selects the
explicit cold molecular comparison. Keep the hot-mode atomic bank above
and set `SNRT_CHIMES_MOLECULAR_TABLE` to the molecular bank with SHA256
`c313497b68633b5f2909eb02aa14160baf98dbdde1784b73f5dd50f3d9fdad0b`.
CHIMES main/nine-group data, FS2010 tables and D03 optics are still required.
Build with SNRT=1, CHIMES=1, DUST_LIVE=1 and HDF5=1; the tested hydro/CR
profile uses NENER=1, NVAR=187, MPI2 and OMP2. Use `SNRT_BACKEND=openmp`;
this spectral transport has no CUDA implementation or implicit model fallback.

Select `dust_cooling='chimes_neq_v1'`, `cooling=.true.`, gamma=5/3,
`dust_mass_model='carbon_olivine_2size_v1'`, DL01 composition material and
`dust_optics_model='d03_transport_v1'`. Grains may have nonzero masses but
are co-advected and fixed against local mass processing: condensation all
zero, growth/sputtering/SN shocks/coagulation/shattering off, no Fe, PAH,
relative motion or sublimation. Gas must remain in **10--10^4.98 K**
(approximately 95,499 K upper limit), including internal dark thermal
trials; require reduced-c times dt no greater than the cell length.
An invalid state rejects the collective transaction, not its molecules.

Ordering is transport/scattering without absorption, joint atomic/molecular/
grain photo absorption, dark CHIMES chemistry/cooling, then dust material/IR.
Gas and grains compete for the same finite directional N/E. Accepted grain
group energies reach material/IR once; chemistry is not repeated after that
receiver. This is an explicit first-order split, not a temperature/opacity
fixed-point iteration. Shielding is frozen locally; the coarse H2/CO line
projection, translational heat capacity and unresolved fluorescence remain
declared approximations. No new photoelectric gas heating is added.

The existing directional state width is unchanged. Native restart version13
and chemical identity6 bind this model; HDF5 uses its base version +50
(56 for the tested IR+stellar profile). D03, atomic and molecular bank
identities are checked before restore. Both dark and irradiated MPI2/OMP2
restarts reproduced all 402 datasets bitwise in the bounded live test.

For `mkrun.py` comparison generation, opt in using
`SNRT_CHIMES_SPECTRAL_MODEL=chimes_cold_d03_maxent128_fs2010_v1` plus the
two spectral table paths and the existing CHIMES binary/data settings.
The generator exports the runtime selector, disables the mass-processing
flags with an explicit notice, and records the domain in the generated
README. The generic namelist generator/GUI also describes the selector;
no new namelist variable was added. See the
[implementation and integrated evidence](../../provenance/snrt_chimes_molecular_coupling_implementation_2026-09-11.md#live-completion).

### Energy-aware high-temperature transition (2026-09-11)

`SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1` selects
the operator-approved rapid-dissociation comparison (kind7). Use the same
atomic/molecular banks, FS2010 data, D03 grains and build options as
the cold comparison above, but rebuild CHIMES with receiver **ABI6** from
`data/chimes_native_receiver.patch`. Existing ABI4/5 libraries remain usable
for their original models; they cannot run kind7. The build must include
the matching library, not merely set a new environment selector.

Gas is admitted on **10--1e9 K**, not an unbounded/relativistic domain.
At entry, after photo heating, or at an upward dark-solve temperature root
at the loaded molecular ceiling (10^4.98 K), remaining molecules are
converted to atomic fragments. ATcT v1.130 0 K formation enthalpies specify
each of the twenty dissociation costs. Nuclei and charge are conserved,
electrons unchanged, and the binding cost is debited from thermal energy
before temperature is recomputed using translational gamma=5/3 capacity.
For a +1 molecule the charge goes to its lowest-ionization-cost H/C/O atom.
This fragment prescription is a declared approximation, not branching data.

ABI6 locates dark crossings using CVODE's accepted solution and an upward
root, with the requested integration endpoint as stop time. Only numerical
RHS trials beyond the boundary use boundary rates; this does not extend
the physical molecular trajectory above the data ceiling. The remainder
of that step uses atomic NEQ, even if it cools below the ceiling; molecular
formation can resume next step. The operator split is first order. This
is not a finite-time high-temperature molecular shock network.

Joint gas/grain absorption and subsequent material/IR wiring remain. Kind7
also admits existing C/silicate `dust_growth`, `dust_sputtering`,
`dust_coagulation` and `dust_shattering` switches. The mass step runs before
RT and reconstructs gas carriers, grain material energy and current opacity.
Molecular heavy nuclei are reserved for accretion; atomic ions deplete in
their existing proportions and electrons reconcile to charge. Erosion
returns neutral atoms. This retains the existing comparison convention,
not a new charge-dependent sticking/recombination heat model. Size exchange
uses TOTAL H nuclei, including H2; this density correction also applies to
the pre-existing grey CHIMES mass path. C/olivine do not accrete hydrogen.
No new cell carriers, material binding-energy model or numerical ledger.

The subsequent [process connection](../../provenance/kind7_dust_process_bundle_plan_2026-09-11.md)
admits existing source condensation, unresolved SN destruction, sublimation
and C/silicate relative dynamics in kind7. Fe/PAH remain excluded. Process
questions default off; fixed grains and the existing effective erosion laws
remain comparison assumptions, not general dust survival predictions.
CPU/OpenMP remains required. Kind5/kind6 restrictions stay unchanged.

With relative motion, competitive node-resolved photo absorption accumulates
four phase energies and directional moments on the accepted trajectory.
Absorption impulse uses physical c; midpoint mechanical work is subtracted
once from solid heating. Cold chemistry reads the staged phase-aware kinetic
energy. Stationary primary scattering is disabled; after photo, reconstructed
nodes supply photon-weighted group-grey phase transport opacity to the existing
nine-group work-conserving moving-scatter solver. This is not node-resolved
moving scattering or cross-group Doppler transport. Split sublimation is
allowed with drift, but coupled-IR sublimation remains coadvected-only.
Grain opacity is frozen during photo: subsequent material/IR sublimation
updates the next call's opacity, not an implicit photo/sublimation solution.
Transient photo counters add no checkpoint carriers.
Relative C/silicate material also uses the existing DL01 low-temperature
continuation and per-band Planck continuation for absolute IR, not the
net-bath energy floor. Phase enthalpy reconstruction and mass exchange use
that same mixture energy. Below the sublimation table, a cold no-op is
accepted only when the warmer first-knot evaluation gives exactly unchanged
represented masses. No bath heating or finite-rate extrapolation is added.
The separate CHIMES grain-temperature input bounds still apply; this is not
permission to extrapolate its reaction tables to zero temperature.

Photo integration uses accumulated optical depth with survival `exp(-tau)`:
the old survival-fraction ODE could produce tiny negative terminal values
despite CVODE constraints. Nuclear/charge/photon budget acceptance thresholds
are unchanged; local solver accuracy is tightened for the nonlinear survival
and independently integrated absorption counters. At cold/transition live
storage only, subnormal FP32 counts and their FP64 energies round together
at unchanged E/N. This numerical energy change must be no greater than
64 FP64 epsilons of the cell's incoming primary energy, otherwise it rejects.
No physical abundance floor or heat source is introduced. The rounded tail
is storage roundoff, not claimed gas/grain absorption. Normal FP32 counts
retain the existing FP64 energy-rebasing convention.

Native restart version14, chemical identity7, and HDF5 base version +58
(64 for the tested IR+stellar profile) distinguish this model. The ATcT
data identity is additionally bound to SHA256
`0bfa808eda56b7ea41b1a9e083e5379106caa2960f4b8a242f6e2e902606a045`.
Existing checkpoint species/carrier widths are unchanged; older model
checkpoints must not be relabelled as kind7.

For wizard/GUI generation use
`SNRT_CHIMES_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1`
and the existing binary/data settings. Generated environment and README
record this opt-in and its restrictions. No new namelist field is needed.
Existing dust mass/size restart identities bind the selected process flags;
changing them across restart still rejects. See the
[grain mass connection and evaluation](../../provenance/snrt_transition_dust_evolution_2026-09-11.md).
See [implementation and live evaluation](../../provenance/snrt_hot_transition_plan_2026-09-11.md#approved-implementation-and-live-evaluation).

Kind7 additionally admits the explicit existing-sink, non-MAD Bondi
`partition_reference_v1` comparison with **NENER=0 hydro** and active SNRT.
Use `sink`, `sink_AGN`, `agn`, `bondi`, `create_sinks=.false.` and
`accretion_scheme='bondi'`; retain the C/silicate coadvection restrictions
above. Optional condensation is admitted; CR, MHD, new sink formation,
SN shocks, sublimation and relative motion remain excluded with sinks.
The checked AGN scalar map transports all157
chemical densities plus eight dust mass/solid-energy descriptors through
gross accretion and jet loading/deposition. Solid thermal energy is separate
from gas mechanical energy; swallowed solid energy leaves the modeled grid.
The setup wizard exposes this default-off option using an explicit
`SNRT_CHIMES_SINK_BINARY`. It changes neither the default stellar population
nor the physical approval status of its reference AGN model. For multilevel
live IR, leave `SNRT_RT_LEVEL` unset. See the
[actual-source connection and its evaluation status](../../provenance/snrt_actual_source_transition_2026-09-11.md).
The subsequent [multilevel AMR repair](../../provenance/snrt_multilevel_amr_fix_2026-09-11.md)
allows empty leaf batches without bypassing MPI collectives and uses the
current grid's parent-based ownership for halo faces. Fully covered and
mixed/derefining MPI2/OMP2 cases pass; this is not a general subcycling,
cosmological or GPU convergence claim.
The new wizard profile retains uniform levelmin3 but allows levelmax4 with
refinement explicitly disabled: the sink cloud radius4*dx_min is thenL/4,
strictly inside the periodic half-box boundary. Do not interpret changing
this allowed level as a mesh-convergence experiment; it also sets minimum
stellar particle mass and AGN support scales.

### Fixed-H charged PAH comparison (2026-09-10)

`dust_pah_model='pah_charge_fixed_h_v1'` is a separate opt-in model, not a
relaxation of `pah_neutral_absolute_v1` or general PAH survival admission.
Build with `DUST_PAH=1 DUST_PAH_CHARGE=1 CHIMES=1 DUST_LIVE=1 NENER=1`:
hydro NVAR=443, with 128 neutral plus 128 monocation excitation mass states.
Use receiver ABI5 and both `SNRT_PAH_NEUTRAL_TABLE` and `SNRT_PAH_ION_TABLE`
original optical tables. The wizard additionally requires an explicit
`SNRT_DUST_PAH_BINARY` and `SNRT_DUST_PAH_CONTRACT`.

The carriers co-advect with gas and bind HDF5 restart to both charge optics,
the energy grid and source convention. Cations carry their ionization
energy through mass times (vibrational energy + 7.02 eV); injection is neutral.
Primary and IR absorption use the charge-resolved opacity. Already-debited
primary captures retain the full-step initial charge mixture through IR
subcycling. Gas photoelectron/recombination heat and electron count commit
together; heat capacity tracks electron count. CHIMES charge closure includes
solid cations, converting molecule number to its own carrier-mass convention.
The gas chemistry and PAH charge solves are first-order split, not a fully
implicit unified network. Empirical CHIMES grain terms exclude these explicit
PAHs to avoid duplicating their charging contribution.

This fixed-H C24H12 model admits photons <=13.6 eV and gas 10--10000 K in
PAH-occupied cells. It rejects cosmology, Fe and relative dust motion. It
does **not** model H loss/addition, carbon-skeleton destruction, anions or
dications; the two charge states share the DL01 vibrational-mode prescription.
These are declared comparison approximations, not demonstrated PDR survival
or publication/production qualification. The original neutral option and its
4 eV bound remain available unchanged.

### Graphite sublimation comparison (2026-09-09)

`PHYSICS_PARAMS dust_sublimation='gd89_graphite_bulk_v1'` explicitly opts in;
the default is `none`. It requires the two-size carbon/olivine model and
`dl01_composition_v1`. Silicate evaporation is **not** enabled. The graphite
vacuum bulk law follows [Waxman & Draine, section 3.1, equation 10](https://arxiv.org/html/astro-ph/9909020):
`da/dt = -nu*(m_C/rho_s)^(1/3)*exp(-81200/Tdust)`, with `nu=2e14 s^-1`.
The two fixed-radius bins use `lambda=3*abs(da/dt)/a_bin` and backward Euler
mass loss; this is not a resolved shrinking-grain size distribution.

The native `dust_sublimation_step` solves the remaining grain temperature
and masses together using the same tabulated `U(log T)` interpolation as
the IR receiver. Lost graphite takes `L=k*81200/m_C` per gram as latent
energy and `2*k*T/m_C` as an explicit effusive-atom kinetic-energy closure;
the latter thermalizes immediately in gas. Total elemental carbon stays
unchanged; the native chemistry reconciler returns evaporated carbon to
neutral gas. Silicate masses remain unchanged in this operator.

The energy convention is `Egas + Edust + L*(C_total-C_solid)`. The phase
term is derived from existing conservative mass carriers, not an additional
passive scalar. The same binding reference applies to graphite growth,
sputtering and SN destruction when this option is on. Incoming ejecta carry
phase energy according to their incoming gas/solid carbon masses; this is
not added SN thermal energy. HDF5 `dust_sublimation_values` binds the option
and constants and rejects mismatched restart selections.

The operation is split before RT, not jointly integrated with photon
absorption. A common grain temperature, frozen D03 dielectric functions,
vacuum evaporation, bulk graphite and fixed bin radii remain explicit
approximations. Vapor backpressure, finite-size PAH dissociation, silicate
evaporation and radiation-heated time-convergence validation remain absent.
Do not interpret the 3000 K material endpoint as a phase-stability limit.

`mkrun.py` offers this selection and requires existing files supplied through
`SNRT_DUST_SUBLIMATION_BINARY` and `SNRT_DUST_SUBLIMATION_CONTRACT`; the run
environment exports that contract as `SNRT_DUST_CONTRACT`. It does not
silently reuse a binary without the new native implementation. Generator
validation and the existing GUI tests cover the same option.
Live verification status is recorded in
`provenance/dust_remaining_implementation_2026-09-09.md`.

### Explicit olivine sublimation extension (2026-09-09)

`dust_sublimation='gd89_xu25_olivine_v1'` retains the graphite operator and
also evaporates MgFeSiO4 from the two silicate bins. The graphite-only
mode and `none` are unchanged. Use the same explicit binary/contract
environment variables above, with a binary built for this extension.

The Xu table-2 crystalline face rates are averaged over equal exposed
areas. They are **not** used as latent heats. A separate ideal
forsterite/fayalite formation-enthalpy reference, neutral atomic gas EOS
and DL01 sensible energy define `LS=2.21423908e11 erg/g`. Congruent atomic
Mg/Fe/Si/O vapor returns to the existing element/CHIMES fields; gas receives
the explicitly modeled effusive heat. Ordinary silicate growth/destruction
uses the same phase reference. The extra HDF5 `dust_olivine_phase` attribute
binds the 15 source/model constants and rejects incompatible restarts.

This is a crystalline-rate/ideal-phase comparison, not validation of
amorphous optical properties, melting, vapor chemistry or back-pressure.
Radiation and evaporation are still operator-split. Source links, assumptions,
the closed native energy test and the MPI2/OMP2 fresh/restart evidence are
recorded in the remaining-implementation provenance. All 565 final physical
datasets matched bitwise in that test; this does not qualify the whole bundle.

`dust_drag.f90` additionally contains neutral Epstein stopping times and a
conservative implicit gas/multibin-dust momentum/heat primitive. It is not
yet wired into persistent live relative-velocity transport and has no
namelist selector. Do not interpret its existing native smoke test as
completion of relative dynamics. Separate Fe and stochastic PAH remain in
the approved, unfinished bundle.

The unconnected `dust_stochastic.f90` primitive now provides DL01 generic
PAH vibrational modes and canonical U/Cv, thermal-continuous stationary
populations, and an energy-moment-preserving representative-photon adapter
with explicit above-grid photon/energy accounting. These are native routines,
not a renamed equilibrium-temperature curve. They remain linked only into
the existing dust smoke: actual PAH optical inputs, independent mass/source
carriers, radiation/transport/restart and namelist activation are unfinished.
Neither a normalized probability nor a solved truncated matrix proves
spectral coverage or validity of a stationary approximation.

### Radiation-coupled sublimation (2026-09-09)

`dust_sublimation='gd89_xu25_olivine_rt_v1'` moves graphite/olivine
evaporation into the live IR material receiver. It requires the CHIMES+D03
profile, `SNRT_RT_ENABLE=1`, a version-4 hot material contract with gas
exchange enabled, and an explicitly selected new binary/contract through
the mkrun environment variables above. The earlier split selections and
the default `none` retain their behavior. The ordinary pre-RT mass operator
does not also apply sublimation in this selection.

At fixed radiation/heating coefficients, the cell solves grain sensible
energy, BE mass loss, latent energy and gas accommodation/effusive vapor
heat together. The root variable is net emitted power above the background;
this avoids subtracting nearly equal Planck powers at long timesteps.
The IR transaction accounts for latent energy without storing it as sensible
grain heat. The accepted chemical trial receives the neutral evaporated
atoms. Mass, energy and chemistry are committed only after all-rank acceptance.

A single stiff endpoint can miss an initially hot evaporation pulse even
when energy closes. The current receiver therefore uses native adaptive BE
step doubling with relative tolerance `1e-4`, a 25% state-change pre-limit
and at most 4096 trials. It accepts the two positive half-step states and
returns the time-averaged IR spectrum; it does not extrapolate mass or
renormalize emission to repair a failed energy check. The error estimate
checks sensible energy, individual bin masses, gas transfer and emitted
band energy. A failed local integration rolls back rather than committing
part of the requested time interval. Material/emission tables are prepared
once per cell call, not regenerated for every short local interval.

Opacity and geometric coefficients stay at the start of each outer IR
substep; gas collision speed follows the evolving gas energy while gas Cv
stays fixed. This is not a fully implicit opacity/chemistry solver or a
substitute for outer-timestep/mesh convergence. The same common-temperature,
fixed-radius, vacuum/crystalline/ideal-phase and frozen-optical-data assumptions
listed above still apply. No metallic Fe or stochastic PAH is activated.

The new nonlinear material integration is native CPU/OpenMP. `auto` can
still dispatch primary/IR transport and scattering through their existing
hybrid paths; an explicit CUDA-only dust-material request is rejected rather
than silently executing the old CUDA material equation. Sublimation restart
identity version 3 separates this adaptive selection from the preserved
version-2 un-subcycled comparison. Old comparisons require their old binaries.

Only this adaptive path requests undamped IR fixed-point updates in cells
whose reabsorbed-emission fraction is below 0.25 in every band. Other cells
keep the old 0.5 damping. This reduces repeated hot-transient integrations;
it does not relax conservation/convergence tolerances or change opacities.

Both the preserved un-subcycled comparison and the current adaptive
MPI2/OMP2 fresh/restart pair match all 565 final physical datasets. The
adaptive pair has nucleus error 6.04e-16, IR balance error at most 1.02e-10
and 0.6075% initial silicate loss in the deliberately 3000 K test. Native
adaptive tests resolve the missing early evaporation and show tightening-
tolerance convergence with energy closure around 1e-11. Remaining
restrictions and artifact paths are recorded separately in
`provenance/dust_remaining_implementation_2026-09-09.md`;
the approved four-item bundle remains unfinished, not production-qualified.

### Metallic Fe material building blocks (not a live selection)

`dust_iron_material` uses the pure-Fe condensed-phase table on printed p1225
of [NIST-JANAF Fourth Edition](https://janaf.nist.gov/pdf/JANAF-FourthEd-1998-Iron.pdf),
with its source molar mass 55.847 g/mol. It stores bulk enthalpy relative to
zero kelvin, neglecting condensed pV. Linear enthalpy interpolation preserves
the alpha/gamma, gamma/delta and melting jumps at 1184, 1665 and 1809 K;
1042 K is a heat-capacity anomaly, not an extra latent jump. The inverse
returns temperature and four phase fractions, including plateaus in mixtures
of graphite, olivine and metallic Fe. Supported maximum is 3000 K; no gas
branch, supercooling or nanoparticle melting-point correction is implied.

Below 298.15 K it uses the Debye-plus-electronic bulk energy shape of
[Hensley & Draine (2017), eq. 1](https://arxiv.org/html/1611.08607), explicitly
scaled by 1.00576510 to meet JANAF's zero-to-298.15 K enthalpy. This is an
approximation, not interpolation through every low-temperature datum: energy
at 100/200 K is 6.8122%/2.3013% above the respective JANAF entries. The source
table, normalization and conventions are carried in the material identity;
that identity is not yet a live restart attribute. Hydro chemistry's existing
integer mass-number convention is not changed by this material table.

`iron_radiative_cell` supplies a native six-bin C/silicate/Fe material solve
for fixed masses and supplied optical band powers. It uses backward Euler,
net emitted power and the existing gas-collision helper. At each phase
transition it solves the enthalpy interval explicitly. Fe's latent energy is
already in returned material energy and must not be added again to the IR
ledger. This callback does not implement Fe evaporation or provide Fe optical
constants. Phase-interval tests use declared synthetic spectra; the additional
six-component IR test uses the provisional electric/eddy Fe optical base
described below, not a completed Fe optical model.

`iron_radiative_batch` provides a native OpenMP batch implementation with
explicit reference mass, six mass densities and optional gas exchange. It
only publishes material energy, temperature, phase fractions, band power and
gas transfer after **every** cell succeeds. The existing four-basis CUDA
material ABI is unchanged and must reject six bases. Generic IR tables now
retain the supplied number of optical bases instead of allocating four;
transport, absorption and scattering still use the existing native backends.

`dust_fe_electric_base_data.inc` contains a provisional causal dielectric/Mie
calculation for 10/100 nm Fe spheres, using actual Werner (2009) DFT/REELS
and Henke Fe inputs with the DH13 finite-size Drude prescription. Its nine
source groups and 136 IR samples match the D03 group contract. The offline
builder is `simulation/snrt/tools/build_fe_grain_optics.py`; pinned sources,
generation hashes and approximations are in
`simulation/snrt/data/dust_fe_electric_base_generation_v1.json`. There is no
runtime Python dependency and no rescaling of Fe into silicate opacity.

This is **electric/eddy response only**, not an admitted full Fe model.
Spin-magnetic absorption, temperature/phase-dependent optics and X-ray
photoelectron energy partition are absent. The causal composite differs
from Werner DFT complex epsilon by 11.56% median and 61.87% maximum under
the declared metric; Kramers--Kronig consistency and a 26.45-electron sum
do not establish agreement with measured grain opacity. In particular, no
single-domain magnetic sphere model is silently applied to 100 nm grains.
`fe_full_optics_admitted` remains false. The six-component radiation test
heats grains with sub-photoelectric-threshold 2.366 eV photons, not keV
photons incorrectly assigned wholly to dust heat.

The elemental partition and chemical-state adapter accept optional separate
metallic Fe. The reserved-Fe wrapper around the existing four-bin growth
operator subtracts that inventory from both elemental Fe and total available
metal; it does not invent a Fe accretion/sputtering law. No current live caller
provides a Fe carrier, and no namelist/GUI option advertises one yet. The old
four-bin paths, selectors, defaults and restart identities remain unchanged.

### Unified Fe/PAH/relative-motion kernels (initial implementation)

`dust_stochastic_evolve` advances the finite-time thermal-continuous master
equation, accepting probabilities or transported energy-bin grain number
densities. It preserves their total, including empty cells, and returns the
integrated excitation/radiation energy exchange. It does not supply a PAH
optical spectrum, photon-overflow/destruction model or live population field.

`dust_mixture_split`, `dust_mixture_drag` and `dust_mixture_radiation_kick`
provide barycentric-to-component momentum conversion, conservative internal
drag and an energy-funded radiation impulse. Relative momentum is signed,
and its kinetic energy is additional to barycentric kinetic energy; neither
may be treated as a normal positive dust-mass passive. Source photon energy
must pay the kinetic work before grain heat is assigned. No hydro flux,
pressure, field layout or restart activation follows from these routines.

The gas-element/CHIMES adapter and reserved-solid growth receiver accept
explicit PAH H/C masses alongside metallic Fe. Only PAH carbon counts
against metallicity; both atoms are withheld from gas/other grains. The live
PAH connection is described below; this earlier kernel work alone did not
activate it. Relative motion remains a separate item.

### Explicit cold Fe electric-only comparison

The operator approved a limited live comparison, not full-band physical
admission. Select `dust_iron_model='fe_electric_compare_v1'` together with
CHIMES/DL01 two-size/D03 and `dust_sublimation='none'`. Build with
`DUST_IRON=1 CHIMES=1 DUST_LIVE=1 SNRT=1 HDF5=1 NENER=1` (NVAR189 in the
standard virial profile). Fe small/large masses follow the 157 chemical
fields. They enter the same element-conserving hydro fluxes, source
transaction, IR receiver and HDF5 restart as the C/silicate carriers.

This mode uses the electric/eddy optical bank, common grain temperature
and native CPU/OpenMP material callback. It rejects T>300 K and absorbed
primary group representative energies >4 eV before transaction commit.
No hard-spectrum extrapolation or silicate substitute is permitted. The
restriction is on grey representative energies, not the resolved spectrum.
Magnetic absorption, charge/cascade physics and Fe size exchange/sublimation
are omitted. Optional Fe mass kinetics is described below. Gas collision exchange remains geometric
hydrogen accommodation; Fe is excluded from C/silicate catalytic H2 and
grain-recombination area. Fe co-advection is not relative dust motion.

`dust_fe_condensation=0` is the default. An explicit fraction [0,1] condenses
non-Ia Fe remaining after olivine into large Fe grains, partitioning thermal
energy from the existing source; separate SNIa ejecta remain gas. This is
an uncalibrated comparison parameter. Restart binds exact Fe optical and
material arrays, source fraction, limits and field offset, including absence
checks when the mode is disabled. The wizard requires
`SNRT_DUST_IRON_BINARY`/`SNRT_DUST_IRON_CONTRACT` and omits the incompatible
hard BPASS SED. Defaults and full-Fe admission remain unchanged.

As of 2026-09-10, `dust_fe_kinetics=.true.` optionally adds geometric seed
accretion with explicit `dust_fe_sticking` in [0,1], and Fe-specific thermal
sputtering (Choban 2026 Table 3 fit to Nozawa 2006). The existing
`dust_growth` and `dust_sputtering` switches apply. Both new parameters
default off/zero. The implementation uses the same Fe radii/density as the
optics, reserves the updated olivine Fe, reconciles gas Fe in CHIMES, stages
the aggregate and two Fe carriers together, and conserves donor momentum
when relative dynamics is selected. Sensible-energy remapping uses the same
analytic Fe mixture enthalpy as its native IR callback. It does not add adsorption
latent heat. The existing Fe optical temperature/photon bounds still apply.

Thermal erosion uses the source's low-Z projectile mixture and resolved
density, neglects erosion below 1e4 K and rejects gas T>1e9 K when active
with nonzero Fe. No nonthermal/Coulomb/finite-size correction is implied.
`dust_sn_shocks` is rejected with Fe kinetics, since there is no selected Fe
unresolved-shock efficiency. Kinetic parameters and exact sputtering law
bind restart; the off profile retains the old Fe attribute exactly. This is
a bounded optional comparison, not a calibrated complete Fe dust model.
See `provenance/medium_fe_sources_2026-09-10.md` and
`provenance/medium_physics_implementation_2026-09-10.md` for sources/status.

Detailed implementation and measured evidence:
`provenance/dust_remaining_implementation_2026-09-09.md`, section
“Approved bounded comparison: live Fe connection”. Relative momentum still
needs its own live connection in the same bundle; PAH is described below.

### Live neutral PAH stochastic population (item 1)

Select `dust_pah_model='pah_neutral_absolute_v1'` with the existing
CHIMES/DL01/two-size/D03 noncosmological comparison and no sublimation.
Build `DUST_PAH=1 CHIMES=1 DUST_LIVE=1 SNRT=1 HDF5=1 NENER=1`:
NVAR315, or NVAR317 when `DUST_IRON=1`. VPATH is unchanged. Forced CUDA
material is rejected; CPU material is supported alongside existing SNRT
transport/scattering dispatch. This is not a charging/destruction model
or arbitrary-spectrum/cosmological production qualification.

`idust_pah:idust_pah+127` stores **molecular mass density per excitation
state**, not number fractions. It follows the CHIMES species and optional
two Fe fields (188:315 without Fe, 190:317 with Fe in this profile).
The live mass convention matches CHIMES: C24H12 has mass `300*m_p`, with
H/C fractions 12/300 and 288/300. Grain optics and modes use Nc=24, Nh=12.
The stand-alone atomic-mass inventory helper is not the live mass convention.
`idust` and `idust_energy` remain C/silicate/Fe only. PAH mass is already
included in total rho/H/C, and excitation energy is derived from these
128 masses and the common level grid; never add its mass a second time.

One shared IR field transports absolute energy, initially zero in this
isolated comparison. There is **no untracked thermal bath subtraction**;
the older modes retain their old excess-above-bath convention. Mixed
absorption uses the sum of bulk and PAH physical opacities and partitions
each captured spectral node once. The primary angular photon debit is also
single, including PAH-only cells. The finite-state backward-Euler solve
conserves PAH number and excitation/absorbed/emitted energy; overflow rejects
the whole trial. The bulk receiver, gas collision exchange and PAH population
share the existing transactional IR/halo/reflux machinery. PAH scattering,
surface H2/recombination and gas accommodation are not invented from bulk area.

`SNRT_PAH_NEUTRAL_TABLE` selects the original uncompressed Draine neutral
table. Cross sections are projected onto the existing 136-node IR grid.
Beyond its 1000-micron boundary an explicitly anchored E^2 Rayleigh/Drude
continuation is used, not claimed measured data. Unsupported occupied
high-energy IR and absorbed primary representative energies >4 eV reject;
the original table's X-ray range does not authorize neutral-grain heating
without electron-loss/charge/survival physics. Spectral features remain
limited by this IR quadrature. The absolute PAH receiver now extends bulk
enthalpy below the first 5 K knot using the DL01 cold Debye limit and the
existing Fe electronic/Debye law, with per-band Planck emission. It does
not impose a temperature floor or an untracked CMB bath. Upper bounds
(including <=300 K with Fe) remain. Material mass exchange and relative
advection use the same inverse; Fe material and PAH restart identities are
version 2, rejecting incompatible old physics. No new namelist flag.
See the [cold integration and remaining group-7 scope](../../provenance/medium_physics_implementation_2026-09-10.md).

`dust_pah_condensation=0` is the default; an explicit [0,1] fraction of
non-Ia carbon remaining **after graphite** forms PAHs, limited by hydrogen
from the same ejecta. This is an uncalibrated source parameter, not an
inferred stellar PAH yield. The existing injection temperature determines
mean vibrational energy; adjacent energy-state populations preserve that
mean and molecular mass. The same source transaction pays this energy from
gas, including its existing CR partition; separately coupled Ia stays gas.

Hydro flux carrier normalization, AMR passive transfer, stellar mass removal,
gas-phase CHIMES H/C and reserved-solid grain growth include the PAH carriers.
They co-advect with gas; this does not implement item-2 relative dust motion.
HDF5 stores every carrier and binds actual optical coefficients, level grid,
source convention/fraction and field index. A disabled or changed PAH model
cannot silently read a PAH restart as the old IR convention.

Both `mkrun.py` CLI/GUI and the namelist editor expose these selectors.
The wizard requires `SNRT_DUST_PAH_BINARY`, `SNRT_DUST_PAH_CONTRACT`, and
`SNRT_PAH_NEUTRAL_TABLE`; with Fe+PAH it needs one combined binary, not two.
It removes the incompatible hard BPASS SED and writes the model restrictions
into the run README. Native evidence is recorded in
`provenance/dust_remaining_implementation_2026-09-09.md`.

### Relative-motion spatial operator (runtime connection still pending)

`dust_multifluid` supplies a conservative gas + pressureless-grain spatial
operator and its transactional transport/implicit-drag composition. Total
rho/p/E and individual grain rho/p are conservative; barycentric relative
momenta are derived, not advected as positive scalars. Gas pressure removes
all component kinetic energy. CR advection and pressure work, gas-owned
chemical tracers and grain-owned energy/population carriers are supported.
Grain thermal/excitation energies remain outside mixture hydrodynamic E.

The first-order Rusanov face and periodic driver are checked using seven
components and 128 PAH carriers, including opposite flows with zero net mass
flux. Stiff local drag stability does not establish asymptotic-preserving
spatial accuracy or multi-stream dust physics. This is not a selectable
RAMSES relative-motion model yet: it must replace the co-advection face
interpretation together with pressure/CFL/source/chemistry/radiation impulse
and AMR/MPI/restart wiring. No dormant namelist switch has been added.

The primary transport now has an optional `dust_moment(leaf,group,3)`
output: the first angular moment of accepted dust photon **number**, after
the H/He inventory cap and returned-photon subtraction. Existing scalar C
entry points remain valid. Moment-enabled OpenMP, CUDA and stream-lease
hybrid entry points use the same cell partition; the prepared RAMSES
transport accumulates all substeps. Ordinary callers that omit the output
do not allocate or compute it. Photon bins already include angular weights.

`dust_fv_absorption_kick` consumes this moment (transposed to `(3,group)`
for one cell), accepted photon counts, group energies in erg, grain
absorption shares and the number/energy/momentum code-unit scales. It uses
physical c, not reduced c. Hydro E receives mechanical work only; the
remaining energy is returned for the separate solid/PAH heat receiver.
This absorption-only adapter does not implement scattering/IR recoil or
photoelectron energy partition. Its C/Fortran receiver and actual hybrid
GPU/CPU paths are tested, but the hydro driver does **not yet request or
apply** this optional moment. This is not live relative-motion admission.

### Optional PAH hydrogen-state comparison

`dust_pah_model='pah_hydrogen_m13_dl01_v1'` extends the fixed-H charge
comparison with H=0..13, neutral/monocation and 128 excitation states:
3584 independently advected mass carriers. Build with `DUST_PAH_H=1
DUST_PAH_CHARGE=1 DUST_PAH=1 CHIMES=1 DUST_LIVE=1 SNRT=1 HDF5=1 NENER=1`
(hydro NVAR=3771). Both original neutral and ion optical tables and CHIMES
receiver ABI5 are required. The wizard and GUI generate this choice and its
environment; neither existing PAH choice nor the default changes.

The rate prescription follows [Montillaud et al. 2013](https://arxiv.org/html/1301.6507v1):
H-loss prefactor 6.8e17/s, even/odd-H thresholds 4.8/3.2 eV, harmonic
Beyer--Swinehart density of states, and cation H attachment 1.4e-10/5e-11
cm3/s for even/odd H. H13 loss uses H11 modes; attachment stops at H13.
Neutral H attachment and H2 pathways are absent, not established zero in
nature. Generic DL01 modes replace molecule-specific spectra, so this is
an explicitly approximate comparison, not a reproduction of that paper.

Photon excitation, IR cooling and H loss use one backward-Euler system.
Attachment consumes finite gas HI and moves binding/kinetic energy to
the daughter excitation state. Recombination competes over all H states
using a single electron inventory. Gas HI/electrons, PAH states, IR and
thermal energy commit together or roll back together. Per-state molecular
mass is `(288+H)*1.66e-24 g`, matching the CHIMES nuclear-mass convention;
the hydro elemental fields retain gas+solid nuclei. H0 remains solid C24,
not 24 free carbon atoms. Source injection remains neutral H12, with the
existing declared condensation comparison fraction.

Normal-H charge-dependent optical/cooling coefficients and IP=7.02 eV
are shared across H states; H-dependent bands/IP are not provided. This
mode admits occupied photons <=13.6 eV and PAH-bearing gas 10--10000 K;
cosmology, Fe and relative dust motion are rejected. It does not establish
carbon fragmentation, PAH survival under hard sources, or a general
multi-charge network. Restart identity binds H rates, DOS spacing,
binding energies, optical data and carrier layout. The existing neutral
and fixed-H restart conventions are preserved. See the active
[medium-term implementation record](../../provenance/medium_physics_implementation_2026-09-10.md)
for bounded native/live evidence, not universal production qualification.

### Optional PARSEC pair-instability feedback

History version4 supports five explicit fates (CCSN, failed SN, PPISN, PISN,
direct BH) on a declared domain no wider than 8--600 Msun. Select
`source_consistent`, enable wind/SNII/PISN with identical source windows
nested inside the IMF, and supply the native table/history together. Both
namelist generators expose this opt-in; defaults and legacy v1--v3 stay
unchanged. Channel3 owns all remnants; channel5 owns pair ejecta/energy and
releases them at the source event age. Pair energy shares the declared SN
CR/shock prescription; direct pair-event dust condensation is excluded.

`tools/build_parsec_pair_feedback.py` constructs the traceable nonrotating
Z=.008/.014, 14--600 Msun comparison from the original PARSEC tracks/yields,
HW02 and Woosley2017 energy sources. Its manifest declares the constant wind
speed, endpoint-calibrated phase composition, unresolved pair pulses and
62--64 Msun helium-core energetic bridge. No Python is used in the runtime.
Changing IMF upper mass changes normalization of all channels. This is not
a matched SED package, radioactive network or universal production approval;
see [implementation and evidence](../../provenance/parsec_pair_feedback_implementation_2026-09-10.md).
