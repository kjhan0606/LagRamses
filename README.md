[1]: https://ramses-organisation.readthedocs.io/en/latest
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music
[7]: https://github.com/osyris-project/osyris
[8]: https://github.com/pynbody/pynbody
[9]: https://yt-project.org

# cuRAMSES

**cuRAMSES** is an in-house, performance-optimised fork of the
[RAMSES](https://bitbucket.org/rteyssie/ramses) cosmological AMR code,
developed at the Korea Institute for Advanced Study (KIAS) for the
Horizon Run 5 and HR6 simulations. The cuRAMSES contributions
(recursive *k*-section domain decomposition, Morton-key hash-table
neighbour lookup, FFTW3 direct Poisson solver, spatial hash-binning
for SN/AGN feedback, variable-rank HDF5 restart, and CUDA-accelerated
multigrid/Godunov kernels) are bundled under [`patch/cuRamses/`](./patch/cuRamses/).

### Upstream lineage

cuRAMSES is derived from
[`bitbucket.org/rteyssie/ramses`](https://bitbucket.org/rteyssie/ramses)
master branch at commit
[`d92cd7d9`](https://bitbucket.org/rteyssie/ramses/commits/d92cd7d9)
(2016-11-08, RAMSES version 3.10), as identified by `git`-blob-hash
intersection of all unmodified base-directory source files. The
sub-grid physics (sink particles, AGN feedback, supernova feedback,
metal yields) trace back to the Horizon-AGN version of RAMSES
maintained by Yohan Dubois and collaborators, which is itself
not publicly distributed.

### Companion paper

J. Kim, *cuRAMSES: Scalable AMR Optimizations for Large-Scale
Cosmological Simulations*, MNRAS (2026, in revision).
Cite this paper when using cuRAMSES.

### Acknowledgements

The cuRAMSES contributions in [`patch/cuRamses/`](./patch/cuRamses/)
rest on contributions from many people:

- **Romain Teyssier** — public release and continued development of RAMSES.
- **Yohan Dubois** — Horizon-AGN version of RAMSES, on which the present
  in-house lineage is based. We thank him for sharing his code.
- **Taysun Kimm** — chemical enrichment model (alpha-element and stellar-wind yields).
- **Christophe Pichon, Sébastien Peirani, and Julien Devriendt** — design
  discussions on the Horizon-AGN code and simulation.
- The upstream RAMSES developer community for continued maintenance and
  improvements (see [`CONTRIBUTORS.md`](./CONTRIBUTORS.md) for the
  full list).

cuRAMSES is supported by the Korea Institute for Advanced Study and the
KIAS Center for Advanced Computation. J.~Kim is supported by KIAS
Individual Grants (KG039603).

### Build

### Equilibrium table cooling (`exact`)

In `PHYSICS_PARAMS`, set `cooling=.true.`, `cooling_method='exact'`, and
`grackle_table` to an absolute path to the existing Grackle binary table.
This uses precomputed equilibrium net cooling/heating and molecular weights,
not runtime Grackle or a non-equilibrium species integrator. Build without
`-Dgrackle`; `gamma=5/3` and `neq_chem=.false.` are required. `compare` retains
the original RAMSES cooling update and compares it against this solver.

The exact solver bilinearly interpolates the density/metallicity corner rates
and molecular weights first, constructs an energy/rate curve, and analytically
integrates its piecewise-linear segments. Equilibrium roots use the same
curve and cannot be crossed. “Exact” applies to that discretized curve at
fixed density, metallicity and radiation background, not the underlying
continuous cooling physics. Inputs outside the tabulated domain still clamp
to the table boundary, including zero metallicity if no primordial plane is
present; the table coverage and UV/dust assumptions must be checked before
production. Individual metal abundance ratios are not cooling-table axes.

Run the bounded analytic regression with
`bash simulation/snrt/tools/test_eunha_exact.sh` using GNU Fortran.

### Build options

### Nine-species non-equilibrium Grackle

Build a separate hydro executable with `GRACKLE_NEQ=1 HDF5=1 DUST_LIVE=0
CHIMES=0` and no legacy `-Dgrackle`. `GRACKLE_DIR` points to a matching
double-precision installation (default `$HOME/local`). The default field
profile gains nine passive densities, normally NVAR=27 with NENER=0.

Select `cooling_method='grackle_neq'`, `cooling=.true.`, `metal=.true.`,
`gamma=5/3`, `neq_chem=.false.`, `radiation_transport='none'`, and
`dust_mass_enabled=.false.`. Set `grackle_neq_data_file` to the **HDF5**
Grackle rate/UV file, not the equilibrium exact solver's binary table.
`haardt_madau` enables the uniform UV background; use a matching noUVB table
when disabled. In cosmological runs, the UV background follows RAMSES'
instantaneous `z_reion` switch (off for `z >= z_reion`, on below it), with the
same starting-redshift cap used by the legacy cooling initialization.
`self_shielding` selects Grackle's method 3 and requires a
compatible shielded metal table. `grackle_neq_h2_shielding=0|2|3` selects no
H2 shielding, a cell-width column approximation, or the local Jeans length.
The cell-width option is not a resolved line-of-sight RT calculation.

Species order: HI, HII, H-, HeI, HeII, HeIII, H2, H2+, electron. Fields are
Grackle hydrogen-mass-equivalent densities for electrons and physical
species mass densities for nuclei/molecules. Initial HII and H2 nucleus
fractions are `grackle_neq_initial_hii=2e-4` and
`grackle_neq_initial_h2=2e-6`. Existing element fields supply H/He inventories;
otherwise fresh ICs initialize the primordial H/He mixture after subtracting
the total metal mass. Restarts require matching HDF5 species-field metadata
and NVAR; missing chemistry is not silently reconstructed from an old dump.

H/He/H2 chemistry and internal energy evolve together; the equilibrium exact
and CHIMES sinks are not additionally applied. Metals remain equilibrium
and total-metallicity based. The hydro EOS remains fixed gamma=5/3, an
effective-model approximation for molecular gas, not a variable-gamma
molecular EOS. Stellar-return composition changes add neutral H/He; removal
preserves existing ion/molecule fractions. Delayed-cooling cells evolve
chemistry at fixed thermal energy. No dust formation, grain cooling,
photoelectric heating, or local RT rates are enabled in this route.

Private Grackle contexts isolate OpenMP owners. Cells with identical H/He
ratios and cooling masks are batched without rounding the abundance ratio.
Hydrogen, helium, charge and positivity checks precede commit; failed tiles
do not overwrite their input state. MPI, multi-level AMR and restart
qualification must be judged from integrated tests, not the bridge alone.

AMR caution (2026-10-08): fixed-level Grackle qualification does not cover
cold-gas prolongation. The 128³ adaptive baseline reproduced negative
internal energy on level 8 even with cooling disabled. A scoped conservative
prolongation limiter for `grackle_neq`, `interpol_var=0`, hydro/NENER=0
passed the original-IC two-step MPI4/OMP2 check (job 415996), including
level-8 hydro/cooling. This is not late-time production qualification.
The limiter contracts child deviations while preserving parent averages;
its roundoff reserve is measured in energy ULPs, not a temperature floor.
`mkrun.py` selects this conservative path for Grackle; other existing
model defaults are unchanged. `interpol_var` belongs to
`REFINE_PARAMS`, not `HYDRO_PARAMS`. Internal-energy interpolation
(`interpol_var=1`) is a diagnostic alternative, not proof of total-energy
conservation. See the [active AMR repair and calibration record](./provenance/galaxy128_grackle_calibration_baseline_2026-10-08.md).
The [native prolongation regression](./simulation/snrt/tests/test_grackle_conservative_prolongation.f90)
retains the four thermal cases executed in job 415996. Link it against the
matching hydro/NENER=0 build objects (excluding the RAMSES main object)
and that build's libraries; it calls the actual `interpol_hydro` routine.
Its optional `--reproduce-mass-closure` mode, linked to the retained uncorrected
build, uses a closed, variable-composition
stencil to reproduce the independent-minmod mass defect while verifying the
conservation of each parent field. An expected-defect result is not a repair
or production pass, nor an attribution of the observed cosmological failure.
The corrected Grackle path uses the final child density with parent metal,
11-element and nine-species fractions, preserving each carrier's parent mean.
`--verify-mass-closure` checks this on the same stencil and verifies that an
inconsistent parent is preserved rather than silently normalized. Subcell
composition slopes are first order; density/momentum/energy reconstruction
and unrelated passive fields keep their existing numerical treatment.
Grackle's CPU MUSCL path also closes the temporary face composition before
the Riemann solve, restricted to the face states actually consumed by that
solve. `--verify-face-closure` checks both flow directions and zero-mass-flux
LLF interfaces along each of x, y and z, positive periodic updates, independent
element inventories, nuclear species budgets and electron charge. GPU hydro does not
yet implement this correction and is rejected for `grackle_neq`; GPU Poisson
remains supported. Both namelist generators enforce the same distinction and
write the CPU hydro choice explicitly: CUDA builds otherwise default to GPU
hydro, even though the setup database's default is False.

A reference build is provided as [`bin/Makefile`](./bin/Makefile) /
[`patch/cuRamses/Makefile`](./patch/cuRamses/Makefile). The default
configuration uses the Intel `ifx` compiler with HDF5 parallel I/O and
optional CUDA acceleration (`make HDF5=1 USE_CUDA=1`). 128-bit Morton
keys for deep-AMR runs are enabled with `MORTON128=1`.

CHIMES users choose the cell solver before launching RAMSES:

```bash
export SNRT_CHIMES_COMPUTE_BACKEND=cpu     # default; CPU CVODE
# or
export SNRT_CHIMES_COMPUTE_BACKEND=hybrid  # CPU CVODE/LU + asynchronous GPU dark RHS/cooling
```

With no backend variables set, CPU is the default. `hybrid` selects the
level-queue CPU/GPU scheduler automatically and requires `CHIMES_CUDA=1`
and a GPU allocation. It does not implement a device-resident CVODE controller.
The device chemistry RHS supports dark chemistry. A nonzero `N_spectra` is
eligible only when every retained photon-group density is exactly zero; this
keeps the full CVODE state dimension for hot/atomic transition solves while
still allowing their dark RHS and cooling work on the GPU. Any nonzero photon
field stays on the CPU. The separate M5 spectral photo operator also remains
a CPU path. The earlier blanket `N_spectra != 0` rejection forced zero-photon
transition solves with retained spectrum slots onto the CPU; commit `ef01f5f`
removes that restriction without changing nonzero-radiation behavior. The
M5+dust+CHIMES restart result from job `440742` was numerically restart
consistent, but its resumed leg recorded zero CHIMES GPU calls. The corrected
backend and stricter per-rank runner check are awaiting a repeat run, so that
earlier result does not qualify hybrid restart execution. The cooperative
tile/level queue is currently called by the radiation-off dark-material
runtime; M5's `material_tile` instead uses its OpenMP dynamic cell loop and the
per-RHS cohort dispatcher.
On the paired manual `lageunha` Ada smoke (same binary/inputs, MPI2/OMP4),
CHIMES CPU vs hybrid elapsed times differed by less than 1% in both continuous
and resumed legs; this shows no end-to-end speedup and does not qualify
production scaling.
Do not combine this selector with contradictory legacy
`SNRT_CHIMES_CVODE_BACKEND`, `SNRT_CHIMES_RHS_BACKEND`, or scheduler settings.
Legacy lower-level selectors remain available for explicit experiments.
The qualified 512-cell dark-network fixture measured approximately 1.05x
speedup for the default mixture, but 0.79x for the mixed composition: hybrid
is opt-in, not a universally faster default. These figures do not measure the
photo-coupled M5 production path.

The independent dust/CHIMES path can also use the experimental
`SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated` in a matching CHIMES/CUDA build
even with `radiation_transport='none'` and all other GPU modules disabled.
CHIMES then requests the shared stream pool at serial startup. This backend
is experimental, not a demonstrated production speedup; CPU remains the
default. See the [native runtime guide](./simulation/snrt/NATIVE_RUNTIME.md)
for requirements and current numerical/performance limits.

For experimental CHIMES worker assignment comparisons,
`SNRT_CHIMES_GPU_WORKER=thread0` permits only OpenMP worker 0 to submit
reaction batches to CUDA; other workers compute chemistry on the CPU.
The default `any` allows whichever worker acquires an available stream.
Both use the existing dynamic cohort scheduling. The selected GPU worker
still waits for each current-state RHS batch before continuing host CVODE;
fixed worker assignment alone does not guarantee acceleration.

`SNRT_CHIMES_CELL_SCHEDULER=level_queue` selects an experimental level-wide
queue for the RT-off cold CHIMES/dust path. CPU-only uses all OpenMP workers
on native cell solves. `cuda_integrated` makes every worker a broker unless
`SNRT_CHIMES_GPU_BROKERS` sets a smaller count. Each broker keeps up to 128
CVODE continuations and its own stream. A broker launches once enough of its
own fibers are blocked, so the rest of its host integration overlaps that
batch. One shared launcher was slower on H200: batches grew to about 64
cells, but a single stream cut the 512-cell hybrid to about 0.34× the
8-thread CPU wall. A level with fewer than 16 cells per worker stays on the
CPU, because a short GPU batch loses to the parallel CPU team. Chemistry
completes before the existing dust stage and collective level commit. Each
continuation resumes only after its own current candidate's GPU result is
ready; stale RHS values are not substituted. The matched one-broker 128³
test passed numerical regression but did not accelerate CPU-only. That
one-broker schedule is no longer the default.

The dark GPU RHS normally returns a compact FP64 result: species creation
and destruction rates plus the reaction fields read by host cooling. Full
reaction scratch remains on device. `SNRT_CHIMES_GPU_RESULT=full` retains
the previous full-result transfer for comparison; paired rate verification
and flux callbacks automatically use full results. Compact results do not
change the chemical network, CVODE tolerances, or thermal formulas.

### Binary particle restart compatibility

The binary particle reader supports both legacy RAMSES part files without a
compact `ptypep` record and LagRamses files that store `ptypep` immediately
after `levelp`. It uses `part_file_descriptor.txt` when available; for older
descriptor-free files it checks the complete expected record layout before
reading and stops if the layout is ambiguous. Keep any descriptor alongside
its particle output files. A legacy file without `ptypep` can reconstruct only
the standard DM/star/sink classes from IDs and birth times; SIDM/ADM restarts
require the explicit particle-type record and are rejected otherwise. This
compatibility rule applies to binary particle restarts, not the HDF5 reader.

Refined-level multigrid phase profiling can be enabled independently with
`make MG_PROFILE=1 HDF5=1 USE_FFTW=1`. This retains the production `-O3`
flags and adds cumulative `MGPROF` reports; it does not enable the broader
debug instrumentation.

The CPU refined-level multigrid solver caches the centre and six face-neighbour
grid indices for the duration of each solve. Coarse MG levels also cache the
owner rank. The caches are rebuilt after the MG topology is constructed and
released during MG cleanup; their storage cost is 28 bytes per active fine grid
and 56 bytes per coarse MG grid when default four-byte integers are used. This
removes repeated Morton/hash lookups without changing the stencil or MPI
exchange backend. The fine-grid cache is shared with the optional CUDA upload
path.

Refined-MG halo exchanges honor the global `exchange_method` policy. With
`exchange_method='auto'`, the solver measures each level's fixed communication
graph once after constructing the MG hierarchy, caches separate forward and
reverse choices, and reuses them throughout the V-cycle. Sparse graphs use
nonblocking P2P. MG also keeps P2P for aggregate rank payloads up to 1 MiB,
because its existing rank buffers avoid the extra collective packing cost.
Larger intermediate graphs use the sectioned k-section schedule, while dense
bounded messages may use `MPI_Alltoallv`. The explicit values `p2p`, `ksection`,
and `alltoallv` force the corresponding backend. Backend selection is not
repeated inside the MG iteration.

### Sink and AGN JSONL diagnostics

The active [`patch/lagRamses/`](./patch/lagRamses/) sink implementation writes
two single-writer JSONL products. Both are enabled by default and can be
configured in `&physics_params`:

```fortran
smbh_capture_ledger = .true.
smbh_capture_ledger_file = 'smbh_capture_ledger_v1.jsonl'
agn_coarse_dump = .true.
agn_coarse_dump_file = 'agn_coarse_state_v1.jsonl'
```

The [SMBH capture ledger](./patch/lagRamses/SMBH_CAPTURE_LEDGER.md) records
each sink group immediately before irreversible merge compaction. Its JSONL
transaction contains the original global `sink_id`, surviving
`primary_sink_id`, primary flag, member state, and all unordered pair
diagnostics. A capture row is a numerical merge event, not by itself a claim
of physical binary formation or coalescence. Validate complete transactions
and restart duplicates with:

```bash
python3 patch/lagRamses/aux/validate_smbh_capture_ledger.py \
  smbh_capture_ledger_v1.jsonl
```

The [AGN coarse-step ledger](./patch/lagRamses/AGN_COARSE_STATE.md) writes one
row per `(nstep_coarse, sink_id)` before feedback deposition and accumulator
reset. It includes sink mass, gas and BH angular momenta, their alignment
angle, Bondi/Eddington/accreted rates and masses, radiative efficiency,
bolometric luminosity, feedback state, saved energy, and derived coherent-disk
episode quantities. Undefined disk or angle fields are JSON `null`. Give each
independent run a distinct output path, and deduplicate identical restart keys
while treating conflicting duplicates as provenance errors.

---

![GitHub logo dark-mode-only](./doc/img/full_project_logo_dark.svg#gh-dark-mode-only)
![GitHub logo light-mode-only](./doc/img/full_project_logo.svg#gh-light-mode-only)

## This is the [ramses](https://github.com/ramses-organisation/ramses/) GitHub repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible, radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree.
[ramses](https://github.com/ramses-organisation/ramses/) is written in  Fortran 90 and is making intensive use of the Message Passing Interface (MPI) library.

### ⬇️ Get the code

Download the code by cloning the git repository using
```
$ git clone https://github.com/ramses-organisation/ramses
```

### ℹ️ Get Support

You can go to the user's guide using [Read The Docs][1].
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).
You can get support on Ramses on the [Github Discussion page](https://github.com/ramses-organisation/ramses/discussions)

### 🛠️ Tools and ressources

To generate idealised initial conditions of galaxies, check out [DICE][5].
To generate cosmological initial conditions, check out [MUSIC][6].

To visualize RAMSES data, we encourage you to use [YT][9], [OSYRIS][7] or [PYNBODY][8].

You'll find a lot of useful ressources, links and news about the Ramses community on https://ramses.cnrs.fr/, a website edited by the French "SNO" Ramses.

## LRD/JWST simulation programme

Production launch caution: `simulation/snrt/runs/simulation_ready_preflight.sh`
supports `--mode physical-production --production-manifest PATH`. It performs
the existing G0 audit, not a simulation or job submission. An approved
`release_record` must bind `executable_asset_id`, `namelist_sha256`,
`approval_id` and `environment_variables` to this exact checkout and launch;
`production_rules.production_status` must be `ready`. Clean source identity,
selected-asset license/provenance and a locked dependency recipe remain
required. Undeclared SNRT/PHASE0/CHIMES environment overrides are rejected.

Native builds declare `backend.runtime="native_ramses"`, compiler and MPI
identities in the environment record and retain the existing hashed
`dependency_files.resolved_lock` recipe. They do not require the P0 JAX
runtime; the JAX profile still checks its pinned versions. The historical
manifest bound to `/gpfs/kjhan/LRD_JWST` is **not** a release record for another
checkout. A mismatched launch is refused before examining that worktree's
Git or assets. Reference-control tests, or a successful Grackle RT/dust-off
run, do not approve physical sources or grant a production release.

Bounded M5 coupled checks use
`simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch` with
`M5_CHIMES_BINARY` naming the verified CHIMES-enabled executable.
`M5_SOURCE_PROFILE=effective_population` selects the already evaluated
Kroupa/effective-SSP mixed-v5 source and empirical Ia contract, not the global
default IMF. Its declared truncated source/dust contracts remain reference
controls; this is not an unrestricted scientific production approval.
The normal check writes no raw dump. Optional `M5_VERIFY_RESTART=1` performs
the same MPI2/static-layout continuous-versus-restarted comparison using
three small checkpoints, retained until evaluation. Submit on Slurm with
two ranks/two CPUs per rank and explicit memory/GPU allocation; the runner
uses the existing explicit local-rank CPU-slice recipe via
`slurm_rank_bind.sh`, including per-rank `OMP_PLACES`. Generic `cores` places
alone can be rebuilt over the shared allocation during MPI/OMP initialization.
Use the existing execution record
[`provenance/snrt_m5_radiation_momentum_2026-10-11.md`](provenance/snrt_m5_radiation_momentum_2026-10-11.md)
for actual executable/input hashes, external dependency locations, output
policy, and qualification limits. Do not reuse an RT/dust-off calibration
input or old-root release manifest to claim this coupled profile qualified.
The trace-free 2026-10-12 selected-profile restart check passes 449 datasets
at the unchanged `rtol=1e-12, atol=1e-25`. M_N keeps six physical face-flux
channels separate through reverse MPI exchange, then combines them in fixed
cell-face order; this avoids remap-dependent partial-sum grouping. Shared
faces still supply equal/opposite conservative increments. This bounded
static MPI2 result is not universal AMR/rank-count reproducibility or physical
asset approval. Temporary diagnostic tracing is not part of the live code.

This repository also hosts the HR5 high-redshift SMBH, dual-AGN, and JWST
little-red-dot simulation programme. Its active scientific scope is recorded
in [`docs/lrd_zoom_simulation_design.md`](./docs/lrd_zoom_simulation_design.md).

- `simulation/` contains reproducible configurations, source ledgers,
  analysis tools, test fixtures, and Slurm submission scripts.
- `manifests/` records immutable locations and SHA256 hashes for large assets;
  the assets themselves are not committed.
- `provenance/` records code, build, migration, and external-data provenance.
- `docs/archive/legacy_dualagn_paper3/` preserves the superseded Paper-III
  planning record.
- `docs/lrd_manuscript/` is the placeholder for the future LRD/JWST
  manuscript. The existing `paper/` submodule remains the cuRAMSES companion
  paper workspace.
