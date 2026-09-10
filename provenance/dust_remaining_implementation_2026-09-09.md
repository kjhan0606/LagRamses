# Explicitly approved remaining dust physics

## Active continuation and sequencing — 2026-09-10

The operator explicitly requires completion of this existing implementation
and its integrated RAMSES evolution/restart **before magnetic-field code work**.
MHD literature/code inspection is not permission to interleave an MHD build or
change the hydro field layout now. No MHD implementation has begun.

Relative-motion live wiring and the bounded integrated/restart test are now
complete. This section supersedes the dated intermediate pending lists below;
it does not extend the supported physics domain. Directed
growth/destruction/bin exchange now carries gas + bin absolute momenta through
the rate substeps, including zero-net-mass cycles and fixed Fe/PAH reservations.
The native kinetics tests pass GNU and Intel. Independent radiation energy is
stored as FP64 `reference_energy*N + shift`; photon count remains FP32. Raw
checkpoint v7 retains the signed shift, v6 migrates to zero shift. The primary
paired transport, phase force/work, hydro EOS/source coupling, AMR/HDF5 layout
and IR absorption/emission coupling are connected as one bundle.

PAH absorption accepts actual captured photon counts separately from deposited
excitation energy; recoil work must not change the event count. The existing
GNU mass smoke and receiver smoke pass with these focused cases in
`.dust-relative.MHlIN8/`. Full profile compilation succeeded in the
fresh `.dust-relative-build.TjwEQt/` (NVAR338, NENER1, Fe+PAH+CHIMES,
SNRT/DUST_LIVE/DUST_DYNAMICS, HDF5, CPU/OpenMP hydro). The separate live
evolution/restart evidence follows below. VPATH remains unchanged; old binaries
and outputs are preserved. The default remains off: bounded integration is not
unrestricted scientific production qualification.

Live connections now compile/link in that profile. Primary transport retains
signed energy corrections through MPI/coarse reflux and FP32 source/commit
rounding. Primary grain absorption and moving scattering stage absolute
momenta and mechanical work; PAH gets separate captures/excitation heat.
IR absorption/emission/scattering solve phase work with the material iteration.
CHIMES receives the staged full-component kinetic state; its primary counts no
longer receive a second angular weight. A dry commit validates all ranks before
any publishes or destroys its rollback snapshot. Existing default scalar paths
remain selected when relative dynamics are off.

The hydro implementation includes conservative phase faces, EOS/CFL, CR gas
advection/work, gravity, drag, SF and source donors. Ephemeral channel-condensed
mass provenance retains C/S release segments and Fe/H/C donor contributions;
opposing donor velocities are supported without changing the 164 old real
source outputs (native bitwise comparison in `.dust-source-donors.Dhb89W/`).
CLI/GUI setup carries the same default-off namelist choice and requires an
explicit positive neutral gas collision cross section. This is a bounded
first-order noncosmological CPU/OpenMP profile, not a charged-grain/AP/GPU or
unrestricted sublimation model.

First integrated attempt `.dust-relative.MHlIN8/live.ubLzAC/`:
the initial namelist had `interpol_var` in HYDRO rather than REFINE; corrected
before the actual evolution. Step 1 completed with IR balance residual
2.2313e-13. Its 60 MiB dump has all 338 fields finite, gas+stars mass exactly
10 code units, 512 stars, positive gas internal energy, and nonzero phase drift.
Step 2 exposed a zero-rate saturation-roundoff failure: an already admitted
elementwise representational excess was clipped before donor tracking. The
coupled path now preserves exact zero-rate identity and the actual donor mass
at admitted boundaries, using the existing element tolerance, not a new floor.
GNU/Intel zero-rate and active-rate regressions pass; 512 saved-row replays pass.

A fresh full build after those fixes is preserved as
`.dust-relative.MHlIN8/ramses_integration_fixed3d`, SHA256
`a2666189eee800fb02fbfaacc9cda4e45752ac1170d5aaf499aa79a2f3173600`.
The MPI2/OMP2 two-step rerun at `.dust-relative.MHlIN8/live-fixed.rN91Oi/`
completed with exit 0 in 171.657 s. Restart from a private copy of its first
dump at `.dust-relative.MHlIN8/restart.DYyiJO/` completed with exit 0 in
128.452 s. All **1,068 datasets** in the final `output_00002/data_00002.h5`
are byte-identical between continuous and restarted execution. The failed
attempt and its outputs remain untouched.

The final 512-cell state has all 338 hydro fields finite, 1,024 stars,
gas+star mass exactly 10 code units, minimum gas density 9.930134587863062,
and minimum gas internal-energy density 1.4914709557712008e-4 (code units).
Actual cumulative stellar returned mass is 2.1719309558997498e-5 code units.
Nonzero grain/gas relative velocities are present (largest bin maximum
2.57494607233677e-10 code velocity). IR balance residuals are 2.2313e-13
and 2.5124e-12. Each snapshot is approximately 60 MiB.

This is an actual noncosmological periodic 8^3 RAMSES evolution with SF,
physical stellar return, CR gas support/advection, CHIMES, C/S/Fe dust,
neutral PAH, IR exchange and relative dynamics. It is not an AMR convergence
or long-time astrophysical calibration test. Primary source photons are zero
in this particular integration; nonzero primary absorption/moving scattering
and paired transport are covered by native backend/receiver tests, not claimed
as a positive-primary-source RAMSES evolution. Growth/destruction rates are
off in this integration; active-rate kinetics have separate GNU/Intel tests.
The explicit neutral gas collision cross section 2e-15 cm^2 is a test input,
not a calibrated universal choice.

Native moving scattering passes 473 assertions (maximum relative number,
energy and momentum residuals 4.15e-16, 1.97e-16 and 1.39e-15); moving IR
passes 50 assertions. The separate six-bulk-plus-PAH moving IR callback also
passes. Evidence resides in `.dust-relative.MHlIN8/`,
`.dust-relative-build.TjwEQt/`, `.paired-energy.0rWe07/` and the donor report
`.dust-source-donors.Dhb89W/REPORT.md`.

The existing transaction runner passes Intel/GNU serial dry-commit checks
and MPI2 rank-local invalid-coarse-energy rejection: no rank publishes,
snapshots remain available for collective rollback, and a valid retry passes
dry validation followed by collective acceptance and commit
(`MPI_DRY_COMMIT_PASS`, `MPI_SMOKE_PASS ranks=2`, `SMOKE_RUN_PASS`).
The CLI/GUI status text now reports the bounded evidence rather than the
superseded pending-run status; 44 existing frontend tests run, 43 pass and
one display-dependent test skips. `git diff --check` is clean.

Remaining model limits are explicit, not additional gates: first-order phase
transport/splitting; nonrelativistic moving grey-group scattering (|v|/c <=
0.01), no cross-group spectral redistribution; CHIMES reference-spectrum
reaction yields with a signed actual-minus-reference absorbed-energy
correction; neutral Epstein drag without Coulomb/Lorentz forces; OpenMP-only
paired primary energy transport; no cosmology, sinks/AGN or IR-coupled
sublimation in this profile. Native AMR/reflux/restart wiring checks are not
a claim of full evolving-AMR qualification. Magnetic-field implementation
remains the next task after this closure, not part of this completed bundle.

## Continuation — directional absorption and mechanical work, 2026-09-10

Continued the same relative-motion item; no new stage, selector or audit.
The shared OpenMP/CUDA primary-photon cell now optionally returns the signed
first angular moment of **accepted dust absorption**, excluding H/He use and
returned photons. The existing scalar C entry points remain available and
their results are unchanged. New moment entry points cover serial/OpenMP,
CUDA and the real stream-lease hybrid dispatcher. Hybrid batches stage the
moment alongside their photon/atom outputs; a late failure publishes none.
The Fortran runtime accepts an optional `(leaf,group,3)` output, and the
RAMSES prepared transport path accumulates it across substeps before return.
No extra calculation/allocation is requested by the existing hydro driver.

`dust_fv_absorption_kick` connects these moments to absolute grain momenta
using physical h*nu/c and explicit code-unit scales/absorption shares.
`dust_fv_radiation_kick` adds only mechanical work to hydro E and returns
the remaining absorbed energy for the separate solid/PAH heat receiver.
Gas momentum/thermal energy remain unchanged during the dust-only kick.
Underfunded work, invalid shares or unphysical moments reject atomically;
there is no heat clipping or second full-energy addition. This receiver is
absorption-only: scattering/re-emission and photoelectron channels are not
silently counted as absorption heat.

Evidence: `.dust-momentum.YDPqWT/`. `backend-intel.log`, `hybrid-intel.log`,
`mass-intel.log`, `mass-gnu.log` pass (GNU bounds/FP traps). Actual GPU-2
`backend-gpu.log` passes scalar parity and accepted moments, maximum scalar
relative difference 1.09567e-7; `hybrid-gpu.log` passes mixed CPU/GPU, busy
stream fallback and late-batch moment rollback, scalar difference 1.19209e-7.
`fortran-openmp.log` exercises the actual C/Fortran runtime -> multifluid
absorption receiver, not just synthetic impulses; beam/isotropic cases,
physical-c conversion, gas invariance and work+heat closure pass. Existing
Fe/PAH/sublimation tests also pass. The first GPU attempt failed allocation
on nearly full GPU-0; retry used free GPU-2 without altering another process.
The initial manual link omitted cuBLAS/cuFFT and was corrected.

Full RAMSES relink also passes (`build-ramses.log`, NVAR317/NENER1,
SNRT+DUST_LIVE+Fe+PAH+CHIMES+HDF5, OpenMP). New binary:
`.dust-momentum.YDPqWT/ramses_absorption3d`, SHA256
`05fb91b18bdf145e16e508c342387bbb92b1cfbdddf8fe473c64aea7ef1584e5`.
The old PAH live binary still has SHA256
`a65d6a9f8e3f6aaa8f5f6186047e7b27d1247a61a1a2f71c92aba702f1bb5b74`;
it was not overwritten. This relink is not a new evolution/restart test.

**Remaining, not completion:** the hydro driver does not yet request/apply
this optional moment, and the native receiver is not yet a live grain-fluid
update. Persistent grain momenta, their EOS/CFL/face and AMR/restart handling,
donor-velocity source/phase exchanges, and the full radiation transaction
(including scattering/IR work) must still be connected coherently. Do not
enable a partial namelist option or claim relative-motion production readiness.
VPATH/defaults and existing outputs remain unchanged. No simulation launched.

## Continuation — relative-motion spatial operator, 2026-09-10

The operator requested continued implementation after the bounded PAH live
connection. The next unfinished work is **relative dust/gas dynamics in the
same approved bundle**, not a new audit or approval stage. The PAH completion
record below remains valid. Relative dynamics are not yet selectable in
RAMSES; this continuation implements spatial fluxes and their coupled drag
step, not a claim that its live connections are finished.

`patch/lagRamses/dust_multifluid.f90` now supplies a native conservative
Eulerian gas + arbitrary pressureless-grain face operator. It evolves total
rho/p/E plus each grain's **absolute component momentum**, not the nonlinear
relative quantity J_b=rho_b*(v_b-v_bary) as if J_b were a conservative scalar.
Existing local barycentric drag helpers remain valid. E includes kinetic
energy of every component and gas thermal/nonthermal energy; grain thermal
and vibrational energy remain separate, as in DUST_LIVE. Gas pressure is
computed after removing all component kinetic energies, not just barycentric
kinetic energy. This realizes the component conservation equations discussed
by [Laibe & Price 2014](https://arxiv.org/html/1402.5248), not their SPH
discretization or an asserted asymptotic-preserving method.

The face flux transports gas carriers at gas velocity and grain carriers at
their own component velocities. Positive state carriers and signed momentum
are distinguished. A shared Rusanov wave bound includes gas/CR sound speed
and every dust velocity. CR pressure work is retained in its energy carrier;
total E already contains that energy and must not receive it again. The
periodic finite-volume driver commits all cells together, checks reciprocal
faces/CFL/positive states, and rejects invalid states without clipping.

`dust_fv_periodic_advance` connects spatial transport to the existing implicit
drag receiver in a single transaction: a failed last-cell drag discards the
entire update. Drag conserves total momentum/E; reduced drift kinetic energy
appears in gas thermal energy without a second heat addition. Stopping times
are supplied and frozen during the step. The scheme is first-order, and
stable local stiff drag is **not** proof of small-grain/stiff spatial accuracy.
Pressureless single-velocity grains also do not resolve multi-stream caustics.
These limits must remain explicit in the eventual runtime admission.

The existing `dust_mass_smoke` now exercises seven grain components, CR, a
gas nuclear tracer, four grain-energy carriers and 128 PAH excitation-mass
carriers. No Python framework, new external audit or RAMSES run was added.
MPI/AMR/restart or GPU execution is not claimed by this native periodic test.

Actual Makefile/mpiifx checks in `.dust-pah-live.kBEoZy/`:
`multifluid-coupled-build.log`, `multifluid-coupled.log` pass. Counterflow
L1 errors for 32/64/128 cells are 5.884736899e-5, 2.982480664e-5,
1.501092713e-5; refinement ratios approximately 1.97 and 1.99. Conservation
residuals are <=6.94e-16. Also pass: gas-only pressure acceleration, correct
gas CR advection/PdV, PAH population sum matching its independently moving
grain mass, zero total mass flux with nonzero phase counterflows, positive
stiff-drag heat, CFL rejection and late-cell atomic rollback. Existing
Fe/PAH/source/growth/sublimation tests remain passing.

GNU bounds/FP-trap builds in `.dust-multifluid.cBvGQj/` also pass
(`build-coupled.log`, `coupled.log`, exit0). The same refinement errors are
obtained; maximum conservation residual is 1.041e-15. The initial GNU command ordered the yield
table module before its dust dependency; it was corrected, not treated as a
physics failure. No existing output was removed. VPATH and inactive runtime
defaults are unchanged; no namelist/mkrun selector was exposed prematurely.

The next live edit must connect this as **one coherent RAMSES path**:
component momentum and energy interpretation through face fluxes/CFL,
source/phase changes, chemistry, actual radiation impulses, AMR/MPI and
restart. In particular, the current co-advection carrier renormalizer cannot
be reused on independently signed phase fluxes: zero total mass flux may
still contain nonzero gas and dust counterflows. This is remaining physical
wiring within the existing item, not another list of approval gates.

## Latest checkpoint — 2026-09-10, requested item 1 (PAH)

The latest request's **item 1 is PAH stochastic heating**, not the older
numbered chemistry inventory below. Its remaining live connections are now
implemented: 128 mass carriers, mixed C/S/Fe/PAH radiation, source H/C and
excitation energy, gas chemistry/growth reservations, co-advection, HDF5
restart and namelist/mkrun/GUI. This closes the **bounded neutral PAH live
implementation**, not the whole four-item bundle or arbitrary-environment
production qualification. Detailed evidence and limits are appended below.

## Current operator instruction: one implementation bundle

The operator explicitly grouped all four remaining items into **one bundle**:
sublimation completion, separate metallic Fe, PAH stochastic heating, and
relative dust/gas transport with drag and radiation momentum. This supersedes
any interpretation of the numbered inventory below as four sequential
approval or audit stages. All four stay in scope; none has been deferred.

Work across their shared dependencies together: independent solid mass and
element carriers, material/phase and kinetic energy, radiation exchange,
source injection, AMR/MPI transport, restart and namelist/mkrun selections.
In particular, the current two-material/four-bin state and common-temperature
IR receiver cannot simply be relabelled as Fe/PAH or relative dynamics.
Do not introduce user-selectable modes until their corresponding live paths
exist, and do not mark a downloaded table or isolated kernel as completion.

Native builds and focused checks remain ordinary implementation work, not
new gates. No per-item approval pause or external audit. At the bundle's end,
report the combined implementation and verification result, limitations and
production-readiness status for operator confirmation. A genuine missing
physical input or required model choice must be reported explicitly rather
than invented to make the bundle pass. The previous graphite comparison is
completed evidence within this bundle, not proof that the bundle is finished.

Bundle status: **in progress, not production-qualified**. No new completion
date or percentage is inferred from the four-item inventory.

The operator now explicitly requests implementation of all previously listed
extensions, followed by operator confirmation. This is new scope approval,
not a claim that those extensions were previously missing mandatory gates.
Workspace: /gpfs/kjhan/LRD_JWST; origin: kjhan0606/LagRamses; starting commit
80265d5. Preserve the unrelated 88-line generator deletion and existing runs.

Implementation scope (no added items):

1. H2/CO chemistry, grain surface formation, photodissociation/shielding,
   molecular cooling and competition with dust for C/O.
2. Metal-ion non-equilibrium chemistry/cooling under the actual radiation
   field, replacing CIE only in an explicit new selection.
3. Separate Fe grains, PAH stochastic heating and sublimation.
4. Dust/gas relative motion, conservative drag and radiation momentum.

Keep defaults and existing comparison results intact. Use native runtime
receivers, persistent/transported physical state, existing MPI transactions
and restart identities, and matching namelist/mkrun/GUI selections. A local
kernel, source download or unconnected switch is not a completed item.
No per-item external audits or new generic test framework. Confirmation is
requested only after integrated verification, or if an actual dependency
cannot be resolved within existing authority. Calibration is not included.

## Source/dependency admission in progress

Official CHIMES source (LGPL3+) was cloned at
`a58e5c0311993b51abc63d84fff2958a0104f6d6`; its data repository at
`1a837f823eeaa2489d8b11d5a6e7967a4d07ab47`. The 16.9 MB main reaction table
was downloaded using selective LFS, not all UVB/eq tables. Full source and
data remain in `.dust-extension.AOz7mU/`. Native C library integration is
being evaluated; no runtime Python chemistry or CIE relabelling.

Upstream CHIMES RT coupling omits depletion of dissociating photons. The
adapter below adds a matched molecular photon sink and uses tables built
specifically for the actual nine SNRT groups. Local dust area/temperature
replace the fixed-Milky-Way surface normalization; the adapter removes the
duplicate CHIMES gas/dust accommodation term because SNRT owns that exchange.

Sources: [CHIMES downloads](https://richings.bitbucket.io/chimes/download.html),
[RT coupling](https://richings.bitbucket.io/chimes/user_guide/AddingChimes/Coupling_to_RT.html),
[GD89 sublimation](https://adsabs.harvard.edu/pdf/1989ApJ...345..230G),
[DL01 stochastic heating](https://arxiv.org/abs/astro-ph/0011318).

## Implementation checkpoint (not a completion/approval gate)

Native cell receiver implemented in `patch/lagRamses/snrt_chimes_bridge.c`,
with ISO C binding in `snrt_chimes.f90`. No Python interpreter is called by
the solver. The existing thermochemistry smoke now has 52 passing assertions
with CHIMES enabled, including four concurrent OpenMP cells; the separately
compiled no-CHIMES baseline retains its 33 passing assertions. Logs are in
`.dust-extension.AOz7mU/native-chimes-rt-unit.log` and `native-legacy-unit.log`.
The latest ABI3 cell log is `.dust-chimes-live.GYestM/native-unit-abi3.log`.
Live integration and its separate status are described below.

The receiver currently covers 157 species, source-backed H2 dust formation,
H2/CO line-data photodissociation, local molecular self-shielding, molecular
cooling and non-equilibrium metal ions. A reconciliation API reserves atoms
locked in molecules against grain growth and adds neutral ejecta without
silently destroying CO. The live dust/feedback path now calls this API.
All eleven nuclei and charge are checked before accepting a trial. SHA256
of the main reaction data is pinned; all nine group edge pairs and mean
energies match exactly, and the receiver exposes actual file SHA256 bytes
for the live restart identity. Initialization requires patched ABI3.

Two substantial errors were independently reproduced and corrected:

* Upstream current-rate storage used `malloc`; zero-CR updates skip the CR
  array but cooling reads it. With RT table allocation history this produced
  a spurious negative heating term, dropping 100 K to 10 K in one second.
  Zero-initialized storage restores 99.9999999997615 K and the dark
  1e10-second result agrees with the non-RT calculation (97.64660 K).
* Broad-spectrum outer-shell-only heating integrals can vanish in a hard RT
  bin while inner-shell absorption remains nonzero. Dividing two floored
  zero integrals then gave a spurious 1 erg per absorbed photon. The offline
  builder now integrates primary photoelectron energy against the SAME
  Verner shell cross sections and Kaastra branching probabilities as photon
  absorption, including Auger-producing channels. The prior 43604 K trial is
  now 10004.03 K, within its photon supply. Quadrature uses bin support and
  an absolute tolerance below the small physical integrals, not 1e-30.

Both upstream modifications are preserved as reviewable patches:
`simulation/snrt/data/chimes_native_receiver.patch` and
`simulation/snrt/data/chimes_tools_rt.patch`. Source/tools versions are
`a58e5c0311993b51abc63d84fff2958a0104f6d6` and
`411fd88eb38ae05bd47adf598f1db5dc08fb4579`. SUNDIALS is double-precision
5.8.0 (`e8a3e67e3883bc316c48bc534ee08319a5e8c620`), installed locally in
`.dust-extension.AOz7mU/sundials-install`. No system installation was changed.

Latest complete nine-group tables: `.dust-extension.AOz7mU/groups-v7/`.
Earlier `groups*` attempts remain diagnostic artifacts; do not select them.
The builder is `simulation/snrt/tools/build_chimes_group_tables.py`; invoke
with `--chimes-tools`, `--main-data`, `--leiden-data` and a NEW `--output`
directory. Group moments use an explicitly declared mean-matched within-bin
power law, not a claim to resolve a stellar/AGN spectrum inside a bin.

Leiden input files are retained in
`.dust-extension.AOz7mU/leiden/new_cross_sections_hdf5/`:

* H2 SHA256: `925bff4d51097c93679845e54210b55a3ba81cfdea4990d1be6b8d20db5c57c7`.
* CO SHA256: `18cf140b5332e15244cff2714ab880969337aae3ed8f69e60d9f77622669eced`.

Source: [Leiden cross-section downloads](https://home.strw.leidenuniv.nl/~ewine/photo/downloads.html).
The source tables' non-ionizing absorption and dissociation are distinct;
the native photon sink counts fluorescent absorption, not just destruction.
The local shielding closure uses the CHIMES tables with a 1 km/s unresolved
turbulent width. Other molecular channels retain CHIMES' FUV-shape
approximation, but their coefficients now follow the evolving photons and
have a matched photon sink. Unspecified Auger-cascade/fluorescence energy
is NOT silently converted into heat. Re-emission transport and its energy
ledger must be addressed in the live receiver; this is not a completed
conservative radiation/material-energy implementation.

Build on lageunha from a top-level private build directory using
`make -f ../bin/Makefile snrt_thermochemistry_smoke SNRT=1 DUST_LIVE=1
NENER=1 HDF5=1 CHIMES=1 USE_FFTW=0 CHIMES_DIR=<patched source>
SUNDIALS_DIR=<local installation>` (one shell line). The native test reads
`SNRT_CHIMES_MAIN_DATA`, optional `SNRT_CHIMES_GROUP_DIR`, and the existing
`SNRT_SECONDARY_TABLE_CONTRACT`. Its link now resolves SUNDIALS/HDF5 without
an LD_LIBRARY_PATH workaround. CHIMES=1 alone does not enable live physics.

## Live chemistry integration

Explicit selection: `dust_cooling='chimes_neq_v1'`, native CHIMES build,
`carbon_olivine_2size_v1`, `dl01_composition_v1`, and gamma exactly 5/3.
Existing defaults remain unchanged. `mkrun.py` and its shared namelist/GUI
generator expose and validate this selection. The chemistry runs as native
C/Fortran, not Python; Python builds data tables offline and inspects dumps.

* Species occupy 31--187 in the standard NENER=1 layout, as density-like
  `atomic_mh*n_i` passives. They are not additional baryonic mass.
* The CPU MUSCL path normalizes independent chemical/dust face carriers
  together and derives total-element fluxes from them. Refinement uses the
  parent chemical composition with the conservatively refined density.
* Feedback and grain exchange reconcile the chemical state with the gas
  element ledger. Every molecular species reserves its constituent nuclei
  against dust growth. Thermal energy uses the actual 157-species EOS.
* The old H/He opacity and thermal receiver are disabled in this selection.
  Actual angular photon densities feed native chemistry; its absorbed
  photons are removed from those same groups before accepted state commit.
  Actual source edges AND means must exactly match the chemical tables.
* HDF5 binds all 320 table hash bytes plus closure version, field offset,
  nuclear mass unit, and the library's actual Boltzmann constant. Restart
  initialization validates the stored state without another chemical update.
* RT disables CHIMES column extinction, but molecular line cooling still
  needs a finite escape column. ABI3 evaluates current CO/H2O/OH columns
  as n_species*cell_size with thermal widths (`StaticMolCooling=1`), instead
  of the singular zero-velocity-gradient Sobolev column. The callback guard
  leaves other CHIMES callers unchanged.

The formation-area reference is not an arbitrary new fit: Richings+2014a,
section 2.4, eq. 2.12 uses 1e-21 cm2/H. The receiver replaces this with the
actual geometric two-size area. Its dust-ratio times boost equals
area/(1e-21 cm2/H), so the intermediate 0.01 dust/H mass normalization
cancels. Applying that common area scaling to grain recombination remains
an explicit approximation, not a size-resolved grain-charge model.
[Richings+2014a](https://arxiv.org/html/1401.4719v2).

Radiative limitations: CHIMES photoelectric gas heating is disabled until
it can be debited from the same absorption energy currently deposited in
SNRT grains. Unresolved fluorescent/cascade radiation is an escaping-cooling
approximation, not transported re-emission. `SNRT_CHIMES_PRIMARY_ABSORBED_EV`
is an absorbed-energy diagnostic, NOT proof of complete energy closure.
Additional CR ionization is zero; a CR energy density is not silently
converted to an ionization rate. Translational gamma=5/3 omits molecular
rotational/vibrational heat capacity. These limits preclude claiming an
unrestricted production chemistry/RT model.
At the first live checkpoint the receiver used primary photoelectron
heating without FS2010 secondaries. This connection is now implemented as
described in the continuation below; the earlier 52 assertions alone did
not validate it. Broader hard-spectrum qualification remains a separate
physical claim, not implied by the added receiver.

## Live execution checkpoint

The ABI3 MPI2/OMP2 fresh run in `.dust-chimes-live.GYestM/live-v2` completed
four coarse steps (132.91 s). Restart in `restart-v2` then exposed an OpenMP
stack overflow: disassembly shows the NVECTOR=500 `unsplit` frame reserves
0x6e81e300 bytes, exceeding the configured 512 MiB worker stack. A fresh
run happened to schedule its one batch onto the unlimited-stack master.
The CHIMES-specific default is now NVECTOR=32; the frame is 0x071294c0 bytes.
Legacy NVECTOR=500 is unchanged. A clean reduced-batch build and fresh/restart
comparison are in `.dust-chimes-small.Lw2LZT`. Both completed (fresh about
129 s, restart about 87 s). The final 565 hydro/SNRT datasets are ALL bitwise
identical: 565 equal, zero different, including all 157 species. This fixes
both the worker-stack crash and the restart-only extra-update discrepancy.
Binary: `ramses_chimes3d3d`, SHA256
`c1189aaf2ce5a105f6b5be49e7f89cd10f33c38076e45d6ab4ab606aa11ee720`.
The doubled suffix is only an artifact name from passing EXEC with `3d`;
the documented build command uses EXEC=ramses_chimes. Each data HDF5 file
is 56,068,368 bytes (the directory also contains a separate cooling dump).
The same build passed all 52 native assertions; the GUI suite passed 36
tests with one no-display skip. Logs: `native-unit.log`, `gui-tests.log`.
A separate clean no-CHIMES build in `.dust-chemistry-baseline.T3J78m`
also passes the original 33 native assertions with NVAR=30/NVECTOR=500.
Independent final-dump sums give maximum relative nucleus error 7.22e-16,
charge/H 1.05e-16, nonnegative species, and finite state throughout. The
sum of the eleven tracked total element masses differs from rho by at most
3.17e-8 relative; this is reported separately, not called exact baryonic
closure. The molecular abundances in this Fe-rich short IC are traces:
max H2/nH=4.11e-8, CO/nH=5.50e-23. This is a wiring/restart test, not a
molecular-cloud formation or galaxy-calibration result.

An earlier ABI2 run (`live-fixed`) measured maximum relative nucleus error
7.61e-16, charge/H error 3.75e-17, and nonnegative species in its first dump.
Its restart completed but was not bitwise identical; it is NOT accepted as
an exact restart reproduction. An initial dump failure additionally exposed
an uninitialized `levelp` occupancy array in the non-FDM particle path;
only that occupancy array is now initialized for all modes.

At that checkpoint, still unimplemented within the approved scope: separate metallic Fe grains,
PAH stochastic heating, hot-material/sublimation, and conservative live
dust relative transport/drag/radiation force. Downloaded PAH data cover only
to 0.001 micron, not the 2--10 keV source bin. The existing DL01 material
table ends at 300 K and cannot support a claimed sublimation calculation.
No substitute silicate opacity is assigned to metallic Fe; no dust-only
velocity kick is called conservative relative transport. These items have
not been removed from scope or marked complete. No commit/push was made.

Additional primary-source check: DustEM 4.3 public archive downloaded from
https://www.ias.u-psud.fr/DUSTEM/dustem4.3_web.tar.gz, SHA256
`85c772a3478fe0ec7b10c120342e5ee19eae495487780f3b648b5c7a4e2e78d2`.
Its DL07 PAH heat-capacity data are available, but its supplied optical
wavelength grid begins at 0.04 micron (about 31 eV) and the package has no
separate metallic Fe table. It does not remove the full-nine-group optical
data gap. The archive is retained privately, not redistributed or called a
live dust implementation. No external author contact was attempted.

## Continuation: photoelectron secondary ionization

Implemented within the actual native CHIMES RHS, not as a post-run analysis.
The same pinned FS2010 interpolation used by the original H/He receiver
partitions the kinetic energy of photoelectrons from atomic base and
Auger-producing channels. H/He ionizations and electrons enter creation and
destruction rates; the identical energy share is subtracted from primary
heating. Secondary reactions do not change the primary photon opacity.
This also explicitly uses the admitted base-channel photoelectron energies
for both heat and secondary rates, avoiding inconsistent energy averages.

FS2010 assumes a primordial atomic medium and matched HII/HeII fractions.
The native extension limits HI/HeI/HeII channels when those actual atomic
populations fall below the table's reference state. HI excitation is limited
the same way; unused channel energy becomes heat. This continuous cap stops
negative/absent targets without breaking molecules to create atoms. It is
an explicit comparison approximation, not a recalculated molecular cascade
or arbitrary-metallicity FS table. Molecular photoionization channels keep
their previous heating treatment. Excitation energy is escaping radiation,
not silently deposited in dust IR. Source:
[Furlanetto & Stoever 2010](https://arxiv.org/html/0910.4410v1).

The default external library is preserved at `chimes/` (ABI3); the patched
ABI4 copy is `.dust-extension.AOz7mU/chimes-fs2010/`. The maintained patch
now reproduces ABI4 and its reverse-application check passes. Chemistry
restart identity is version 3, so old physical closures cannot silently
resume under the new energy partition. FS table loading occurs before
fresh hierarchy initialization or restart identity validation, not lazily
inside OpenMP callbacks. Existing SNRT HDF5 identity also binds its FS table
manifest. No new namelist switch or added default activation is introduced;
the mkrun README has been updated to describe the actual receiver.

`.dust-chimes-secondary.H4xQhx/native-unit.log`: 58 assertions pass.
In a 100 K pure-H one-zone test, the ninth photon group lost
9.436063237444614e-7 photons/cm3, HII increased by 6.283875055131066e-5/cm3,
and thermal energy increased by 6.450779764298048e-16 erg/cm3. Multiple
ionizations per absorbed photon are therefore present in the integrated
reaction solver; thermal plus ionization energy remains below the removed
photon energy. Other new checks exercise the atomic FS reference state,
99.9/100.1 eV continuity, unavailable atomic targets, and nuclei/charge.
The integrated MPI/restart test is pending; this section is not its result.

### Cached secondary receiver execution

The native suite now passes 60 assertions (`native-cached-unit.log`), including
vanishing He targets and four concurrent hard-photon cells. The table's
primordial He/H normalization is fixed at 0.24/(4*0.76); dividing by the
actual He abundance would incorrectly retain finite He ionization as He
vanishes. For the hard pure-H test, cached HII gain is
6.283875055131061e-5 and thermal gain 6.450779764298018e-16 in the units above:
the cached and direct calculations agree to roundoff.

The first ABI4 MPI attempt in `live/` was intentionally terminated during
its first illuminated chemical step after several minutes. It did not
complete and is not accepted evidence. Repeated fixed-energy FS table
searches in each RHS evaluation were replaced by read-only initialization
of the six raw FS quantities on their 14 xi nodes. Runtime interpolates
these raw values before the original partition normalization and target
limits; it does not interpolate pre-normalized fractions or change the
physical closure. Even a spectrum without ionizing channels initializes
the xi axis. The old binary is retained as `ramses_chimes_fs_uncached3d`.

The rebuilt executable `ramses_chimes_fs3d` has SHA256
`76a7e67214ce9f242f51eddf384e245514d33d310f8089001f4f377819adc02f`.
Fresh/restart tests use new `live-cached/` and `restart-cached/` directories
under `.dust-chimes-secondary.H4xQhx`. MPI2/OMP2, NVECTOR32, NVAR187,
NENER1, 8^3 periodic cells, four steps; no AGN/sinks/cosmology. Output policy
in both effective namelists: noutput=1, aout=[2], tout=[1e30], foutput=2,
fbackup=1000000. Three new dumps total, approximately 201 MiB, versus
155 TiB free before launch. Both executions completed successfully: fresh
about 188 s and restart about 119 s. The first three illuminated chemical
steps took 53.005, 61.331 and 60.872 s. Their final 565 hydro/SNRT datasets
are all bitwise equal (zero differences), with finite values throughout.
Final active-leaf species are nonnegative; independent nuclear sums have
maximum relative error 7.27e-16 and charge/H error 9.23e-17. Total element
mass versus rho differs by 3.175e-8 relative, reported separately from the
chemical nucleus closure. H2/nH <=4.04e-8 and CO/nH <=1.67e-23: this is still
a deliberately Fe-rich wiring test, not molecular-cloud calibration.

Separate Fe-data check: the authors' [Semenov optical-constant page](https://semenov.www3.mpia.de/Opacities/RI/new_ri.html)
provides a pure-Fe n/k table (`ironk.lnk`, density 7.87 g/cm3), but its
wavelength range starts at 0.1 micron. It cannot alone cover the nine SNRT
groups through 10 keV. This is a candidate long-wavelength input, not an
admitted replacement for silicate or a completed Fe implementation.

### Hot harmonic material (not complete sublimation)

Native `dust_composition_material.f90` now evaluates the same DL01 normalized
Debye mode integrals above 300 K through 3000 K. Existing cold interpolation
is unchanged. Sixteen-point Gauss-Legendre integration in the smooth warm
domain agrees with independent 256-node quadrature at 300, 1000 and 3000 K;
the existing dust test checks continuity, mixtures, invalid-domain rejection
and byte-identical cold identity. No new numerical test framework or runtime
Python was introduced. This extends harmonic internal energy, not solid
phase stability: latent heat and grain evaporation are still unimplemented.
D03 dielectric functions are still frozen at 20 K, not temperature-dependent
optical data. Separate grain temperatures and PAH modes remain absent.

Injection/domain validation and live per-cell material curves use the new
domain. Existing common Fortran/OpenMP/CUDA material ABI accepts the curve
unchanged; the hot path has been tested on Fortran/OpenMP, not newly on CUDA.
The existing backend test's sixth D03 case starts at 1000 K, remains above
300 K, conserves material-plus-IR energy, and has 1.0205e-16 relative backend
difference. Logs: `hot-unit.log`, `hot-backend.log`, `hot-contract-backend.log`
under `.dust-chimes-secondary.H4xQhx`. GUI tests: 36 pass, one display skip.

Cold checkpoint material identity remains version 1 and bitwise unchanged.
For contracts extending above 300 K, the same-sized version-2 identity binds
the curves sampled across 5--3000 K. The normal material/IR contract identity
also binds its actual knots. Default namelist model names and generated
5--300 K runs are unchanged; mkrun documentation now states the distinction.
No hot contract is silently injected into a run using an old binary.

Using the existing offline builders, generated `hot-material.json` (100
knots) and `hot-contract.nml` (162 combined knots, full original IR grid).
The contract's SHA256 is
`8bc8c51d07118b2cd63121e4e379f8c73ba01c836f2cefe02a1d975e6fdf706c`.
It retains the nine-group ledger, D03 receiver, and gas/dust accommodation.
The material builder used numpy-only `.dust-optics-venv`; thermal export uses
the already installed `simulation/snrt/.venv` JAX environment. No package
installation or system environment change was needed.

Full native binary `ramses_chimes_hot3d`, SHA256
`90d66a9d75d0843a5430d518fc3ee880a13dd1b339f9d5b2ff7a9ff2b6f907c7`,
was built successfully on lageunha. Tests in `live-hot/` and `restart-hot/`
use MPI2/OMP2, 8^3 cells, two steps, SNRT/DUST/CHIMES CPU, NENER1/NVAR187,
NVECTOR32. Both effective namelists: noutput=1, aout=[2], tout=[1e30],
foutput=1, fbackup=1000000; three new dumps (about 201 MiB), 154 TiB free at
launch. This run checks the hot table's live connection/restart; the above
one-zone IR test supplies the actually-hot grain state. Both MPI runs
completed; all 565 final hydro/SNRT datasets are bitwise equal and finite.
The 487-value `dust_material_composition` identity has version 2 and also
matches exactly between fresh/restart. The maintained CHIMES ABI4 patch
passes reverse-application validation; `git diff --check` is clean.

Still remaining: metallic Fe grains; PAH stochastic heating; actual
sublimation/latent-energy exchange; conservative dust/gas relative
transport, drag and radiation momentum. None is marked complete by the
hot internal-energy extension. No commit/push or external audit this turn.

### Graphite sublimation implementation and coupled-cell fix

Implemented opt-in `dust_sublimation='gd89_graphite_bulk_v1'`, default `none`.
This is graphite only, not MgFeSiO4 sublimation. Native Fortran uses the
vacuum bulk graphite law in [Waxman & Draine section 3.1](https://arxiv.org/html/astro-ph/9909020),
with `nu=2e14 s^-1`, `B/k=81200 K`, and carbon atomic mass 12.011 amu.
Two fixed-radius bins use backward Euler for `lambda=3*abs(da/dt)/a_bin`;
mass and common grain temperature are solved together. This closure does
not conserve resolved grain number/radius during evaporation. The source's
Mg2SiO4 coefficients were not assigned to the different live MgFeSiO4 solid.

Native phase accounting conserves `Egas+Edust+L*(C_total-C_solid)`, where
`L=5.62096715e11 erg/g`. Evaporated carbon returns to neutral gas, and its
explicit effusive `2*k*T/m_C` kinetic energy heats gas. The same latent
reference applies to ordinary graphite growth, sputtering and SN
destruction. Phase energy is derived from existing advected mass carriers;
no extra passive scalar or runtime Python is involved. The material exchange
also now preserves unchanged composition energy exactly and uses the IR
receiver's `U(log T)` interpolation when composition changes.

Namelist reader, generator, mkrun choice, existing GUI tests and restart
identity are connected. mkrun requires an explicit new binary and matching
contract through `SNRT_DUST_SUBLIMATION_BINARY` and
`SNRT_DUST_SUBLIMATION_CONTRACT`. The six-value HDF5
`dust_sublimation_values` attribute binds the model and rejects switching
it on/off at restart. The selected hot contract still uses frozen D03
dielectric functions and a common bulk grain temperature.

Existing native dust smoke test: a one-gram mixture initially at 3000 K,
dt=1000 s, loses 0.0250399249 g graphite and reaches 2072.31398 K.
Gas heat=7.18412026e8 erg; latent energy=1.40748596e10 erg.
It checks mass/size selectivity, unchanged silicate, energy closure,
small-dt analytic-rate agreement, zero dt and invalid-input rollback.
Initial logs: `.dust-chimes-secondary.H4xQhx/sublimation-unit*.log`.
GUI: 38 tests, 37 passed and one display-dependent skip.

The first live run is **not a pass**: `live-sublimation/` and its restart
reach the next chemistry trial and reject it. The first dump nevertheless
shows actual graphite mass/rho falling from 0.005 to 0.004686744059910256
(6.2651% loss), with total carbon/rho unchanged at 0.01. Diagnostic
restarts (`restart-sublimation-diag`, `-codes`, `-nuclei`) identify a trace
oxygen nucleus error after stellar injection:
expected O/H=9.6442110780925215e-9, actual=9.6456856502476338e-9;
relative discrepancy 1.529e-4 exceeds the unchanged 1e-5 acceptance gate.
This is a coupled chemical integration failure, not a successful full run.
Original dumps/binaries/logs are retained.

The receiver now retries only nucleus/charge-failing cells, at most twice,
from their untouched original input, reducing relative, absolute and
explicit solver tolerances by 100 each time. Conservation acceptance is
not relaxed; exhausted retries still reject without publishing the cell.
Successful first attempts retain their original calculation. The chemical
restart identity advances to version 4, so the final verification starts
fresh rather than modifying a version-3 checkpoint. Final live result follows.

This remains a split graphite comparison, not completed dust physics:
simultaneous radiation-heating/evaporation convergence, vapor backpressure,
silicate sublimation, separate metallic Fe, PAH stochastic heating and
conservative relative dust/gas dynamics are not claimed complete.

Final graphite-coupling verification: **passed**. Binary
`.dust-chimes-secondary.H4xQhx/ramses_chimes_sublimation_retry3d`, SHA256
`3311f3b12782e737a9ce8e2db4978b25c237be7a621e1aa9be487386abc20d92`.
Effective namelists are under that directory in
`live-sublimation-retry/physical.nml` and
`restart-sublimation-retry/physical.nml`. Fresh/restart both use an 8^3
uniform grid, two steps, MPI2/OMP2, CPU SNRT/DUST/CHIMES, NENER1/NVAR187,
NVECTOR32; hydro, gravity, stellar formation/channel feedback and CR are
enabled, AGN/sinks/cosmology disabled. Initial carbon/rho=0.01 and graphite
mass/rho=0.005 at 3000 K; ordinary growth/sputtering/coagulation/shattering
and SN dust destruction are disabled to isolate evaporation (stellar
sources remain active). Each output schedule: noutput=1, aout=[2],
tout=[1e30], foutput=1, fbackup=1000000; three new dumps, approximately
201 MiB total, 153 TiB free at launch. No original output was overwritten.

Fresh/restart both report `Run completed`, with no rejection or
nonconvergence; measured program timers are 59.898/50.601 s. All 565 final
hydro/SNRT datasets are bitwise equal and finite. The sublimation,
487-value hot material and version-4 chemical identities also match exactly.
Initial graphite loss is 6.2651188%; the final chemical nucleus closure
error is at most 1.2987e-15 relative and charge/H error 1.9562e-17, with
minimum species density zero. Total element mass/rho differs by 1.3588e-8
relative, explicitly separate from the chemistry conservation result.
The full live run has stellar sources, gravity and radiation escape: its
ordinary RAMSES `econs` is not a closed-box sublimation energy test.

The existing native thermochemistry test passes all 60 assertions; the
existing dust test passes including the sublimation/latent-energy case.
Logs: `sublimation-final-chem-unit.log`, `sublimation-final-dust-unit.log`,
`sublimation-retry-build.log`, and both live `run.log` files. Generator/GUI
tests: 37 pass, one display skip; `git diff --check` clean. No commit/push,
new external audit, or additional simulation left running by this bundle.

## Next implementation and bounded production-readiness accounting

Operator instruction: proceed to the next implementation and state what
remains before production readiness. Do not count the historical comparison
closeout as qualification of this newly approved extended model.

The already approved unfinished physics groups are, without new scope:

1. Complete sublimation: silicate mass/element/phase-energy coupling and the
   already recorded radiation-heating/evaporation splitting limitation.
2. Separate metallic Fe: real optical/material inputs and conservative live
   species/source/transport/restart connections, not silicate opacity reuse.
3. PAH stochastic heating: source-supported absorption/emission and its
   native radiation/energy connection, not a renamed equilibrium temperature.
4. Relative dust/gas dynamics: transported state, conservative drag and
   radiation momentum, with matching source/restart/backend behavior.

After these, combine the necessary qualification of the selected physical
profile into one closeout effort: source/input applicability, coupled budgets
and relevant convergence, then its actual MPI/AMR/backend/restart operating
conditions. This is not a fifth physics feature or a per-step audit ladder.
Some individual checks already exist; use them, not duplicate them. The
remaining four groups differ substantially in size; there is no defensible
percentage or calendar completion estimate yet. Relative dynamics is a major
remaining change, not a minor finishing patch. Advanced dust extensions are
not universal prerequisites for every restricted production model, but remain
in the user's explicitly approved expanded scope. No deferral is inferred.

Historical checkpoint before the implementation below: silicate source/model selection; no new silicate
runtime option has been enabled. [Xu et al., preprint table 2 and section II.3.1](https://arxiv.org/html/2509.11036v1)
provide MgFeSiO4 surface mass-flux fits, matching the live solid's elemental
formula: (ln S0, E/k) = (20.06,84780 K), (16.29,74420 K),
(22.23,90590 K), for (100), (010), (001), with S in g cm^-2 s^-1.
These are crystalline kinetic-rate fits, not automatically latent heats for
the live amorphous/bulk optical approximation. The rendered preprint eq. 1
contains a double minus; check the published equation/figures before code
admission. Surface weighting, evaporation products and phase-energy closure
must be explicit. Do not silently transplant Mg2SiO4 or graphite constants.
This is source/model work on existing item 1, not completed implementation or
a new standalone research programme. No simulation or external audit launched.

## Unified-bundle progress: live olivine and neutral drag (2026-09-09)

The operator requests reporting each result followed by continued work, not
an approval pause. The four-item bundle remains in progress. No completion
claim for radiation-coupled evaporation, separate Fe/PAH, or relative transport.

`gd89_xu25_olivine_v1` is now an explicit native sublimation selection, retaining
the graphite path and adding MgFeSiO4 loss from both silicate bins. Namelist,
mkrun/GUI, element return, gas heat, phase energy, ordinary mass-process phase
accounting and restart identity are connected. Defaults remain unchanged.

Physical assumptions are deliberately bounded:

* Xu preprint table-2 rates use equal exposed crystalline faces and the physical
  negative Boltzmann exponent. The rendered equation's double minus is not
  reproduced. This is a kinetic approximation, not amorphous optical validation.
* Thermodynamics uses an ideal equal-molar forsterite/fayalite mixture:
  Hf(298.15 K)=-2173.0/-1478.2 kJ/mol from RH95 (printed pp. 33/32).
  Neutral atomic Mg/Fe/Si/O Hf=147.10/416.3/450.00/249.18 kJ/mol. Subtract
  ideal monatomic gas enthalpy 2.5*7*R*T and add the DL01 silicate sensible
  energy at 298.15 K to define the phase reference: 2.21423908e11 erg/g.
  This is an explicit ideal-mixture/EOS approximation, not measured amorphous
  MgFeSiO4 latent heat. Excess mixing enthalpy, melting, electronic excitation,
  gas SiO/O2 products and vapor back-pressure are not modeled here.
* Congruent Mg+Fe+Si+4O vapor carries effusive heat 2*k*T per atom; existing
  integer atomic-mass element-carrier conventions remain unchanged. The
  material molar mass is 172.231 g/mol. These conventions are not identical.
* Closed-cell energy is Egas+Edust-LC*M_Csolid-LS*M_Ssolid, up to the constant
  reference fixed by total nuclei. Native pure-silicate test loses 0.0483568427 g
  from 1 g at initially 3000 K in 1000 s and closes to 6.38e-15 relative.

Sources: [Xu et al. table 2](https://arxiv.org/html/2509.11036v1),
[RH95](https://pubs.usgs.gov/publication/b2131),
[NIST Mg](https://webbook.nist.gov/cgi/cbook.cgi?Mask=35&Source=1965RIS381&Units=SI),
[NIST Si](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440213&Mask=20F&Units=SI),
[NIST O](https://webbook.nist.gov/cgi/cbook.cgi?ID=C17778802&Mask=29),
[NBS tables, Fe table 41](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf).
RH95 PDF retained at `.dust-extension.AOz7mU/robie_hemingway_1995.pdf`.

Real MPI evidence: `.dust-chimes-secondary.H4xQhx/live-olivine/physical.nml`
and `restart-olivine/physical.nml`; same-directory `run.log` and HDF5 dumps.
Binary `ramses_olivine3d` SHA256:
`f6044d883f77d0a07e98293ceb7b63ec7e6ad2e4a3d44cc5e84a58a222c36391`.
CPU SNRT/DUST_LIVE/CHIMES, NENER=1, NVAR=187, NVECTOR=32, MPI2/OMP2,
8^3 grid, stellar feedback/SF/CR enabled, no AGN/cosmology, two steps.
Initial silicate/rho=0.005, 3000 K, total metals/rho=0.01 in matching
Mg/Fe/Si/O proportions. Other grain mass processes off for this comparison.
Hot contract `.dust-chimes-secondary.H4xQhx/hot-contract.nml` unchanged.
For both: noutput=1, aout=[2], tout=[1e30], foutput=1, fbackup=1000000;
three new dumps, about 201 MiB expected, 152 TiB free before launch.
Neither run overwrote earlier comparison output. Timers: 58.703/47.597 s.

Both runs completed without rejection. All 565 final hydro/SNRT datasets
are finite and bitwise equal. Chemical (324), material (487), graphite (6)
and new olivine (15) identities match. Initial-dump silicate/rho becomes
0.004577632489166732 (8.4473502% loss). Final nucleus relative error
1.04954e-15, charge/H error 2.17155e-17, minimum species density zero.
Element mass/rho error 1.35340e-8 is recorded separately. Gravity, stellar
sources and escaping radiation make ordinary run econs unsuitable as a
closed sublimation energy budget. Radiation/evaporation are still split.

The native `dust_drag.f90` now supplies a transactional multicomponent BE
momentum/heat update and cgs neutral Epstein times with supersonic correction.
It distinguishes single-grain stopping time from mixture relaxation time,
uses isothermal sound speed, and rejects radii beyond 9*lambda/4. Caller
supplies the collision mean free path; no unverified collision law is hidden.
Existing dust smoke tests cover momentum/kinetic+heat closure, stiff limit,
rollback, size scaling, linear/supersonic drag and Stokes-domain rejection.
Logs `unified-olivine-unit*.log`, `unified-olivine-build.log`, and
`unified-drag-unit*.log` are in the build directory. Existing generator/GUI
tests: 37 passed, one display skip.

Source: [Laibe & Price 2012b, Epstein law](https://users.monash.edu.au/~dprice/pubs/dust/LaibePrice-Dust-II.pdf).
Old `patch/mrn-scratch` and `patch/tsc-scratch` contain particle drag, but
hard-code dust fields beginning at 9 and use a different particle ABI.
They are not active-runtime drop-ins; VPATH was not changed. The new primitive
is currently linked into the existing native test, NOT live dust transport.
Do not count it as completion or expose a live selector yet.

## Continued work after the olivine result: native stochastic excitation

`patch/lagRamses/dust_stochastic.f90` now implements the thermal-continuous
stationary recursion (DL01 eqs. 53-54), with arbitrary upward jumps and
nearest-lower cooling. Log-space populations AND energy fluxes retain the
strong-heating limit even when a tiny probability itself underflows. This
is not the exact-statistical/discrete-emission solver or a transient solver.
Its probability normalization is not an emission-energy renormalization.

The photon-group adapter accepts per-grain absorption rates and physical
photon energies, preserves the mean deposited energy on a nonuniform grain
energy grid, and explicitly returns above-grid photon count/energy for every
initial state. Self-jumps are omitted from the generator, not from the
incident photon accounting. This center-interpolation is a declared numerical
closure, not the paper's finite-bin transition integrals. A finite-grid
stationary solution must not silently discard its probability-weighted tail.

DL01 generic PAH normal modes and canonical harmonic U(T)/Cv(T) are native
as well: C-C skeleton and C-H bending/stretching are separate; NC and NH are
explicit inputs, not a hidden change to the hydrogen inventory. The model
uses the source's spectroscopic wavenumbers. The C24H12 case has 102 modes;
the zero-T limit, derivative and high-T mode-count limit pass. This is an
excitation model, not proof of molecular survival or a microcanonical
temperature prescription at the lowest excitations.

Source: [DL01 sections II, VIII.2](https://arxiv.org/html/astro-ph/0011318).
The existing dust smoke also passes the analytic three-state distribution,
dark/strong-field cases, absorbed/emitted energy equality, nonuniform-grid
photon moment and explicit overflow checks. Log:
`.dust-chimes-secondary.H4xQhx/unified-stochastic-unit.log`; build log has
the same prefix with `-build`. No separate test framework or external audit.

This is **not a live PAH feature yet**. Physical cross-section/charge and
high-energy input coverage, independent H/C mass and source carriers,
relaxation validity, IR radiation, transport/restart and namelist selection
are not connected. No PAH switch was exposed and no restricted runtime
profile was silently broadened.

Input reuse investigation: DustEM 4.3's temperature-distribution code uses
module-global allocatables and file-supplied radiation; its comments explicitly
refer to subsequent emitted-power normalization. It is not a thread-safe
transactional runtime drop-in. The new native recursion was implemented from
the published method, not copied from that archive. The DL21 Harvard API is
reachable, but lists precomputed spectra for selected stellar fields, not a
transition table for arbitrary evolving SNRT radiation. No multi-GB spectra
collection was downloaded. The neutral/ion PAH tables already downloaded
extend only to about 1.24 keV; do not extrapolate these to 10 keV silently.

Artifact SHA256s:
* DustEM archive: `85c772a3478fe0ec7b10c120342e5ee19eae495487780f3b648b5c7a4e2e78d2`.
* PAHneu_30.gz: `d1d0f8371b730ab473ab90cb98def2c6367507ae7a34c5a1214eab33a24fe09c`.
* PAHion_30.gz: `0a7b7880ad5916daec1a76d14310933d09b46b79b30d84624f85c9322587b836`.

Fe data are still incomplete: the available Henning/Semenov iron index starts
at 0.1 micron and cannot cover the nine groups. DH13 appendix B describes a
microwave-to-X-ray reconstruction, but its underlying full Fe table has not
been obtained. Astrodust tables with Fe inclusions are NOT pure Fe tables.
HD17 describes a low-temperature Fe heat model but uses laboratory data for
higher temperatures and phase transitions; do not extend that low-T Debye
formula through melting and call it physical data. Sources:
[DH13 appendix B](https://arxiv.org/html/1205.7021v2),
[HD17 thermodynamics](https://arxiv.org/html/1611.08607).

Current handoff status: the two olivine simulations and native builds/tests
have finished; no simulation/audit was left running. The approved bundle
is still unfinished. Results were reported and further native implementation
followed without an operator-approval pause. No commit or push in this turn.

## Coupled radiation/material/evaporation continuation

`dust_sublimation='gd89_xu25_olivine_rt_v1'` now selects a native coupled
IR material/phase solve, rather than applying evaporation before primary RT.
It requires D03 optics, CHIMES, SNRT_RT_ENABLE=1 and the exchange-enabled
version-4 hot IR contract. Namelist validation, mkrun/GUI, conservative mass
carriers, CHIMES neutral-vapor return, gas heat, halo updates and restart are
connected. The new material closure is CPU/OpenMP only; auto retains the
existing hybrid transport/scattering kernels, explicit CUDA-only material
selection rejects. The default and the older split modes are unchanged.

The scalar residual is
`Ed'-Ed + phase + effusive_heat + dt*(net_IR_power-heating) - Qcollision=0`.
Final BE mass loss and sensible U(T) occur in that same solve. The gas
transfer returned to the IR ledger is `Qcollision-effusive_heat`; the IR
trial temporarily accounts for `Ed'+phase`, then stores only sensible Ed'.
CHIMES reconciliation uses the accepted chemistry trial, not stale uold.
All persistent changes wait for the existing all-rank RT transaction commit.

Two actual MPI failures were preserved in `live-coupled-sublimation/` and
`live-coupled-sublimation-2/` within `.dust-chimes-secondary.H4xQhx/`.
The first solve used a temperature root, ill-conditioned near the background
at very long timesteps. It was replaced by a NET-emission-power root with
incremental band powers, matching the existing material solver strategy.
The second failure exposed a residual tolerance scaled only by tiny stored
dust energy when the energy transferred from gas was much larger. The root
now scales with the actual terms cancelled, not the full unused gas reservoir.
The new gas-dominated native case has Q/U=3.17960e5 and closes to 1.86e-12.
The isolated irradiated olivine test loses 0.369285589 g from 1 g and closes
to 4.78e-13; its 8/16/32 interval results approach one another. These are
bounded native tests with a synthetic positive two-band emissivity, not
observational calibration or full D03 timestep qualification.

The resulting **un-subcycled comparison** completed fresh and restart runs:
`live-coupled-sublimation-3/physical.nml` and
`restart-coupled-sublimation-3/physical.nml` in that build directory.
Binary `ramses_coupled_sublimation33d` SHA256:
`6c5094eae8237e962d0c6d01831f458e77f51aeeedead240c8f8a32e7aa0d913`.
CPU SNRT/DUST_LIVE/CHIMES, NENER1/NVAR187, MPI2/OMP2, 8^3 grid;
same 3000 K silicate initial state, SF/stellar feedback/CR, no AGN/cosmology.
Both effective namelists: noutput=1, aout=[2], tout=[1e30], foutput=1,
fbackup=1000000; three new ~67 MiB dumps, ~201 MiB total, 147 TiB free
before launch. Previous outputs were preserved. Timers: 68.975/50.322 s.
565 final hydro/SNRT datasets are finite and bitwise equal. Chemical (324),
material (487), sublimation (6, version 2) and olivine (15) identities match.
Final nuclei relative error 7.19426e-16, charge/H 2.45570e-17, minimum
species zero. Total element-mass/rho error 1.35340e-8 is separate.
IR balance errors: 2.0610e-12 and 6.9517e-10, within the existing 1e-9 bound.

This comparison also revealed the accuracy limitation of a single stiff
endpoint: initial hot silicate/rho stays exactly 0.005 after the first
step, missing the brief evaporation pulse while the grain cools. It is
NOT qualified merely because its energy/restart checks pass. Native local
adaptive integration is being added next; version-2 results remain preserved
as evidence of the un-subcycled method, not the current solver's validation.

Physical scope stays explicit: common grain temperature, vacuum crystalline
kinetics, ideal-phase thermodynamics, frozen D03 dielectric data, fixed bin
radii, lagged primary/IR opacity and gas Cv. Coupling does not remove those
assumptions, provide PAH survival chemistry, or complete Fe/relative dynamics.

### Resolving the early hot transient

The new `dust_radiative_sublimation_evolve` performs native adaptive BE step
doubling inside the existing IR callback. Default relative tolerance is 1e-4;
a 25% state-change pre-limit prevents both coarse/half-step endpoints from
skipping the same initial pulse. The accepted answer is the two positive
half steps, not Richardson-extrapolated mass. All four bin masses, sensible
energy, gas transfer and time-integrated spectral band energy enter the
error estimate. At most 4096 trials are allowed; incomplete integrations
return the original state and zero transfers. The IR receiver sees time-
averaged emission, integrated gas/phase exchange and the final sensible U.
The sublimation identity is now version 3. Defaults/split modes are unchanged.

An initial fractional-state-change-only implementation failed the native
tightening-resolution assertion (`coupled-sublimation-subcycle-unit.log`,
STOP 162). It was replaced by actual step doubling, not by relaxing the
assertion. For a 1-g silicate grain ensemble initially at 3000 K, the positive
synthetic two-band test spans 1e15 s and resolves the early evaporation:

| Local tolerance | Accepted intervals | Remaining mass [g] | Relative energy residual |
|---|---:|---:|---:|
| 1e-3 | 388 | 0.950175536 | 1.75e-11 |
| 1e-4 | 1172 | 0.950052084 | -4.13e-12 |
| 1e-5 | 3654 | 0.950006246 | -1.25e-11 |

The differences decrease; this is a bounded native time-convergence test,
not a universal accuracy estimate. Logs `coupled-sublimation-adaptive-unit.log`
and `coupled-sublimation-adaptive-opt-unit.log` preserve the result. The
latter also covers failed-integration rollback. The zero-contribution high-T
tail of band summation is skipped; nonzero arithmetic/results are unchanged.

The first adaptive live comparison uses `ramses_coupled_sublimation43d`,
SHA256 `04ced9df43a10fe0d6ea10a0fb40d2e12c691fbd3349532a89a18fa6475729d2`,
in `live-coupled-sublimation-4/` and its planned `restart-coupled-sublimation-4/`.
It retains the old IR 0.5 fixed-point damping for comparison. The run is
expensive: it repeats the adaptive hot transient while converging IR
reabsorption. These directories must not be confused with the earlier
un-subcycled, already completed `*-3` pair.

To reduce those repeated solves, the new material path enables direct
fixed-point updates only where every spectral band's local reabsorbed
emission fraction `1-response` is below 0.1. Other cells retain 0.5 damping;
the old runtime selections do not request this option. Existing local and
global energy/residual acceptance thresholds are unchanged. The native
backend smoke confirms fewer iterations and agreement with the old damped
reference at 9.14e-12/9.10e-12 for the two material modes. Log:
`coupled-sublimation-thin-unit.log`.

Optimized live binary `ramses_coupled_sublimation63d` SHA256:
`8eac244aceaadaca8477c21c666d6237d65b2d3c4d6e87f37f4fb12e7c215908`.
Its fresh/restart directories are `live-coupled-sublimation-6/` and
`restart-coupled-sublimation-6/`, within the same build directory. MPI2/OMP2
with `I_MPI_PIN_DOMAIN=omp`; the earlier run's rank affinity was one logical
CPU despite OMP_NUM_THREADS=2. Both new pairs use the same audited physics
and output schedule: noutput=1, aout=[2], tout=[1e30], foutput=1,
fbackup=1000000, three ~67 MiB dumps per pair. Free space at their respective
launches was 146 and 145 TiB. All earlier output remains preserved.
The `*-4` and `*-6` fresh runs were explicitly stopped (SIGTERM to their
verified mpiexec processes) before completion; neither produced a snapshot
or started its restart run. Their logs/binaries remain. Inspection of actual
optical coefficients gave max tau=0.23250046 and max reabsorbed fraction
0.10774105: the initial 0.1 acceleration threshold did not apply to these
cells, so both were repeating the old damped iteration. The current solver
uses the still conservative 0.25 threshold, with unchanged final acceptance
tests. Do not mark these interrupted runs as passed or physics failures.
The completed adaptive pair is **`live-coupled-sublimation-7/` and
`restart-coupled-sublimation-7/`**. Both use binary
`ramses_coupled_sublimation73d`, SHA256
`aa46dd8b787edefb9656e8f731a9a8bd13960c4994a48ae242c27b514183c82c`.
The effective namelists retain the same physics/schedule above; free space
was 145 TiB before this launch. MPI2/OMP2 with I_MPI_PIN_DOMAIN=omp.
All previous files were preserved. Native backend/adaptive/rollback log:
`coupled-sublimation-quarter-unit.log`; build log:
`coupled-sublimation-quarter-build.log`. No VPATH ordering change.

Both runs completed without rejection: 144.482/52.534 s. The hot first
coarse step took 86.66 s and the next 52.91 s. All 565 final hydro/SNRT
datasets are finite and bitwise equal across fresh/restart. Chemical (324),
material (487), sublimation (6, version 3), olivine (15) attributes match.
First-dump silicate/rho=0.004969624381738415, corresponding to 0.607512365%
evaporation from 0.005. Final silicate/rho range:
0.004969606175538519--0.004969613617090501. Final nuclei relative error
6.03394e-16; charge/H error 2.09930e-17; minimum species zero. Total element
mass/rho error remains 1.35340e-8 and is not hidden in the nucleus test.
IR balance errors: 4.6216e-12 and 1.0116e-10, below the unchanged 1e-9
acceptance threshold. These are not closed whole-run econs measurements.

The raw split-evaporation comparison had 8.44735% loss, whereas the long
single coupled endpoint missed the pulse entirely. Neither is a reference
truth for the adaptive radiation-coupled method. The differences motivate
resolving competition between hot radiation cooling and evaporation; the
0.6075% result is for this deliberately hot test, not a universal dust loss.
Physical optical-coefficient lag, common T, fixed radii, ideal-phase/vacuum
kinetics and outer timestep/mesh qualification remain explicit limits.
Existing generator/GUI tests: 37 pass, one no-display skip. mkrun README
now describes adaptive coupling and no longer says 'No latent heat' when
sublimation is selected. The unrelated generator deletion stays untouched.

No coupled-sublimation verification jobs remain running. No commit/push or
external audit was performed. This completes the live radiation/evaporation
connection and its bounded integrated check, **not the entire approved
four-item bundle or blanket production qualification**. Fe, PAH and relative
motion still require their live connections; no new operator approval is
being requested for their already-approved implementation.

### Fe input continuation in parallel

Obtained the actual [Werner et al. 2009 NIST-hosted paper](https://www.nist.gov/document/jpcrd3820091013ppdf),
not only an index database. Its table 13 has both DFT and REELS Fe dielectric
data from 0.5 through 70.5 eV; the experimental real dielectric part below
10 eV has large stated uncertainty. This is useful overlap data, not a
standalone complete dust dielectric function or a reason to extrapolate
the REELS fit to microwave/10 keV. Local PDF:
`.dust-extension.AOz7mU/werner_2009_metals.pdf`, SHA256
`f86156bab5f885fb3c1099fb95352092a0e73ae0574c2e0f724e31e387508eb8`.

Also obtained the official [LBNL Henke Fe atomic scattering factors](https://henke.lbl.gov/optical_constants/sf/fe.nff),
10--30000 eV, with the tabulated Fe edges. Below 29 eV, f1=-9999 is an
undefined sentinel, not a usable optical constant. The [official data notes](https://henke.lbl.gov/optical_constants/asf.html)
identify the Fe 600--800 eV edge revision. Local `henke_fe.nff` in the same
directory has SHA256
`cdcd0f4babbc2f07b6268a698a7fd1b61c920afdc7c5ca331efeb80dd1d0e7f5`.
These inputs narrow the optical coverage gap. They still require justified
joins/dispersion and finite-grain electromagnetic treatment, independent Fe
heat/phase data and the live carrier connections. No Fe switch or fabricated
full-band data was introduced by downloading them.

The [NIST-JANAF iron collection](https://janaf.nist.gov/pdf/JANAF-FourthEd-1998-Iron.pdf)
is also retained as `.dust-extension.AOz7mU/janaf_1998_iron.pdf`, SHA256
`96ec56da8d92685cf963a7382435064ef108d981cb36525b09c78fbf06213bbb`.
It is a scanned 33-page PDF; text extraction is empty, so the relevant
table was visually inspected (PDF page 7, printed 1225), not parsed as zeros.
This is the pure Fe crystal(alpha/gamma/delta)-liquid table, NOT the FeO
table at the start of the collection. It identifies the 1042 K heat-capacity
anomaly and distinct transitions at 1184, 1665 and 1809 K, with enthalpy
jumps rather than a smooth single solid curve. Its molar mass is 55.847,
whereas the current WebBook value is 55.845; a future input binding must
declare its convention. This now supplies actual phase data to investigate,
but no smooth interpolation across a latent jump, low-T extrapolation or
Fe live feature is claimed here.

### Approved continuation: native Fe phase material and reserved-element coupling

Continued in `/gpfs/kjhan/LRD_JWST`, verified origin
`git@github.com:kjhan0606/LagRamses.git`. No separate stage audit, new gate
framework, Python test generator, namelist selector or operator approval
pause was added. Existing VPATH order and unrelated generator deletion were
preserved. The four-item bundle is still the scope; no extra tasks were added.

Implemented `dust_iron_janaf_data.inc`, `dust_iron_material.f90` and
`dust_iron_radiation.f90`. The JANAF combined condensed-phase page was visually
rechecked against the individual alpha/delta, gamma and liquid pages (printed
1222--1225). In particular H(500)-H298 is **5.527**, not 5.277 kJ/mol from the
earlier provisional note. The native table preserves the three distinct
latent intervals, 0.900/0.837/13.807 kJ/mol, and no latent jump at 1042 K.
It uses the source's 55.847 g/mol convention, not the chemistry ledger's
integer mass number. Condensed pV, hysteresis and nanoscale melting shifts
are not modelled. Max admitted material temperature is 3000 K.

Low-T shape is HD17 eq.1 (Debye theta=415 K plus electronic gamma=6e-4/K),
with explicit 1.00576510 normalization to JANAF H298-H0. The exponential-tail
Debye integral is evaluated natively; no runtime quadrature/Python dependency.
U100/U_JANAF100=1.06812197; U200/U_JANAF200=1.02301335. These are admitted
approximation discrepancies, not a claim to interpolate all raw low-T points.
The material identity has 90 values binding table and approximation; it is
not yet bound to a live restart because Fe carriers are not yet activated.

The common-temperature inverse returns alpha/gamma/delta/liquid fractions
and conserves enthalpy within a plateau in C/olivine/Fe mixtures. The new
six-bin radiative callback accepts supplied band-power bases, performs
backward Euler in net emitted power with the existing native gas-exchange
helper, and explicitly resolves phase intervals. Material energy already
includes Fe transition enthalpy (no second latent-energy ledger). Optical
coefficients and masses are fixed during this solve; no Fe evaporation law
or real Fe optical table is supplied by this callback.

`dust_gas_elements` and `chimes_cell_state` now accept optional metallic Fe.
`dust_size_step_reserved_iron` reserves it in both the elemental and total
metal inventories before calling the existing four-bin C/olivine operator;
it validates the final budget and only publishes an accepted result. With
zero Fe it is bitwise equal to the old growth call. Current live callers
still omit Fe: no persistent field/source/transport/restart wiring is claimed.

Verification reused **only existing native smoke programs**, Intel ifx through
the current Makefile on lageunha, OMP2. Logs in
`.dust-chimes-secondary.H4xQhx/`:

* `iron-reservation-build.log`, `iron-reservation-unit.log`: passed source
  values, 301 temperature samples per mixture, all three latent plateaus
  and endpoints, ordinary/dilute mixtures, invalid-state rejection, elemental
  reservation, growth and zero-Fe bitwise parity; old dust tests passed too.
* `iron-reservation-backend.log`: existing backend/sublimation tests pass,
  plus new Fe radiative heating AND cooling through each phase interval,
  constant-state balance, stiff gas/trace dust and rollback. Radiation uses a
  **synthetic two-band input**, not purported Fe optical data. Trace case
  T=5.00500750 K, gas transfer=9.99456623e-13 erg/cm3,
  relative material+gas+radiation residual=4.20632703e-13.
* An initial monotonicity test rejected its generated last value slightly
  above 3000 K under optimized arithmetic (`iron-material-unit2.log`). Fixed
  the test's prescribed endpoints to exact values; material domain bounds
  were not relaxed. `iron-material-unit3.log` then passed.
* `iron-prereq-build.log`: full CPU SNRT+DUST_LIVE+CHIMES, NENER1/NVAR187
  build/link passed. Binary `ramses_iron_prereq3d` SHA256
  `e413e341b17f20a6a3f7b2748f59e804eb866f74e5b3819f2039316c75a4f076`.
  This is **compile compatibility**, not an Fe-enabled run. The dormant Fe
  material modules are linked into the smoke receivers, not into live RT.

No RAMSES calculation/dump/restart was launched in this continuation; the
previous completed coupled-sublimation MPI evidence stays preserved.
No commit/push or external audit performed; no verification job left running.
Remaining next work is still full-band Fe optical admission and live Fe
source/carrier/chemistry/transport/restart binding, then the already-approved
PAH and relative-dynamics connections. Fe material progress is not bundle
completion or a production-readiness claim.

### Continuation: physical Fe electric base and six-component native IR

The same approved bundle continued in the GPFS LagRamses worktree. No new
selector, additional audit gate or Python test framework was added. The
offline optical-data builder is preparation for native Fortran, not runtime
simulation code. No persistent Fe carrier or complete Fe physics is claimed.

`build_fe_grain_optics.py` reads the pinned Werner PDF table 13 and official
Henke Fe file retained above. All **150** Werner optical rows, 0.5--70.5 eV,
are recovered: the final surface-loss column becomes a dash above 49.5 eV,
but its missing values do not invalidate the separately tabulated optical
columns. Requiring all ten fields to be numeric would silently lose 42 rows.

The provisional composite uses DFT epsilon2 below 40 eV, DFT/REELS blending
40--45 eV, REELS through 70.5 eV, a bridge to 75 eV, then actual Henke f2
through 30 keV, retaining the Fe edge knots. Undefined Henke f1 is never used.
Bulk DH13 three-Drude imaginary response is subtracted before the bound
response is integrated. The low-energy linear epsilon2 and >30 keV E^-3
tail are explicitly **dispersion assumptions**, not measured grain data.
An analytic principal-value Kramers--Kronig integral collects adjacent-knot
logarithms, with the continuous x log|x| limit at a pole knot. Extended
precision reduces cancellation. An independent ad hoc QUADPACK Cauchy
integration at eight energies from 0.01 eV to 10 keV differed by at most
4.3e-9 in bound epsilon1; this verifies the integral, not the physical model.

Finite-size Drude damping uses vF/a, then ordinary complex-index Mie theory
supplies Qabs, Qsca and g for 10/100 nm spheres at rho=7.87 g/cm3. There is
no clipping, f-sum normalization or replacement by silicate cross sections.
The resulting electric/eddy bank covers the actual nine primary groups
(including 2--10 keV) and all 136 IR samples. Group edges, representative
energies, IR energies, radii and density have exact native binding checks.

The generated input is `patch/lagRamses/dust_fe_electric_base_data.inc`, SHA256
`9bf5ca47769050ef4c7a76ec07872b27e14db925ade58c40c55bde428cd96a13`.
The offline generator SHA256 is
`a95999118ce4d5c286fdba62454e1eee76c5161d6ea12d41e332f678509d6db1`;
the input contract `.dust-chimes-secondary.H4xQhx/hot-contract.nml` SHA256 is
`8bc8c51d07118b2cd63121e4e379f8c73ba01c836f2cefe02a1d975e6fdf706c`.
The checked-in JSON manifest records these hashes, source URLs, miepython
3.3.0 and the construction assumptions. Earlier provisional generated files
in `.dust-extension.AOz7mU/` are preserved but not included by native code.

**Admission remains false.** The composite effective-electron integral is
26.45372 versus Fe's 26, and its complex dielectric difference relative to
Werner DFT, normalized by max(|epsilon_DFT|,1), is 11.56% median / 61.87%
maximum. For example, at 0.5 eV reconstructed epsilon1 is about -34.99
versus DFT -121.135 despite matching epsilon2=68.667. Agreement with the
assumed epsilon2 and the dispersion identity is not observational validation.
This is not DH13's original optical/VUV composite. Spin magnetic absorption,
phase/temperature-dependent optical response and photoelectron energy
partition are missing. DH13's single-domain Fe-sphere diameter bound
17 nm does not justify silently applying that model to the current 20/200 nm
diameters. No manufactured magnetic table or whole-keV dust heat was added.

Native implementation:

* `dust_iron_optics` combines the existing four D03 bases and two explicit Fe
  electric bases. `fe_full_optics_admitted` is a constant false.
* `snrt_dust_ir` allocates/interpolates the supplied number of optical bases,
  instead of hard-coding four; the old four-basis material ABI still rejects
  six rather than truncating or indexing beyond its layout.
* `iron_radiative_batch` runs the phase-aware material solver on OpenMP.
  Its actual six mass densities and reference mass are explicit. Gas is an
  input reservoir; the surrounding IR transaction applies the returned -Q.
  All output arrays commit only after every cell succeeds. No old CUDA
  material kernel is claimed to implement this phase solver.

Reused `snrt_dust_backend_smoke`, compiled natively with Intel ifx on
lageunha. Logs in `.dust-chimes-secondary.H4xQhx/`:

* `fe-six-optics-build.log`, `fe-six-optics-unit.log`: actual C-only,
  silicate-only, Fe-only and mixed cells through the native six-basis IR
  transport/absorption/scattering path. Relative energy closure 7.29453e-13.
  Heating is from 2.366 eV optical photons below the Fe photoelectric
  threshold, explicitly not a hard-photon heating-partition test.
* `fe-six-batch-build.log`, `fe-six-batch-unit.log`: OpenMP1/2 bitwise parity,
  all three Fe latent intervals, empty cell, and a later-cell root failure
  preserving **all** earlier outputs. Six-basis IR additionally includes
  gas exchange: relative closure 6.58402e-13, scaled by exchanged/dust
  energy rather than hiding error in the larger gas reservoir. Existing
  four-basis backend, gas and adaptive-sublimation cases also pass.
* `fe-six-full-build.log`, `fe-six-mass-unit.log`: full SNRT+DUST_LIVE+
  CHIMES CPU build/link (NENER1/NVAR187) and existing native dust-mass tests
  pass. The new binary `ramses_fe_six_ir3d` SHA256 is
  `5f0f0b34c1331c72924bdb1595440b0bc36c4a3c34f6edaab0d8bc9c1a5cf09b`.
  This checks compatibility of the shared IR change: Fe optical/phase
  modules still are not linked into or activated in the RAMSES driver.

These checks exercise the native IR receiver through a test-host callback
that supplies the optical reference mass; they are **not Fe-enabled RAMSES
evolution**, source injection, AMR transport of Fe mass or Fe restart checks.
Remaining live Fe connections, PAH and relative dynamics are still within
the already-approved bundle, not deferred or counted as completed.

No RAMSES evolution or new dump/restart was launched in this continuation;
all of its build/test commands finished. No commit, push or external audit
was performed. VPATH order and the unrelated generator deletion remain
untouched. `git diff --check` passes.

### Unified continuation: population time evolution and common material receivers

The operator reiterated that Fe, PAH, relative motion and their common live
connections are **one bundle**, not separate approval/completion stages.
No external audit or additional selector was introduced. Work remains in
the GPFS LagRamses repository. No existing production run was touched.

Implemented native `dust_stochastic_evolve` for finite-time backward Euler
of the thermal-continuous master equation. It accepts either probabilities
or **number density in each energy bin**, preserving the incoming total
without normalization; an empty population stays empty. This is the form
needed for conservative transport and source injection, rather than forcing
every advected state to a unit-normalized stationary distribution. Positive
rate elimination keeps the identity/escape term separate in the stiff
lower-Hessenberg system, O(N^2). Returned absorption/emission are integrated
energies, consistent with the population units. The routine does not supply
physical rates, source photons, PAH destruction/charge, overflow treatment
or an emission spectrum. A solved truncated system is not a coverage claim.

`dust_drag_physics` now provides the exact barycentric split
J_b=rho_b*(v_b-v_bary), gas/solid momenta, and internal drift kinetic energy.
It supports an arbitrary component count, not a fixed four-bin layout.
Drag is solved in the barycentric rest frame, preserving total momentum
and converting lost relative kinetic energy to gas heat without subtracting
large bulk energies. A radiation kick receives separately accounted photon
impulses/energy, updates total and relative momenta, and debits kinetic work
before returning grain heat. A scattering force without the necessary
radiation energy debit is rejected; RSLA does not replace h*nu/c by h*nu/chat.
All outputs commit together on success. These are receiving operators, not
a photon-force calculator or a live dust spatial-transport implementation.

The shared elemental receiver and native CHIMES state adapter now accept
optional PAH H/C **mass** inventories as well as separate metallic Fe.
The reserved-solid growth call withholds PAH C from carbon-grain growth and
the metal capacity, and PAH H from the available hydrogen inventory. PAH H
is never subtracted from metallicity. Existing callers with neither PAH nor
Fe remain unchanged. This does not impose an H/C ratio or invent a PAH
stellar-source yield. Live callers still lack the new persistent carriers.

Existing `dust_mass_smoke` was extended, not replaced by a new framework.
`bundled-material-time-unit.log` / `bundled-population-unit.log` in the
existing private build directory pass the old mass/Fe/sublimation checks,
the mixture momentum/kinetic/radiation-energy check with seven carriers,
rollback on unfunded impulses, transient master-equation residuals,
zero and dilute population, and the stiff dt*rate=1e200 case. Two-level
time-refinement errors for 8/16/32 steps are 5.531662e-3, 2.462698e-3,
1.145466e-3 against the analytic solution. No stationary approximation was
silently used to obtain these finite-time results.

Rechecked [HD17 sections III--IV](https://arxiv.org/html/1611.08607): it
provides an Fe-specific UV charging/yield prescription, but its ISM charge
calculation uses a 13.6 eV upper photon energy and does not establish a
complete 10 keV deposition/cascade model. [WDB06](https://arxiv.org/abs/astro-ph/0601296)
addresses carbonaceous/silicate high-energy emission; it is not a supplied
Fe calibration. These are input limits to handle in the existing bundle,
not new review gates or an excuse to pass incomplete physics as complete.

The actual RAMSES field layout is still NVAR187: no Fe mass, PAH population
or signed relative-momentum fields have been activated. Spatial fluxes,
pressure/relative-kinetic accounting, source injection and restart must be
connected coherently; appending ordinary positive passives would be wrong.
No namelist/mkrun option was exposed for an unconnected path. This entry
records implementation progress only, not bundle completion or production
qualification. No additional operator approval is being requested.

Combined native compatibility check finished: `bundled-receivers-build.log`,
`bundled-receivers-mass.log`, `bundled-receivers-ir.log` in
`.dust-chimes-secondary.H4xQhx/`. CPU SNRT+DUST_LIVE+CHIMES (NENER1/NVAR187)
build/link passes, as do the new common receiver checks, PAH H/C reservation
and existing six-basis Fe/gas/IR tests. Those IR closure values remain
7.29453e-13 (without gas) and 6.58402e-13 (with gas). Binary
`ramses_bundle_receivers3d` SHA256:
`a7970549c4d10ae7f4552c00617ba7971b90746bef0b4340123cfd6be35510b8`.
This is compatibility evidence, not a new-physics live run; PAH time and
relative-motion routines still have no RAMSES-driver caller. No simulation,
snapshot or restart was launched. All commands from this continuation have
finished; no commit/push or external audit was performed. VPATH and the
unrelated generator deletion were preserved; `git diff --check` passes.

## Approved bounded comparison: live Fe connection (2026-09-09)

The operator's latest `허용함` authorizes explicitly limited comparison
modes, with missing physics/applicability recorded, rather than claiming
full-band production admission. This supersedes the earlier prohibition on
exposing *any* live Fe comparison, not `fe_full_optics_admitted=.false.`.
It does not close PAH or relative-motion work by relabelling co-advection.

`dust_iron_model='fe_electric_compare_v1'` uses the actual electric/eddy bank,
neutral-grain/common-T material, and two conservative Fe density carriers
after CHIMES (188:189 in NENER1/virial/NVAR189). Fe is already included in
rho, total elemental Fe and Z; `idust` becomes C+olivine+metallic Fe. Fe
atoms are reserved from chemical reconciliation and C/olivine accretion.
Hydro face flux normalization and AMR prolongation include the Fe carriers;
mass restriction/halo synchronization and HDF5 carry the same fields.

The actual six-component primary absorption/transport-scattering and IR
absorption/emission/scattering use these fields. The material callback is
native CPU/OpenMP, not the incompatible four-component CUDA material ABI.
Shared radiation transport still uses its existing runtime backend. The
comparison rejects dust T>300 K and positive primary absorption at group
representative energies >4 eV, collectively before state commit. This is
a representative-energy/grey restriction, not proof of a resolved spectrum
below the neutral-grain threshold. Thermal IR remains its complete original
quadrature; no opacity is clipped or swapped for silicate.

Explicit approximations: room-temperature electric/eddy optics at all
allowed T; no magnetic absorption, charging, photoelectron cascade, Fe
growth/sputtering/shock loss/size exchange/sublimation or relative drift.
Gas collision exchange uses the existing geometric hydrogen-accommodation
comparison, without Coulomb focusing. Fe contributes no adopted C/silicate
H2 formation or grain-recombination surface; CHIMES uses C/silicate mass/area
only for those reactions. These omissions are NOT production-qualified.

`dust_fe_condensation` defaults to zero, with an explicit [0,1] comparison
fraction of non-Ia ejecta Fe remaining *after* source-node olivine formation.
Injected Fe is all-large; thermal energy at `dust_injection_temperature` is
partitioned from the same source energy, with the existing lock/MPI/progress
transaction. The separate SNIa budget remains gas. This is not a measured
Fe condensation yield or a new population/DTD approval.

Namelist requires CHIMES, D03 two-size DL01 and no sublimation. Build profile
`DUST_IRON=1 CHIMES=1 DUST_LIVE=1 NENER=1` reserves NVAR189; existing default
profiles remain unchanged. `mkrun.py` and the shared GUI/namelist generator
expose and validate the same choices. The bounded wizard intentionally
omits the incompatible BPASS hard SED and requires explicitly selected Fe
binary/IR contract paths. HDF5 restart binds exact optical and material
arrays, field offset, limits and condensation fraction; missing, altered,
or disabled comparison identity is rejected.

Native work and evidence are in `.dust-fe-live.J71qGw/`. A clean parallel
build encountered an existing missing Morton/AMR module dependency; the
serial build succeeded without changing VPATH or unrelated AMR code.
Full NVAR189 binary plus existing mass/IR smoke builds pass. The existing
six-bank IR checks retain energy residuals 7.29453e-13 and 6.58402e-13.
Setup/GUI checks: 39 run, 38 pass, one display-dependent skip. Native IC
preparation checks mixed U(T), six weights and the 301 K rejection.
Actual MPI evolution/restart results are recorded below when completed.

### Live Fe evidence and resolved cold-cell failure

The first startup attempt (`live/`) accidentally selected the old-format
fallback yield input and stopped before producing any dump. The environment
was corrected to the existing `.dust-mass.w1Jo9A/live/yields.dat`; neither
the physical yield data nor the fallback file was changed.

`live2/` then passed step 1 but rejected step 2. This was NOT a 300 K domain
violation. Increasing the bisection budget from 100 to 256 alone did not
fix it (`restart/`). Stage-specific diagnostics (`restart-trace/`) isolated
the final spectral energy closure, native status16: a converged material
root was stopping too close to the final 2e-12 energy tolerance, so spectral
summation rounding could move the final residual across the threshold.

The material root now stops at one quarter of the final tolerance, reserving
rounding headroom; the final 2e-12 closure criterion is UNCHANGED. The actual
cell and adjacent representable conductances reproduce the original
failure in a native standalone baseline, without a RAMSES rerun or Python
solver. The fixed cell residual is 1.54977e-33 erg/cm3 for Q=1.636286e-20,
about 9.47e-14 relative to exchanged energy. All 33 neighboring inputs pass.
This case is retained in the existing `snrt_dust_backend_smoke`, not a new
test framework. Final native `regression.log` passes old closures, Fe phase
and six-component radiation/transport/scattering tests; their energy
residuals are now -1.49206e-13 (without gas) and 5.52616e-14 (with gas).

Final binary `.dust-fe-live.J71qGw/ramses_fe_fixed3d` SHA256:
`7b5be377810f39b6c5da2bef560420511131bfecd0e856e0548d7e0a862c94ac`.
`fixed/physical.nml`: noncosmological periodic 8^3 cells, NENER1/NVAR189,
CHIMES, C/silicate plus six-bank Fe optics, CPU MPI2/OMP2, no sinks/AGN or
stellar SED. Initial metallic Fe fraction is .001 (half small/large),
silicate fraction .005, common T=20 K; gas moves at .01 code velocity.
C/silicate growth/erosion/size exchange are disabled for this bounded check;
stellar mass/energy/element sources, all-large Fe condensation fraction .2,
chemical reconciliation, gas exchange and IR transport are active.

Two coarse steps completed. IR energy residuals: 5.0114e-13 and 1.7068e-14.
Dust aggregate/species relative mismatch <=3.5295e-16. Step 2 has positive
injected large-Fe excess over small Fe (3.894316945e-10 summed code density),
while solid Fe+olivine Fe stays below the total elemental Fe reservoir.
This exercises actual source injection, hydro/chemistry and the IR caller,
not just a manually invoked material kernel. It does NOT exercise a hard
primary source or establish spatial-resolution convergence.

`restarted-fixed/physical.nml` reads the completed first dump through a
read-only-use symlink and evolves to step 2 with the same binary/environment.
All **587** hydro, SNRT and particle datasets are bitwise identical to the
continuous run, as are the exact Fe and CHIMES restart identities.
The two fresh-run dumps and the resumed final dump each occupy about56 MiB.
Every launch was reported with its namelist, physics/MPI/OMP and output
schedule; all trial logs and existing outputs are preserved.

This completes bounded Fe live integration evidence, not the entire bundle.
PAH populations and relative-motion spatial transport remain unconnected;
full magnetic/charging/high-energy Fe physics remains unqualified.

Final negative restart check: `reject-identity/physical.nml` changed only
the Fe condensation comparison fraction (.2 -> .3), plus suppressed output
clocks for the reader-only rejection test. The same NVAR189 binary rejected
the checkpoint collectively with `dust mass-evolution restart identity
mismatch` (MPI abort11), before evolution and with no new dump. The linked
original completed output remains untouched. No commit/push, external
audit, or global shared-context modification was performed in this turn.

## PAH spectral radiation connection (2026-09-10, same approved bundle)

Operator requested implementation of item 1 in the latest two-item remaining
summary: PAH stochastic heating. This is NOT a new bundle or an approval
gate. The separate Q-GOAL / Q-LEAN directive is recorded in the existing
audit cadence/governance documents and project AGENTS.md; it adds no review.

Implemented native `patch/lagRamses/dust_pah_radiation.f90` and connected it to
the existing `snrt_dust_ir_advance`, rather than an independent transport code:

* Actual primary angular photons are attenuated using the same per-grain
  optical coefficients as emission. Primary debit, IR radiation and energy
  populations commit together only on success.
* The IR nonlinear iteration now optionally carries the absorbed **spectrum**
  to a stochastic material callback. Every iterate starts from the original
  number distribution; both spectral and total-energy residuals must converge.
  Existing common-temperature callbacks and their default behavior remain.
* Finite-time populations use the existing positive backward-Euler solver.
  DL01 eq32 supplies low-excitation emission temperatures, with the canonical
  harmonic inverse above the twentieth mode. Spontaneous thermal-continuous
  cooling integrates only photon energies below the excitation. No emitted
  spectrum is renormalized to an imposed luminosity. Number-weighted T is a
  diagnostic, not a thermodynamic/common grain temperature.
* H/C solid masses and excitation are derived from the same population.
  A significant above-grid absorption tail rejects the transaction; no
  high-energy tail is deleted, clipped or silently converted to heat.
* The original neutral PAH table is read natively, including all 30 radius
  blocks and 1,201 wavelength nodes. Fixed-width negative `g` entries can
  adjoin the preceding column: whitespace parsing was the initial reader
  failure, NOT unequal wavelength grids. The corrected reader verifies the
  grids match. Radius interpolation uses the source's carbon-count mapping.
  No wavelength or radius extrapolation is made.

Evidence is in `.dust-pah-native.LWUJVc/` (GNU bounds/FP checks) and
`.dust-pah-ifx.1qgTOo/` (actual Makefile/mpiifx). The existing dust mass smoke
was extended, not a new Python framework. Both builds pass two-cell spectral
transport/reabsorption, equal-energy/different-spectrum response, dark
cooling, H/C conservation, unsupported-photon rejection and late-overflow
rollback including the staged primary debit. Original neutral C24H12 optics
also pass the full primary -> population -> IR transaction:

* GNU IR energy/spectral residuals: `6.135629e-14`, `6.121559e-14`.
* Intel IR energy/spectral residuals: `6.193125e-14`, `6.121559e-14`.
* Total primary + IR + excitation balance rounds to zero in both cases.
* Existing C/S/Fe backend parity/rollback smoke passes with two OpenMP
  threads; the actual `snrt_dust_live.o` caller recompiles against the extended
  optional interface. No VPATH order changes or RAMSES evolution launches.

The original compressed table remains at `.dust-extension.AOz7mU/PAHneu_30.gz`
with SHA256 `d1d0f8371b730ab473ab90cb98def2c6367507ae7a34c5a1214eab33a24fe09c`.
The decompressed input used by `SNRT_PAH_NEUTRAL_TABLE` is
`.dust-pah-native.LWUJVc/PAHneu_30.dat`, SHA256
`ea5ecc3ef7f51a00925eff099c2efaac58cab3d58fdd062550ff088890fbee2a`.
Sources: [DL01](https://arxiv.org/html/astro-ph/0011318),
[Draine optical tables](https://www.astro.princeton.edu/~draine/dust/dust.diel.html).

**Item 1 is still incomplete.** This is a PAH-only native spectral receiver,
not yet selected by RAMSES: live population fields, source H/C bookkeeping,
mixed C/S/Fe/PAH opacity competition, hydro/AMR co-advection, restart identity
and matching namelist/mkrun remain to connect. No disconnected selector was
added. Neutral charge/hydrogenation/survival and spectral coverage must stay
explicit: the original PAH optical interval is about .00124--1240 eV, whereas
the existing D03 live grid extends beyond both ends. Neither truncating that
grid nor treating equilibrium excess IR as absolute PAH illumination is
silently admitted. No claim of charging/destruction, general production
qualification, external audit approval or completion of the whole dust bundle.

## Item 1 live-connection closure — 2026-09-10

The immediately preceding incomplete-status paragraph is historical and is
superseded for the live connections by this section. Work stayed in
`/gpfs/kjhan/LRD_JWST`, verified origin `kjhan0606/LagRamses`. No new bundle,
external audit, generic Python test framework, commit or push was added.
Existing runs and the unrelated generator deletion were preserved.

### Connected implementation and physical conventions

* `dust_pah_model='pah_neutral_absolute_v1'` selects neutral C24H12; defaults
  remain off. `DUST_PAH=1 CHIMES=1 DUST_LIVE=1 SNRT=1 HDF5=1 NENER=1`
  builds NVAR315, or NVAR317 with `DUST_IRON=1`. NVECTOR32; VPATH unchanged.
* `dust_pah_live_model` defines the 128 excitation levels and injection
  distribution, original-table optical projection, molecular mass and
  restart identity. RAMSES fields hold mass density per state. CHIMES's
  integer nuclear-mass convention is used consistently: molecule=300*m_p,
  H/C mass fractions=12/300 and 288/300, not two independently rounded
  atomic-mass definitions. Bulk `idust`/`idust_energy` stay C/S/Fe only.
* `dust_pah_mixed` couples the stochastic receiver to the existing native
  IR transaction. Absolute IR starts empty, with no hidden bath subtraction.
  Bulk and PAH physical opacities partition each absorbed spectral node
  once. The primary driver adds physical PAH optical depth to the same dust
  sink, including bulk-free PAH cells. Spectral and total-energy iteration
  starts every trial from the original populations, not the prior iterate.
  The first iterate includes known absorption of transported old radiation.
* Source deposition calls the same native `dust_pah_condense` tested by the
  mass smoke: an explicit [0,1] fraction of non-Ia C remaining after graphite,
  limited by H from the same ejecta. The injection temperature supplies mean
  excitation via two adjacent levels; its energy is deducted in the same
  source/CR/gas transaction. Separate SNIa ejecta remain gas. This fraction
  is an uncalibrated comparison parameter, NOT a physical stellar PAH table.
* Consistent hydro carriers include PAH mass and reconstruct total H/C;
  chemistry withholds both atoms from gas, and grain growth reserves them.
  Existing passive prolongation/restriction and ghost synchronization carry
  all 128 fields. This is gas co-advection, not relative dust dynamics.
* HDF5 carries the fields and an 822-value PAH identity (actual projected
  optics, axes, molecular/injection conventions, index and source fraction).
  Presence/absence and changed data are checked. No restart can silently
  reinterpret PAH absolute IR as the old excess-above-bath convention.
  The hydro descriptor now prints wide scalar indices without `**`, and
  labels PAH states explicitly. The descriptor-only formatting change does
  not change stored HDF5 arrays or their numerical meaning.
* Both namelist generation frontends and mkrun CLI/GUI expose and validate
  the selection. The wizard uses one combined PAH/optional-Fe binary and
  `SNRT_DUST_PAH_CONTRACT`/`SNRT_PAH_NEUTRAL_TABLE`, omits incompatible hard
  BPASS, and records restrictions in its README. Forced CUDA material is
  rejected; existing CPU transport/scattering support is retained.

The original optical input/hash are unchanged from the previous section.
The live projection uses the existing 136 IR nodes. At wavelengths beyond
1000 microns it explicitly anchors an E^2 long-wave Rayleigh/Drude tail;
unsupported high-energy bins are rejection masks, not assumed transparency.
Occupied primary photons above a 4 eV representative energy are rejected
before commit even in PAH-only cells. The finite excitation grid rejects
significant overflow rather than losing or thermalizing energy.
The fixed neutral/H-C choice and long-wave behavior follow the model form
in [Li & Draine 2001, equations 4 and 12](https://arxiv.org/html/astro-ph/0011319);
they do not supply charging, photoelectron escape or survival physics.

### Integrated native evidence

Private build and all effective inputs/logs: `.dust-pah-live.kBEoZy/`.
Actual mpiifx/Makefile build succeeded. A fresh parallel build exposed the
existing missing AMR module dependencies; sequential build succeeded without
changing VPATH, AMR algorithms or unrelated production processes.

Final artifact: `.dust-pah-live.kBEoZy/ramses_final3d`, SHA256
`a65d6a9f8e3f6aaa8f5f6186047e7b27d1247a61a1a2f71c92aba702f1bb5b74`.
`build-complete.log` records the final successful link. A separate no-DUST
SNRT syntax compile initially exposed an unguarded PAH gas-reservation block;
it is now guarded and `nodust-syntax-fixed.log` is empty with exit0. This is
a compile check, not a full alternate-profile evolution qualification.
The post-evolution changes were descriptor text formatting and conditional
compilation guards; the exercised DUST_LIVE numerical path is unchanged.

| Evidence | Result |
|---|---|
| Existing mass smoke | PAH post-graphite H/C source, H-limited injection, invalid carbon budget, off-default and reserved-growth checks pass; original-table primary/IR/population rollback checks still pass |
| Existing backend smoke, with physical contract argument | C/S/Fe regressions pass; mixed bulk+PAH and PAH-only primary partition / IR / excitation balance residual `5.01759764e-13`; unsupported hard source rolls back radiation, emission, temperatures, energies and populations |
| Frontend CLI/GUI tests | 40 run: 39 pass, 1 display-dependent skip; PAH+Fe needs one binary; invalid mode/fraction/cosmology/sublimation rejected |
| `live/`, MPI2/OMP2, NVAR317 | 2 actual updates, 512 cells, nonzero PAH excitation. IR balance residuals `2.4549e-16`, `6.5192e-17` |
| `restart/` from live step 1 | Step-2 HDF5: **1,005 datasets byte-identical** to continuous execution; PAH model identity equal |
| `restart-mismatch/` | Condensation changed .05→.06: collective identity rejection (MPI abort11), no new dump |
| `source-advection/`, MPI2/OMP2 | Two PAH fractions across a moving periodic gas. Stars 512→1024, nonzero physical stellar return and PAH injection; IR residuals `2.9680e-13`, `5.8296e-14` |
| `source-restart/` | Physical source + moving PAH restart: **1,005 datasets byte-identical**, including all 128 PAH fields and particle source progress; no double injection, no NaN/negative PAH |

For the source/advection case the gas+star mass is **10.0 code units in
both snapshots**. Step 2 has gas 9.990278150361132 and stars
0.009721849638868569; actual cumulative returned stellar mass is
2.1719298672602546e-5. Initial PAH fractions were 1e-5 and 2e-5; after step2
their range is 1.000163406458414e-5 to 2.0001613716217637e-5, with 512
distinct advected/source-processed cell fractions. Excited-state mass sum
increases from 3.23349824739962e-9 to 6.5164706676577515e-9 (sum of cell
code densities, not volume-integrated solar masses). Gas-only `mcons` in
the log excludes converted stellar mass; it is not a baryon-loss signal.

The first `preflight/` used the old very dilute/cold-gas Fe input and an
initially empty absolute IR field; it rejected the IR update, no dump.
It is **not** included among passing simulations. Bulk tables still have
a lower temperature bound of 5 K, and with Fe an upper bound of 300 K;
a long-step isolated radiative solution outside these bounds is not admitted.
The passing warm-gas and high-density/source cases are separately named,
preserved tests, not edits or relabellings of the rejected case. No tolerance
was loosened and no artificial heat bath or floor energy was added to pass.

Every simulation launch reported its effective absolute namelist, complete
output clocks, build/physics/MPI/OMP and storage. Passing continuous runs
each wrote two ~57 MiB dumps; restarts each one. Reader-negative runs wrote
none. Existing step-1 outputs were linked for read-only restart use. Free
GPFS space was 114--116 TiB during these checks. No large production job was
launched; all finished test processes have been collected.

### Closure scope

All previously listed **item-1 live wiring** is now connected and exercised;
the result is no longer a detached stochastic kernel or unconnected selector.
The source-bearing integration/restart proof is fixed-level 8^3, not a new
AMR refinement-convergence or calibration claim. No extra acceptance gates
were created. The neutral comparison still excludes charge/hydrogenation
evolution, PAH destruction/photoelectric heating, arbitrary hard spectra,
cosmological radiation backgrounds and dust relative motion. Its 136-node
spectral resolution and existing bulk temperature bounds are explicit model
limits, not general production admission. The other approved bundle items
are not declared finished by this result.
