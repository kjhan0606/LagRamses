# Cosmic-ray pressure / SF implementation and dust evolution inspection

Worktree: `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`,
base commit `3d51b88`. This is implementation evidence, not a production
or publication approval. No external audit or parameter-calibration job
was launched. The operator preapproved continued implementation, but
requires a COLIBRE-style calibration proposal AFTER completion and a new
approval before any calibration execution.

## Approved closeout result (2026-09-08)

Following explicit preapproval of the entire closeout, the current restricted
CR/bulk-dust implementation has completed source review and focused regression
checks and is being handed off with the code changes. This closes this
implementation bundle, not the full unresolved physics programme. The advanced
dust scope below stays medium-term; source-data/population limitations in the
parent plan remain open. No new simulation, calibration or external audit was
needed for this closeout.

- Native CR partition/support, dust condensation/growth/sputtering/energy,
  and OpenMP photon closure/atomic rejection smoke executables pass again.
  Their build targets were already up to date. CUDA was unavailable in this
  CPU-linked smoke, so its result is explicitly SKIP; the earlier actual
  GPU result below remains the CUDA evidence, not a fresh GPU qualification.
- The preserved final combined and NENER=0 compatibility binaries still
  match their recorded SHA256 values. Final live and restart logs both
  report step 4 and `Run completed`; existing HDF5 comparison evidence is
  reused rather than regenerating the same snapshots.
- Shared setup/GUI tests: worktree 34 run, 33 pass, display test skipped.
  Against the exact staged generator with unrelated deletions excluded:
  34 run, 32 pass, display and missing-sector-only tests skipped. The latter
  skip is expected because the staged generator retains those sectors.
- Corrected the new Makefile dependency guard from `ifdef SNRT` to
  `ifeq ($(SNRT),1)`: explicit `SNRT=0` no longer pulls in the new dust runtime.
  Make dependency inspection verifies both 0 and 1 configurations. VPATH
  remains lagRamses first. This does not change the tested physics binary.
- Runtime documentation now includes the CR setup path, restrictions and
  local-input/build dependencies. A fresh clone alone is not a packaged
  physical-data release. New switches remain off by default.
- Commit scope excludes the pre-existing 88-line deletion of unrelated
  generator options and all scratch inputs, executables, logs and snapshots.
  Nothing was deleted or restored. Diff whitespace checks pass.

## Current scope decision: close out first, defer advanced dust (2026-09-08)

The operator approved finishing the current implementation before starting
the newly researched advanced dust programme. This decision supersedes any
reading of earlier continuation instructions as a requirement to implement
all newly identified dust physics in the current bundle.

Current closeout is limited to essential RT/feedback/CR/bulk-dust wiring,
defects that invalidate the admitted model (including conservation, restart
and runtime defects), usable run settings, and an honest scope/limitations
handoff. Reuse the native evidence below; do not introduce another suite of
audit gates or repeat successful runs without a relevant code change.
Pending code review and repository handoff are not completed merely by this
scope decision. Current bulk dust remains a restricted comparison model,
not full physical production/publication approval. Existing stellar source,
population-matching and metallicity-data limitations remain open in the
parent plan; this decision neither resolves nor silently defers them all.

Advanced dust is the highest-priority medium-term extension, with the
following retained scope. These are planning groups, not new per-item gates:

1. **Element/composition physics:** choose explicit carbon/silicate/iron
   carriers; enforce shared gas/dust elemental budgets and gas-phase cooling;
   connect stellar C/O and limiting-element condensation, distinguish fresh
   ejecta processing from ambient SN destruction, and retain the factor-three
   radius-to-mass sputtering conversion. Do not invent an oxygen carrier.
2. **Size solver selection:** compare the bulk control and a two-size model
   with a bounded 8/16-bin piecewise-linear reference; consider moments only
   where their closure can meet the needed observables. Include coagulation,
   shattering and positive conservative updates. Assess unresolved growth
   against the SF assumptions and actual CPU/GPU/MPI cost. Reuse free-stream
   GPU dispatch with OpenMP fallback; do not add a separate autotuner.
3. **Coupling:** connect size/composition to opacity, scattering, material
   energy and IR; connect depletion, grain-surface H2, shielding and CO
   competition without treating CR pressure as CR ionization chemistry.
   Assess relative dust motion/drag and AGN radiation pressure separately;
   independent dust particles are not a prerequisite for the other work.

Select the least costly adequate model by comparing numerics at fixed
physics, then physics at fixed numerics. Implementing every alternative is
not the objective. PAH/porosity/extinction-bump detail remains research scope.
The literature informs candidates, not an already validated winner; see
[RAMSES two-size dust](https://arxiv.org/html/2402.18515v2),
[NewCluster](https://arxiv.org/html/2508.18374v1),
[COLIBRE dust](https://arxiv.org/html/2505.13056v2),
[FIRE size-resolved dust](https://arxiv.org/html/2603.08504v1), and
[the moments approach](https://arxiv.org/abs/1606.02272).

Order: current approved-scope closeout -> advanced dust/model selection ->
final COLIBRE-style tuning. A calibration proposal may identify dependencies
earlier, but final tuning must not precede the model choice. Calibration
execution still requires separate operator approval. No advanced-dust or
calibration job is launched by this planning update.

## Implemented CR reference

`patch/lagRamses/cosmic_ray_physics.f90` selects a relativistic trapped
fluid, gamma=4/3. Existing RAMSES NENER pressure, enthalpy flux, CFL,
advection and PdV work are used; CR energy is NOT a passive metal scalar.
The cold-gas thermal energy is `Etotal - Ekinetic - Ecr`. SNRT chemistry
and joint gas/dust/IR exchange now subtract NENER energy in both thermal
receivers, leaving CR energy unchanged by thermal radiative cooling.

The channel-resolved stellar transaction assigns
`Ecr_source = cr_sn_fraction * E_SNII + cr_snia_fraction * E_Ia_coupled`.
The latter is the existing Ia thermal-coupling budget. No bulk ejecta
kinetic energy is charged to CR. Total deposited energy is unchanged:
CR is a component of it, not an additional source. The CR increment is
published by the same cell lock, reverse MPI exchange, stellar mass/clock
commit as the other source fields. AGB/wind energies are not diverted.

Optional `cr_sf_support` adds `gamma_cr * Pcr / rho` to the effective
sound-speed squared in virial SF models 1,2,4. The gas temperature cut
still sees only thermal energy. This is an explicitly chosen effective
compressibility closure, NOT a universal CR suppression law or a
calibration of the turbulent-density PDF. The tested/default comparison
uses model 4. CR remains in the gas reservoir when gas forms stars.

Namelist fields in PHYSICS_PARAMS: `cr_enabled`, `cr_transport`,
`cr_sn_fraction`, `cr_snia_fraction`, `cr_sf_support`; all new physics is
off by default. `mkrun.py` comparison modes expose the same choices and
select the NENER=1 CPU executable, turn off sinks/AGN, insert CR initial
pressure, and shift the virial/chemical/dust IC columns. The shared GUI
database also exposes them. `sf_virial` is correctly typed logical.

Build from a new one-level-below-root directory (do not reuse NENER=0
objects):

```
make -f ../bin/Makefile -j1 HDF5=1 SNRT=1 DUST_LIVE=1 NENER=1 USE_CUDA=0 USE_FFTW=0 EXEC=ramses_cr
make -f ../bin/Makefile HDF5=1 SNRT=1 DUST_LIVE=1 NENER=1 USE_CUDA=0 USE_FFTW=0 cosmic_ray_smoke
```

Makefile VPATH retains lagRamses before cuRamses. No CUDA hydro admission:
its CR fluxes have not been qualified. NENER=1 is required even when the
CR initial pressure is zero. Only hllc/hll/llf and noncosmological periodic
domains are admitted. HDF5 binds six CR model/fraction/SF identity values
and refuses missing/changed identity on CR restart or disabling CR on a
CR checkpoint. Old non-CR checkpoints with CR disabled are unaffected.

Not implemented by this reference: diffusion, streaming, collisional
losses/heating, shock acceleration, CR ionization/spectral bins, sink/AGN
CR partition, or cosmological super-comoving CR expansion source. These
must not be inferred from the existing pressure/advection implementation.
The general cosmological wizard does not silently turn this mode on.

Physical context: the distinction between advection and anisotropic
diffusion follows [Dubois & Commercon](https://arxiv.org/abs/1509.07037).
CR support depends on transport/trapping, as studied by
[Commercon et al.](https://arxiv.org/abs/1811.11509); that work does not
validate the particular SF sound-speed substitution used here.

## Direct native evidence

Scratch root `.cosmic-ray.kyySgK`. Physics-tested executable SHA256
`44ecf8a616253f5a7b90627ba7598c352c02ff83b31dbb7fdbfa20ba3529a5bd`.
Final executable after additional admission-only guards SHA256
`583fddf888cc35607ea7d1fc85dc435fef48b1704c9c81087691eaf3c4682d1e`.
The intermediate guard build `588c244c...` incorrectly treated an empty
BOUNDARY_PARAMS block as a physical boundary: the native parser clears
that derived flag later. It rejected the generated nml before any output.
The corrected guard checks the actual boundary count. Its failed log is
preserved as `generated/run.log`; the corrected run is `run-final.log`.
All builds CPU-only, NDIM=3/NVECTOR=500/NVAR=30/NENER=1/HDF5/SNRT/DUST_LIVE.

The default NENER=0 profile also rebuilt successfully as
`.snrt-cpu.OKoz9T/ramses_cr_compat3d`, SHA256
`85f5d4a7b070852dd7c9a04e37de5a9c5e8c1baf463af5c0bba252c173be45fb`.
This is compile compatibility evidence, not another NENER=0 runtime run.
Shared wizard/GUI suite:33 tests,32 passed and one display-dependent skip.
No external reviewer was invoked for these small connection fixes.

- `cosmic_ray_smoke`: partition, total-energy budget, pressure derivative,
  effective support, and invalid fractions/build/gamma checks passed.
- `sf-on/physical.nml`, `sf-off/physical.nml`: 1 MPI / 2 OpenMP / 512M
  worker stacks, noncosmo level3 fixed 8^3 mesh, four steps, ordinary LC18
  + AGB pulse/effective Ia source, independent BPASS radiation, no AGN,
  DL01 scattering and joint nonlinear gas/dust/IR exchange. About39s each.
  Initial PCR=1e-8, Ecr=3e-8. Final Ecr is positive and increased by SN
  sources. All hydro arrays finite; thermal energy positive.
  Joint gas/dust/IR balance residual at most9.724e-10 in sf-on.
- The on/off final stellar masses are 2.9379572446e-8 and5.1727765905e-8
  in code units, BUT final times differ (.10975597 vs .10175563), because
  the feedback histories change subsequent CFL steps. These numbers show
  a live connection; they are NOT a matched-time SF suppression measurement.
- `gradient/physical.nml`, `zero/physical.nml`: 2 MPI / 2 OpenMP,
  four-step periodic hydro, no RT/cooling/SF/gravity. Gas initially at
  uniform pressure; only PCR differs between two half boxes. Gradient
  run develops momentum extrema +/-6.94472554e-9, with sum exactly0,
  kinetic energy mean9.41121831e-15. Mean total energy remains
  2.999999999895e-8 and density.001. Zero-CR control remains exactly at
  zero momentum/kinetic energy. Thermal/CR energies stay nonnegative.
- `restart/physical.nml`: COPY sf-on output_00001, resume to step4.
  All hydro and SNRT datasets exactly match continuous sf-on. All2047
  distinct-ID stellar masses, initial masses, birth/progress clocks and
  metallicities also exactly match. About13.2s.
- `reject-restart/physical.nml`: same COPY, change SN fraction .1->.2;
  MPI abort10 during HDF5 restore, before any new output_00002.
- `generated/physical.nml`: actual shared mkrun preview, real pulse source
  files copied into its new directory; Z=.01, 1 MPI/2 OMP, same four-step
  schedule and scattering/joint exchange, no sink IC file. Final binary
  completes in31.784s (`run-final.log`). All129 floating datasets finite,
  gas thermal energy positive, CR density3.05945e-8--3.08377e-8; maximum
  joint gas/dust/IR balance residual9.414e-10. HDF5 CR identity matches
  the generated .1/.1 injection fractions and SF support flag.

Effective paths are the above directories beneath the absolute scratch
root. Each positive run used noutput1, aout2/tout1e30 unreached,
foutput2/fbackup1000000; expected <=56MB/dump,161TB GPFS free at launch.
No old output was overwritten/moved. The complete scratch root occupies
392MB before the final generated run (two further dumps of order56MB).
No MG nonconvergence was accepted.

## Dust inspection before the subsequent mass-evolution implementation

Reading the native source map, deposition scratch and the final SNRT
commit confirms that `idust` is an advected passive mass reservoir. SF
removes it with consumed gas. The RT transaction reads it to scale opacity
and material energy, but writes `idust_energy`, not a mass-growth source.
The stellar source map has no dust-condensation deposition field.

Connected: dust mass/thermal-energy advection, astration, absorption,
isotropic elastic scattering, bulk material emission, IR transport and
absorption, gas/dust collisional thermal exchange, restart identity.

Not connected: stellar dust condensation/injection, accretion of gas-phase
metals onto grains, thermal sputtering or SN shock destruction, shattering,
coagulation, element-resolved depletion, evolving size/composition,
grain/gas drag or dust radiation pressure. Fixed material and collision
size parameters are NOT a grain-size evolution model.

The pure-hydro pressure test conserves mean dust mass at
2.409638554216868e-6 (code density), while the star-forming run decreases
it to2.409559416778099e-6 and changes thermal energy. This is consistent
with advection/astration and thermal response; it does not test growth.

For comparison, [McKinnon et al.2018](https://arxiv.org/abs/1805.04521)
explicitly evolve grain sizes under accretion, sputtering, shattering and
coagulation. Those processes must not be attributed to this code merely
because a dust thermal solver and optical tables are present.

## Completion boundary

The CR reference connection is implemented and directly exercised; the
entire preapproved physics programme is NOT complete. The parent plan's
remaining physical-data dependencies (same-population SED, low-mass CCSN
and edge AGB coverage, microscopic binary/SNIa etc.) remain. This check
also identifies substantive dust mass-evolution work, not another audit
gate. CR cosmic/AGN/transport limitations above remain explicit.

Do not replace these dependencies by scalar normalizations, fake source
clocks, zero yields or calibration coefficients. No COLIBRE calibration
run is authorized yet; completion and a bounded, resolution-dependent
proposal must precede the user's final calibration approval.

## Subsequent approved implementation: native bulk dust mass evolution

The user's further approval was applied to continued implementation, not
to a COLIBRE calibration launch. Added `dust_mass_physics.f90` (condensation,
cold geometric accretion, thermal sputtering, exact frozen-rate bounded
integration and conservative material/gas energy exchange) and
`dust_mass_runtime.f90` (native leaf update, all-rank validation before
commit, restriction and halo exchange). The stellar deposition transaction
now includes dust condensed from returned wind/AGB/SNII metals, with no
SNIa dust. Dust is already part of returned rho/total metals; no double
mass addition. Its injection energy is partitioned out of the source
energy remaining after CR allocation. Evolution transfers heat only
between material and gas thermal reservoirs, not CR/kinetic energy.

This supersedes the earlier statement that condensation/growth/sputtering
were disconnected. It does NOT supersede the explicit limitations for
element depletion, evolving composition/size, shock destruction,
shattering/coagulation, drift or dust radiation pressure. The model is
one-moment with fixed characteristic radius and composition; it is not a
full grain distribution. External metal cooling is rejected because
depleted gas-phase metal cooling is not connected. No latent heat model.
Noncosmological periodic scope and no-sink/AGN restrictions are explicit.

All new namelist parameters are also in the shared GUI generator and
mkrun comparison setup. Evolution defaults off. HDF5 `dust_mass_values`
stores version/enabled/growth/sputtering, three condensation fractions,
radius/density/sticking/temperature cap/effective atom mass/injection T.
Restart rejects changed or missing identity rather than changing physics.
Default coefficients and physical equations are documented in NATIVE_RUNTIME.
The grain erosion fit follows
[McKinnon et al. 2018, section 3.4](https://arxiv.org/html/1805.04521v1);
the cold growth coefficient is the geometric collision rate with a
gas-phase availability factor. Comparison coefficients are not calibration.

### Focused native execution evidence

Scratch root `/gpfs/kjhan/LRD_JWST/.dust-mass.w1Jo9A`.
Build: HDF5/SNRT/DUST_LIVE/NENER1/NVAR30, CPU, lagRamses VPATH first.
Existing CR executable preserved; mass executable is
`.cosmic-ray.kyySgK/ramses_dust_mass3d`.
The first successful build after the photon cap correction had SHA256
`259dc46388a5dfb0a45f7c96406ee40db1baa80f1a4eefcf839f8a36c3a0f473`.

Effective run files were explicitly announced before execution. All are
periodic level3 8^3 references, one MPI/two OpenMP unless noted later,
512M thread stacks, noutput1/aout2/tout1e30 unreached/fbackup1000000.
Four-step live/control use foutput2; two-step cold/hot use foutput1.
Two ~56MB dumps per positive run; restart copies step2 and adds one dump.
159TB free before launch. No old output was overwritten or removed.

- `live/physical.nml`: real LC18/AGB7-pulse/effective Ia sources plus
  independent BPASS radiation, CR pressure/SF and scattering/thermal
  exchange. Initial dust mass AND energy zero. Completes step4 in5.175s.
  Mean dust density4.282749639924887e-11 code units, range
  [3.513600744850457e-11,5.0309896329057433e-11]; max D/Z5.0321e-6.
  Mean material energy5.224765313601319e-22; minimum gas thermal
  energy1.0639570661233599e-7. All hydro finite. Total2047 stars.
- `control/physical.nml`: same zero-dust initial condition, evolution off,
  completes4 steps in5.546s. Dust mass/energy remain exactly zero.
- `cold/physical.nml`: no stars/gravity, rho1 code, initial P1e-6,
  initial dust density.0060240963855, Zdensity.01. Two steps finish47.564s.
  Dust density.006129571369711066, +1.7508847%; total metal unchanged.
- `hot/physical.nml`: no stars/gravity, rho.001, initial P1e-4,
  initial dust density6.0240963855e-6, Zdensity1e-5. Two steps11.248s.
  Dust density5.8434173524677935e-6, -2.9992719%; total metal unchanged.
  Both thermal controls keep finite hydro, positive thermal/material
  energy and 0<=D<=Z. Their `thermal.env.sh` explicitly unsets stellar SED
  because star formation is disabled; the first incompatible setup was
  correctly rejected and its original log is retained.
- `restart/physical.nml`: COPY of live step2, continues through step4 in
  2.723s. All90 hydro datasets and4 SNRT datasets exactly equal continuous
  execution. All15 per-particle fields exactly equal for2047 matched IDs.
  Dust identity equal. `reject-restart` changes sticking .30->.31 and
  aborts with identity mismatch before any output2 (MPI_Abort11).

Fortran `dust_mass_smoke` passes condensation, analytic growth/decay,
balanced rates, split-step identity, stiff bound, zero seed, invalid mass,
both thermal-transfer signs, insufficient energy, cold/hot rates.
GUI/model tests:34 run, one display-dependent skip, all others pass.

### RT defect exposed by the zero-dust source run

The first live run rolled back at step1, transport status103. A standalone
cell case in the existing C++ backend smoke reproduced a negative photon
value -2.98023e-8 at trial3: FP32 H/He weighted shares summed above the raw
absorbed photon budget, yielding a directional cap >1 in a fully absorbed
cell. This was not fixed by inserting an artificial dust seed or clipping
negative photons afterwards. The shared CPU/CUDA cell now reduces the
overbudget species shares before BOTH atom debit and directional return;
already admissible shares are unchanged. The same4096 zero/trace-dust
cases pass, along with the existing closure/atomic reject CPU checks.
The direct-cell test uses precise FP arithmetic, matching the production
CPU kernel; fast reassociation is not a meaningful test of rounding caps.
CUDA equivalence is recorded separately below when actually exercised.
First failed `run.log` and diagnostic `run-debug.log` remain preserved;
successful live results are in `run-fixed.log`. Failure-code diagnostics
were retained in the native driver so transport/scatter errors are visible.

This completes the stated bulk mass connection and its targeted CPU
reference checks, not all physical programme items listed above. No
calibration simulation, new audit bundle, or source-table fabrication was
used to close unresolved physical data requirements.

### Final arithmetic, MPI and compatibility follow-up

An initial proportional rescaling of all three overbudget H/He shares
passed CPU execution but failed the existing CPU/CUDA inventory comparison
(absolute1.76326e-6 in a nearly exhausted atom inventory). Changing CUDA
FMA contraction alone did not fix it. Replaced that correction with the
minimal downward-ULP adjustment to the largest share only, until the
FP32 sum is within its original target. This avoids perturbing all three
residual inventories and their existing threshold branches. No comparison
tolerance was loosened. CUDA primary compilation now explicitly disables
FMA contraction, consistent with the other shared material/IR kernels.
Both CPU backend checks and the actual CUDA comparison on visible GPU1
pass; maximum normalized difference1.09567e-7, including the original
zero/nonzero-dust backend cases. The4096 fully absorbed zero/trace-dust
direct-cell cases pass as well. This is a kernel comparison, NOT a claim
that CR hydro or the new mass update has gained CUDA dispatch.

`hot-mpi2/physical.nml`, same thermal-only2-step case and output policy,
finishes with two MPI/two OpenMP in6.737s. Global512 leaves, split32 grids
per rank. Density, gas total/CR energies, total metal, dust mass and dust
energy are exactly equal to the one-rank uniform result; dust identity
matches. No new MPI framework or domain-decomposition gate was introduced.

Final NENER1 mass executable SHA256:
`ae65c23ae49ffa07dc3e5eff776ee562ef881f3da821dd67c5d483a78c199943`.
Final NENER0 CPU/SNRT/DUST_LIVE compatibility build succeeds, executable
`.snrt-cpu.OKoz9T/ramses_dust_compat3d`, SHA256
`ed91c93a0d9f89a36d16d41711323ad5437a0911bec97209c6de343e48538fff`.
Compatibility here is a build check, not a new NENER0 evolution run.

The final shared-cell version was rerun from scratch in `live-final` and
from a COPY of its step2 in `restart-final`; old evidence remains intact.
Same declared4-step1-MPI/2-OpenMP namelist/output policy as live; initial
dust zero. Live4-step completion5.129s. Detailed final dataset comparison
follows below. The entire scratch root is about842MB, including preserved
first attempts and checkpoint copies, not a production simulation archive.

Final restart completes in2.547s: all90 hydro datasets,4 SNRT datasets,
and15 particle fields for2047 ID-matched stars exactly equal continuous
execution; dust model identity equal. Final mean dust density
4.282749639924932e-11, mean dust energy5.224763797576013e-22, maximum
D/Z5.032097723348812e-6. All hydro finite and dust within metal bounds.
`git diff --check` passes. No commit/push was requested in this turn.
