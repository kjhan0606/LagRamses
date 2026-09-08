# Medium-term dust: ordered implementation plan

**Current implementation status (2026-09-09 closeout):** native composition,
two-size evolution, DL01 material, and local D03 primary/IR optics are now
connected and tested together. See the final section below; earlier
"prepared/not live-selected" entries are historical. This closes the live
comparison-model connection, not the separate full-NEQ/multibin physical
qualification or cosmological production/calibration work. The subsequent
size-shift reference section records the first measured 8/16/32-bin comparison;
the earlier "unperformed" description is historical, not its current status.

Operator approved proceeding in this order after closeout commit 5be56ae.
Workspace / repository: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses.
Final objective: native, scientifically defensible RT/feedback/dust galaxy
simulation physics, followed by model selection and separately approved
COLIBRE-style calibration. Do not turn this into instrumentation work.

## Bundle 1: composition, shared element budgets, gas-phase cooling

Retain bulk_v1 unchanged and default-off. Add an explicit comparison model
with carbon and MgFeSiO4 olivine (and evaluate a separate Fe carrier without
double-consuming silicate Fe). Condensation uses actual gross channel ejecta,
AGB number-ratio C/O availability and limiting stoichiometric elements, not
total-Z fractions relabelled as species. Coefficients are named assumptions,
not a claim to reproduce any paper's calibrated dust yield.

Transport composition as native passive densities after existing dust fields;
total element fields continue to mean gas+dust. Derive gas-phase element
abundances by subtracting the stoichiometric locked masses. Growth and
sputtering preserve those budgets; source injection participates in existing
cell lock / reverse MPI / particle-progress transaction. Keep aggregate dust
mass and material energy consistent with constituent masses. Check SF removal
and restart identity; no new generic checkpoint framework. No silent bulk-to-
species initial/restart conversion. Fixed optical/material mixture remains a
labelled comparison approximation until bundle 3.

The actual current cooling consumers (cooling_fine.kjhan and Eunha tables)
take scalar total Z, not element-resolved rates. Locate source-backed individual
element rates before calling gas-phase cooling complete. Exporting abundances
or scaling a solar-mixture curve by residual total Z is NOT equivalent. Keep
unsupported external cooling combinations rejected until a valid receiver is
connected; report this dependency honestly rather than inventing rates.

Use existing dust native smoke plus one small native coupled/restart case
for changes that need runtime evidence; no per-item audits or new test engine.
Namelist changes update mkrun and shared generator together. Preserve unrelated
generator deletions and all existing outputs. Free CUDA stream -> GPU, otherwise
OpenMP remains the dispatch design; no new autotuner.

## Bundle 2: size distribution / growth / destruction

Start with two-size composition bins and jointly implement coagulation,
shattering and ambient SN-shock processing, distinguished from fresh-ejecta
survival. Compare against a bounded piecewise-linear 8/16-bin reference only
when choosing accuracy/cost; do not implement every solver alternative fully.
Moments remain a candidate if they can reproduce needed observables. Preserve
element/mass budgets and assess SF-compatible unresolved growth assumptions.

## Bundle 3: RT / thermochemistry / model selection

Connect composition and size to opacity, scattering and material/IR response;
H2 grain surfaces/shielding and CO competition need their own source-backed
chemistry, not an inference from CR pressure. Evaluate dust drift/drag and AGN
force separately, not as a prerequisite for all other work. Compare numerical
methods at fixed physics, then physics at fixed method; select adequate lowest
total cost. Final calibration awaits model choice and separate approval.

## Planning review request

Fable: read-only review focused on feasibility, necessity for the final
objective and overinstrumentation. Assess shared budgets, active-runtime
connections, cooling-data dependency and the fixed-opacity intermediate scope.
Give a verdict and only actionable relevant conditions. Do not edit files,
run jobs, contact anyone, add extra gates or demand unrelated AMR infrastructure.
The driver independently evaluates findings. End-of-bundle review is by the
driver under the operator's latest instruction, not a new external audit chain.

Sources already inspected: [RAMSES](https://arxiv.org/html/2402.18515v2),
[COLIBRE](https://arxiv.org/html/2505.13056v2). Their prescriptions motivate
comparisons; neither validates the existing local coefficients or supplies
an element-resolved cooling table to this worktree by itself.

## Initial inspection / implementation choices

- Start with two noncompeting carriers, pure C and MgFeSiO4. Fe is explicitly
  tracked within olivine; a separate metallic-Fe population is not silently
  assumed. Both carriers share the existing radius/density/sticking reference
  until composition-specific rates are selected. No unknown oxygen carrier.
- Condense source-node release segments before age/Z/IMF mixtures. AGB CO
  formation is based on number abundances (C/12 versus O/16); CO-bound excess
  is not dust and remains in the total element/gas budget. No SSP-averaged
  C/O switch. Existing AGB jump/left-limit semantics must also apply to dust.
- Existing idust/idust_energy stay in place. Append carrier densities, retain
  total elements, and reject inconsistent old scalar seeds/restarts rather
  than inventing a composition. The current NVAR=30 comparison has headroom.
- SF already converts all imetal:nvar passives to fractions and back with
  remaining gas density. Thus the appended carrier mass should follow the
  same astration path; no new star-particle dust reservoir is required.
- Scalar cooling is confirmed in cooling_fine.kjhan and Eunha; public
  [Oppenheimer/Schaye element tables](https://noneq.strw.leidenuniv.nl/) were
  located but the HM12 tar endpoint returned HTTP403, including direct curl.
  [Hybrid-CHIMES](https://www.sylviaploeckinger.com/hybridchimes) links its
  public table-generation pipeline but its download section still says
  "Coming soon". Neither is claimed available/qualified in this worktree.
  Other official CHIMES code/data routes remain candidates. No contact/email
  or large table-generation calculation was initiated.
- Fable planning call failed before review with HTTP429 usage limit, session
  edde6cd4-5083-4875-bc87-101c52d85071. This is not a negative verdict. The
  existing backup reviewer, claude-opus-5, was invoked read-only instead.

## Planning review disposition

Opus returned **proceed with conditions** (read-only session
c119b6a7-b13a-4513-8da1-8c1518f16f0c, actual model claude-opus-5;
approximately 321 seconds / 46 turns). No further per-item audit was added.

- Accepted: keep returned-minus-H/He aggregate Z; avoid reusing olivine Fe;
  one level collective for both carriers; extend the existing dust smoke,
  not a separate testing framework.
- Accepted with explicit scope: an opted-in depleted scalar solar-mixture
  cooling receiver is useful intermediate work. It does NOT finish the
  individual-element cooling requirement or validate a unified NEQ closure.
- Applied condensation on physical source segments, earlier than the
  review's aggregate-deposition suggestion: mixing AGB C/O first is nonlinear
  and can erase a carbon-rich source contribution.
- Corrected a review assumption: increasing used fields does not necessarily
  change NVAR (our NVAR30 already has spare fields), so a separate small
  composition/cooling identity is necessary. No generic new restart framework.
- Retained the operator's 1 -> 2 -> 3 ordering. Did not bring a speculative
  opacity layer forward just for planning approval. Fixed mixture optics
  still respond to total dust, but not directly to the new constituent ratio.

## Bundle 1 native implementation and driver review

Implemented `carbon_olivine_v1` alongside unchanged/default `bulk_v1`:
source-node segment condensation (including AGB terminal jump/left limits),
age/Z/IMF propagation, deposition in the existing source transaction, native
passive C/olivine fields, element-limited growth/sputtering, material/gas
energy exchange, and composition/cooling restart identity. All tracked
elements remain total gas+dust; no untracked metal mass is discarded.
No separate Fe carrier, rate calibration or evolving size is claimed.

The actual `cooling_fine.kjhan` consumer now has opt-in `depleted_scalar`,
original solver only, without HM/J21/self-shielding backgrounds. Inspection
found J21=0 still retained the original electron/UV floor: this comparison
now zeroes J0min/J0min_ref after legacy model initialization, before table
creation. A native trial then exposed log10(0) in the original table's zero
photoheating and derivative. The table now uses the existing metal-zero
log convention 1e-100 for zero rates/species; not a new physical heating term.
Its H/He equilibrium thermal closure remains an explicit approximation
beside SNRT's photoionization/recombination, not element-resolved NEQ cooling.

An actual step-4 failure also identified independently limited aggregate
dust transport: D=8.229557581863915e-12 versus constituent sum differing by
7.788341675712907e-22. The fix shares the constituent face fluxes and
prolongation with their dependent aggregate; tolerance was not widened and
there is no post-step mass clipping/creation. All trials, including failures,
remain in `.dust-composition.BVPvN3/`.

### Bounded execution evidence

- Existing `dust_mass_smoke`: PASS, including prior bulk checks, source
  condensation, C/O premixing counterexample, telescoping, element budgets,
  constrained growth and rejection of overdrawn carbon.
- Existing GUI suite: 34 tests, 33 passed and one display-only skip. mkrun
  and shared generator updated together; defaults retain bulk/none.
- Final CPU native binary:
  `.cosmic-ray.kyySgK/ramses_dust_composition_zero_uv3d`, SHA256
  `ffb033cae9f399dbacb6d0a62b506a5218efdf9aab72109d01472f9c87c3f312`.
  Build command from `.cosmic-ray.kyySgK`:
  `make -f ../bin/Makefile -j4 SNRT=1 DUST_LIVE=1 HDF5=1 USE_CUDA=0 USE_FFTW=0 NENER=1 EXEC=ramses_dust_composition_zero_uv ramses`.
  NDIM3/NVAR30/NVECTOR500; lagRamses precedes cuRamses in VPATH. Objects
  and link identify the edited cuRamses cooling/hydro implementations.
- Effective nmls `live-zero-uv/physical.nml` and
  `restart-zero-uv/physical.nml` below the evidence root: fixed level3 8^3,
  one MPI rank / two OpenMP threads, four steps, CR SN fraction .1 and SF
  model4, zero dust seed, 1408-row physical yield input and independently
  labelled BPASS radiation comparison. No AGN/sinks/cosmological expansion.
  Environment sourced from `.dust-mass.w1Jo9A/live/physical.env.sh`.
- Both logs end `Run completed`; step2/4 COMPLETE outputs preserved.
  Final leaf hydro fields finite, gas-phase elements positive (minimum
  7.425104090420321e-13 code density), thermal energy excluding CR/kinetic
  positive (minimum 9.076620817839741e-08). Carbon density range
  [5.5791480199648e-12,7.892724050737332e-12]; olivine
  [8.092380872721322e-13,1.2023826832574808e-12]. Maximum aggregate/sum
  relative discrepancy 3.0366034771209613e-16; max dust/total-metal
  9.09698028532559e-07. These are short, low-condensation tests, NOT a
  saturated-depletion or production-abundance result.
- Saved original cooling table: all 12 rate/derivative/species arrays finite,
  161 density x 101 temperature nodes. Log photoheating exactly -100 and
  its log derivative exactly 0 everywhere. The actual cooling call uses
  depleted Z, but this trial does not quantify a significant dust-induced
  thermal change or validate individual-element metal rates.
- Restart: all 90 hydro datasets and all four SNRT datasets bitwise equal;
  15 particle fields bitwise equal after ID alignment for 2048 particles.
  Header composition identity `[1,2,1,22]` present. No broad AMR/distributed
  infrastructure review inferred from this fixed-grid comparison.
- Negative restart `restart-wrong-model/physical.nml`: unchanged NVAR30,
  only carbon_olivine_v1 -> bulk_v1. Rejected during HDF5 restore with
  `dust mass-evolution restart identity mismatch`, MPI abort 11, before
  time evolution or any new output_00002. Original checkpoint untouched.
- Existing noncosmological SFRD text reports NaN/Infinity from its volume
  diagnostic; this is not a NaN in stored hydro, particles or cooling fields.
  No unrelated reporting-system rewrite was added.

### Remaining ordered work (not closed by this evidence)

Bundle 1's physical individual-element cooling rates/receiver are still open;
the scalar option is explicitly intermediate. Fixed mixed optics and shared
fixed-size growth constants likewise remain assumptions, not calibrated
composition physics. The dependent aggregate flux fix is compiled for AMR
prolongation, but these native trials qualify only the stated fixed-grid
CPU comparison, not arbitrary AMR/GPU hydro or strong-depletion transport.
Bundle 2 (size distribution) and bundle 3 (composition/size RT chemistry and
model selection) have not started. Do not relabel this as completing all
three bundles or as production-ready medium-term dust physics.

No calibration, external contact/email or large data-generation job launched.
Pre-existing unrelated generator deletions and prior binaries/outputs were
preserved. This implementation is currently uncommitted.

## Continuation: individual-element CIE receiver

Operator requested continued implementation, then explicitly asked about
LTE/NLTE cooling/heating differences. Keep three distinctions separate:
LTE/NLTE level populations, CIE/photoionization equilibrium, and equilibrium/
time-dependent ionization. A low-density CIE table is NOT an LTE model and
does not provide the final radiation-dependent NEQ closure just because
SNRT operates alongside it. Shared RT/ion/electron/rate state remains a
bundle-3 requirement; no user approval of CIE as the final model was inferred.

The blocked previous HM12 download was not the only data route. Retrieved
the author-published [WSS09 CIE ASCII table](https://local.strw.leidenuniv.nl/WSS08/z_collis.txt)
from the [paper's data page](https://local.strw.leidenuniv.nl/WSS08/).
Local numerical source has 352 x 25 entries, T=100--959070000 K; file
SHA256 `acbea8c6faf8b316debada411b0d433c05036056b35293748e3b5eda267b47b0`.
Comments/numerical tokens retained, final newline normalized. Embedded
Fortran import SHA256
`3449524be2e849849aa63b5d185e245e9dff1638607d953750e5820e2f720f2f`.
This is a CLOUDY 07.02 comparison source, not a claim to use the latest
atomic database or to have selected the best final cooling model.

Implemented `dust_cooling='wss09_cie'`, composition mode only. Native
cooling calls `dust_gas_elements` and passes actual gas elemental fractions,
mixture density, T/mu, gamma and physical dt to `wss09_step`. It replaces
the original cooling call; there is no duplicate scalar metal cooling.
H/He rates/electron density interpolate in the supplied He/H nodes; nine
metal contributions are weighted by actual gas-phase number abundances as
in [WSS09 equation (3)](https://arxiv.org/html/0807.3748v2).
Equation (4)'s solar electron density is not in this ASCII table: do not
invent it from a H/He-only column. Metal electrons are a labelled trace-
element approximation in the energy/temperature mapping. Untracked metals
remain in total mass; no spurious cooling species is assigned to them.

Numerical update: exact integration of a signed, piecewise-linear net-rate
curve in T/mu at fixed composition/density. It traverses table intervals,
retains heating and equilibrium zeros and rejects departure from published
T/He support. No new timestep tuning, benchmark framework, per-item audit,
external job or automatic model switch. The source table has 343 negative
individual-metal entries; those were not clipped away. No optically thin
cooling energy is secretly injected into the dust IR ledger.

### Continuation evidence

- Built on lageunha in the existing NENER1 CPU build directory (lagRamses
  VPATH first), EXEC=ramses_dust_cie. Same HDF5/SNRT/DUST_LIVE/NVAR30/
  NVECTOR500 options. Binary SHA256
  `c98c32c4ab7f23f5a9144a022906accf35ead187895d3b4b3316bdc4c17372f9`.
- Extended the existing dust_mass_smoke only: equal-Z C/Fe contrast,
  linear removal of locked C from its elemental rate, signed low-T heating,
  out-of-domain rejection and full/half timestep agreement PASS. Independent
  DOP853 integration of the source-based same interpolant returns
  T/mu=146424.0602205506 versus native 146424.06022055328,
  relative difference 1.83e-14 (initial 151913.77815298695, rho=1.66e-27,
  dt=1e13 s, XH=.74/XHe=.25/XC=.001/XFe=.009). This establishes the
  integrator, not CIE's astrophysical validity under a radiation field.
- GUI suite remains 34 tests: 33 PASS, one display-only skip. Defaults
  remain bulk/none; invalid CIE+bulk is rejected. mkrun creates the selected
  namelist and selects the matching native binary, without launching jobs.
- Native `.dust-cie.5yWa6S/live/physical.nml` and `restart/physical.nml`:
  periodic noncosmo level3/8^3, 1 MPI x 2 OMP, four steps, physical source
  input and independently labelled BPASS RT, CR/SF and zero-seed composition
  dust. Both Run completed; continuous/restart step4 agree bitwise for
  all 90 hydro arrays, four SNRT arrays, 15 ID-aligned particle fields
  (2047 particles), and all three composition/CIE identity attributes.
- Final fields finite; gas-element minimum 7.113364822600874e-13 code
  density, minimum thermal energy excluding CR/kinetic 1.06262817119515e-07;
  D=sum(C,silicate) maximum relative discrepancy 3.47351469292606e-16,
  maximum D/Z=8.784143677570757e-07. This is still a short weak-depletion
  source test, not full grain-growth saturation or dense-gas qualification.
- Compiled/checkpoint H/He and used individual-metal numbers match the
  downloaded numeric source bitwise. One unused total-metal column entry
  (1.5932e-22) differs by one ULP in Intel's decimal literal conversion;
  no physical table correction or widened mass-budget tolerance was made.
  Restart binds the exact compiled data (5310+3520 doubles in two attributes
  below the individual HDF5 attribute-size limit), not just the filename.
- Original successful and failed comparison trials, binaries and pre-existing
  unrelated generator deletions remain untouched. Current changes uncommitted.

The previously open **individual-element data-to-native-receiver connection
is now implemented for this CIE comparison**. It is not full NLTE/NEQ or
publication approval. Bundle 2 size distribution and bundle 3 consistent
radiation/thermochemistry remain in the approved order; do not quietly move
the latter into "done" based on this short CIE run.

## Bundle 2 implementation and driver review (2026-09-08)

Implemented native `carbon_olivine_2size_v1`: four mass carriers, common
composition reservoirs, size-dependent growth/thermal erosion, coagulation,
shattering, and optional ambient SN destruction protected against fresh
same-step ejecta. `dust_mass_physics`, native runtime, stellar deposition,
CPU MUSCL face fluxes, AMR prolongation, read_hydro_params, HDF5 identity,
mkrun, shared generator and its existing GUI test are connected. No new
Python validation framework, unrelated AMR audit or per-item external audit.
The unrelated pre-existing generator deletions were preserved.

Implementation/admission details and numerical limits are in
`simulation/snrt/NATIVE_RUNTIME.md` under "Two-size carbon/olivine comparison".
The primary model reference is
[Dubois et al. 2024](https://arxiv.org/html/2402.18515v2). This is an explicit
hybrid of that two-size collision/SN reference and our retained condensation,
effective-atom accretion and Tsai-Mathews erosion assumptions, not a claim to
replicate every prescription or calibration in the paper.

### Actual native wiring

- NENER1/virial layout: D20, Ed21, C22, silicate23, four bins24--27.
  SN-on only: energy28, fresh C29/silicate30. NVAR30 suffices. Native
  admission rejects insufficient NVAR and nonzero transient uniform ICs.
- Source-node composition still precedes age/Z/IMF mixing. Only the size
  split is new at native deposition. All rows enter the same source-progress
  transaction/cell lock and reverse MPI sum. No nonlinear destruction is
  performed separately on MPI ranks or separately for each star.
- Owning leaves process ambient SN exposure after hydro/cooling, then
  size evolution, then existing SNRT. The mass stage exchanges Ed with gas
  thermal energy, derives aggregates, clears transient rows and refreshes
  restriction/halos. Its existing one-per-level MPI validation is reused.
- Godunov flux and prolongation derive composition from the four size bins
  and total D from composition. Existing SF removes all passives with gas;
  it sees zero transient fields. This is wiring review, not a renewed generic
  AMR infrastructure certification or a dust-drift implementation.

### Evidence

- Built on lageunha in `.cosmic-ray.kyySgK` using existing Makefile/VPATH
  (lagRamses first): SNRT=1, DUST_LIVE=1, HDF5=1, USE_CUDA=0, USE_FFTW=0,
  NENER=1, NVAR30, NVECTOR500. Final `ramses_dust_sizes3d` SHA256:
  `76b286c10702fe0fb07ce201bfb57c17b71fac13f180651b4f88e5205eebe316`.
  Build log `dust_sizes_final_build.log`; native smoke `dust_sizes_smoke.log`.
- Existing Fortran smoke PASS: two-size source split, exact coagulation and
  shattering donor solution/composition conservation, zero-transfer identity,
  common-reservoir growth bounds and timestep split, analytic radius-dependent
  erosion, invalid parameter rejection, analytic SN survival, fresh-only
  protection, fixed-environment exposure composition and excess-fresh rejection.
  Previous bulk/composition/CIE tests also PASS. GUI suite 34 tests:
  33 PASS, one display-only skip; no new test suite.
- `.dust-sizes.cn4sDJ/live` and `shock` are preserved initial one-MPI/two-OMP
  trials. Both completed four steps. Final dust differed by -1.783677265e-5
  fraction (-0.00178368%) with SN destruction, a weak-processing smoke, not
  a physical destruction-efficiency measurement. Initial binaries preserved
  as `ramses_dust_sizes_noshock_initial3d` (SHA256
  `8bac1f110fb4c0920099d4a216e7db18333191d5792d5861a77505de00cc4c65`)
  and `ramses_dust_sizes_sn_initial3d` (SHA256
  `eb09aa451748342fecbfaeffc08b0df63bc7d3104329972e142b34a461a8858e`).
  Final integrator avoids changing an untouched tiny receiver by subtracting
  it from a much larger total when the transfer is zero.
- Final native pair: `/gpfs/kjhan/LRD_JWST/.dust-sizes.cn4sDJ/final/physical.nml`
  and `final-restart/physical.nml`, 2 MPI x 2 OMP, periodic noncosmo 8^3,
  levelmin=max=3, four steps, SN-on/CR/SF/physical KL16+LC18+effective Ia,
  independently labelled BPASS radiation, WSS09 CIE, zero initial dust.
  Both Run completed (about 2.99 s continuous, 1.46 s restarted).
  Native leaf grid ownership is 32+32. Both source ranks participate.
- Output policy verified before launch: noutput1/aout2/tout1e30/foutput2/
  fbackup1000000; observed two snapshots 55/67 MiB, free space 159 TiB.
  Original snapshot copied to a new restart directory; none deleted.
- Continuous/restart final state is **bitwise equal** for all 90 hydro
  datasets, four SNRT datasets, 15 ID-aligned particle fields (2047 stars),
  and five dust/composition/size/CIE identity attributes. Final material/size
  identity includes all radii/densities/injection/flags and field indices.
- Final minimum gas-phase element density 7.113364822600874e-13,
  minimum thermal energy excluding CR/kinetic 1.0626281711951004e-7,
  maximum D/Z 8.783979677376048e-7. Maximum aggregate discrepancies:
  D vs C+sil 2.3675e-16; C vs its bins 2.7371e-16; sil vs its bins
  4.3358e-16. All four bins nonnegative, all fields finite, transient fields
  28--30 identically zero. Small bins are produced from initially all-large
  injection; final summed bin densities are approximately
  [2.32686e-18, 3.32745e-9, 2.46651e-20, 4.99426e-10].
- git diff --check PASS. No commit or push requested in this continuation;
  source changes remain uncommitted together with earlier composition/CIE work.

### Driver verdict / remaining model work

PASS for native two-size **comparison wiring and the tested short MPI/restart
case**, not blanket production/publication approval. No extra approval gate
is being inserted. Bundle 2's listed native operators are implemented; the
ordered next implementation is bundle 3's size/composition RT/material and
consistent thermochemistry connection. The accuracy/cost 8/16-bin comparison
is still open for that model-selection decision, not falsely marked complete.

Retain explicit limitations: unresolved density/PDF and SF-Mach coupling are
not assumed; CR is excluded from thermal growth temperature; collision model
is within-composition and dense-gas coagulation is resolved-only. Effective
metal mass and shared sputtering fit are retained approximations. SN energy
equivalent exposure is default-off, not physical event counting or calibrated
overlap handling with thermal sputtering. Four bins still use the old fixed
mixed optical/material receiver. CIE versus time-dependent RT/NEQ is not
resolved here. These are physical model-selection/connection tasks, not a
reason to proliferate generic tooling, per-step audits or unrelated gates.

## Bundle 3 continuation: composition material / collision geometry (2026-09-09)

Current project and origin reconfirmed: `/gpfs/kjhan/LRD_JWST`,
`git@github.com:kjhan0606/LagRamses.git`, main at 5be56ae, previous uncommitted
work preserved. Initial login-shell commands stalled; non-login shells and
lageunha access worked. No workspace move, rsync, shared-context rewrite,
external audit, commit or push was performed.

### Physical implementation

Added explicit `dust_material_model='dl01_composition_v1'` for the two-size
model. Default fixed_mix and its arithmetic remain. Existing DL01 bulk
material builder emits separate graphite/silicate curves (162 temperatures,
5--300 K); native code mixes U(T) by local C/silicate dust mass. Bulk energy
per gram has no radius dependence; do not fake separate small-grain heat
capacities. This is one common temperature, not grain-resolved stochastic
heating. Source data/derivation reuse the existing
[Draine & Li 2001 bulk limit](https://arxiv.org/abs/astro-ph/0011318).

Same U(T) now connects source injection, mass-evolution energy exchange,
old-temperature decoding and the implicit gas/dust/IR material solve.
Mass change preserves the old common T but accounts for the changed
composition's specific energy; Ed is no longer simply scaled by D_new/D_old
for this option. Gas receives the opposite energy transfer. Gas collision
area uses all four masses, radii and solid densities; no optical Q is used
as a geometric area and CR energy is not treated as gas thermal energy.

Per-cell U(T) was added as an optional final callback input to the existing
IR/material path, with no new transport solver. Shared C++/CUDA cell code,
device byte counts and hybrid batch packing were updated together. Mode bit
2 selects appended cell-specific U data; lower bits preserve previous
material/exchange behavior. The CUDA lease/busy-CPU policy is unchanged.
No new MPI collective or audit gate. Namelist admission, mkrun and shared
GUI generator are updated; restart binds the complete compiled material
tables in one small attribute.

### Verification

- Native build: `.cosmic-ray.kyySgK/ramses_dust_composition_material3d`,
  CPU/NENER1/NVAR30/SNRT/DUST_LIVE/HDF5, lagRamses VPATH first. SHA256
  `221e8a5f869e7f2d5fbdc8c011c811de29d91da60d1cd0a7947f1898f62883c6`.
  `dust_dl01_composition_data.inc` SHA256
  `37562da8135aa24b2cf5ecfefb4a05f730737e178daecb4db88f145a3779c54d`.
  Build log `dust_composition_material_build_final.log`. The initial build
  caught two missing field imports; those were fixed before the successful
  final native build/launch.
- Existing dust_mass_smoke PASS for pure-material contrast, linear mixture,
  small/large area ratio 20, domain rejection and all previous mass/CIE tests.
  Existing backend smoke extended for heterogeneous cell U across hybrid
  batch boundaries, with/without gas exchange, against per-cell ordinary
  material calls. All prior conservation/rollback/parity checks also PASS.
  An initial synthetic test had gas at 1000 K but an artificial material/
  emissivity table bounded at 100 K: correctly rejected, not a backend defect.
  The accepted bounded test uses 80 K gas; no physical table limits widened.
- Backend smoke: CPU-only hybrid PASS. Also built the existing GPU smoke
  in `.ir-hybrid.dYEXir` (preserved prior executable as
  `snrt_dust_backend_smoke_pre_material`) and ran on the lageunha RTX 5000
  Ada. Auto mode actually used CPU=3/GPU=2 material batches; forced CUDA
  PASS. Baseline backend/reference discrepancies <1.9e-15; heterogeneous
  material batch/scalar comparisons pass the declared 1e-11 tolerance.
  Logs `composition_material_gpu_build.log`, `composition_material_gpu_auto.log`
  and `composition_material_gpu_forced.log`. This is material-kernel evidence,
  not a full GPU hydro or multi-GPU scaling qualification.
- Effective native namelists:
  `/gpfs/kjhan/LRD_JWST/.dust-composition-material.jR6vJD/live/physical.nml`
  and `restart/physical.nml`. Two MPI x two OMP, periodic noncosmo 8^3,
  four steps, real KL16/LC18/effective-Ia feedback, CR/SF, SN destruction,
  WSS09 comparison and independently labelled BPASS RT. Zero dust ICs.
  Both Run completed; continuous wall time about 3.29 s.
- Launch policy checked: noutput1/aout2/tout1e30/foutput2/fbackup1000000,
  free 158 TiB; actual snapshots 55+55 MiB. Snapshot 1 copied to a new
  restart directory. Original outputs/binaries preserved.
- Continuous/restart state bitwise equal: 90 hydro datasets, four SNRT
  datasets, 15 ID-aligned particle fields for 2047 stars, six material/dust/
  CIE identity arrays (plus two unchanged field-index attributes).
  Final common T_d=10.250163828157014--10.293097264990259 K;
  minimum thermal gas energy 1.0626281710885389e-7, minimum gas-phase
  element density 7.113364822614062e-13, max D/Z=8.78397967737597e-7.
  Aggregate relative differences <4.30e-16, all masses nonnegative,
  all fields finite, transient source rows exactly zero.
- GUI suite 34 tests, 33 PASS/one display skip; git diff --check PASS.
  Source changes and generated material data remain uncommitted.

### Optical data limitation / next work

Inspected the official
[optical grain data page](https://www.astro.princeton.edu/~draine/dust/dust.diel.html)
and fetched/decompressed Gra_81.gz and suvSil_81.gz in memory. Actual files:
81 radii, 241 wavelengths, 0.001--1000 microns; about 1.24 keV upper energy.
This does not cover our 2--10 keV group (representative energy 4023.59 eV),
nor all the long-wavelength IR nodes. Existing WD01/D03 whole-mixture optics
remain active. No spectral group removed, unmeasured tail zeroed, hidden
mixture switched, or extrapolated cross section invented.

Driver review: PASS for the implemented common-T material/geometry option
and tested MPI/restart/backend paths. **Bundle 3 is not complete.** Next:
source-backed full-band size/composition optics and consistent emissivity,
then the radiation/thermochemistry closure and accuracy/cost model selection.
Obtaining only a UV-sized table is not enough to close the hard-X/IR path.
Stochastic small-grain heating, common-T adequacy and CIE versus evolving RT
ionization remain explicit physics/model-selection limitations. No additional
approval wait has been inserted and no broad production approval is claimed.

## Bundle 3 continuation: D03 native optical coefficients (2026-09-09)

### Corrected source assessment and implementation

The earlier 2 keV limitation was inferred from the old Draine web-page text.
Reading the actual `callindex.out_*D03*` files shows usable graphite data
through 18 keV and silicate through 19 keV, before anomalous final rows.
All nine primary groups including the 10 keV upper edge and all 136 current
IR nodes are covered without spectral extrapolation. The separate old
Gra_81/suvSil_81 cross-section files still have their reported narrower range.

Five originals are retained in `simulation/snrt/data/draine_d03/`; exact
SHA256s, byte counts and usable ranges are recorded in
`simulation/snrt/data/dust_d03_optics_generation_v1.json`. Primary sources:
[Draine dielectric data](https://www.astro.princeton.edu/~draine/dust/dust.diel.html),
[D03 X-ray/material model](https://arxiv.org/html/astro-ph/0308251).
The model densities are **2.2/3.8 g/cm3**, not the current evolution defaults
2.2/3.3. The supplied graphite files are for **0.01/0.1 micron at 20 K**,
not the default 0.005/0.1 micron bins. Existing defaults were not changed.

`build_d03_grain_optics.py` performs offline sphere Mie calculations using
miepython 3.3.0, installed in a local `.dust-optics-venv`. There is no new
simulation-time Python dependency. Re(n)-1 is linear in log E and Im(n)
log-log interpolated, with source edge knots as anchors; no extrapolation,
negative-coefficient clipping or fabricated tail. The graphite 1/3--2/3
approximation is explicit. Its angular moment is scattering-weighted, not
an unweighted average of g. Size parameters remain below tested x=10000.

`dust_composition_optics.f90` plus a generated 64 KiB include are in the
native Makefile object list. Four bases produce mass-weighted absorption,
scattering/g and primary absorption optical depth. Units use Q*3/(4*rho*a),
either cm2/g or cm2/reference-H. Binding rejects different edges, primary
representative energies, IR nodes, sizes or solid densities. Zero dust gives
zero opacity; scattering is not charged as heating.

Each of four pure populations and a mixed population was passed to the
existing native IR table initializer and implicit advance: **the same IR
absorption coefficient** enters attenuation and Kirchhoff emission, with
matching DL01 composition U(T). No separate solver/test framework added.

**Not yet live-selected:** the driver and persistent live IR still use the
old fixed-mixture optical contract. This is a built native coefficient API
and tested receiver integration, not a heterogeneous AMR evolution. No
namelist/mkrun option advertises an unconnected model, and no checkpoint
identity changed. The previous material/geometry path remains intact.

### Evidence

- Existing `draine_dust_opacity.py --d03` PASS: published Mie amplitudes,
  Rayleigh limits, raw hashes, full-group support, invalid/domain rejection,
  graphite angular weighting and generator/compiled identities. Its original
  sidecar entry remains; unrelated JAX is imported only for that old entry.
- Regeneration from cached originals gives bit-identical include and manifest.
  Include SHA256 `37c462cec715f0e3009e277249558ab364246c78e373d4049d6e88605f284b54`.
- Existing `snrt_dust_backend_smoke` extended for five mixtures, mass-unit and
  binding checks, closed-cell energy conservation and stationary 10 K bath.
  Reference/backend differences <=6.09e-16, tolerance 1e-10. Old material,
  gas exchange, transport and rollback checks also PASS.
- CPU: `.cosmic-ray.kyySgK/d03_optics_cpu.log`. GPU:
  `.ir-hybrid.dYEXir/d03_optics_gpu_auto.log` and `d03_optics_gpu_forced.log`.
  Auto actually used CPU=3/GPU=2 on the existing 1031-cell material case;
  forced CUDA passed all new physical cases. This is kernel evidence, not
  full MPI evolution or multi-GPU qualification.
- Native CPU link `.cosmic-ray.kyySgK/ramses_dust_optics_prepared3d`,
  SNRT/DUST_LIVE/HDF5/NENER1/NVAR30, lagRamses VPATH first, SHA256
  `6977b9ea2cd6559bb79d585ddecea4da5882476f61578768d0c217875ad3af53`.
  No lagRamses evolution launch or new snapshot this turn.
- `git diff --check` PASS; all previous dirty work/unrelated generator
  deletions preserved. No commit/push, external audit or shared-context edit.

### Next implementation boundary

Switch primary absorption, scattering and live IR absorption/emission together
using cell-local four-bin weights and shared optical bases in existing hybrid
paths. Bind optical data and physical size/density selection on restart;
expose namelist/mkrun choices only when the path is complete. Avoid a
cell x temperature x frequency emissivity cube just to mix four components.
Keep angular information: computed group-9 g=0.99986--0.999998, so full
Qsca used isotropically is not an accurate X-ray transport closure. Any
isotropic comparison must remain explicitly labelled; physical angular
transport remains necessary. Common temperature, frozen 20 K optics,
graphite versus PAHs, and radiation-dependent metal NEQ remain explicit
limitations. Bundle 3 is open; no new approval wait or production claim.

## Bundle 3 live optical connection / driver closeout (2026-09-09)

### Implemented native path

`dust_optics_model='d03_transport_v1'` now connects all four evolving mass
populations (C-small/large, silicate-small/large) to primary absorption,
primary transport scattering, IR absorption/emission and IR scattering.
Requires `carbon_olivine_2size_v1`, `dl01_composition_v1`, radii 0.01/0.1 micron
and densities 2.2/3.8 g/cm3. Defaults remain the existing fixed mixture.
The nine edges/representative energies, 136 IR nodes and material dimensions
are checked at existing startup preflight; there is no silent extrapolation.

Scattering uses the explicit **delta-isotropic transport approximation**:
each basis contributes Qsca*(1-g), never Qsca as isotropic opacity and never
scattering energy as heat. This preserves photon number/weighted IR energy
and the first angular moment for the symmetric native quadrature, not the
full Mie phase function, resolved X-ray halo, or higher angular moments.
IR uses the existing FP64 shared scattering kernel through a layout-aware
wrapper: intensities -> weighted bins -> scatter -> intensities. It is
Lie-split after IR transport/absorption/material exchange within existing
substeps and the existing level transaction. No extra MPI gate was added.

Kirchhoff absorption/emission uses the same local absorption coefficient.
The joint gas/dust/IR solve uses the same local DL01 U(T) and geometric
collision area as injection and mass evolution. Four shared optical bases
plus four weights/cell feed the existing free-stream-or-OpenMP scheduler;
no separate solver or cell x temperature x frequency cube. The 222-temperature
bases occupy 0.928 MiB instead of a 237.5 MiB 1031-cell emissivity cube.

HDF5 `dust_optics_d03` contains 1900 doubles, identity version **2**, binding
both primary and IR transport scattering plus all optical values. Missing,
added or changed optical identity rejects restart. The temporary absorption-
only IR version-1 binary/output were preserved, not silently reinterpreted.
Namelist, mkrun and its common GUI generator were updated together. D03
selects the base contract's scattering flag automatically; absorption-only
cannot accidentally be emitted for this model. The existing fixed path is
unchanged. Model-specific README text no longer claims fixed 30/70 material
or omitted IR scattering when D03 is selected.

### Evidence (native tests, not extra audit infrastructure)

- Existing `dust_mass_smoke` PASS on the final objects: source-segment
  condensation, shared element budgets, growth/erosion, coagulation/
  shattering, SN/fresh protection, WSS09 signs/timesteps and DL01 area/domain.
  Log `.cosmic-ray.kyySgK/dust_d03_mass.log`.
- Existing `snrt_dust_backend_smoke`: 1031 heterogeneous cells, zero dust,
  scalar reference comparisons, local U(T), implicit nonlinear gas exchange,
  energy balance, primary number/flux moment, unequal-weight FP64 IR energy
  and angular relaxation, and invalid-input transaction rollback PASS.
  Reference/backend error <=6.09e-16 on five pure/mixed optical cases;
  heterogeneous comparisons use 1e-10 and IR scattering 2e-14 tolerances.
- CPU log `.cosmic-ray.kyySgK/dust_d03_final_cpu.log`; GPU auto and forced
  logs `.ir-hybrid.dYEXir/dust_d03_final_gpu_{auto,forced}.log` PASS.
  Existing mixed/busy GPU lease test `dust_d03_final_hybrid_cpp.log` PASS
  (free slot actually GPU, held slots CPU, rejected work unpublished).
  Full MPI evolution below uses CPU hydro; this is not multi-GPU hydro proof.
- Offline D03 source/Mie/Rayleigh/spectral/angular test PASS; data include
  hash remains `37c462cec715f0e3009e277249558ab364246c78e373d4049d6e88605f284b54`.
  Generator/manifest status now describes the connected explicit comparison.
- Native binary `.cosmic-ray.kyySgK/ramses_dust_d03_live3d`, SHA256
  `695950be845bc6b34c62416732a1ecbb2447cfc8681dd65058235d7301b0d454`.
  lagRamses VPATH first; SNRT/DUST_LIVE/HDF5, NENER=1/NVAR=30, USE_CUDA=0.
- `.dust-d03-final.pLPhPw/{live,restart}/physical.nml`: noncosmo periodic
  8^3, four steps, MPI2 x OMP2, KL16/LC18/effective SNIa, CR/SF, independent
  BPASS RT, WSS09 CIE, two-size dust with SN exposure and D03/DL01 exchange.
  Both completed. HDF5 final states agree bitwise: 90 hydro datasets,
  4 SNRT datasets, 14 physical particle fields after ID sort (2047 stars),
  and dust/CR identities. Dust mass closure max 4.29e-16; all four bins and
  gas-phase element budgets nonnegative; thermal/CR energies positive;
  source transient fields 28--30 zero. Td=10.2286--10.2674 K.
- Output policy: noutput=1, aout=2, tout=1e30, foutput=2, fbackup=1000000;
  step2/4 dumps, measured 55--67 MiB each, free space 157 TiB. All output
  and checkpoint copies retained. The existing noncosmo SFRD text NaN is
  the previously recorded volume diagnostic, not nonfinite stored state.
- Same binary, same physics/radii/densities, `fixed/physical.nml` changes
  only optics to fixed_mix and output frequency to 4. Final dust/gas total
  energy unchanged at this weak-dust test's precision; Ed decreases 0.639%,
  mean Td 0.219%, IR energy 31.483% with D03. These are model sensitivity,
  not accuracy or galaxy calibration results. Single-run wall times 3.345 s
  D03 vs 3.098 s fixed are indicative only (different dump schedules, shared
  machine); no performance/8-bin equivalence claim.
- `reject/physical.nml`: optical model disabled while reading a D03
  checkpoint. MPI exits 11 with dust restart identity mismatch, before
  evolution and without new output. Dummy aout/tout as above and both
  output periods 1000000; only the copied input checkpoint is present.
- GUI/namelist: 35 tests, 34 PASS / 1 unavailable-display SKIP. Offline
  optical test and `git diff --check` PASS. No commit/push or external audit.
  Previous uncommitted work and unrelated generator deletions preserved.

### Model choice and honest remaining scope

Use the explicit D03 option for size/composition-dependent **comparison**
work; retain fixed_mix for backwards-compatible baselines. Local bases
provide the required native connection without a new general multibin
solver. This closes the live optical/material implementation item, including
MPI restart and OpenMP/CUDA kernel paths; it does not make every previously
listed physical research item complete.

The 8/16-bin accuracy/cost comparison remains **unperformed**, not waived or
passed. Current data bind two graphite radii; an interpolated many-radius
reference would add material assumptions, and cannot by itself validate
coagulation/shattering dynamics. No "best possible dust solver" claim follows
from this smoke case. Existing size-transfer/SN coefficients are uncalibrated.

WSS09 remains low-density gas-phase **CIE** cooling, not LTE or a consistent
radiation-dependent metal NEQ network. Native RT adds photoheating, not a
second H/He cooling sink; this does not fix the CIE/RT ion-population mismatch.
Common grain temperature, frozen-20-K dielectric, no PAH stochastic heating,
H2/CO surface chemistry, separate Fe carrier, sublimation or dust drag/AGN
force remain explicit physical extensions. The noncosmo/no-sink/CPU-hydro
comparison guards are retained. Cosmological galaxy production and
COLIBRE-style parameter calibration are not released by this closeout.

## Continued without approval wait: native size-shift reference (2026-09-09)

Operator clarified that the standing approval remains active. Completed
composition/size/material/optics work was committed and pushed as **93d1719**
to `kjhan0606/LagRamses/main`. Only the unrelated 88-line generator deletion
was left unstaged; binary/output directories were not committed. Original
Draine/WSS09 source whitespace was retained to preserve their byte hashes
(source-code diff checks exclude those raw-data whitespace warnings).

Added a bounded multi-bin radius-shift operator to **existing native
`dust_mass_physics.f90`**, exercised through the existing `dust_mass_smoke`.
This is the planned numerical reference, not another Python framework or a
new live many-bin namelist selector. Fixed/two-size production code is unchanged.

The number distribution is linear in radius within each bin. An imposed
radius shift is integrated analytically over source/destination overlaps;
three-point Gauss integration is exact for the degree-four mass integrand.
Number and mass determine the new linear reconstruction. Its positivity
limiter preserves mass but may change number, following
[McKinnon et al. 2018, section 3.2](https://academic.oup.com/mnras/article/478/3/2851/4995927).
The limiter correction is returned separately, not hidden as physical
destruction. Grains below the explicit minimum size are destroyed, including
return of their residual mass to gas; upper-bound overflow is rejected rather
than silently redistributed. Growth must fit the available gas reservoir.
All bins publish together after validity/budget checks. This assumes frozen
da/dt for the transaction; it is not a new multibin collision prescription.

Measured a common truncated n(a) proportional to a^-3.5, 0.001--1 micron,
with delta-a=-0.0005 micron, normalized initial mass 0.001 in the test units.
Analytic shifted moments and the same minimum-size destruction boundary
provide the independent reference. Three resolutions give:

| Bins | Mass error | Geometric area error | Number error | CPU per call |
| --- | --- | --- | --- | --- |
| 8 | 0.156641% | 3.22084% | 8.80473% | 1.03 microseconds |
| 16 | 0.0102393% | 0.00384499% | 1.86336% | 1.94 microseconds |
| 32 | 0.00113329% | 0.0216182% | 0.301715% | 4.73 microseconds |

The 16-bin area error happens to cancel more strongly than at 32; the
initial test incorrectly required every area error to decrease monotonically.
That assertion was corrected, not the computed values. Both refinements
improve on 8, while mass converges monotonically. Timings are 100 repeated
single-cell calls on the shared CPU, not an end-to-end simulation benchmark.
The 8-bin limiter changes number by -8.55391% of its discrete initial number.
Conservation checks explicitly include minimum-radius destruction and this
numerical correction. A same-total-mass two-size fixed-radius closure gives
0.935742% mass error here; this also changes the physical size representation,
so it is not a formal convergence point for the multi-bin discretization.

Native checks PASS: old dust physics tests, total dust+gas mass to 3e-14
relative, transported number with the boundary/limiter terms to 1e-13,
positive reconstruction, allowed growth, complete erosion, gas exhaustion,
and upper overflow with unchanged trial outputs. Evidence:
`.cosmic-ray.kyySgK/dust_multibin{_build,}.log`. No new RAMSES calculation or
snapshot was needed. Existing tested live D03 binary/output are preserved.

Disposition: this completes the **first bounded size-shift accuracy/cost
comparison**, not the combined coagulation/shattering or optical-SED study.
Do not promote 8/16 bins solely from mass accuracy. The importance of bin-edge
handling and number errors is consistent with the independent discussion in
[Sumpter & Van Loo 2020](https://academic.oup.com/mnras/article/494/2/2147/5813268).
Use 32 as the better-resolved reference in this test; a globally optimal live
bin count is not determined. Continue the physical work under standing
approval; no audit or user-approval wait was inserted. Full NEQ/molecular
chemistry still requires actual reaction/rate data and native integration,
not relabelled CIE tables.
