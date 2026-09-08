# Population-matched SED / physical Z coverage / microscopic SNIa bundle

User request: implement medium-term items 2, 3 and 5 after commit 357bc2f.
Project: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses. Final objective is usable
native RT/feedback/dust simulation physics with a defensible publication
model, not an expanding collection of Python gates. Preserve current model
defaults, approved empirical DTD, effective-SSP comparison, and old restart
identities. Existing ramses_nml_generator.py deletions are unrelated.

Continuation instruction (2026-09-07): report when the current stage is
complete and enter the next in-scope task without an approval pause. This
does not permit declaring unresolved physics complete or silently selecting
a materially different population model.

Current authorization (2026-09-08): the operator additionally preapproved
**all of #2 (additional RT/AGN/dust physics) and #3 (backend/operations
expansion)** and requested continued implementation. This supersedes the
historical #2/#3 approval holds below; no further substep approval pause.
Data dependencies and source/model limitations are still reported honestly.

## Source facts and proposed scope

1. Population-matched radiation: current BPASS HDF5 is already integrated
   for bin-imf135_300 (0.1--300 Msun, slopes -1.3/-2.35). The feedback control
   is Kroupa 0.08--120 Msun, effective binary fraction 0.5. An IMF header edit
   or scalar normalization cannot make these the same population. Look for
   track-level spectra and population weights with a compatible binary model.
   Keep the current independent-radiation path explicitly independent. Do
   not change feedback IMF/binary assumptions without a user model decision.
2. Physical Z coverage: existing KL16+LC18 common hull is .007--.01345.
   Generic native interpolation already exists. The local COLIBRE source
   snapshot also contains yield_z001.txt, attributed to Fishlock et al., not
   KL16, but that file alone supplies no source-matched total lifetime/core
   classification. Verify its provenance and original evolution data before
   making a distinct low-Z adapter. High-Z AGB data cannot extend the common
   hull while LC18 still stops at .01345. No endpoint clamping/extrapolation,
   invented lifetime fit, relabelled Z, or silent cross-source stitching.
3. Microscopic SNIa: current strict_wd and effective_ssp modes are empirical
   DTD accounting choices, not binary evolution. A real alternative requires
   a chosen binary population and source-linked histories of WD formation,
   accretion/merger, exploding mass, donor mass and event rates. A tabulated
   aggregate DTD alone cannot supply those quantities or replace mass debit.
   Examine public BPASS/BPS products for the required common-population data.
   Do not call a reservoir/normalization change microscopic BPS.

## Implementation boundaries

Implement source-backed native consumers only once the necessary numerical
inputs and model conventions are identified. Reuse existing age/Z integration,
particle transaction and HDF5 identities; add no generic approval engine.
Any new namelist choice must update mkrun.py and its shared generator together.
Offline conversion is permitted for physical input preparation, not a Python
runtime replacement. A material choice of a different stellar/binary model
must be returned to the operator with its impact before activation.

Testing, once executable scope exists: extend existing native source tests
for actual-table integration, shared population normalization, WD/event mass
closure and interval splitting. One bounded coupled run/restart is sufficient
for wiring; no generic AMR/restart gate expansion. End review by the operator.

## Planning review request

Fable: read-only review of feasibility, scientific sufficiency, fit to final
goal and overinstrumentation. Identify which subset can honestly be implemented
with existing inputs and which needs a model decision/additional numerical data.
Do not demand extra registries, hashes or synthetic harnesses to substitute
for missing physical data. Do not run jobs or edit files.

Primary source entry points inspected:
- https://bpass.auckland.ac.nz/9.html (BPASS release)
- https://zenodo.org/records/6338460 (integrated BPASS spectra)
- https://zenodo.org/records/7340797 (BPASS starter kit)
- https://arxiv.org/abs/1410.7457 (Fishlock et al. low-Z AGB models)
- https://academic.oup.com/mnras/article/482/1/870/5123725 (BPASS transient rates)

## Operator decision and focused implementation progress

The operator subsequently allowed a separately selected binary-population
comparison option. This permits developing an alternative, not changing the
default or relabelling the current BPASS/Kroupa mixture as one population.
The parallel COLIBRE calibration/interpolation discussion did not authorize
or trigger a new resolution-tuning campaign in this bundle.

Fable planning review is retained in
`fable_population_sed_z_snia_plan_audit_2026-09-07.json`. Its warning that a
scalar binary fraction or aggregate DTD is not a microscopic binary model is
accepted. Its proposed transcription of Fishlock total lifetimes is not:
Table 1 does not contain total stellar lifetimes. Its numbering objection
referred to a different historical list; the present scope remains SED/Z/SNIa.

### Low-Z native connection completed as an explicit comparison

`build_kl16_lc18_native.py --low-z-agb fishlock2014_raiteri96` adds 15
Z=.001, 1--6 Msun Fishlock gross-yield nodes to the existing 58 KL16 nodes.
Source SHA256: `b7379ba6eda1018bfeab527ff9bd17fa57fc4da45b2aa5ff845014c31f7eee00`.
The 7-Msun ONe model is excluded. Gross mass summed over all source elements
sets the return, and its complement sets remnant mass, without renormalization.
Fishlock endpoint-envelope assumptions are retained. All 15 inferred remnants
also agree with independently published Table-1 core masses within their
0.0005-Msun rounding intervals (maximum difference 0.00045562 Msun). The
reader checks this without altering the tabulated ejecta. Low-Z ages use the
published Raiteri96 Padova fit (Valiante 2009 eqs. 3--6), explicitly NOT
Fishlock/Monash lifetimes. Existing KL16 ages and normalization are unchanged.
The common active-channel hull is now .001--.01345 for this option only.
Native linear-Z source mixtures and terminal events are reused. This does not
extend LC18 to high Z or to primordial composition.

The existing native source fixture was extended, not replaced by a Python
runtime. Optimized production Intel objects were reused without changing
Makefile VPATH, native source ABI, or main namelist fields. Consequently no
generator schema change is required; existing generator deletions are untouched.
Tests pass for both default and optional tables: 73 node budgets/event ages,
interval splitting at Z=.004, full 0--13.7-Gyr effective-SSP DTD mass closure
at Z=.001,.004,.007,.01,.012,.01345, and explicit rejection at .00099/.01346.
The strict 40-Myr DTD/WD-supply incompatibility remains detected; no rate was
clipped or delayed. For a 10000-Msun SSP at Z=.004, full return including SNIa
is 2283.84391724636 Msun, remaining mass 7716.15608275364 Msun.

Scratch/evidence root: `/gpfs/kjhan/LRD_JWST/.population-extension.jE6JDf`.
`native-default.log`, `native-low-z.log`; exported optional table/history:
- yields.dat: `b69d478a2a7470a21115a8b2a29ccdf69dc4d018beb5eb94038a77d5fba3db75`
- history.nml: `75e21909bbb948db983f6029ed02753d9d5b59f1c2e3ed08ac0fd620373ce860`
Default outputs remain byte-identical to the prior .agb-physical comparison:
- yields.dat: `de15e4e62b757acaac581b4bf536c733cc866d9db49f88a2d4393991cb7c87d2`
- history.nml: `894339fd9e1b1238b1a658cb5829b20c25cdb94b992c308e405012792091690d`

Bounded native wiring run: `live/run.nml`, initial gas/star Z=.004,
four fixed level-3 steps, MPI1/OMP4, CPU backends, live independent BPASS RT,
AGN, DL01 dust and effective-SSP SNIa. Qualified binary remains
`.ir-hybrid.dYEXir/ramses_ir3d`, SHA256
`d6d1f8f3e5848250572ae1fc068843d08afb49cdfa2acbf80b856d2583285209`.
Run exits 0 in 22.41 s: four RT and four IR commits; peak IR balance
5.0293e-10; all 161 floating HDF5 datasets finite; positive leaf gas/dust mass
and dust energy; no MG nonconvergence. This short run is not evidence of
late AGB events or galaxy-scale calibration; the native source test covers
long-time AGB return. The independent BPASS population is still labelled as
such, and the large numerical BH seed remains a stress-control assumption.

`restart/run.nml` resumes a preserved COPY of the step-2 dump to step 4.
94 hydro/SNRT arrays compared: 91 identical, only three leaf momenta differ,
maximum array-normalized difference 1.4068416179391722e-18. Both source-model
and SED identities are identical. The first ad-hoc inspection assumed
`snrt/level_3` was a group; it is a dataset. That diagnostic AttributeError
occurred after the successful parity calculation; the corrected inspection
passes. No simulation code or tolerance was changed for it.

Launch policy was reported before each run: noutput=1, aout=2, tout=1e30
(unreached), foutput=2, fbackup=1000000, nstepmax=4. Live produces two
55,090,544-byte output directories, restart uses a copy of one and adds one
new dump, about 221 MB total. Free storage before launch: 166 TiB. Existing
outputs remain preserved; no long or production calculation was launched.

Remaining: source-linked common-population SED/SNIa implementation, and
additional physical Z coverage beyond this explicitly mixed low-Z option.
No claim that the entire three-item medium-term bundle is complete.

### BPASS track and population recipe reconnaissance

Actual numerical source members, not only aggregate SN rates, were inspected
via verified HTTP byte ranges from [Zenodo 7340797](https://zenodo.org/records/7340797).
The 8,151,863,444-byte model ZIP has 256098 entries, ZIP64 central-directory
offset 8111316989 and size 40546357. Only metadata and selected members were
read; the full 8-GB archive was not installed. ZIP member names, decompressed
lengths and CRCs were verified. The following are **two independent examples**,
not a claimed consecutive genealogy:

1. `bpass-v2.2-newmodels/NEWBINMODS/NEWBINMODS/z010/sneplot-z010-5-0.5-1`
   - local-header offset 3616814451; compressed size 37378; raw size 195199;
     CRC32 1840578054; SHA256
     `86e266db9945545a93755ee2d9dbe9a9d68788ab89bd40959418d18e46685c42`.
   - 127x96 numeric values; age 0--125178000 yr; primary 5 -> .8313 Msun,
     companion 2.5 -> 2.670092 Msun; final primary CO core .81037 Msun.
2. `bpass-v2.2-newmodels/NEWBINMODS/NEWSECMODS/z010_2/sneplot_2-z010-3-1.00000-2.40000`
   - local-header offset 6880198090; compressed size 41262; raw size 298178;
     CRC32 2271612651; SHA256
     `165c4ee518dd042c782b891badf3a52524fab439773e9975359e5f4f22f89289`.
   - 194x96 values; own track age 0--97442970000 yr; evolved star 3 ->
     .69724 Msun, remnant companion 1 -> 1.417531 Msun. This is NOT by itself
     an SNIa event/WD debit: the stored companion mass persists and the
     population clock, event classifier and post-event states need resolution.

In `bpass_v2.2.1_imf135_300.zip` (1537193948 bytes), the Z=.010 binary recipe
`bpass_v2.2.1_imf135_300/input_bpass_z010_bin_imf135_300` is at local-header
offset 16034957, compressed size 251667, raw size 2008517, CRC32 4229237729,
SHA256 `4296683815fc42cbec9d0c0d366b5bcb93210c3bc693afbed30f6c625399af29`.
It has header normalization 1, and 21298 entries with type counts
{-1:241, 0:2150, 1:10513, 2:8332, 3:62}. The examples above occur with recipe
weights 4.8647494/type1 and 5.5797787/type2 respectively. Secondary recipe
mixed-weight/age fields for the selected example are both zero. These weights
must NOT all be interpreted as independent primordial systems and normalized
by summing initial masses: the recipe includes later/rejuvenated components.

Source interpretation references:
- [BPASS manual, stellar models and synthesis inputs](https://warwick.ac.uk/fac/sci/physics/research/astro/research/catalogues/bpass/v2p2/bpassv2.2_manual-arial.pdf)
- [Provider Hoki compilation recipe](https://heloises.github.io/hoki/ModelDataCompiler.html)

Next implementation target: resolve the synthesis recipe's component/clock
semantics and event ownership before exporting native SED/feedback histories.
Keep population normalization, stellar mass loss, surviving companions and
WD explosion debit coherent. Do not add another empirical DTD or activate an
unconnected track reader merely to claim a microscopic population is done.

### Post-push component-clock findings (2026-09-07)

Low-Z implementation committed and pushed to `kjhan0606/LagRamses`, main,
as `d1505550aaf2b61d3f53eff3f3a4d01f8ac05917`. The pre-existing generator
deletions and scratch outputs were excluded.

The full Z=.010 recipe identified above was scanned, not just the two sample
tracks. Of 8332 type-2 components, 5095 have zero rejuvenation age and 3237
have finite positive ages. Of 62 type-3 components, 24 have zero age and
**38 have Infinity with nonzero rejuvenated weight**. These are secondary
components referencing single-star tracks, including masses in the WD
progenitor range. They cannot be silently discarded or assigned age zero
in a physical WD/SNIa population. This is a finding at Z=.010, not a claim
that all metallicities have been scanned.

Provider [Hoki CMD code](https://github.com/HeloiseS/hoki/blob/1a2dce0d5907b3c14a176b0aae2cfd74d6d9752a/hoki/cmd.py#L143)
confirms that total component weight includes the rejuvenated part. For a
finite shift, contributions are unshifted `(total - mixed)` and shifted
`mixed`, not two full-weight independent systems. It rounds weights and
maps Infinity to zero for CMD construction. Its primary-WD omission is
also a light-counting convention, not physical remnant destruction.
Neither convention is adopted for feedback/event accounting.

Printed weight precision also differs. One actual type-2 record has total
8.6395101547241211 and mixed 8.6395102 (shift 300174580 yr); naive subtraction
produces a small negative weight. This is distinct from the nonfinite
clock problem; any future reconciliation must preserve total weight and
be justified by source precision, not a general clipping policy.

The [BPASS v2.2 paper](https://academic.oup.com/mnras/article/479/1/75/5003394)
ties rejuvenation to mass-transfer onset and the initial/final secondary
hydrogen-burning lifetimes. The primary terminal age is not a general
replacement. Next required physical connection is therefore the actual
primary-to-secondary mapping and those timing inputs, followed by an
explicit explosion/remnant debit. The public aggregate recipe alone has
not yet established that connection. No new runtime mode, empirical DTD,
standalone parser framework, or simulation launch was added for this
reconnaissance. Existing independent BPASS SED and effective-SSP SNIa
remain unchanged; common-population SED and microscopic SNIa remain open.

### Genealogy and post-explosion feasibility follow-up (2026-09-07)

The BPASS v2.2.1 manual, pp. 19 and 31, was read from the starter-kit PDF.
Primary filenames contain initial M1, mass ratio and period; secondary
filenames contain post-primary M2, remnant mass and period. The input recipe
contains marginal model weights and rejuvenation fields, not parent IDs or
conditional transition weights. The [BPASS stochastic-population paper](https://academic.oup.com/mnras/article-abstract/522/3/4430/7140546)
explicitly describes many-to-many primary/secondary mapping after sampled
supernova kicks. Thus matching masses/filenames cannot uniquely invert this
recipe. This establishes a limitation of the inspected products, not proof
that no additional provider data exists.

Next, the single-degenerate event interpretation was checked against
[Eldridge, Stanway & Tang 2019, section 2.1.3](https://academic.oup.com/mnras/article/482/1/870/5123725).
Their model selects a WD initially below 1.2 Msun which accretes to 1.4 Msun;
this is a model prescription, not a universal explosion criterion. The
previously identified secondary track (SHA256 `165c4ee5...89289`, full hash
above) was fetched again by byte range and verified against its full hash.
Its first 1.4-Msun companion crossing is bracketed by these **1-based rows**:

| Row | Track age (yr) | Donor mass (Msun) | Compact companion (Msun) |
| --- | ---: | ---: | ---: |
| 79 | 431262800 | 2.769640 | 1.224356 |
| 80 | 431273700 | 2.364330 | 1.415459 |
| 81 | 431278000 | 1.999380 | 1.415459 |
| 194 | 97442970000 | 0.697240 | 1.417531 |

Linear interpolation gives a diagnostic crossing at 431272818.26 track yr
with donor mass 2.39711696 Msun. It is NOT an exported SSP event age or an
assertion of a resolved physical explosion time. The age column is monotone;
114 rows follow the first sampled crossing. Keeping those rows while also
ejecting the WD would retain an already-exploded object and continue donor
evolution in its presence. Truncating the whole system would instead lose
the surviving donor. A threshold detector alone therefore cannot complete
conserved population feedback or its SED.

Implementation remains dependent on parent/branch timing and weight data,
and source-consistent post-explosion donor evolution and ejecta. No WD debit,
runtime option, or source-table normalization was altered to bypass this.
The existing empirical SNIa and independent SED remain usable comparisons.
As a next-source feasibility check, public COSMIC documentation identifies
linked initial/bpp/bcm histories through bin_num, unlike the inspected BPASS
marginal recipe. That is only a candidate: it uses a different evolution
model and is not population-matched to BPASS spectra by sharing an IMF.
No alternative model was installed or activated.

The local `pdftotext` utility was absent; PDF inspection was completed with
the already-installed pypdf reader. No new diagnostic framework or job was
created. This follow-up narrows missing physical inputs; it does not mark
the common-population implementation or this bundle complete.

### Implemented separate birth-to-remnant source path (2026-09-07)

After the operator requested resolution, an actual alternative was implemented
within the previously permitted **separate binary comparison** scope. It is
not a repair of BPASS v2.2.1, and the BPASS missing genealogy is not relabelled
as recovered. No default, native namelist, or existing yield input was changed.

`simulation/snrt/tools/build_cosmic_binary_histories.py` invokes COSMIC's
compiled Fortran BSE engine for explicit ZAMS systems and retains initial
conditions, all phase records, final states and kicks with common system IDs.
It requires COSMIC 4.2.0 and the pinned release examples/Params.ini (full hash
in the script). Source commit: `f9b90f451bca014e9e0adb3a410bee8752e30e53`.
The CLI is bounded to 1--256 low/intermediate-mass initial binaries, a fixed
explicit seed and existing SSE metallicity range .0001--.03. It does not
sample an IMF, infer population weights or extrapolate a BPASS age.

Actual execution root: `/gpfs/kjhan/LRD_JWST/.binary-source.AIstqO`.
COSMIC was installed in its isolated `venv` (608 MB); the four-system physical
output `reference/` is 99 KB. No existing Python environment was modified and
no lagRamses job was launched. The input fixture is four 5+3-Msun ZAMS binaries
at Z=.01, eccentricity zero, periods 1,10,100,1000 days, followed to 13700 Myr.
These are unweighted examples, **not** a population-normalized SN rate.

The 10-day system forms a CO WD and then records its disappearance at
156.51229100340674 Myr. Its mass remains zero for the rest of the calculation.
The donor has 1.267163010011896 Msun at that record, continues evolving and
ends as a 0.6255253387250216-Msun CO WD. The 60 phase records across all four
systems preserve simultaneous transitions. Merger-created massless components
in the control systems are NOT counted as CO-WD no-remnant events. This is
not an independent reproduction of the BPASS 3+1-Msun secondary sample.

The focused real-engine test `simulation/snrt/tests/cosmic_binary_histories.py`
passes birth-clock/ID checks, persistent WD removal, donor survival, merger
exclusion and rejection of nonfinite/missing/resurrected histories. Two fresh
executions reproduce bpp/bcm/initC/kick CSV bytes exactly with seed 20260907.
No separate mock evolution engine or generic approval infrastructure was added.

Important physics boundary: the pinned BSE source `evolv2.f`, lines 3200--3246,
includes a CO-WD helium-accretion destruction prescription after 0.15 Msun
accumulation as well as Chandrasekhar destruction. These are not interchangeable
with N100. A BPP previous phase row is not the instantaneous pre-event state;
the adapter therefore leaves event mass and elemental yields **unset**, instead
of using a phase-to-phase mass difference as an explosion yield. The output
explicitly has `runtime_ready=false`. This source path resolves shared clocks,
component ownership and post-removal survival for the comparison examples;
it does not validate thermonuclear physics or finish common-population coupling.

Next implementation dependency: obtain channel-resolved instantaneous explosion
mass/ejecta consistent with the selected evolution prescription, then apply
explicit population weights and a same-evolution atmosphere/SED calculation.
Retain the current BPASS SED and empirical SNIa defaults until those choices
are physically consistent. Do not attach BPASS spectra or N100 through a
normalization-only adapter. No claim of publication-ready binary physics.

### Scope reset and ordinary CCSN connection (operator instruction, 2026-09-07)

The operator selected unfinished-work **section #1 (physical inputs and their
native coupling)** for continued implementation. Sections **#2 (additional
RT/AGN/dust physics)** and **#3 (backend/operations expansion)** require further
approval. Do not turn verification work into new implementation gates.

Priority correction: the previous KL16+LC18 comparison enabled SNII but supplied
only 40/60/80/120-Msun wind-only endpoints, with zero terminal ejecta/energy.
It did not include ordinary CCSN. This dependency now takes priority over more
binary-evolution feasibility scripts. Existing empirical SNIa is retained;
microscopic binary population work remains separate and incomplete.

Implemented explicit `lc18_set_r` alternative in the existing LC18 and combined
KL16/LC18 offline builders, consumed by the **production Fortran modules**:

- LC18 nodes 13,15,20,25,30,40,60,80,120 at all four source metallicities.
  The first four nodes have terminal ejecta = table8 minus table9; higher
  nodes follow the Set R wind-only assumption. All 15,984 isotope differences
  across 48 low-mass/rotation/Z models are nonnegative; no clipping was needed.
- Wind and terminal isotope sums define returned mass. Remaining initial mass
  is the baryonic residual, not a gravitational remnant mass or rounded iron
  core mass. No new isotope decay/net-yield inference is introduced.
- History v2 explicitly spans 13--120; v1 retains 40--120. Native terminal
  events have no pre-lifetime leakage. The >=40 endpoint preset cannot
  suppress ordinary SN or mix across the 40-Msun seam. Nearest source-cell
  mass-fraction scaling is unchanged; 25/30 source-cell boundary is 27.5,
  **not a physical individual-star explodability threshold**.
- SN energy is a required explicit builder argument, never inferred from the
  pre-SN binding energy. The local comparison selects 1e51 erg at each exploding
  source node, effective massive wind 1000 km/s, AGB wind 15 km/s, rotation 0,
  Kroupa/effective-SSP SNIa and the optional Fishlock low-Z AGB source.
- `mkrun.py` offers `comparison_ccsn`; old comparisons/defaults remain unchanged.
  The shared namelist generator documents v1/v2 history support. Its pre-existing
  unrelated 88-line deletion was preserved, not attributed to this work.

Evidence: `/gpfs/kjhan/LRD_JWST/.ccsn-source.lobKc9/`. Actual input table SHA256
`48ecb3a36450a93834e4454d30c140330a5d41e20c962823156e9b27384e98e4`,
history `92a9fb7e9ff93ee4bbf16b702ca9ef5103f184c7f600c6642c2dc13dfb5ca263`.
Existing builder default output was byte-identical to HEAD; all three LC18
rotations generate complete 36-node sources. Extended existing
`kl16_lc18_native_test.f90` passes GNU bounds/FPE checks and Intel production
objects: 16 exploding endpoints, source-hull rejection, timestep/IMF-bin
independence, AGB/effective-SNIa closure at six Z values through 13.7 Gyr.
At initial SSP mass 10000 Msun, Z=.004, by .2 Gyr the ordinary SN channel
returns 527.075650954212 Msun and 3.77999123914469e52 erg. Old default and
low-Z-only fixtures and high-mass endpoint/history regressions pass.

RAMSES `ramses_ccsn3d` was rebuilt with the existing Makefile, lagRamses-first
VPATH, HDF5/SNRT/DUST_LIVE/CUDA linkage and NVAR=30. SHA256:
`c76dd129a091e16b858e9c18ab930d19455d343b86f7cb0e1689c5b3fa3c996b`.
The old object cache and binaries were preserved. Fresh `live/ccsn.nml` uses
1 rank/2 OpenMP threads, 512M worker stacks, CPU hydro/RT/dust, level 3, four
steps; `restart/ccsn.nml` restarts from step 2 and stops at 4. Both finish;
the short live run is an initialization/wind-coupling check and does **not**
claim terminal SN/AGB ages were reached. See native endpoint/SSP tests for
actual terminal release, and `build.log`, `native-{intel,gnu}.log`, live/restart
logs for evidence. No production-scale or external-auditor jobs were launched.
Final HDF5 comparison: source identity and every SNRT dataset are exact;
hydro maximum absolute difference is 4.14e-25 in momentum. The 512 stars,
matched by positive unique ID, have exact masses, birth clocks, metallicities,
progress clocks and positions; velocity differences are at most 2.12e-22.
Sink-cloud particles share negative IDs, so sorting the entire particle array
by ID alone is not a valid correspondence (the initially reported apparent
position mismatch was from those clouds, not stellar-source state loss).

Remaining within #1: first fill the genuinely unsupported **8--13 Msun**
CCSN/ECSN transition with an explicit compatible source/fate choice; do not
extrapolate LC18's 13-Msun node down to 8. The same-population SED problem,
microscopic SNIa, AGB edge masses, wider metallicity support, phase-resolved
winds, net/isotope treatment and alternate high-mass/PISN source packages
remain on the previously stated #1 list. This closes the LC18-supported
ordinary-source connection, **not all of section #1 or publication readiness**.

### Full #1 preapproval and physical extensions (2026-09-08)

Operator: “#1 전체를 지속 진행. 세부항목은 모두 사전승인.” All subitems
within #1 are approved for implementation; do not pause at tiny substeps.
#2 and #3 are still approval-held. Scientific data dependencies do not become
validated physics merely because implementation is preapproved. No external
audit, production-scale calculation, commit or push was performed in this
implementation turn. Work remains in /gpfs/kjhan/LRD_JWST, origin
kjhan0606/LagRamses, main. Shared global context was not changed.

Implemented in the production patch/lagRamses modules, not just Python:

1. **Solar 9--13 CCSN comparison**: 17 Sukhbold nodes, Z9.6 through 12,
   W18 above 12. Stable-segment wind/terminal gross ejecta and actual PHOTB
   energies feed native endpoints/SSP. Source-specific history v3 bounds
   [9,13] reject out-of-range mass/Z. Residual remnant, solar coordinate,
   Raiteri96 ages and effective wind speed are named comparison choices.
   Maximum remnant difference from PHOTB .09434268 Msun (12.5-Msun node).
   8--9, non-solar extension and complete isotope decay remain unsupported.
2. **LC18 phase timing**: Table5 cumulative-loss shape sets wind timing;
   integrated composition/energy endpoints stay unchanged. Explicit uniform
   fallback at six rotation-zero nodes with zero printed phase loss. Sparse
   per-star age knots reduce the AGB7 input from about 3.6 MB to 176 KB;
   they are the same piecewise model, not a timing approximation. Native
   validation no longer treats physically distinct hour-separated phases as
   duplicates under the 0.1-year budget tolerance. At M=20/Z=.01345, H-end
   wind at 9,158,200 yr is .499961981828186 Msun; final 12.4590525871584.
3. **AGB envelopes through 7 Msun**: 61 KL16 nodes. Z=.007/M=7 CO(Ne)
   contributes its envelope normally but cannot supply a strict CO-WD Ia.
   CO classification is applied before Z mixing, and participates in restart
   identity. Other selected nodes retain normal remnant accounting. This
   does not manufacture a 7--8-Msun common-Z grid or microscopic Ia events.
4. **Optional signed AGB net**: match source initial M/Z/Y, require identical
   vectors for duplicate initial entries, normalize initial fractions, then
   subtract initial composition times returned mass from normalized gross.
   Initial raw fractions do not sum exactly to one: e.g. Z=.007/M=1 sum
   1.006061391734592. The normalization is explicitly named and recorded,
   not presented as exact raw-source abundance recovery. Native channel
   mask F,T,F,F,F distinguishes supplied AGB net from absent table channels;
   unavailable channels cannot carry nonzero net columns. Gross gas return
   is unchanged. The channel mask participates in restart identity.

All are explicit alternatives. Existing builder defaults and the preceding
ordinary-CCSN input remain byte-identical. mkrun terminal/GUI selections and
shared namelist-generator history help were updated together. No Makefile
or VPATH order changes. The unrelated pre-existing generator deletion was
preserved. No source-derived large table or binary was added to Git.

Evidence root: `.physical-extension.7rcxv4/` (local, retained).

| Input directory | yields.dat SHA256 | history.nml SHA256 |
| --- | --- | --- |
| sukhbold | f69d1690d76f3bcd60936f43049924afcf807bce5ae2bc875412e8cf24959a78 | 5be5e46b0a19799e531cdc89a2f2ff334024c971fdd79362ca27afe30b02f0f7 |
| lc18-phase-sparse | b71c1ac741924aa98a2b3ff63d119a247329d8e85669cae28981171bc3d16f6c | d2a62e64bba38c66eeba9d3957debba54675b9b3812fc1b4165ef890ef7d91f0 |
| agb7-sparse | 2d61de4b3381ed72ceed472f9d45c08173d06376068efc5616d959c2a3ad2e49 | 63dff962d20081279f86108d1100b1deb0d3723560fa25c5195306280c0ff0ed |
| agb7-net | 3c3af500c33162cea184e89f072e7841e612a86042c3a140254915961fab817f | 5afb3011d8284836c08b9b2e843fc9d7fc71b59443093dc6023202642c8589e2 |

GNU bounds/FPE and Intel production-object fixtures pass endpoint budgets,
age splitting, Z-hull rejection, hybrid exclusion and signed net. The
9--13 comparison at SSP initial 10000 Msun/Z=.02/.1 Gyr returns
357.497985956766 Msun (SN 335.862657285151), SN energy
2.11024969527311e52 erg. The AGB7 effective-SSP model at Z=.01/13.7 Gyr
returns 3343.3696147592855 Msun. These are IMF-integrated native consumer
results, not a claim that short live runs reach all terminal events.

RAMSES builds used existing Intel -O3 objects plus rebuilt changed modules,
HDF5=1 USE_CUDA=1 SNRT=1 DUST_LIVE=1 USE_FFTW=0, NVAR=30, lagRamses-first
VPATH. Latest `ramses_physical_net3d` SHA256:
`a5a8f269525a0e05f30028dca632529315adbf0cd4e1e681a9c77e582fb79121`.
GNU/Intel fixtures cover the real terminal release; full RAMSES tests use
1 rank, 2 threads, OMP/KMP stack512M, CPU hydro/RT/dust, level3, nstepmax4,
noutput1, aout2/tout1e30 unreachable, foutput2/fbackup1000000. New live run
`live-net/physical.nml` completed in 17.54s, step2-to4 restart
`restart-net/physical.nml` in 4.02s. Output size about53MB each, free space
165TB at launch. No MG nonconvergence or runtime error in positive runs.
`restart-bad-net-mask` changed only net availability (not source model ID):
exit1 in init_part before any new output_00002. Earlier `live-sparse` and
`restart-sparse` pass; `restart-bad-core` likewise rejects a changed CO(Ne)
classification before any new output. Existing outputs were preserved.
GUI/model suite: 29 tests, one display-dependent skip, no failures.
Final live-net/restart-net HDF5 comparison: native source identity and all
SNRT datasets exact; hydro maximum absolute difference 6.204e-25 in
momentum. Matching 512 stars by positive unique ID gives exact masses,
birth/progress clocks, metallicities and positions, velocity differences
at most4.236e-22. Negative duplicate-ID sink clouds are not a stellar-ID
correspondence. Previous ordinary CCSN input was regenerated in memory and
again matched both table/history checksums above without changes.

#### Remaining #1: physical dependencies, not additional test gates

| Original #1 item | Current extent and actual remaining dependency |
| --- | --- |
| Ordinary CCSN 8--40 | LC18 13--40 coupled; solar-only Sukhbold 9--13 comparison coupled. 8--9 ECSN fate and non-solar 9--13 still need compatible evolution/yields. |
| Same-population SED | BPASS independent radiation remains; require spectra/atmospheres and feedback histories from the same IMF/evolution population, not a normalization-only join. |
| Microscopic SNIa | Effective SSP retained. COSMIC genealogy exists for comparison, but transfer/retention/outflow/orbit re-evolution, population weights and channel-matched ejecta are not complete. |
| BPASS clock/genealogy | 38 nonfinite source clocks and missing binary genealogy are not repaired by substituting an unrelated lifetime fit. |
| Wider Z | Optional Fishlock extends CO<=6 to .001; AGB7/net common hull .007--.01345. Outside each source hull requires new physical data. |
| AGB edge masses | 1--7 envelopes now coupled; 0.8--1 and common-Z 7--8/superAGB fate remain unsupported. |
| Wind phases | Massive loss timing coupled; phase composition/terminal speeds and detailed AGB time profiles remain approximations. |
| Net/isotope treatment | AGB normalized-initial net connected. Other channels' initial composition/complete radioactive inventories and decay projection remain absent. |
| Alternate high-mass source | Solar Sukhbold tail reviewed; cannot blindly join PHOTB post-fallback remnants to KEPLER integrated yields. See findings below. |
| PPISN/PISN | No selected core-mass/metallicity/rotation-resolved yield and population source. Switch remains disabled; no universal ZAMS threshold invented. |

Read-only full Sukhbold check: stable yield segments vary in length (e.g.
14-Msun table), so the historical 283-row review contract is not a full-grid
parser. More importantly, W18 27.3 gives stable wind+terminal minus
(ZAMS-PHOTB remnant) = +4.34056860167437 Msun, PHOTB fallback4.3871;
N20 29.6 gives +5.464652019701305 Msun, fallback5.4944. These differences
are comparable to fallback, not tiny rounding. This suggests a mismatch of
ejecta/mass-cut stages, but is **not** sufficient evidence to subtract uniform
element fractions or change the source tables. Tail120 also differs by about
.37--.39 Msun. No full-grid promotion or increased remnant tolerance was made.
Consult the [source archive](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/index.html)
and [paper, sections4/5 and tables7/8](https://arxiv.org/pdf/1510.04643).
The archive describes stable vs selected-radioactive segments and PHOTB versus
KEPLER postprocessing. Source-derived data remain local; redistribution
permission is not inferred from research access.

Continuation order within the already approved #1 scope: resolve source/fate
compatibility for the remaining mass/Z cells and fallback/decay stages;
connect the corresponding native inputs; then complete shared-evolution
SED and channel-matched binary population coupling. Reuse existing fixtures
and native integration checks. Do not grow verification gates to replace
missing physics, and do not relabel this partial progress as all-#1 complete.

### Continued #1: low-Z AGB7 envelopes and net yields (2026-09-08)

Corrected an overly restrictive earlier source assessment: the pinned
Fishlock element table `yield_z001.txt` already contains each model's X0(i).
No separate external initial-composition file is needed for this source.
The existing reader now optionally computes signed net from the same-model
normalized X0 and gross return. All 78 initial element fractions sum to
.9999998846994237 before normalization. Fishlock gross return is unchanged,
not normalized using the KL16 prescription. This extends the selected
net-yield path to Z=.001, including intermediate-Z cumulative mixtures.

The same source contains a genuine 7-Msun ONe model. It is now included when
both `kl16_envelopes_to7` and `fishlock2014_raiteri96` are explicitly selected.
Gross return5.855482487283754 Msun, residual remnant1.144517512716246 Msun;
the latter is inside the .001-Msun printed rounding interval of Table1's
1.145-Msun core. No source-value correction was needed. Its 53.248752-Myr
Raiteri96 lifetime remains a declared approximation, not a source-matched
Monash age. [Fishlock2014, Table1 and section3](https://arxiv.org/pdf/1410.7457)
identify the evolution/envelope assumptions; no extra thermal pulses or ECSN
event were invented.

New selected grid: 61 KL16 + 16 Fishlock =77 AGB nodes, envelope1--7 Msun,
common active Z=.001--.01345. The history names two non-CO nodes:
ONe at (7,.001), CO(Ne) at (7,.007). The existing production Fortran
classification map excludes both from strict CO-WD SNIa inventory while
retaining their envelope and total-remnant contributions. This filtering
also works before interpolation: at 7 Msun/Z=.004 neither mixed component
can create CO-WD supply. Effective SSP SNIa and the original DTD are unchanged.

Generated input: `.physical-extension.7rcxv4/agb7-lowz-net/`.
Table SHA256 `3b693dcfc12fafc735b8bd3f03a15bf55c2d2edcda2acd4c0b753d7d09dfc1b8`;
history SHA256 `2554cc6f6fdb5c7c41e3bf3bab0985b7ca6e0a6a5a85a42217ce558e46b70a09`.
The prior `ramses_physical_net3d` already implements this multi-node map and
channel mask; no new simulation module or Makefile change was required.
Native test links its actual Intel production objects; GNU bounds/FPE test
also passes. Reused `kl16_lc18_native_test.f90`, option
`lc18_set_r_agb7_lowz_net`, rather than creating another test framework.
Checks include exact endpoints, non-CO exclusion at exact/mixed Z, signed
H/He net, rejection outside mass/Z coverage, timestep telescoping and
effective-SSP closure at six Z values through13.7Gyr. For SSP initial10000
Msun/Z=.004/.2Gyr, native AGB return662.150277439322 Msun, total return
1580.23727066374 Msun before effective SNIa. Fishlock7 signed netH/He:
-.617030706434305 and +.606194676556949 Msun. The strict40-Myr DTD still
fails causal CO supply, as expected; it was not modified to force a pass.

mkrun terminal/GUI adds explicit `agb7_lowz_net`, changes channel2 upper
mass to7 and binds these inputs to the tested native binary. No new RAMSES
namelist keyword. GUI/model suite29 tests passes (one display skip). Old
ordinary CCSN and `agb7-net` inputs were regenerated and are byte-identical.
No full RAMSES calculation, external audit, commit or push in this follow-up.

High-mass fallback follow-up: Sukhbold section3 explicitly separates PHOTB
remnant/energy from KEPLER isotope postprocessing. The available full-zonal
samples are only12.0,15.2,20.1,25.2 Msun, not the large-fallback examples.
Griffith2021/VICE W18F is a **forced explosion alternative**, not a documented
repair of those fallback yields. Neither was used to normalize away the
4--5-Msun discrepancy. See [Sukhbold2016](https://arxiv.org/pdf/1510.04643)
and [VICE W18F](https://vice-astro.readthedocs.io/en/v1.2.x/users_guide/pkgcontents/vice.yields.ccsne.S16.W18F.html).
Doherty's staged super-AGB tables were also checked: missing Ca and matched
total ages remain actual inputs to resolve, not reasons to pad missing
elements with physical zeros. #1 remains partial; #2/#3 remain held.

### Continued #1: source-timed AGB winds and terminal WD formation (2026-09-08)

The earlier terminal-envelope option is preserved. New explicit selection
`--agb-release fishlock_tp_mass_loss` connects actual Fishlock thermal-pulse
mass histories to native RAMSES. This is a physics/timing change, not another
Python verification framework. Source:
[CDS J/ApJ/797/44 ReadMe](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/797/44/ReadMe)
and [table2.dat](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/797/44/table2.dat).
The CDS table number differs from the article's printed table numbering.
Local files are under `external/g2_candidates/fishlock2014_pulses/`:

- ReadMe SHA256: `6fea0aa376e1677e92da47edb21f8a34a2f3e7399049134c98e9b4ef721ad57c`.
- table2.dat SHA256: `4711e3662df69c490a697fba5725bdbba1a10271076b1b6383041f0077fa9d93`.

772 rows, 17 columns. Exactly repeated (3.5 Msun, TP25) and (4 Msun, TP35)
rows are identical in every field and collapsed; 770 unique pulses remain
across 16 models at Z=.001. The initial mass minus tabulated TP total mass
sets cumulative wind; preceding interpulse periods set relative timing.
No raw source files are rewritten or presumed redistributable.

Explicit approximations: last-TP left limit at the already selected
Raiteri96 Padova terminal age; uniform pre-TP loss; unresolved loss after the
last reported TP included in the terminal jump, not invented later pulses.
Integrated mean composition and fixed 15 km/s remain. The source is not a
Monash total-lifetime table. KL16 nodes retain terminal-envelope timing.
At 7 Msun there are 135 pulses over75687.63yr; first/last TP total mass is
6.9846/2.1117 Msun. Native return at the first TP is .0154 Msun, at the last
left limit4.8883 Msun, with **zero remnant at both times**. At the terminal
event gross return becomes5.855482487283754 and ONe remnant1.144517512716246.
Its exclusion from strict SNIa CO inventory is unchanged.

Native table admission/interpolation now accepts `wind_history_terminal_remnant`
and source-specific `agb_wind_jump_*` fields. Ordinary rows interpolate to
the left limit; the terminal row contains the right limit and creates the WD.
Admission rejects a wind row exceeding that left limit in total mass,
elements, untracked mass or energy. Restart identity binds the consumed jump
fractions as well as the full cumulative table. Old all-terminal identities
remain unchanged. No new main RAMSES namelist keyword or VPATH change.
mkrun CLI/GUI adds `agb7_pulses`; shared generator history help is updated.

Generated `.physical-extension.7rcxv4/agb7-pulses/`:

- yields.dat SHA256: `4a42f7d2fd089b0f8fabc3127f3ef27019fcb9eda418f2342b453f1d2b6c3c50`.
- history.nml SHA256: `67271d116d1c791ffd92c787dbe5f09bf16d2842e0d7e81a71ba0547f6831865`.
- New binary `ramses_physical_pulses3d` SHA256:
  `fa879828536954987597e2267a29a8210ef87a51086ca42876c6b2dc2598a897`.
- Build: `make -C .physical-extension.7rcxv4 -f ../bin/Makefile -j4 HDF5=1
  USE_CUDA=1 SNRT=1 DUST_LIVE=1 USE_FFTW=0 EXEC=ramses_physical_pulses`.
  Build log `build-pulses.log`; earlier binaries/outputs preserved.

Reused `kl16_lc18_native_test.f90`, mode `lc18_set_r_agb7_pulses`:
GNU bounds/FPE and Intel production-object tests pass for all77 endpoints,
terminal left limits/no early WD, signed net, non-CO filtering, Z mixture,
SSP interval telescoping, ordinary CCSN and full0--13.7Gyr effective-SSP mass
closure at six metallicities. Existing40-Myr DTD still fails strict CO-WD
causality; no DTD change. Old `agb7-net` and `agb7-lowz-net` generation remains
byte-identical. GUI/model29 tests pass, one display-only skip.

Bounded live check: `live-pulses/physical.nml`, MPI1/OMP2, worker stacks512M,
level3,4steps, Z=.004, H=.746/He=.25; initial metal/Fe/dust/thermal dust energy
adjusted together to retain dust-to-metal ratio and20K initial temperature.
RT, real-source feedback, independent BPASS, AGN and live DL01 dust remain
the selected comparison, not galaxy calibration. noutput1/aout2/tout1e30
(unreached), foutput2/fbackup1000000; about2x56MB,165TB free before launch.
Completes in20.386s. `restart-pulses/physical.nml` resumes a COPY of output1
(step2) to step4 in4.197s. Both logs complete without MG nonconvergence.
All161 floating datasets are finite. Hydro/SNRT and48054 source-identity
values are exact; matching512 positive unique star IDs gives exact masses,
birth/progress clocks, metallicities, positions and velocities. Negative
duplicate sink-cloud IDs are not used as unique stellar correspondence.
Short live evolution tests wind wiring; native long-age tests cover TP/WD
events. This closes this AGB timing subset, not all of section #1.

Negative restart: `restart-bad-pulse-jump/physical.history.nml` decreases
only the first terminal-jump fraction by .002, leaving all table values,
source labels and hashes unchanged. `run-short-path.log` exits1 via
MPI_Abort at stellar source identity comparison in init_part, before a new
output2; the consumed timing change cannot hide behind an unchanged label.
The first invocation's `run.log` did not test this: the legacy80-character
command-line filename buffer truncated its long absolute path and exited0.
It was retried with `physical.nml` relative to the same audited workdir.
The existing files are preserved. No generic path/error-exit rewrite was
added to this physics task.

Authorization update: the subsequent user instruction preapproves #2 and #3
as well. Earlier holds are historical. Continue source/native work plus
additional radiation/dust physics and backend expansion without substep
approval waits; do not revive unrelated generic HDF5/AMR infrastructure gates.

### #3 execution portability: toolkit-free native CPU build (2026-09-08)

Implemented `SNRT=1 DUST_LIVE=1 HDF5=1 USE_CUDA=0 USE_FFTW=0` in the existing
Makefile, preserving VPATH. Omitted USE_CUDA is also CPU-only. Primary RT,
material emission and IR transport/absorption retain their existing OpenMP
kernels; feedback remains the same native Fortran. `auto` uses only CPU
batches when the device implementation is absent. Explicit `cuda` rejects at
startup; unavailable GPU diagnostic ABI calls return errors without output
writes. No physics, source identity or HDF5 format change.

New local build `.snrt-cpu.OKoz9T/ramses_cpu3d`, SHA256
`7f875ab3a7cbfb6bb755d91de1f880cb03564d0e6e51563fca4c7aa276187bef`.
`build-ordered.log` records successful clean CPU compilation/link with no
nvcc, HYDRO_CUDA or CUDA libraries. `ldd` likewise shows no CUDA dependency.
Intel MPI/Fortran/C++/OpenMP and HDF5 remain required. Initial `-j4` exposed
missing base-module dependencies; amr_parameters->amr_commons->random is now
explicit. The first clean build used `-j1`; no general build-graph rewrite.
Never reuse CUDA-compiled objects for a CPU build or vice versa.

Also fixed a concrete pre-existing diagnostic ABI mismatch:
`snrt_cuda_weighted_sum_fp32_c` returns its error as a C function result,
not a sixth output-pointer argument. The Fortran interface now matches it.
Reused CUDA smoke: real GPU angular/ledger absolute errors1.1921e-7/7.7486e-7;
CPU unavailable-device mode rejects both and preserves sentinel outputs.

Reused hybrid smoke passes all primary/dust/IR batches and rollback with
GPU=0 (primary exact to OpenMP). Dust native reference parity5.063e-16 for
constant capacity,6.521e-16 for U(T), in auto and explicit OpenMP. Two-rank
backend initialization/operator smoke passes; this is not a new multi-rank
domain/AMR qualification. Forced CUDA test exits2 at backend initialization.
mkrun terminal/GUI now selects `cuda_linked` (unchanged default executable)
or `cpu_only`, and rejects forced CUDA with the latter before file creation.
The existing main namelist schema is unchanged. GUI/model30 tests pass,
one display-dependent skip. No additional Python runtime or gate framework.

Full CPU live `live/physical.nml`: one rank/two OMP, auto primary/dust,
level3/Z=.004, same AGB-pulse/LC18/SNIa/BPASS/AGN/DL01 reference physics;
four steps complete in24.351s. Before launch: noutput1/aout2/tout1e30,
foutput2/fbackup1000000, two expected~56MB dumps,162TB free. Source flags and
worker stacks512M unchanged. All161 floating HDF5 datasets finite; maximum
array-normalized hydro/SNRT difference from previous CUDA-linked/OpenMP run
2.6642e-15. Source identity, birth mass/clocks/Z/positions exact; largest
stellar mass/velocity absolute differences1.034e-25/1.006e-20.

`restart/physical.nml` uses a preserved COPY of the previous CUDA-linked
step2 dump, then CPU auto to step4 in7.080s, one new~56MB dump. SNRT and
source identity exact, hydro momentum maximum absolute difference6.204e-25
(maximum array-normalized4.221e-19); star masses/clocks/Z/positions exact,
velocity difference<=4.236e-22. Both complete without MG nonconvergence.
No production calculation, external audit, commit or push in this continuation.
All #1--#3 remain authorized; closing this portability subset is not a claim
that all remaining physics or backend expansion is complete.

### #2 primary dust scattering, native CPU/CUDA/hybrid (2026-09-08)

Resumed actual implementation after the operator queried the repeated stops;
no new approval request. Implemented an explicit `isotropic_elastic` primary
scattering comparison in patch/lagRamses, not a Python runtime prototype.
The old default is `none`. Draine C_ext*albedo is sampled at the existing
source-ledger representative energies; positive intervals use log-log
interpolation and exact-zero endpoints avoid log(0). Same original optics,
thermal grid and DL01 U(T), with a separate v4 reference sidecar:
`simulation/snrt/config/dust_dl01_bulk_030_scattering_reference_v4.nml`, SHA
`88d7efbd290f5cfc40775d36660cad3aaa47b11781d71b09c6f107464d9e76c3`.

Algorithm: exact local isotropic elastic mixing of direction-integrated
photon bins, tau=nH*relative_dust*sigma_sca*reduced_c*c*dt, first-order split
after transport/absorption on owned leaves. Each nonlinear trial starts from
the same incoming state, and the existing transaction commits only once.
No absorption/heating contribution from scattering; group photon number and
energy conserved to FP32 rounding. Shared CPU/CUDA cell kernel, nonblocking
stream lease per cell batch; busy stream means the OpenMP worker executes
that batch. Forced CUDA fails instead of silently switching; whole-call
outputs remain untouched on any late-batch failure. New native layout is
(direction,group,cell), independently checked at the Fortran boundary.

This is explicitly not the measured anisotropic Draine phase function; no
HG scattering, gas/grain recoil or radiation-pressure coupling, Doppler,
IR scattering, stochastic heating, sublimation, or unresolved diffusion
limit/AP claim. It remains a selectable reference comparison, not production
approval or completion of #2. Existing material compatibility limits remain.
Source definitions: https://www.astro.princeton.edu/~draine/dust/dustmix.html

HDF5 binds actual scattering coefficients and model marker in the existing
dust attribute (1640 values for old none; 1674 when enabled). Old default
identities remain unchanged. mkrun GUI/terminal selects the sidecar and a
capable CPU/CUDA-linked executable; generated README states limitations.
No new main RAMSES namelist keyword, no unrelated schema/AMR gate expansion.

Preserved old binaries; new CPU executable SHA
`8d75b5251e6910387d14a843789cf07e97e21a6741e127c54979afd60b8207bf`
at `.snrt-cpu.OKoz9T/ramses_scatter_cpu3d`; CUDA-linked SHA
`c832f47e895119d37f95a2495d47bc56f8cc5f4535f56417d7a7753cb2ef35ed`
at `.physical-extension.7rcxv4/ramses_scatter3d`. Both builds pass. Existing
hybrid smoke covers zero opacity, unequal angular weights, analytic solution,
two half steps, positivity/conservation and atomic rollback. On visible A10
physical device1, auto scattering used CPU16/GPU1 batches; occupied stream
CPU17/GPU0; forced CUDA CPU0/GPU17. Fortran analytic/layout test passes on
CPU/forced CUDA; existing material/IR parity <=1.871e-15 on GPU. CPU two-rank
operator test passes; this is not a full multi-domain MPI qualification.
Existing Intel/GNU contract tests pass opt-in/reset, thermal/export tests
pass, mkrun31 tests pass with one display-only skip.

Live evidence `.dust-scatter.CWEvh3/`: audited before launch, one rank/two OMP,
auto primary/dust, noncosmo level3/Z=.004, same physical source/feedback/AGN/
BPASS/DL01 control and 512M worker stacks. live/off both finish four steps,
foutput2 (two ~56MB dumps each), noutput1/aout2/tout1e30 unreached,
fbackup1000000. GPFS161TB free. `live` takes48.196s while running alongside
off. The parent `cooling` timer includes SNRT work; it cannot separately
attribute cooling versus scattering cost. No MG nonconvergence.
All161 floating datasets finite. Default-off regression versus the earlier
CPU run is not bit-identical: max array-normalized hydro2.6642e-15,
SNRT1.9984e-15; near-zero sink angular momenta differ <=4.66e-24 absolute.
An initial all-dataset exact-equality diagnostic correctly failed; no exact
regression claim is made. On/off SNRT level3 differs2.2766e-2 in normalized
payload, dust material energy1.011e-3 relative; these measure short-control
model response, not a convergence/accuracy bound.

`restart` copies (does not move) live step2, advances to step4 in7.23s.
SNRT arrays and dust identity exactly match live; hydro momentum differs
<=6.204e-25 absolute (4.221e-19 array-normalized). Particle slot order changes
on restart, so physical stars must be aligned by unique positive ID, not
raw array position. All512 unique positive-ID/PTYPE_STAR=1 stars have exact
mass, birth mass/clocks, Z and positions; velocities differ <=4.236e-22.
`bad-restart` uses the same copy but the old none sidecar:
MPI_Abort10, SNRT checkpoint rejection13 at init_hydro, before new output2.
All input/output evidence is preserved; no production calculation, external
audit, commit or push. #1--#3 remain preapproved, not approval-held.

### #2 native gas/dust/IR thermal exchange (2026-09-08)

Continued implementation without requesting approval. Added explicit
`hydrogen_accommodation` to the existing v4 reference dust contract and mkrun
GUI/terminal, with `none` default. No main RAMSES namelist keyword, new Python
runtime, external audit or generic gate framework. The physical collision
law is the geometric hydrogen-equivalent accommodation of McKinnon et al.
2021 eqs20--21 (https://academic.oup.com/mnras/article/502/1/1344/6067372).
This implementation uses finite existing U(T), not that paper's instantaneous
grain-equilibrium approximation. The control explicitly adopts area/H
3.495e-22 cm2 and alpha=.5: effective0.1-micron spheres of density3g/cm3 with
the existing reference dust mass1.398e-26g/H, not a recovered WD01 size law.
No electron/ion Coulomb, H2 colliders, grain evolution or radiation pressure.

First implemented the conservative implicit isolated gas/dust operator on
both CPU and CUDA with nonblocking leases. Its analytic BE residual, heating,
cooling, zero-transfer, stiff and late-error rollback tests pass. However,
inserting it AFTER IR failed in the actual live control at step0: the dust
reservoir would leave the material table because it could not radiate during
that separate kick. The failure returned0 through legacy clean_stop but was
recognized from the pre-commit rejection log and absence of output. Preserved
`.dust-exchange.WKZl1j/live/run.log`, initial binaries
`ramses_exchange_cpu3d` / `ramses_exchange3d`; they are NOT mkrun's selection.

Corrected the physics integration: eliminate the BE gas equation analytically,
Q=r*(Eg-Cg*Td), r=Kdt/(Cg+Kdt), and solve material+IR emission+Q together in
the existing monotone material root/fixed-point IR absorption solve. Reuse
CPU/CUDA/hybrid material dispatch; material mode2 appends gas E/Cv/K and Q,
while old modes0/1 retain their ABI and independent Fortran parity. Gas Cv
and thermal speed freeze per existing IR substep; conductance updates between
substeps. Gas Q enters the existing total-energy residual. No accepted gas,
dust or IR state changes on failure. Only the accepted sumQ is passed back
to RAMSES total gas energy, with kinetic energy unchanged. Chemistry remains
operator split, and the existing bath/material domain is not extrapolated.
The isolated operator remains available as a no-radiation reference, but the
live driver uses ONLY the joint material/IR path.

Sidecars: `dust_dl01_bulk_030_exchange_reference_v4.nml` SHA
`05670f96f7f461db95e46d3b3b232da9d02c49297b566a5d48cb6d799500100a`;
`dust_dl01_bulk_030_scattering_exchange_reference_v4.nml` SHA
`b1f27df3a6fcd2f76b7f6ca8a14c794b5811540905cb1c04480586bf7614909d`.
The latter has1678 HDF5 dust-identity values (old scattering1674 + model6,
joint algorithm2, area, alpha). All old no-exchange identities remain intact.

Native CPU build `.snrt-cpu.OKoz9T/ramses_exchange_coupled_cpu3d`, SHA
`ba6bcc99843a69b199693b30dd376c6428b87f0507e99bb3403e5bb01cff61c7`;
CUDA-linked `.physical-extension.7rcxv4/ramses_exchange_coupled3d`, SHA
`23b95ab47aa0b32fdc58151b090f35aa84ab2e195af4a089f46b107d5eb24459`.
Existing joint Fortran tests pass on CPU, forced CUDA and two MPI ranks:
stiff gas cooling into IR, total energy, material domain and rollback of gas,
dust, IR and photon outputs. Extended existing hybrid material test mode2:
auto CPU15/GPU2 free, CPU17/GPU0 occupied, full output parity with CPU and
energy closure pass. This does not constitute a full multi-domain AMR run.
Intel/GNU contract opt-in/reset, thermal exporter and mkrun32 tests pass
(one display-dependent skip). No extra external auditor was invoked.

Actual short evolution: `.dust-exchange.WKZl1j/live-coupled/physical.nml`,
one rank/two OMP, CPU-only auto, noncosmo level3/Z=.004, same LC18/AGB pulse/
effective SNIa/BPASS/AGN/DL01+isotropic-scattering reference. Native hydro,
sink and feedback flags rechecked; MG nonconvergence abort enabled. Four
steps complete27.390s. Scheduled noutput1/aout2/tout1e30 unreached; foutput2,
fbackup1000000, two~56MB dumps,161TB GPFS free at launch. The parent cooling
timer includes SNRT and is not an isolated collisional-kernel cost measure.
Joint balance residual <=5.0311e-10 under1e-9 tolerance; no MG failure.

`restart/physical.nml` uses a COPY of live step2, CUDA-linked build, visible
physical GPU1, primary OpenMP/dust auto, then step4 in9.531s. Material, IR
transport and absorption each visibly use CPU1/GPU1 batches. All161 floating
datasets finite. Final hydro arrays, source identity and all512 positive-ID
PTYPE_STAR masses/clocks/Z/positions/velocities exactly match continuous CPU.
SNRT maximum absolute difference3.312e-277 (not an exact-bit claim); dust
identity exact. Reconstructed dust temperatures22.682--79.642K remain inside
the material table. `bad-restart` disables exchange with old scattering
sidecar: MPI_Abort10 / checkpoint rejection13 before new output2. Existing
outputs and failed-trial evidence preserved. No production or commit/push.

This closes the native hydrogen-equivalent gas/dust/IR connection, not all
of #1--#3 or publication readiness. Remaining dust choices include actual
size/composition-matched collision areas, electron/ion charging, anisotropic
scattering, radiation pressure, grain evolution and cosmological-domain
qualification. Source/population gaps previously listed remain open.

### Push and next thermal-coupling bundle (2026-09-08)

User requested commit/push then continued implementation. Committed the prior
native source, CPU-portability, scattering and gas/dust/IR work as `8e7491f`
and pushed `origin/main` to `kjhan0606/LagRamses`. Forty-four code/config/doc
files; no scratch binaries, source archives or simulation dumps. Preserved
the unrelated 88-line deletion in `ramses_nml_generator.py` outside the index;
only this project's history-field help update was included. Remote and local
main matched before the commit. Work stayed under `/gpfs/kjhan/LRD_JWST`.

Next bounded #2 bundle implements thermal-speed-dependent gas/dust exchange
inside the already joint native material/IR solver. The collision law still
follows [McKinnon et al., equations20--21](https://arxiv.org/pdf/1912.02825);
no new physical coefficients, grain species, source range or approval gate.
Instead of freezing K at the start of the IR substep, solve
`Cg*(Tg-Tg0) + dt*K0*sqrt(Tg/Tg0)*(Tg-Td)=0` along with the material equation.
Gas Cv/chemistry and geometric/accommodation factors remain fixed during
each IR substep; this is not a coupled chemistry or multi-species collision
network. The pre-existing no-radiation helper is still explicitly fixed-K.

The gas equation reduces to a positive-root cubic in sqrt(Tg), scaled to
avoid stiff-coefficient overflow. Bracketed Newton solves it between the old
gas and trial dust temperatures. Q is evaluated without subtracting nearly
equal gas energies in the weak limit. The outer material equation includes
this Q at both temperature-domain bounds and each trial. Failure rejects the
whole existing native transaction, not a CPU replay or temperature clamp.
Modes0/1 retain no-exchange layouts; old fixed-K mode2 remains available in
the low-level ABI, new mode3 is selected by the Fortran live callback. CPU,
CUDA and the existing nonblocking stream/OpenMP dispatcher share the kernel.

Existing native tests extended, no new test framework: independent BE gas
residual/bounds for heating/cooling, zero and stiff coefficients; fixed-Tdust
analytic ODE refinement in both directions; material modes0--3 parity,
conservation and late-error rollback; joint gas/material/IR Fortran checks.
The analytic relation is `u=(sqrt(T)-sqrt(Td))/(sqrt(T)+sqrt(Td))`,
`u(t)=u(0)*exp(-alpha*sqrt(Td)*t)` for `dT/dt=-alpha*sqrt(T)*(T-Td)`.
At32/64/128/256 steps, heating T0=1 errors are .0238667/.0119272/.00596190/
.00298051 K; cooling T0=1000 errors7.99556/4.01593/2.01255/1.00743 K.
This establishes first-order convergence for the isolated gas equation,
not arbitrary full-AMR time convergence. An initial T0=5 ratio-only test
was unsuitable across the ODE's curvature change; used monotone-curvature
heating/cooling controls, retaining independent residual/bound checks.
CPU and CUDA native tests pass. Mode3 mixed dispatch uses CPU15/GPU2 free,
CPU17/GPU0 with its sole stream held. Joint Fortran tests pass CPU, forced
CUDA and two MPI ranks. GUI/model32 tests pass with one display skip.

Built separate executables, preserving old ones:
`.snrt-cpu.OKoz9T/ramses_exchange_nonlinear_cpu3d` SHA
`8c7fa6c124dba287c34541b6197c69782610b54e59a43dee947e067c253bd7b4`;
`.physical-extension.7rcxv4/ramses_exchange_nonlinear3d` SHA
`1c5ae8665fbf5c24dbd92abf97826f56471311cb2a4e4802aceeab0d809f3486`.
VPATH and native main-namelist fields unchanged. mkrun now selects these
only when exchange is enabled and documents the solver/restart distinction.
The v4 physical sidecars are unchanged. Exchange HDF5 algorithm marker2->3
is the only changed entry (index1675, zero-based) of the1678-value dust
identity. Old no-exchange identities remain unchanged.

Evidence directory `.dust-nonlinear.TtLIOc/`. Effective live namelist
`live/physical.nml` audited before launch: noncosmo level3/Z=.004,
SNRT/DUST_LIVE/HDF5/NVAR30/NVECTOR500, one MPI/two OMP, 512M thread stacks;
LC18/AGB pulse/effective Ia/BPASS/AGN/scattering+exchange reference unchanged.
Four steps, noutput1/aout2/tout1e30 unreached, foutput2/fbackup1000000;
two~56MB dumps,161TB GPFS free. MG nonconvergence abort enabled, none seen.
CPU live completes26.259s; joint balance<=5.0311e-10 under1e-9 tolerance.
CUDA-linked restart from a COPY of step2 completes step4 in7.425s with
material/IR transport/absorption each using CPU1/GPU1 batches (visible GPU1).
Timings are short-run observations, not a performance speedup measurement;
the cooling parent timer includes SNRT.

All161 floating datasets are finite. Continuous/restart hydro differs only
in momenta, maximum absolute1.8612e-24 (normalized<=1.27e-18), while SNRT
maximum absolute difference2.728e-253. Dust and source identities match;
particle arrays must be matched by unique positive star ID, not file order.
All512 matched stars have exact masses, birth/progress clocks, metallicities
and positions; velocity differences are at most4.2352e-22 in code units.
`reject-old/physical.nml` uses a COPY of the old fixed-speed step2 and is
rejected before output2 (radiation checkpoint status1 / MPI_Abort10).
No simulation outputs were deleted, overwritten or moved.

This closes temperature-dependent thermal speed in the existing native
gas/dust/IR coupling, not all of #1--#3. Previously listed source/population
data gaps and additional dust/radiation physics remain explicitly open;
no new verification bundles or external audits were inserted.
