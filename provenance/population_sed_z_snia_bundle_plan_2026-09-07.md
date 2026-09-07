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
