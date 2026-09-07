# Population-matched SED / physical Z coverage / microscopic SNIa bundle

User request: implement medium-term items 2, 3 and 5 after commit 357bc2f.
Project: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses. Final objective is usable
native RT/feedback/dust simulation physics with a defensible publication
model, not an expanding collection of Python gates. Preserve current model
defaults, approved empirical DTD, effective-SSP comparison, and old restart
identities. Existing ramses_nml_generator.py deletions are unrelated.

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
