# Actual sources with CHIMES transition and evolving grains

Status: completed bounded source-coupling bundle; driver end evaluation PASS.
This is not universal production/scientific approval. Earlier in-progress
and failed-run sections below are historical; final disposition is at EOF.

Operator approved one bundle: existing stellar/AGN radiation into kind7,
joint mechanical feedback, and a small integration plus restart. Actual
project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`, base `0136f39`.
No new physical source model, universal production claim, or new audit gate.

Existing source calls already enter the paired N/E transaction. Stellar
feedback precedes chemical reconciliation and hydro; dust mass evolution
precedes RT. This bundle checks those existing calls jointly and repairs
concrete faults only. Reuse the PARSEC eleven-Z matched high-mass population,
its existing atmosphere closure, and `partition_reference_v1` AGN. Neither
is promoted to a full calibrated population by this test.

AGN currently requires NENER=0. Respect that restriction: no CR or MHD in
this integration. Do not remove the guard without separately addressing
nonthermal energy in AGN material loading and temperature limiting. Existing
stellar thermal/advective-CR evidence remains separate. Fe/PAH, grain drift,
condensation, SN shock processing, and sublimation remain excluded in kind7.

Private build/input/evidence root: `.chimes-sources.r1TaD8`. Fresh serial
compilation reuses the repository Makefile and unchanged VPATH; NVAR187,
NENER0, NVECTOR32, SNRT/DUST_LIVE/CHIMES/HDF5, CPU/OpenMP. Same ABI6 CHIMES
library and physical data banks as the previous two completed bundles.

Final effective input: `.chimes-sources.r1TaD8/coupled-ready2/physical.nml`.
Purpose: bounded short evolution, 4^3 periodic noncosmological hydro,
gravity, actual star formation and initial sink/Bondi accretion, live gas
chemistry and evolving C/silicate grains. Stars and sink start with zero
emitted photons; no photon/restart seeding. PARSEC wind/CCSN/pair returns
enabled, AGB/SNIa deliberately off (not matched by this population).

Output policy before launch: noutput1, aout2/tout1e30 (outside test),
foutput2, fbackup1000000, nstepmax4. Two new full dumps, plus one restart
input copy and one resumed dump. Initial gas-only estimate10--12MiB/dump
was superseded by measured source-bearing HDF5 size about62MiB; diagnostic
and restart working allowance is below1GiB. Free space about100TiB.
Retain evaluated inputs/logs/hashes/results; remove exact raw snapshot
directories only after integration and restart evaluation are complete.

## Focused review and implemented material connection

[Fable](snrt_actual_source_fable_2026-09-11.txt) returned CONDITIONAL with
Q-GOAL yes and Q-LEAN yes. The driver adopted the physical requirements:
extend the checked scalar map, retain gas-vs-solid energy conventions,
test material budgets/rollback and audit actual source receipts. Accretion
scales every material carrier by GROSS removed mass, not net BH growth.
Gas internal energy stays in the donor during existing cold jet loading;
solid thermal energy travels with the loaded material, outside gas energy.
Accreted solid energy leaves the modeled grid with swallowed material; it
is not injected again as AGN heat or claimed to be a new BH thermal ledger.

The current selector appends eight dust densities and all157 chemical
densities; it does not append disabled shock/fresh-source slots. The pure
map retains duplicate/reserved/extent checks. Existing native helper tests
now exercise nonzero ions/molecules/electrons, gas+solid elements, charge,
aggregate/four-bin identities, gross accretion, jet load/receiver balance,
separate solid U and invalid-state rollback. The existing complete AGN
efficiency/source/deposition suite passes with Intel bounds checks.

The suggestion to omit mkrun was not adopted: the operator specifically
requires both namelist generators to track supported runtime settings.
The new existing-sink option is default-off, kind7-only, requires an explicit
NENER0 binary and leaves the original mass/CR profiles unchanged. Shared
GUI/generator tests pass49 total,48 pass plus1 display skip. No new namelist
field or automated audit framework was added. Full157 map testing is retained
to catch trailing-field omissions; no new per-step material diagnostic.

SINK_PARAMS is actually read AFTER read_hydro_params. The first attempted
`coupled` run exposed the initial implementation checking its startup
defaults; it was stopped before evolution and is not evidence. Final scheme
and create_sinks checks now run after SINK_PARAMS, while existing AGN model
consensus/admission remains authoritative. The model lookup itself is lazy
and can safely be checked earlier; no reference-admitted flag is consulted
before admission. `combined/physical.nml` is an unused initial draft with
wrong offsets; launched `stellar`/`coupled*` use H at passive2, ichem7 and
exact5/3 gamma. No existing snapshot was rewritten to repair a layout.

## Integration progress and failures (not passes)

`stellar` (no sink) forms real stars and completes three RT/IR commits. Its
step2 checkpoint has128 stars and returned mass3.63127e-9 code, gas+solid
element residual8.47e-16 and charge5.30e-16. At step4 it rejects photo
status51. Diagnostic continuation from its untouched step2 checkpoint
reproduced HCO+ = -1.6958303447609539e-167 in CVODE's accepted internal state,
not merely interpolated output. This is NOT a completed stellar run or a
reason to waive the nonnegative-state check.

The photo solver now checks each accepted internal step and retries a
negative step from the complete last valid state at half step size. Spectral
nodes, rates, accumulated optical depths and ledgers remain unchanged;
no abundance clipping or physical tolerance relaxation. Nonfinite states,
more than64 retries or10000 attempts still reject. Native thermochemistry
passes256 checks, including a vanishing HCO+ tail. `stellar-positive` tests
the former failure segment with this fix.

`coupled-short` uses a100-times smaller Courant factor for a short source
wiring check, not a long-time/SN-lifetime validation. Its particle exhaustion
was NOT solved by increasing capacity to65536 in `coupled-final`: at
levelmax2 the sink cloud radius equals the periodic box length, so offset
cloud particles alias the centre and survive cloud removal. The next cloud
creation multiplies their count (14763 cloud entries at step2). No generic
tree code was changed. Final `coupled-live` retains4^3 cells with levelmax3,
explicit refinement suppression, radiusL/2 and unchanged4step output schedule.
The allowed finest level also affects stellar particle mass and minimum AGN
deposition radius; this is a corrected numerical fixture, not a resolution
comparison. `coupled-ready` omitted mandatory REFINE_PARAMS and exited before
evolution; it is superseded. Failed attempts are not passed tests.

Measured first stellar HDF5 file61,638,240 bytes, about62MiB with auxiliary
output. Large matched-source identity arrays explain the difference from
the prior10MiB gas-only estimate. Revised total working output allowance
is below1GiB for evaluated runs plus necessary restart copies; available
space still about100TiB. Diagnostic raw input is retained until evaluation ends.

## Completed numerical regression

Final executable `ramses_sources_positive3d` SHA256:
`4ffaef1b0e79ed33f1983cb2e1e658eeb0dfda4ccab2decc55ae7fda4e86be26`.
`stellar-positive` completes the original step2-to-step4 interval after the
photo fix (MPI2/OMP2). At step4 it has256 actual stars, returned stellar mass
1.37787867516e-8 code, gas+solid element error7.05895e-16, charge error
3.91715e-16 and temperature1374.33--2529.96K. Both gas and grains absorb actual
stellar radiation; maximum IR balance residual9.1263e-10. This is a resumed
failure regression, not an independent fresh-vs-restart comparison.

`coupled-live` first rejected a retained `SNRT_RT_LEVEL=2` setting during
preflight, since multilevel live IR requires all levels. Its final launch
explicitly unsets this variable; the rejected log is retained separately as
`preflight-level-filter.log`. No partial evolution was accepted from it.

The subsequent `coupled-live/refined-attempt` input used `m_refine=1d99` to
try to suppress refinement, but with this artificial initial condition the
reference `mass_sph` is zero, so the density threshold remains zero. It
unexpectedly created64 level3 grids and rejected a nonfinite hydro state
after level3 chemistry and coarse upload. This is not a passed AMR test;
multilevel chemistry/material evolution requires separate investigation.
The final bounded fixture uses the actual disabling value `m_refine=-1`.
The failed effective input and log remain beside the final run.

## Superseded actual-source integration (not final evidence)

Final `coupled-live` completes4 coarse steps on64 leaves with MPI2/OMP2.
Its trajectory is nevertheless superseded by the empty-level Bondi bug
below. Completion and material residuals alone did not validate its rates.
Level3 remains empty, and the sink cloud stays at2109 particles whose sum
matches the BH mass (relative FP64 roundoff). Four RT/IR commits pass;
maximum IR balance residual7.5567e-10. Step4 contains256 stars with returned
mass2.14476677e-12 code, gas+solid element residual7.42605e-16 and charge
residual4.84052e-16. Temperature1603.76--1609.59K; grains and gas both absorb
radiation. The coarse AGN ledger records positive accepted accretion and
JET mode at steps2--4; peak bolometric luminosity1.62761e37 erg/s. Both
`average_AGN` and `AGN_blast` execute without rejected material transfers.
Pending radiation is consumed; mechanical queues remain restart state and
are not required to be zero at a checkpoint between operator-split calls.

This artificial source-wiring fixture is not a realistic galaxy/AGN
calibration: the Mpc unit box and initial BH1e-5 code imply2.45e11 Msun.
Its short duration tests actual wind/source and AGN coupling, not the full
SN delay-time population; the longer stellar failure regression and native
mechanical tests provide separate evidence. Live AGN uses JET mode here;
thermal/deferred material deposition is covered by the native AGN suite,
not falsely claimed as a second live AGN branch.

## Restart-discovered Bondi initialization bug

The first `restart` restored its checkpoint but rejected an invalid accretion
event before the first RT step. `init_sink` allocated four per-level weighted
arrays without initializing them, and `grow_bondi` summed even globally empty
levels. The empty level3 never enters `bondi_hoyle` to fill those arrays.
The prior fresh `coupled-live` was also contaminated: its ledger recorded
gas-relative speed22.056 code despite essentially stationary initial gas/BH,
explaining the spuriously small Bondi rate. Its pass status is withdrawn.

Both allocation paths now initialize weighted density/volume/momentum/c2 to
zero, and the sum excludes globally empty levels (also preventing stale
contributions after derefinement). No rank-local ownership filter is used.
This bounded fix needs no new monitoring or restart format. New
`coupled-verified` and `restart-verified` use `ramses_sources_final3d`; the old
midpoint is not reused for their comparison. Same effective output policy;
all raw snapshots including superseded runs remain within1GiB.

Fixed executable SHA256:
`060aa40b6b37afad8ac459eb8002dca69c4b11c301199f27e4a3921cee26011a`.
The corrected fresh step1 ledger has zero initial gas-relative speed, not
22.056. Step2 accepts8.49200431e-10 code BH mass and enters THERMAL mode
with bolometric luminosity2.57023e48 erg/s. This confirms that merely keeping
the old completed run would have hidden a scientifically significant error.

This fixes a real defect but is NOT yet the complete restart fix:
`restart-verified` still rejects accretion before RT, and the corrected fresh
run rejects photo status6 after its step2 checkpoint. The dedicated
`coupled-diag`/`restart-diag` runs record failure-only budget/input details.
The new trace-metal native case passes with the existing tolerances (native
258 total); that synthetic case does not reproduce the live failure and
does not justify waiving it. Final completion remains pending.

## Resolved failure causes and final verification

The diagnostic fresh run identifies sulfur, initial S/H
9.1584670172097568e-15 versus final9.1584789384892034e-15. Absolute error
1.19e-20 is below the old universal species absolute tolerance1e-16, but
relative error1.30e-6 violates the unchanged1e-8 nuclear budget. The solver
now sets species absolute tolerances from their limiting element inventories
using the receiver's own stoichiometry: min(1e-16,max(1e-32,1e-12*N_e/count)).
Electron tolerance and relative tolerance remain unchanged; no elemental
renormalization, clipping or budget relaxation. All258 native tests pass.

The diagnostic restart confirms r2k=0 and NaN weighted volume/rate. After
coordinate synchronization and MPI tree migration, the pre-migration local
canonical particle index is stale; the new owner fails to sample the centre.
`kjhan_refresh_canonical_sink_map` now rebuilds the map from owned active
particle lists after migration, using exactly synchronized coordinates and
requiring one global centre per sink. It precedes Bondi sampling and the
downward tree merge. It does not pick a nearby cloud particle or store a
stale Bondi radius. Refined/subcycled later migration remains outside this
uniform-fixture verification; the unrelated refined failure is not waived.

Final executable `ramses_sources_resolved3d` SHA256:
`34ee8827ac39237583a36d230f8822e37328b5f8bad9a83244e534f9b19e6cc5`.
`coupled-resolved` and `restart-resolved` use this same executable/input
physics and an independently regenerated midpoint. Their verification is
in progress. Evaluated obsolete raw snapshots were removed per the operator
instruction; see [cleanup manifest](snrt_actual_source_cleanup_2026-09-11.md).

`coupled-resolved` now completes all4 steps and passes its snapshot evaluation.
Final state:248 stars, returned stellar mass8.98069386e-13 code, BH mass
1.000104245e-5 code,2109 clouds with matching mass,64 leaf cells. Maximum
gas+solid element residual8.47035e-16, charge residual5.33755e-16, dust
redundancy residual4.23517e-22 absolute. Temperature1603.30--21004.25K;
maximum IR balance residual8.7671e-10. All4 RT commits pass; actual stellar
and AGN radiation is absorbed by both gas and grains. Corrected AGN enters
THERMAL mode at steps2--4; native tests cover jet material transport. The
superseded JET trajectory is not retained as scientific evidence.
`restart-resolved` is being compared from its regenerated step2 checkpoint.

Both resolved runs complete, but they are not bitwise equal. Gas energy
differs4.23e-13 relative to its peak, chemistry around1e-13, final clock
5.42e-20 code; these are small. A separate sink position-sum accumulator
differs by1 code length and cannot be dismissed on that basis. This input's
cloud radiusL/2 puts antipodal points exactly at the periodic minimum-image
branch in `move_fine`: an ulp perturbation selects +L/2 versus -L/2. The
instantaneous BH centres differ only2.22e-16, but their saved accumulator
would affect a subsequent update. No production claim is based on this
degenerate fixture.

Final `coupled-bounded`/`restart-bounded` keep the same4^3 actual mesh and
refinement disabled, but allow levelmax4 so cloud radius is strictlyL/4.
This also changes minimum-cell-derived stellar particle mass and AGN radius;
it is a corrected source-wiring fixture, not a refinement/convergence study.
The wizard's new sink option similarly keeps its uniform levelmin3 but
sets levelmax4 with explicit disabled refinement; the shared validator
explains the periodic support requirement. No particle-tree or periodic
force algorithm is rewritten for this artificial boundary case.

Final comparison records bitwise equality separately from measured floating
differences. It does not relabel unequal arrays as bitwise identical or
silently ignore a saved sink accumulator. Physical budgets/finite-state and
discrete identities remain required; clocks must agree within1e-12 relative.
The driver reviews the measured field differences, including sink statistics,
rather than adding a new general testing framework.

`coupled-bounded` then exposes a separate donor-search defect when it enters
JET mode at step3. In non-ksection ordering, AGN bin width is about rmax;
here rmax=0.0625 but half the actual donor-cell width is0.125. Searching only
one neighbouring bin misses the containing cell centre, despite the correct
half-open donor predicate. `average_AGN` now searches
max(1,ceil(0.5*dx_cell/bin_width)) bins per axis. Physical rmax and all jet/
thermal geometry tests remain unchanged; this only expands candidate lookup.
The current field-transfer binary is `ramses_sources_ready3d`, SHA256
`120256551c35d94c5fc0295d61e1a0c2c546fe1d2c1c2f67c495e5c60ba16e08`.
Final pair is `coupled-ready2`/`restart-ready2`, with the same bounded input.

`coupled-ready2` completes4 steps. It has245 stars,2109 clouds, final BH
mass1.000091912987e-5 code and temperature1604.26--2.87143e6K. All4 RT commits
pass; maximum IR residual5.8196e-10, gas+solid element residual6.35281e-16,
charge residual5.34671e-16. Global gas+stars+BH mass agrees to2.15e-16
relative with initial mass minus the epsilon0.1 radiated rest-mass loss;
cloud pseudo-particles are excluded from that mass sum.

The ledger records THERMAL at step2 and JET at steps3--4. On this deliberately
underresolved injection support, do not infer that every thermal entitlement
was deposited: final heat pending6.0618e59erg and jet pending8.2971e55erg
remain explicit checkpoint state. The high-temperature gas response and
material budgets are real; the queue identities must survive restart.
No claim of resolved jet morphology or a calibrated coupling radius.

## Final driver disposition

PASS for the approved existing-source/material connection, MPI2/OMP2 uniform
noncosmological periodic hydro with NENER0, kind7 evolving C/silicate dust,
actual matched high-mass stellar sources and reference non-MAD Bondi/AGN.
Final `coupled-ready2` and `restart-ready2` both complete and satisfy the same
finite-state, charge, gas+solid element, dust redundancy and global mass checks.

Restart comparison:449 datasets inspected,441 bitwise identical. ALL hydro
fields (including gas energy, dust/solid energy and157 chemical densities),
SNRT fields, gravity, BH mass/position, pending AGN energy/mass queues and
all six time/step attributes are bitwise identical. Eight differing arrays
are three nearly zero sink angular momenta, three velocity-sum statistics,
sink z velocity and particle z velocities. Maximum absolute differences:
angular momentum7.755e-26 code; statistics1.776e-23; velocities6.463e-27.
These are parallel summation roundoff, not the prior unit-length position-sum
ambiguity; all saved position sums now match exactly. Full measurements:
`.chimes-sources.r1TaD8/final-results.json`. This is NOT a449-array bitwise pass.

Native chemistry258 PASS; existing compiled AGN efficiency/source/deposition
suite PASS (including full material-map tests); frontend49 tests,48 PASS
and1 display skip. No new routine end audit: the driver performs this
evaluation under the approved reduced cadence. Fable's indispensable
material conservation/ownership/admission conditions were addressed; optional
extra instrumentation was not made a new gate. Both generators were updated
as explicitly required by the operator.

Limitations remain explicit: no CR/MHD, new sink formation, Fe/PAH, drift,
condensation, SN shocks or sublimation in this new admission. The unexpected
refined run's state-validation failure is retained as a separate unresolved
multilevel issue, not a passed AMR test. This small fixture neither calibrates
AGN coupling nor validates resolved jet morphology or the entire stellar
population. Do not reopen this completed uniform-source bundle as a sequence
of new instrumentation gates.

Evaluated raw dumps are removed under the standing retention instruction;
inputs, logs, binaries, build identities, compact JSON evidence and the
[cleanup manifest](snrt_actual_source_cleanup_2026-09-11.md) remain. The operator
subsequently requested committing/pushing this completed bundle before fixing
the separate multilevel AMR issue; this record accompanies that bundle commit.
