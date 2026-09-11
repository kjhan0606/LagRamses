# Bundle 4: effective-hybrid frozen He-contact comparison

Operator authorizes implementation as a named approximation, not a
self-consistent binary/WD population. Repository: /gpfs/kjhan/LRD_JWST,
origin git@github.com:kjhan0606/LagRamses.git. Main includes this design in
the already planned aggregate selective Fable review; no separate audit.

## Physical definition

Use retained COSMIC4.2.0 Z=.01 system169 (initial7.5+5.274Msun), first
post-common-envelope CO-WD/He-star row: age85.96190548931196Myr,
WD1.034476609229475Msun, donor1.0071270883432488Msun. The pair is detached:
contact is IMPOSED, not inferred from recorded Roche geometry. Donor/net-WD
deltam are not measurements of incident helium transfer. Prescribe
2e-6Msun/yr, eta=1 stable He burning using Wang+2017 eqs1/2 and its
2.05e-6 central/off-centre boundary. Apply this solar-Z reference at Z=.01;
freeze the mass argument at1.35Msun up to the explosion threshold. This
continuation and instantaneous threshold ignition are comparison assumptions,
not resolved ignition. No KH04 unsupported flashes or COSMIC0.15Msun
He-trigger disappearance are turned into N100.

Threshold equals unchanged N100 returned mass/debit1.4004633930489443Msun.
The analytic event age is86.14489888122169Myr. Converter verifies the event
lies within the retained history, donor background mass loss plus prescribed
transfer fits its finite donor budget, and the reference rate lies within
the selected stable/central interval. N100 gross elements, untracked
residual, kinetic energy and zero terminal remnant remain unchanged.

## Source and normalization

One event row (Z,age_yr,weight_per_initial_SSP_Msun) is embedded in an
additional snia_frozen_he_event group inside the existing ordered runtime
sidecar. Population event_model selects it; omission selects the historical
empirical power law. No new main namelist/environment variable.

Weight1.8439455443768248e-5 is copied from retained grid-v2/z0.01-grid.csv
system169. It already includes the illustrative primary-IMF, binary-system
fraction, secondary/period quadrature and FULL initial population mass
denominator. Do not renormalize to selected/exploding binaries, multiply
binary fraction again, or pass it through the single-star IMF integrator.
It is not an observed Ia rate or a matched individual-star SSP IMF.

Only birth Z=.01 is admitted (roundoff tolerance only). All ages use the
original ZAMS clock; age units convert yr -> Gyr exactly once. Evaluate
Ninterval=Minitial*(CDF(tnew)-CDF(told)), with CDF(t)=w if t>=tevent else0.
Thus old<event<=new, equal endpoints give zero, and the same CDF supplies
the prior-return reconstruction. The single-event model has no additional
late events; it does not extrapolate donor evolution beyond100Myr.

## Native accounting and restart

Only effective_ssp is admitted; strict_wd rejects this event model. Ordinary
SSP returns remain unchanged: the comparison replaces empirical Ia counts,
not ordinary returns, and supplies no second binary ordinary-ejecta source.
Existing aggregate remaining-particle accounting rejects overdraw/replay.
This does NOT prove disjoint microscopic progenitor ownership.

The event type stores all numerical row/closure parameters and SHA256
bindings of retained BPP/BCM/grid/selection/Params, converter, retention
reference and N100 input. A module helper returns the complete numerical
identity and string bytes. The empirical helper result is empty, preserving
old identities. Main appends this to phase0_snia_identity so existing MPI
source consensus and HDF restart comparison bind both contents and model.
Changed event age/Z/weight/parameters/source hash cannot reuse old identity.

## Exact main-owned runtime integration

In stellar_ramses_runtime.f90:

1. Both evaluate_snia_interval_events calls (interval and cumulative prior)
   append birth_metallicity=population%birth_metallicity. Old callers remain
   valid for the empirical model; the table rejects absent birth Z.
2. phase0_snia_identity appends snia_event_model_identity output. Empty for
   empirical, complete versioned payload for the frozen model.
3. load_snia_runtime_contract calls validate_snia_event_source after reading
   the physical source, binding threshold to its returned/debit masses,
   zero remnant and N100 source identity. Retain existing unity factor guard;
   table Z selection is a separate actual birth-metallicity argument.

No stellar_yield* changes; those remain Gibbs-owned. Native population
module, converter, fixture and this document are this worker's scope;
general runtime/HDF/setup remain main-owned pending coordination.

Exact stable helper interfaces (all in stellar_snia_population_contract):

```fortran
evaluate_snia_interval_events(realization, initial_mass_msun, age_old_gyr, &
    age_new_gyr, metallicity_factor, expected_events, ierr, birth_metallicity)
! birth_metallicity: optional real(stellar_dp), mandatory for the table.
snia_event_model_identity(realization, values)
! values: allocatable real(stellar_dp), intent(out); empty for empirical.
validate_snia_event_source(realization, returned_mass, wd_debit, &
    terminal_remnant, source_id, source_sha, ierr)
```

Native file representation uses `row%...` members of
`snia_frozen_he_event_t` in the additional group. There is one Z/event row,
not a general table framework. Legacy DTD parameters remain present but do
not multiply or reshape this model's CDF. `events_per_initial_msun` must
match the row weight, with Kroupa ID, binary_fraction=.5 metadata,
baked-in binary normalization and imf_conversion_factor=1.

## Bounded verification

Existing native SNIa tests plus one focused fixture exercise default parity,
missing/wrong Z, strict rejection, half-open endpoints, split/cumulative
agreement, actual N100 event construction, effective overdraw/prior-return
rejection and identity changes. Private Intel/GNU native builds only here;
no MPI launch, production jobs, commits or raw-output deletion.

## Implementation and measured native evidence (2026-09-11)

Main connected the two actual-birth-Z calls, event identity append and
startup N100-source validation in stellar_ramses_runtime. Read-only check
confirmed argument/signature agreement. This worker did not edit shared
runtime, general HDF, setup or stellar_yield*; Volta's isotope scope is
unaffected. Aggregate design review and full integration build/run are main's.

Converter: simulation/snrt/tools/build_snia_frozen_he_contract.py.
Complete ordered sidecar:
simulation/snrt/tests/fixtures/phase0/snia_frozen_he_runtime_v1.nml.
Native test: simulation/snrt/tests/fixtures/phase0/fp2_snia_frozen_he_test.f90.
Private build root: /gpfs/kjhan/LRD_JWST/.snia-frozen-he.RYbxQw.

The converter verifies the retained seed, birth clock, finite donor budget
and stable reference, and produces age86.14489888122169Myr and donor
survivor0.6410843108150832Msun after background loss. Two independent renders
match the checked-in fixture byte-for-byte. Physical N100 and thermal groups
are byte-identical to the existing effective-SSP template.

GNU13 and Intel2025 native bounds-checked coupling tests PASS:
`FP2_SNIA_FROZEN_HE_NATIVE_COUPLING_PASS`. Tests execute the event evaluator,
N100 budget, source-to-cell mass/elements/energy and effective SSP accounting.
They verify half-open endpoints, cumulative/split agreement, no restart
replay from reconstructed prior mass, strict/missing-Z/wrong-Z/factor
rejection, donor and aggregate shortfalls, and changed identity payloads.
All seven existing SNIa native fixtures also PASS on both compilers: DTD,
population, physical, deposition, event ledger, runtime contract, accounting.
No RAMSES/MPI execution or full HDF restart comparison was performed here.

Physical reference: Wang, Podsiadlowski & Han2017,
https://academic.oup.com/mnras/article/472/2/1593/4094897, equations1/2 and
central-growth reference. Its Z=.02 models, point-mass stable-growth
approximation and uncertain ignition mapping are not promoted into measured
Z=.01 binary evolution by this fixture. All prescribed choices remain
explicit, including original COSMIC detachment and the frozen1.35Msun tail.

## Actual RAMSES-runtime native test (2026-09-11, after Fable PASS)

Private directory: `/gpfs/kjhan/LRD_JWST/.snia-native-runtime.c21WIP`.
`runtime_event_test.f90` compiles and executes the actual current
`stellar_ramses_runtime`, all current stellar source modules and the real
unit conversion/cell locator, with Intel2025 `-O0 -g -check all -traceback`
and OpenMP. It reuses the existing synthetic `high_mass_snia_z` ordinary-SSP
fixtures and the physical frozen-He/N100 sidecar. No source/receiver mocks;
only the fixture's fatal `clean_stop` terminates with an error instead of
RAMSES shutdown. The unused cooling dependency `f_de` is copied verbatim
from the current native source into the private link, without modification.

Exit0: `SNIA_ACTUAL_RUNTIME_EVENT_PASS`. For initial SSP mass10000Msun,
the runtime deposits0.184394554437682 events and0.258237823367545Msun.
Assertions cover initial-SSP normalization, particle/cell mass, bulk momentum,
N100 Fe, thermal plus bulk kinetic energy once, actual-birth-Z rejection
(ierr78) with unchanged particle/cell/progress, retry, same-age no-op,
continuous versus split agreement, restored particle-mass/progress no replay,
and inconsistent prior-Ia mass rejection (ierr83) without publication.
The public runtime identity contains5552 values and ends with the exact
823-value event-model payload. This tests both actual-Z call sites and the
successful source-binding/identity integration, not only helper routines.

An initial mixed link using the September7 shared `nbors_utils.kjhan.o`
crashed in the locator; recompiling that actual source privately against the
current common modules resolved it. The passing executable uses that fresh
object. Main's build must not reuse the stale common-module-dependent object.
Remaining Intel406 diagnostics are array-temporary warnings only.
No production source/interface changes were needed. This is a serial native
particle-to-cell test, not a RAMSES/MPI launch or HDF checkpoint/restart test;
the ordinary SSP remains synthetic and no full WD/binary microclosure is claimed.

## Actual MPI event and restart: PASS

Evidence `.snia-live.hUBFAz/metrics.json`, logs, effective inputs and frozen
binary SHA256 `985fd740526206a9b17290010118ea805f7e5853d7bcffdcbd620ff27d1681f7`.
MPI2/OMP1,64 uniform cells, four steps. The first attempt in
`.snia-live.94g7SR` formed no stars because the uniform-grid reader skipped
REFINE_PARAMS/mass_sph. The bounded reader fix consumes the optional group
on uniform grids; existing refined-grid requirements are unchanged.

In the retry,64 SSP particles born at68.645376Myr and Z=.01 cross the
86.144899Myr event age after the step2 checkpoint. Ordinary synthetic sources
are already flat before that checkpoint, isolating the later Ia receipt.
Expected count180918.5924217 is an effective population expectation, NOT
180918 individually resolved binary systems. Ia particle debit is
253369.8658085Msun, matching N100 to1.21e-13 relative. Gas Fe increment
133879.7583874Msun agrees to3.46e-11; energy agrees to2.40e-15. Subtracting
large gas totals gives1.10e-8 relative error in the small return increment,
while total gas+stars closes to2.17e-16. All83 datasets and51 compared
attributes agree after restart; no duplicate Ia debit occurs.

Four evaluated raw HDF files (450848bytes) were removed after evaluation;
inputs, binary, logs, hashes and metadata remain with cleanup.json. No raw
copies are retained; recreate by rerun. This qualifies the named frozen-contact
effective comparison only, not disjoint WD/donor ownership or a calibrated
cosmological binary population.
