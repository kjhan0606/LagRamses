# PARSEC low-Z printed-precision common source

Status: bounded source conversion, printed-budget, GNU native, repaired
MPI2 fresh/restart and driver end evaluation PASS; raw cleanup COMPLETE.
Not a full-SSP or all-medium-physics completion.
Work root `/gpfs/kjhan/LRD_JWST/.parsec-lowz.zbQdfU`, LagRamses origin.

## Reviewed scope and driver decisions

The [existing review](parsec_low_z_precision_fable_2026-09-10.txt) proceeds
with corrections. Under the operator's reduced audit cadence no further
routine review is requested; the driver performs end evaluation.
Reuse the converters/native binder/fixtures, no separate precheck framework,
new runtime control, passive carrier, MPI format or capacity increase.
Makefile VPATH is unchanged. Both namelist generators need no new option:
select the already supported yield/history/SED file paths together.

Explicit `--metallicity-grid precision_eleven --wind-km-s 1000` selects
`parsec2025_w17_hw02_precision_rate_v1`. The fixed speed is a declared
thermalized-wind comparison, not a measured velocity. The phase-escape
model is rejected on this grid rather than extrapolating its metallicity
law over eight decades. Existing `solar_pair`/`metal_rich_five` models and
defaults are preserved.

Actual nodes: Z=1e-11,.0001,.002,.004,.006,.008,.01,.014,.017,.02,.03;
45 masses from14 to600Msun per branch,495 nodes. Source tracks above600
are outside the selected IMF, whose full denominator remains .08--600.
No rotation, no <14Msun radiation/returns closure. Z=0 is outside the grid,
not silently a primordial floor.

Entire author branches omitted (no single-mass hole or energy clamp):

- Z=1e-6,24Msun: exact printed feasible interval is empty by
  `3.197519368631e-8Msun`; reproduced by the converter helper/test.
- Z=.001,150Msun: PPISN He core64.586 exceeds the existing W17/62--64
  bridge; the energy interpolation continues to reject it.

Queries at those Z use the selected neighboring branches' linear-Z
same-age cumulative mixture, NOT the omitted author's predictions. The
.001 bracket is .0001--.002, a factor20 range and an explicit scientific
limitation. The actual integration test is placed at birthZ=.001.

## Baryons and wind shape

Decimal80 feasibility uses each original Mfin/Mbar string's half-last-digit
interval. Preserve gross isotope/element values and fates; prefer unchanged
Mfin, then the closest feasible remnant. Failed/direct collapse requires
Mfin=Mbar; PISN requires exact zero remnant. Select binary representatives
with outward ULP steps only inside the original intervals. An infeasible
Decimal node never reaches that arithmetic correction. This is a bounded
representative, NOT recovery of hidden author values or blanket41-isotope
normalization. Per-node original strings, bounds and selected values are
in the generated manifest.75 selected nodes change a printed baryonic value;
all old five-Z node return/remnant/explosion-energy/age/fate budgets remain
unchanged. Wind energy uses the new model's explicitly selected fixed speed.

All nodes of this NEW model use the positive RATE column28 trapezoid as
wind-loss shape; normalize it to the selected mass budget. Nine primordial
14--30Msun tracks have unresolved printed MASS differences but positive
RATE. Reject negative/accreting RATE or unsupported positive endpoints.
Keep the existing positive interval/element balancing and1e-4 cumulative
compression. Existing models retain their mass-difference shapes/bytes.

Actual surface H/He strings can each have only5e-12 half-last-digit bounds,
while their projected composition can exceed unity by9.821e-12. Contrary
to the preliminary single-largest-fraction proposal, BOTH printed H and He
intervals may be required. Remove excess in descending major abundance
order without crossing either printed bound; preserve all trace species.
This affects the relative timing shape only, not published gross endpoints.
No global metallicity normalization, trace floor or yield rescaling.

## Matched radiation

Use the same495 M/Z/fate/death nodes and exact model identity. Author ZIP,
track and photon primordial prefixes differ (1E-11/1e-11/1D-11); explicit
source naming handles these, without guessing source coordinates.
All consumed archives are pinned in both converter manifests.

Two measured photon tables overhang material death:270Msun/Z=.004 by
31.4861yr and70Msun/Z=.006 by27.3327yr. Restrict the existing positive
piecewise-linear Q/E rate integral to the measured material death,
interpolating between the two bracketing rate vectors. This is a convex
operation preserving group-energy and bolometric constraints; no logQ/Teff
extrapolation or post-death photons. Record discarded overlap duration.
Shorter tables retain the existing distinct full-track Planck tail.

Conversion produced100544 material rows and161392 radiation knots, below
the existing500000 native cap. Maximum bolometric quadrature discrepancy
1.8344902269312498e-14; maximum constrained energy/Lbol .9792262207477095.
These are internal closure checks, NOT atmosphere accuracy estimates.

## Retained evidence

`input/`, `sed/`, converter logs and `source-check.log` under the work root.
The focused `simulation/snrt/tests/parsec_printed_precision.py` checks495
source endpoints and actual printed bounds, unchanged old-node budgets,
and both known excluded-node failures. No raw simulation output was made
by conversion. The first matched Intel CPU binary, before the arithmetic
repair below, had SHA256
`8a9665010a81e910104c0b6d4426620e116cbe1e6e0217fe9436284921d3db40`.

Initial converter attempts correctly refused a surface adjustment larger
than a SINGLE H bound. Those attempts published no input package; final
conversion uses the two separately bounded major fractions described above.

## Actual integration failure and bounded arithmetic repair

The first `live/` invocation ended after354.432s with transport code103,
despite process exit0. It is a failed run, not a pass. It created no raw
dump. The low-Z input exposed an existing numerical incompatibility:
Intel Fortran's optimized MAIN enabled FTZ/DAZ, flushing subnormal FP32
photon stencil increments/survivors while retaining FP64 energy. Resulting
group7 means fell below500eV; no source spectrum or conservation gate was
relaxed to accept them.

Extended the existing C++ backend smoke with normal FP32 photons near
5e-38 and a500.5eV admissible mean. Under `-ftz` it independently fails
with a transported mean464.7500014313eV (`backend-ftz.log`). The SAME
kernel/test under `-no-ftz` passes every existing HHe/D03/secondary/ledger
test (`backend-gradual.log`), including maximum relative ledger errors
N=1.077921896186e-8 and E=1.554312234475e-15. No packet floor, tolerance
increase, spectral clipping or extra carrier was added.

For SNRT builds only, Makefile adds `-no-ftz` to the Fortran flags and
makes cached `ramses.o` depend on the Makefile so MAIN is actually rebuilt.
The selected MAIN remains `patch/cuRamses/ramses.f90`; VPATH is unchanged.
The standalone backend smoke MAIN also uses `-no-ftz`. Per
[Intel's FTZ documentation](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/ftz-qftz.html),
this setting belongs to MAIN; precise C++ kernel compilation alone does not
override the parent Fortran floating-point environment. This enables gradual
underflow, not arbitrary dynamic range beyond FP32 representability.

`live-gradual/physical.nml` is byte-identical to the failed input;
`restart-gradual/physical.nml` changes only nrestart0->1. New binary
`ramses_parsec_lowz_gradual3d` SHA256
`f28ea597da805555e47261c8d8325aa19c4ecd5b696b570f8499f140e0b84fd3`.
Fresh four steps now PASS in389.999s total wall time (reported evolution
21.795s). Input validation dominates startup, not an external review wait.
Both step2/4 dumps are55104576bytes; density/thermal fields are finite and
positive, minimum final density .0009995136440064322 and thermal energy
density6.17490009929032e-7. Gas plus stellar mass .0009999999999999998
(initial .001), actual stellar birthZ=.001 throughout. CR energy minimum
4.917895472830974e-8 and mean5.03660605144466e-8. No transport/MG failure.
The standalone native source identity has3553993 doubles and all495 nodes;
four IMF mixtures, below-grid rejection and mismatched actual-grid rejection
PASS. Old two-Z native behavior and regenerated default phase-yield/history
bytes also PASS. GNU underflow/denormal flags are informational, with
invalid/divide-by-zero/overflow trapped.

Restart2->4 completed in379.447s (reported evolution11.339s). All82 physical
datasets are EXACT, all numeric datasets finite, SNRT identity/attributes
and physical clocks equal. Gas+stellar mass and positive thermal/CR fields
pass; total photon-code N=3.5570265249346886, energy-code eV=23.081461144885672.
Group7 final mean506.37361985eV is inside500--2000eV. No transport, MG,
source staging or MPI failure. Driver accepts this bounded common-source
connection and arithmetic repair. Effective paths are
`/gpfs/kjhan/LRD_JWST/.parsec-lowz.zbQdfU/live-gradual/physical.nml` and
`/gpfs/kjhan/LRD_JWST/.parsec-lowz.zbQdfU/restart-gradual/physical.nml`.
MPI2/OMP2, noncosmological4^3, Intel CPU NVAR19/NENER1/SNRT/HDF5, HHe
node-FS2010 RT, star formation, gravity, wind/CCSN/pair feedback with thermal
and .1 advective-CR injection. Dust, CHIMES, ordinary cooling, AGN/sinks,
CR-SF pressure support are deliberately off. H=.7487368421052631 and
He=.25026315789473685; initial metals all in Fe are a coupling-test IC,
not a realistic abundance pattern or a calibrated galaxy.
Output noutput1/aout2/tout1e30/foutput2/fbackup1e6: two fresh dumps and one
restart dump, measured total165313728bytes; free102TB checked. Fresh and
restart use the same gradual-underflow binary. This is not a closed-box
isolated energy calibration: gravity and physical source injection are on.

All three evaluated raw HDF5 files have now been permanently deleted and
their absence verified; [cleanup manifest](parsec_low_z_raw_cleanup_2026-09-10.md).
Inputs, source archives, logs, compact results and binary identities remain.
Source output hashes: yields
`114ca33dd65d9e93ceed0a589aadbfc7a9de979b7fe5bb51ec2a9da7b880ddf3`,
history `9d9d002d9fb7f3c83d849c51f7d63161a815694f3720ebd94e1326d347317260`,
SED nodes `cb28f81219c40bea79f551863da345d703f92494251997cbb111dd31b335d286`.
