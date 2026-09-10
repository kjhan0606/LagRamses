# PARSEC common high-mass radiation/feedback — bounded PASS

Driver evaluation: actual native implementation, physical input conversion,
Intel/GNU tests and MPI2/OMP2 continuous/exact-restart execution PASS.
This closes the **14--600Msun common-population connection**, NOT the
full SSP/low-mass population, exact NLTE atmosphere or all medium tasks.
Default BPASS/reference choices are unchanged. No approval wait, extra
end audit or new RAMSES namelist key. Workdir `/gpfs/kjhan/LRD_JWST`,
origin `kjhan0606/LagRamses`.

[Plan](parsec_common_radiation_plan_2026-09-10.md),
[Fable disposition](parsec_common_radiation_fable_2026-09-10.md),
[physical source inspection](parsec_common_radiation_sources_2026-09-10.md).

## Physical model and native implementation

Later same-day correction: the [phase-wind bundle](parsec_phase_wind_implementation_2026-09-10.md)
identified PARSEC's luminosity unit3.846e33 versus the initial converter's
nominal3.828e33. New inputs use the source unit; the numerical measurements
below describe the retained original input, not that correction. Native
algorithm/population wiring is unchanged. New corrected inputs pass the
subsequent native and MPI2 exact-restart tests.

The same90 nonrotating PARSEC mass/Z nodes as the v4 feedback package,
Z=.008/.014, use the SAME full .08--600 initial-mass normalization, native
IMF mass-fraction function and nearest-source cells INCLUDING the40Msun
seam. A public `build_source_mass_edges` helper was extracted from the
existing feedback integrator and reused; physics and VPATH order unchanged.
Per-node radiation weight is IMF mass fraction/source-node mass, not an
independently number-normalized or covered-subset-renormalized IMF.

`build_parsec_native_sed.py` is an offline physical input converter, not
the simulation implementation. Published instantaneous QWerner/QHI/QHeI/
QOII/QHeII differences constrain disjoint nominal11.2/13.6/24.6/35.12/
54.4eV intervals. Positive Planck photon shapes at actual tabulated Teff
supply the within-interval spectrum and energy moments. Remaining Lbol
below11.2eV uses a Planck prior. Actual nine SNRT edges are integrated
without nearest-group reassignment. IR/X tails remain escaped, not heat.
Finite/infinite integrals use64-point Legendre/Laguerre quadrature.
No clipping of Q, negative residuals or lower-edge-energy fallback occurred.

Pre-first-photon ages12--236yr hold the first source state. All missing
late photon intervals use the actual full track L/Teff Planck spectrum;
this is a named **unmeasured blackbody tail**, NOT measured atmospheric Q.
The two important Z.008 60/70Msun cases are recorded separately. Each side
of the join is integrated independently, permitting a physical-model rate
jump rather than adding an invented smoothing interval. Radiation ends at
the exact feedback terminal time. No post-death emission. The source is
HIGH-MASS ONLY and supplies no below14Msun radiation. In particular it is
not a qualified full dust-heating SSP after15Myr.

The existing cumulative knot-compression rule bounds error by1e-4 of EACH
positive final component.90nodes retain30444 knots. The first candidate
`input/` was REJECTED: subtraction of two large cumulative endpoints lost
precision in tiny late hard-photon moments. Accepted `input-increments/`
stores positive interval integrals summed from original rate intervals.
Native integration multiplies these by overlap fractions, avoiding endpoint
cancellation. This is mathematically the same piecewise-linear cumulative
interpolation, with no arbitrary spectral cutoff or photon deletion.

`snrt_parsec_source.f90` reads/validates this immutable source, binds prepared
feedback M/Z/age/fate/model coordinates and configured source windows/IMF,
integrates age intervals and same-age linear-Z mixtures. It publishes nothing
before successful binding. `snrt_stellar_source.f90` admits explicit source
v4 through existing `SNRT_STELLAR_SED`, exclusively with other SED versions.
`read_params` invokes `phase0_bind_radiation`/idempotent source preparation
before startup consensus, not at the first radiation transaction. All670707
consumed SED identity values, including interval moments and native weights,
use existing MPI/HDF5 identity. Separate feedback identity retains actual
yield rows. This is value binding, not cryptographic authentication of a
user-supplied model label. No per-cell source-history tensor was added.

Native injection already accepts Q/E and transports its actual mean energy;
no nominal-E replacement or new driver source path. Old v1/v2/v3 identities
and public calling signatures remain unchanged. No RAMSES NML change was
made, so neither NML generator needs a new option for this external asset.

## Evidence

Task root `.parsec-sed.ohf9MX/`. Accepted payload:

- `input-increments/nodes.dat` SHA256
  `979ec5c9e1331cdf844155aa65b907deac6137a23524bbe8f0b13516ecb65f68`.
- `input-increments/source.nml` SHA256
  `caade3ce5f3ce742a846e3b6ee3d900b0b495e24be9816c1b65e77a4d1f09488`.
- `input-increments/manifest.json` records every source hash, tail duration/
  lifetime/bolometric contribution, moment units and compression error.

Full-source64-point bolometric closure max1.85070e-14. Maximum Q-constrained
high-energy/Lbol=.98383073.64/128-point samples of the two problematic
Z.008 tails' parent tracks and Z.014600Msun differ by at most5.15e-13 for
moments above1e-14 of the row peak. This is quadrature evidence, not an
atmosphere-model accuracy claim.

Intel bounds checks and GNU bounds/FPE checks pass: all four IMF choices;
time telescoping; exact metallicity knots/interior mixtures; energy support;
no pre-birth or post-death photons; omission of required E rejected; five
mutated actual feedback bindings reject without subsequent publication.
Identity is deterministic after rebuilding native weights. Independent
analytic Salpeter mass-cell integration agrees with native Q/E to3e-14.
For an initial1e4Msun SSP atZ.011 over0--30Myr, Chabrier yields
Q=1.148692649103698e66photons, E=7.038122177065079e66eV (all nine groups).
These are photons produced over a stellar lifetime, not instantaneous density.
Existing BPASS tests pass: legacy/v2/v3, interval/Z/QE parity,9 old admission
rejections and8 energy rejections. `git diff --check` passes.
Six additional actual v4 file rejections pass (IMF, policy, energy semantics,
legacy-array shape, invalid SHA text, malformed node payload).

Planck-prior sanity comparison atZ.011,0.2Myr intervals around1/3/6/10Myr:
group4--7 mean energies (eV), PARSEC versus independent BPASS:

| age/Myr | PARSEC | BPASS |
| --- | --- | --- |
| 1 | 12.3284,17.4980,29.4627,59.1221 | 12.3051,17.8032,29.9594,57.4129 |
| 3 | 12.3053,17.5768,34.2363,67.4021 | 12.2868,17.5045,30.2363,57.4351 |
| 6 | 12.2914,16.7923,27.5931,57.5128 | 12.3077,17.2992,30.1910,63.9715 |
| 10 | 12.2367,16.4264,27.1570,57.0503 | 12.2415,17.7508,32.5548,63.1738 |

Different stellar/binary populations and atmospheric closures prevent
interpreting these differences as measured model errors.

## Live and restart

Checked binary `.parsec-sed.ohf9MX/ramses_parsec_sed3d`, SHA256
`647d767bb83c7eb0e4c94f93166bf545978d5308c567f93cb18110c05b2360ee`.
Intel MPI/ifx/OpenMP, NVAR19/NENER1/NVECTOR32, SNRT1,DUST_LIVE0,HDF51,
USE_FFTW0,FDMDEBUG1. Subsequent source-only cleanup removes unused helper
variables, and a test-only rejection branch does not change binary physics.
Inputs `live/physical.nml`, `restart/physical.nml`, `environment.sh` retained.

Periodic4^3 noncosmological,4steps, MPI2/OMP2, gravity/star formation,
actual PARSEC v4 wind/CCSN/P(P)ISN thermal+CR.1feedback and matched radiation,
`hhe_maxent64_fs2010_v1`. No dust/CHIMES/AGN/sink/ordinary cooling, CR
advective, CR SF support off deliberately. Legacy RAMSES rt flag remains
false; SNRT is activated by its explicit environment switch. Output schedule
noutput1,aout2,tout1e30,foutput2,fbackup1000000.2fresh+1restart dumps.

All82 physical datasets EXACTLY equal: AMR12/coarse3/domain1/gravity8/
hydro38/particles17/SNRT3. SNRT attributes and physical clocks also agree.
HDFformat23, cellwidth1444, stellar-source identity670707values. All numeric
datasets finite. Gas+star mass=.001 exactly, minimumrho=.000999542141984,
minimum internal-energy density5.911150322368409e-7, CRmean4.557530053634920e-8.
Primary code-number sum3.523442603930908, code-eV energy sum20.58016176559802;
not volume-integrated source energies. All transported group means in support.
No MG nonconvergence/Fortran/MPI abort. Noncosmological SFRD NaN is the
pre-existing diagnostic, not a physical field. This is not an isolated
global-energy box: star formation, feedback and gravity are active.
Wall42.846s fresh,34.238s restart. Compact `evaluate.py`/`evaluation.txt`
and logs retained. Evaluated raw cleanup is recorded separately.

Remaining: below14Msun common population/full SSP, original atmosphere grid
instead of the Planck prior/tails, broader Z/rotation/terminal timing, and
other original medium groups. Do not mark all medium tasks complete.
