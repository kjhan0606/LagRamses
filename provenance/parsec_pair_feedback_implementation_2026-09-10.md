# PARSEC P(P)ISN native feedback comparison

Status: native physical-source tests and MPI2/OMP2 thermal+CR evolution/restart
PASS. The separate dust/CHIMES/dynamics configuration build also passes.
This is a bounded alternative, not completion of all medium-term tasks or a
production-approval token. [Plan](parsec_pair_feedback_plan_2026-09-10.md),
[Fable disposition](parsec_pair_feedback_fable_2026-09-10.md),
[physical sources](parsec_sources_2026-09-10.md).

## Delivered runtime behavior

History v4 adds five explicit source fates: CCSN, failed SN, PPISN, PISN,
direct BH. Wind is channel1; CCSN material and all baryonic NS/BH remnants
are owned by channel3; P(P)ISN material/energy is channel5, which owns no
second remnant. Its release is a step at the declared source age, never the
old generic Cartesian age ramp. Per-node wind+both terminal channels+remnant
must close. Unknown fates, inconsistent channel payloads, missing channel5,
mixed endpoint presets, windows outside the IMF and unrepresented nodes fail
before publishing the prepared history.

v1--v3 preserve their old <=120 Msun domain and reject channel5. v4 is limited
to <=600 Msun, source_consistent, and matching wind/SNII/PISN windows. It
uses the existing mass-cell fractional budgets and same-age linear-Z SSP
mixture without extrapolation. Restart adds a version/fate block to the
existing identity containing every consumed table value. No new persistent
hydro carrier, schema framework or namelist key. Both namelist generators
and shared GUI are updated; defaults remain unchanged.

Channel3+5 now share the declared delayed-cooling shock tracer, CR SN-energy
fraction, and dust shock-energy path. These are shared subgrid prescriptions,
not independently calibrated PISN efficiencies. Direct pair-event graphite,
silicate, Fe and PAH condensation is excluded, including elemental donor
allocation for moving dust. Pair metals remain gas and can subsequently
participate in ambient growth. SNIa remains separately coupled.

## Physical package and approximations

Offline converter: `simulation/snrt/tools/build_parsec_pair_feedback.py`.
Runtime consumes only its ordinary native ASCII table and v4 history, not
Python. Package `.parsec-pair.dnTIS0/input-balanced/` has 90 source nodes,
24,884 rows, nonrotating PARSEC Z=.008/.014, 14--600 Msun. Fate counts:
12 CCSN, 23 FSN, 24 PPISN, 24 PISN, 7 DBH. Default comparison IMF Chabrier,
individual-star domain .08--600 Msun; changing its upper mass renormalizes
all channels. No AGB <14 Msun or microscopic binary-Ia claim.

Terminal material is published total minus wind; gas return uses baryonic
Mbar, not gravitational remnant mass. Untracked elements remain untracked
gas; no normalization of the eleven elements to one. Fe is decay-inclusive,
not a fresh radioactive inventory. No neutrino energy is injected.

Wind phase shapes use track surface composition and lost mass. Heavy surface
Si/S/Ca/Fe follow the initial mixture under the author's prescription. A
positive row/column balance preserves BOTH each interval's lost mass and
the published final elemental masses, including the untracked gas column.
The initial column-only endpoint scaling failed the native monotone
untracked-mass invariant: it was rejected, not admitted or hidden by a
tolerance relaxation. Original `.parsec-pair.dnTIS0/input/` is that failed
candidate; only `input-balanced/` is used. Balancing took <=51 iterations.
Maximum original integrated endpoint difference is 0.1493% of initial mass.
Linear time-knot compression is bounded by 1e-4 of each final component
(measured maximum 9.999989633225664e-5). This is a declared calibrated phase
shape, not a recovered complete isotope network or phase-dependent speed.

Wind energy uses the explicitly supplied comparison speed1000 km/s;
CCSN energy uses the declared 1e51 erg/event. PISN energy uses HW02 versus
He core. PPISN uses Woosley2017 Table1 with an explicit 62--64 energetic
bridge to the full-disruption endpoint for two source nodes (He62.601 and
63.7932). Fate/material/remnant still come from PARSEC, not interpolation of
fate labels. Pulses occur as one unresolved event at PARSEC track end (pair
instability onset); neither pulse chronology nor final-collapse delay is
resolved. Selected pair nodes have Mfin=MHe to printed precision. These
approximations are all recorded in the source manifest/identity.

Photon files are not used: five cumulative Q thresholds cannot supply nine
native Q/E groups. This package does not resolve the common-population SED
problem. Other Z/rotation populations, source precision in other branches,
lower-mass gaps, full decay inventories and microscopic Ia remain open.

## Evidence

`parsec_pair_history_test.f90` is one native fixture, using the actual package.
Intel bounds-checked and GNU bounds/FPE-checked physical tests pass all five
fates, terminal timing, same-age Z mixture, material/energy telescoping, IMF
64/128-cell agreement, invalid admission, no Z extrapolation and terminal
plateau. The final Intel fixture also tests shared CR shock partition and
zero direct pair-event Fe/PAH donors. Existing v1 three-preset history test
passes. GUI suite: 47 tests, 46 pass and one environment-dependent skip.

For an initial1e4 Msun SSP at Z=.011, after20 Myr:

- returned1537.98680082280 Msun; baryonic remnants670.823145723004 Msun;
- wind/CCSN/pair returned994.116338113309 /428.955305545197 /
  114.915157164293 Msun;
- wind/CCSN/pair energy9.883852574240858e51 /2.861600964241552e52 /
  3.258050887085736e52 erg.

Fresh checked CPU binary `.parsec-pair.dnTIS0/ramses_parsec_pair3d`,
SHA256 `a4d9706a398a06cf77eab988bf00f1064876cac313f6379732e9b1bf85ce8d79`,
Intel MPI/ifx OpenMP, SNRT1/DUST_LIVE0/HDF51, NVAR19/NENER1/NVECTOR32.
VPATH unchanged. Effective `live/physical.nml` / `restart/physical.nml` in
that directory: periodic noncosmo4^3, MPI2/OMP2, four steps, actual star
formation and v4 return plus advective CR (SN fraction.1), no active
RT/AGN/dust/CHIMES/ordinary cooling. These inactive modules are intentional.

128 stellar particles, all older than the maximum table lifetime at the
end. Total gas+star mass exactly.001 code units, residual0. Density minimum
9.990902985782017e-4; thermal-energy minimum3.2816297425809743e-7;
volume-mean CR energy6.113456149466493e-8. Returned stellar mass
1.3344027616663722e-7 code units. Source identity contains797,773 values.
All79 datasets match exactly across continuous/restart: AMR12, coarse3,
domain1, gravity8, hydro38, particles17. Physical numeric datasets finite.
Noncosmological SFRD log remains NaN (pre-existing diagnostic, not a field
or usable observable). The gas-only mcons log is not total gas+star mass;
the direct total is the mass check above.

Wall times30.941s fresh /30.337s restart; actual evolution3.65s /2.93s.
Initialization dominates (existing complete table audit); feedback accounts
for ~3.47s /2.80s of evolution. No new performance framework was added.
Inputs, logs, binary, `evaluation.txt`, evaluator and physical manifest are
retained. Each raw dump is6,456,376 bytes; exact evaluated cleanup is recorded
in [the cleanup manifest](parsec_pair_feedback_raw_cleanup_2026-09-10.md).

The independent conditional build `.parsec-dust-build.RcTyuJ/ramses_parsec_dust3d`
passes with SNRT1/DUST_LIVE1/DUST_DYNAMICS1/CHIMES1/HDF51,
NVAR4096/NENER1/NVECTOR32, Intel MPI/ifx, OpenMP and bounds checks.
SHA256 `fb2588fde72f7a27221c265cee7797a95c5c9c1c99d50ecc69a3f43d619b8b69`.
This binary was not evolved: it establishes conditional compilation, not
a full pair-feedback+dust live validation. Driver end evaluation: PASS for
the declared source package, native coupling and thermal+CR restart bundle;
the physical approximations and remaining medium-term scope above remain.
