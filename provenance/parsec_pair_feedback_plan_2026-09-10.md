# PARSEC source-node P(P)ISN feedback bundle

Final objective: physically justified, simulation-ready RT/stellar and AGN
feedback/dust in lagRamses. The operator preapproved all eight medium-term
groups and continuous implementation. This bundle advances group 5 (actual
ZAMS/Z/fate-resolved P(P)ISN), reusing the native source, deposition and restart
paths. It does not claim to finish all medium-term groups.

## First review questions

1. Q-GOAL: does the concrete runtime delivery below advance that objective?
2. Q-LEAN: remove excessive instrumentation/gates or artificial subdivisions.

Then assess physical justification and feasibility. One plan review, ordinary
engineering tests during implementation, driver end evaluation. No new
operator approval wait and no per-helper audit.

## Available inputs and bounded model

Primary source: Costa et al. 2025, https://arxiv.org/abs/2501.12917,
author data https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/ .
Downloaded all_ejecta.zip SHA256
49c0f0dac42ffb9643afc82eefe28793b44556cf6bbf0ed7d6df651a5aaa9b63.
Its separate wind and total yields include M_initial/M_final/He and CO cores,
gravitational and baryonic remnant masses and explicit CCSN/FSN/PPISN/PISN/DBH
fates. Terminal yields are total minus wind, never total plus wind.

First physical package: nonrotating Z=0.008 and 0.014, 14--600 Msun.
These two full branches have nonnegative eleven-element residuals relative
to the baryonic returned mass in every wind/terminal row; other Z branches
need separate source-precision treatment, not an extrapolation or silent fix.
Gross tracked species are H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe. Other species remain
explicitly untracked gas mass. Use Mbar, NOT gravitational Mrem, for baryon
bookkeeping; escaped neutrino energy is not injected as gas mass or heat.
Published Fe is decay-inclusive: do not inject fresh Ni56 again.

Acquire author evolutionary tracks at these Z for terminal ages and wind
histories. Their atmosphere photon files supply five cumulative Q thresholds,
not the native nine Q/E groups: RT source remains separate and no common-SED
claim is made. AGB below 14 Msun is not supplied by this package.

Energetics are an explicitly named composite prescription, not a claim that
PARSEC tabulates explosion energy:

- CCSN: declared fiducial 1e51 erg/event, distinct from the yield source.
- FSN/DBH: zero *explosive* energy/material as specified by their source fate;
  wind mass/energy remains nonzero where supplied.
- PISN: interpolate HW02 public E_expl versus He core, within tabulated
  64--133.3 Msun only: https://2sn.org/DATA/HW01/bulk_yields.txt .
- PPISN: Woosley2017 Table1 summed pulse kinetic energy versus He core
  (34--64 Msun); the same core mapping family used by PARSEC for remnant
  budgets. Explicit bare-core energetic approximation when an H envelope
  remains, not a fully matched explosion calculation. Do not interpolate a
  fate label or cross into a different fate to determine energy.
- Pulse yields/energy are one unresolved terminal event at the PARSEC end
  age, not resolved lightcurves, shell collisions or pulse chronology. Expose
  this approximation in source identity and documentation.
- Wind mechanical energy uses an explicitly supplied comparison speed
  (default package value 1000 km/s), not falsely called tabulated velocity.
  Isotropic stellar-frame ejecta have zero vector momentum; their mechanical
  energy follows the existing deposition channel, without adding it twice.

Preserve defaults and existing user-selected v1--v3 histories. The package
does not claim universal rotation/Z coverage, full radioactive inventories,
phase-dependent wind velocity or common-population radiation.

## Native implementation

Extend the existing source-node history contract to v4 with an explicit fate
per node and domain maximum up to 2000 Msun (actual package stops at 600).
CCSN, failed collapse, PPISN, PISN and direct collapse retain distinct identity.
Admit only source_consistent for this version. No mixed-remnant suppression
of an enabled PISN channel. Existing models retain their old 120 Msun limit.

Wind=channel1; ordinary CCSN return=channel3; P(P)ISN return=channel5. The
existing channel3 terminal ledger owns all NS/BH baryonic remnants, including
post-PPISN remnants; channel5 never owns a second remnant. Full PISN has zero
remnant. Every node supplies explicit zero rows for its inactive terminal
channel, rather than extrapolating across source gaps. These zeros follow
known fates; they do not fill missing yields.

Use the same existing nearest-mass fractional budgets and same-age linear-Z
SSP mixtures. For v4, require matching enabled wind/SNII/PISN mass windows,
and reject a missing channel or unsupported model. Bind actual fate and
payload values to restart identity. Update both namelist generators/GUI.
PISN dust condensation remains off unless independently supported; do not
silently inherit a CCSN dust-condensation efficiency. Gas metal return still
feeds the existing dust growth path.

No new generic audit framework, carrier tensor or automatic production token.
Source conversion produces the existing canonical table plus v4 history and
a compact source/approximation manifest. Wind histories should reuse physical
track mass loss where possible; any compression/endpoint normalization must
be quantified, not conceal unsupported elemental histories.

## Completion evidence for this bundle

- Existing native history/endpoint tests remain passing.
- New fixture: distinct fate/channel assignment, no double remnant,
  same-age Z mixture, endpoint mass/element/energy conservation, timestep
  telescoping and invalid/missing fate rejection.
- Load the actual physical package through the Fortran SSP/runtime path.
- One small fresh MPI/OpenMP hydro feedback evolution and split restart with
  actual P(P)ISN source activity; compare conserved fields and exact source
  identity. Report output schedule/free space before launch. Delete only
  evaluated raw test dumps after driver end evaluation, retaining evidence.

Review should identify material physics defects without turning optional
full stellar evolution, new atmosphere synthesis or supernova lightcurves
into conditions for this explicitly bounded alternative feedback model.
