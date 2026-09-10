# Static Fe spectral connection: driver evaluation

Status: bounded native/live/restart PASS, NOT general stellar/AGN Fe
production approval or completion of medium item8. The [single Fable review
disposition](snrt_fe_spectrum_fable_disposition_2026-09-10.md) records the
CHIMES/layout conflict and retained physical limits. No routine end audit.

## Scope and native changes

Opt-in `hhe_d03_fe_maxent128_fs2010_v1` extends existing D03/FS2010 positive
128-node N/E transport to six actual mass columns: C/S four bins plus Fe two.
The old four-bin C ABI is intact; a bounded four/six-bin entry shares the
same absorption/scattering kernel. Zero Fe is bitwise equal in every output
and ledger to the old four-bin entry. Total actual absorbed energy enters
the existing shared-temperature six-bin material/IR transaction once.
No per-cell spectral tensor, new dust temperature model or CUDA claim.

The existing electric+eddy Fe comparison remains <=300 K and primary
absorption<=4 eV. The node operator rejects hard absorption in an Fe-bearing
cell before publishing any state (status10). In particular a broad group2
packet with mean2.366 eV fails because its actual nodes extend to5.6 eV;
no tail clipping or normalization. Broad stellar/AGN irradiation is therefore
NOT admitted. Missing photoelectron escape, spin response, variable-temperature
dielectric and independent component temperatures remain limitations.

Static Fe may now run without CHIMES only with cooling=none, all grain
condensation/growth/sputtering/size-transfer/SN-shock switches disabled,
no Fe kinetics and no relative motion. Fe fields31:32 follow the reserved
dust window in the NENER1/virial hydro profile (DUST_IRON1,NVAR32).
Existing CHIMES offsets188:189 and physical paths are unchanged. Both mkrun
and generic namelist/GUI selection/validation match this conditional layout;
defaults remain unchanged and no physical stellar SED is silently selected.

The first MPI launch reproduced a true missing hydro connection: the
non-CHIMES dust aggregate face flux omitted Fe, although its two independent
carriers were transported. Step2 rejected aggregate/component disagreement.
`patch/cuRamses/umuscl.kjhan.f90` now adds both Fe fluxes to the dependent
aggregate, as the existing MHD path already does. CHIMES already replaces
these fluxes with its consistent-carrier sums; its behavior is unchanged.
No tolerance widening or post-hoc mass normalization. VPATH order unchanged.

Nodal source: reuse pinned Werner09/Henke/size-Drude mu=1 Mie builder, not
new Fe dielectric assumptions. `dust_fe_band128_manifest.json` records
source/generator hashes. Compiled include SHA256
`7c5859f538419f26d206d7f6a45e62dd5552df1d9dd94aba3901f4499ec461a4`;
content identity
`36c48264e9fedef3c8441b8bae360cefc0a4b0f1da79ca33421b99bae5253dbd`.
Runtime requires the Fe and D03 node energies to match exactly; both content
hashes bind native version11 and HDF format40. Four-bin identities unchanged.
The nodal table does not override the old representative/IR optical identity.

One actual-source resolution measurement compares128 and2048 log nodes with
uniform-dE weighting. In .01--1 eV, Fe absorption differences are0.01305%
(10nm) and0.01108%(100nm); transport scattering0.12996%/0.15707%.
Across all nine groups the largest absorption difference is0.51727% in
2--10keV (100nm). This is quadrature evidence, NOT a universal near-edge
convergence claim or admission of keV Fe heating.

## Evidence retained in `.fe-spectral.MDeWoz/`

Intel/GNU backend regressions pass. Extended existing Fortran dust smoke
covers actual Fe-only/mixed opacity, zero-Fe gas/FS2010, zero dust, opaque
limit, invalid sixth column and both hard endpoint/broad-tail rollback.
Maximum photon relative closure8.8390196273380184e-9; energy2.7755575615628914e-14.
Native checkpoint tests pass new mode and old fixed mode; tampering either
D03 or Fe nodal hash rejects before mutation. Setup/GUI48 tests:47 pass,
one display-dependent skip. Python syntax and `git diff --check` pass.
The old full dust/material/IR smoke also passes. Rebuilding the original Fe
representative/IR include with its MATCHED original hot contract gives
byte-identical output (`legacy-matched.inc`). An initial regeneration with
the newer nodal run's contract differed only in its supplied temperature
axis; that was an input mismatch, not a replaced legacy table. Both logs
are retained. Final Fe smoke explicitly checks status10 for hard absorption
and status2 for invalid sixth-column mass, not just any nonzero failure.

First failed executable SHA256
`b371327d3783db2b41eb039e958a13dad5640c699f5259602c662c22e4b3441f`;
`live/run.log` retained, no raw dump was produced. Final executable
`ramses_fe_spectral_fixed3d` SHA256
`7cad07fa7044cb3d00d392ab46d4492656cb607131d033c04285b710c1fc11a8`.
Intel/ifx/OpenMP/bounds, CPU SNRT1/DUST1/HDF51/CHIMES0,NVAR32,NENER1,NVECTOR32;
the MAIN retains `-no-ftz`. Effective input `live-fixed/physical.nml` is
byte-identical to the failed input; `restart/physical.nml` changes only
nrestart0->1. Fresh8.133s/restart4.377s wall.

MPI2/OMP2, periodic4^3, noncosmo4steps; C and Fe each1e-4 initial mass
fraction (equal small/large), native shared U(20K). Existing stellar mass
return, star formation, gravity and advectiveCR remain as fixture physics.
`subev_control.nml` is explicitly SYNTHETIC: constant group1 photon rate
per initial stellar mass, not BPASS/PARSEC or a physical SED qualification.
AGN/CHIMES/ordinary cooling and grain mass reactions are off. This exercises
nonzero primary Fe absorption and shared material/IR rather than only zeros.

All108 physical datasets exactly match uninterrupted vs restarted evolution
(AMR12,coarse3,domain1,gravity8,hydro64,particles17,SNRT3); all numeric data
finite, SNRT/source/Fe/material/clock identities equal. HDF carrier width
12324 is unchanged from the existing primary+IR state. Gas+star mass.001,
residual0; Fe mass mean9.996649766493149e-8, dust aggregate closure exactly0.
Minimum rho9.99543354340077e-4, gas thermal8.355594177446613e-8,
dust energy1.1757668717229772e-18. Primary N sum.008037053594478039 and
energy sum.0008036175238828132 in code-photon-density*eV; IR energy-density
sum8.612704127918662e-16 erg/cm3. These are stored-value sums, not global
volume-integrated budgets. IR local balance max2.1878e-11. Noncosmo SFRD
NaN is the existing log-only diagnostic, not a physical field/observable.

Driver accepts only this bounded capability. No additional review gate or
new completion requirement is created. Evaluated raw snapshots are removed
under [the cleanup manifest](snrt_fe_spectrum_raw_cleanup_2026-09-10.md):
three files22900904 bytes permanently deleted, absence verified.
Inputs, logs and compact results remain.
