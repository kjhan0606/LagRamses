# Approved pilot pre-enrichment: Z=1e-10

Operator approved a trace-metal initial condition, not primordial Z=0. Applied
to a new IC view at `/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007`.
The original density, velocity and displacement files are unchanged symlink
targets. Twelve new GRAFIC passive files supply total Z and eleven elements;
no initial dust is supplied. This does not authorize cosmic PopIII claims.

Heavy-element pattern: existing pinned PARSEC source archive
`.parsec-sources.9ZIUlJ/all_ejecta.zip`, member
`ejecta/Z0.014_Y0.273_winds_ejecta.dat`, initial-composition header, NOT ejecta.
Use the same eleven-element projection as the source builder; normalize its
nine tracked C--Fe fractions to unit metal mass, then multiply by1e-10.
This is an explicit scaled-source comparison mixture: untracked elements are
not added as hidden mass. H/He retain the existing .76/.24 relative ratio,
rescaled to1-Z. The stored float32 He is the nearest residual after H/metals.

All twelve files have128^3 finite, nonnegative constant fields with48-byte
headers and matching slab markers. New passive payload100,676,256 bytes.
Stored Z=1.000000013351432e-10, sum(C..Fe)=9.999999954620961e-11;
relative metal closure error1.78894e-8 is float32 storage precision, not a
runtime floor. After residual-He rounding, total element mass fraction is
1.0000000001 (absolute error1e-10). Do not modify evolution abundances or
change scientific acceptance tolerances to enforce an artificial Z floor.
The unused first independent-rounding IC view is retained as
`ics-before-closure`, not used by the active namelist fragment.

Read actual native yield rows and SED node records: both bracket the requested
Z between1e-11 and1e-4. No extrapolation. These are still the existing declared
mixed-source/late-SED-truncated effective model, not a complete stellar library.
Age and later enriched-Z coverage are separate from this initial-Z check.

Reproducible preparation script and per-file SHA256/absolute input identities:
`.galaxy-preenriched.IJ8aeX/prepare.py` and `manifest.json`.
Configuration fragment: `simulation/snrt/config/galaxy_pilot_initial_composition_v1.nml`.
Existing z_ave uses a legacy .02 reference, so Z=1e-10 corresponds to z_ave=5e-9.
Explicit pvar files bypass the old missing-file metallicity fallback. No global
default, public namelist schema or runtime abundance-clamping policy changed.
The existing mkrun frontend imports the same namelist generator definitions;
this run-specific fragment does not add a frontend option.

Required layout: nhydro5,NENER0,metal enabled, no delayed cooling, sf_virial,
aton or SGS: imetal6,ichem7. The fragment is not a standalone executable setup.
Native reader confirmation passed in jobs539024 and539031: all12 files read,
clean first coarse step, no NaN counters or raw output. Batch postprocessing
initially failed because rg was absent on the compute node; standard-awk
re-evaluation of the preserved logs passed. See the coupling progress record
for executable/input hashes and CPU timings. Full coupled evolution is still
pending. Cosmological dust guard is unchanged; approval of this Z does NOT
approve omission of radiation redshift or use of a fixed dust bath at high z.
