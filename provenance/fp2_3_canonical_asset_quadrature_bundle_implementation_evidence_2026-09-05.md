# F-P2.3 canonical asset synchronization and source-closure quadrature — implementation evidence — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Bundle: F-P2.3, approved after the F-P2.2 Claude Opus 5 conditional-pass
audit. Work was performed on `/gpfs`; no live RAMSES job was launched.

## Delivered scope

Q1 synchronizes all ten entries in
`simulation/snrt/data/agn_nine_group_external_assets.json`, including the
explicit SED, ledger, metadata, static input, and transport control. The
manifest now records current size/SHA-256 pairs, the explicit static metadata
digest is repaired, and reproduction commands include the required
`--time-averaged-absorption-iterations 32` setting for the passing transport
controls.

Q2 makes the pilot and explicit canonical artifact tests symmetric. Both
recompute the current live-file hashes, current HEAD, and the same
`simulation/snrt` porcelain-status digest. The attestation is evaluated from
the repository root, so it does not silently query the nonexistent nested path
`simulation/snrt/simulation/snrt`. The validator records the dirty development
state as provenance; it does not call that recorded state an independent
physics pass.

Q3 makes the explicit source boundary deterministic and auditable. The
tabulated `energy_fraction_per_ev` is interpolated piecewise-linearly and
converted to photon number as
`q_E = f_E/(E_eV * EV_ERG)`. Photon moments, Verner H I/He I/He II closure,
and source-weighted dust closure use logarithmic refined union grids and a
2048-versus-4096 convergence guard at `5e-6`. The interpolation and
quadrature settings are included in the path-free SED identity.

Q4 derives explicit support and rate expectations independently of the
serialized ledger from the SED,
validates the source status/identity contract, cross-checks duplicated dust
energy-fraction ledgers, and labels the canonical explicit SED as a synthetic
non-physical wiring fixture. It does not approve that SED for astrophysical
AGN use.

## Evidence executed

From `/gpfs/kjhan/LRD_JWST/simulation/snrt` using the project `.venv` and CPU
JAX:

```text
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
P0_SED_CLOSURE_OK synthetic_groups=4 metadata_groups=9 species=3 subthreshold_opacity=zero
AGN_PHOTON_LEDGER_TEST_OK p0_groups=9 legacy_groups=5 hard_xray_group=positive subthreshold_opacity=zero edges=exact
STELLAR_PHOTON_LEDGER_TEST_OK sources=3 groups=9
DRAINE_DUST_OPACITY_TEST_OK rows=812 groups=9 max_consistency=8.078e-04
DUST_OPACITY_TEST_OK groups=1 weighted_energy_ev=7
P4_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
AGN_NINE_GROUP_PASS mode=pilot hard_q=3.52996e+51 hard_to_soft_q=0.219425 hard_supported_sed_fraction=0.0437846 hard_bolometric_fraction=0.00794605
AGN_NINE_GROUP_PASS mode=explicit hard_q=2.87679e+53 hard_to_soft_q=1.16096 hard_supported_sed_fraction=0.800001 hard_bolometric_fraction=0.800001
AGN_NINE_GROUP_ARTIFACT_OK hard_q=3.52996e+51 hard_to_soft_q=0.219425 hard_supported_sed_fraction=0.0437846 hard_bolometric_fraction=0.00794605
AGN_NINE_GROUP_EXPLICIT_ARTIFACT_OK mode=explicit criteria=all_true
```

The regenerated explicit ledger reports source-moment convergence maximum
relative error `6.32e-7` and Verner-closure convergence maximum relative error
`1.57e-6`, both below `5e-6`. These are base-versus-refined convergence
differences, not formal Richardson error bounds. The passing pilot and explicit short transport
controls use S4, float64, `0.001 Myr`, and 32 time-averaged absorption
iterations. The full-CFL one-step failure probe remains preserved and is still
classified as an expected gate failure.

The changed Python modules and focused tests compile with `py_compile`, and
`git diff --check` is clean. Python remains an offline reference/oracle and
provenance boundary; no claim is made that these helpers are the final Fortran
RAMSES runtime implementation.

## Remediation after the first Opus disposition

The first bundled Opus audit is recorded at
`provenance/claude_opus5_fp2_3_canonical_asset_quadrature_bundle_end_audit_2026-09-05.md`
with `CONDITIONAL PASS`. Its two blocking local findings were addressed in the
same approved work package: `AGN_NINE_GROUP_VALIDATION.md` now carries current
pilot/explicit hashes, reproduction commands, and the engineering-only SED
label; and the explicit validator now independently compares serialized group
means and energy fractions with the SED reconstruction. The empty-group
support convention, failed-probe reproduction record, pilot fresh-validator
rerun, and convergence-difference wording were also corrected. Canonical
pilot/explicit validators and both artifact tests were rerun and passed after
these changes. The follow-up Opus audit is recorded at
`provenance/claude_opus5_fp2_3_canonical_asset_quadrature_bundle_followup_audit_2026-09-05.md`
and returned `PASS` for the F-P2.3 scope. Its non-blocking observations are
carried forward to the physical-source/publication bundles.

## Boundary

This bundle closes the F-P2.2 local integrity findings and source-closure
quadrature control. It does not select a physical AGN or stellar SED, solve
the `[40,120] M_sun` yield gap, admit mixed STAR+AGN dust, add scattering/IR
re-emission or grain evolution, or connect the closure to the live Fortran
RAMSES runtime. Those remain later high-level RT/feedback/dust work.

The required end-of-bundle review is the single bundled Claude Opus 5
read-only audit. Its verdict is intentionally not pre-filled here.
