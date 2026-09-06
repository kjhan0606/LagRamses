# DUST-2 grain thermal balance and IR candidate implementation evidence — 2026-09-06

## Scope and status

- Workspace: `/gpfs/kjhan/LRD_JWST`
- Bundle: DUST-2 grain thermal balance and one-pass IR re-emission
- Status: `conditional_candidate`
- Pre-bundle plan audit: Fable `CONDITIONAL APPROVE`; amendments F1–F10
  incorporated in the plan before implementation
- Final-project boundary: static JAX/SNRT/P5 evidence only; no claim of
  astrophysical dust-mixture approval, recursive IR transport, native Fortran
  parity, live RAMSES coupling, MPI/restart qualification, or production run

## Implemented components

1. `snrt_core/dust.py` now validates `snrt_dust_thermal_v1`, including exact
   group edges, strict temperature/power grids, temperature-dependent IR
   energy fractions, temperature-dependent emission-weighted photon energies,
   explicit out-of-band energy, source-table SHA-256, dust mass/H, code
   manifest, and payload hash.
2. The JAX thermal operator adds the same-table CMB power, solves the local
   single-temperature balance with a fixed 32-step log-temperature bisection,
   returns a JIT-safe out-of-range mask, and emits zero new IR with a 0 K
   sentinel for zero local dust absorption.  It uses the thermal photon
   energies rather than source-ledger mean energies.
3. The operator is called inside every existing P5 thermochemical subcycle.
   Tracked IR energy, explicit untracked energy, photon counts, and power
   residuals are accumulated with the same timestep weights as DUST-1
   heating.  The IR source is recorded one-pass and is not recursively
   transported.
4. P5 exposes `--dust-thermal-metadata` only with a provenance-pinned v3
   opacity sidecar whose source-table hash and dust mass/H match.  It records
   grain temperature, tracked/untracked IR energy, photon rates/counts,
   closure residuals, the background temperature, and one-pass semantics.
5. `tools/build_draine_dust_thermal.py` derives the candidate power/fraction/
   photon-energy table explicitly as `4*pi*C_abs(E)*B_E(T)` over the raw Draine
   energy range.  The generated candidate is stored under the existing ignored
   `external/draine_wd01_rv31/` asset area.

## Physical asset provenance

- Draine raw source SHA-256:
  `b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c`
- Current DUST-1 v3 opacity sidecar SHA-256:
  `7521ef988a47b590f375f49cdedf375109f5ee306968e54749b38f5e43a1faa8`
- Current DUST-2 thermal sidecar SHA-256:
  `e0065f2b6de47b43f1b2739ff7fedb4ba90f074d17f775955672e00f42da5259`
- The DUST-1 opacity sidecar was regenerated after the DUST-2 code-manifest
  change; its numeric opacity/scattering arrays, mixture, and mass/H are
  unchanged.  Superseded sidecars remain as explicit `.pre_dust2*` backups
  under the ignored external asset directory.
- Thermal table: 64 temperatures from 5 to 300 K; IR group `[0.01, 1.0]` eV;
  the out-of-band fraction is explicit and large at cold temperatures rather
  than being silently renormalized.
- The physical asset remains a candidate. It assumes one equilibrium
  temperature for the reference mixture and omits stochastic PAH/small-grain
  heating, dust-gas exchange, obscuration, and IR self-absorption.

## Focused verification

All commands ran in `/gpfs/kjhan/LRD_JWST` with the project CPU virtual
environment and passed:

```text
DRAINE_DUST_OPACITY_TEST_OK rows=812 groups=9 scattering_groups=9 max_consistency=8.078e-04
DUST_THERMAL_TEST_OK temperatures=64 ir_groups=1 source_sha256=b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
DUST_OPACITY_TEST_OK groups=1 weighted_energy_ev=7
P4_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata scattering=isotropic
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
P4_INGESTION_OK format=v3 shape=4x3x2 sources=2 groups=2
REFINE_STATIC_INPUT_TEST_OK factor=2 source_luminosity_conserved=true
```

After the Opus bundle-end audit, the focused tests were extended and passed
again.  They now include a genuine two-IR-group photon-number/energy check,
wrong source-table hash and wrong dust-mass binding negatives, a non-v3
thermal-attachment rejection, and a zero-dust thermal-flag P5 control.  The
thermal-on/off P5 comparison now asserts bitwise equality for the DUST-1
absorption, scattering, and heating ledgers in addition to H II and gas
internal energy.  The runner also records the declared
`dust_ir_energy_closure_tolerance=1e-5`.

The regenerated external thermal sidecar was independently loaded against the
current physical v3 opacity sidecar and passed its edge, source-table, mass/H,
code-manifest, and payload checks:

```text
DUST_THERMAL_EXTERNAL_METADATA_OK sha256=e0065f2b6de47b43f1b2739ff7fedb4ba90f074d17f775955672e00f42da5259 groups=9 ir_groups=1
```

The sidecar algorithm block now describes the actual linear-in-power,
log-temperature interpolation and fixed bisection.  It also records the
candidate's conservative CMB excess spectral split: total excess power is
distributed with the emitting-temperature fractions, not with the exact
group-by-group differential subtraction.  This remains a pre-transport
approximation.

The P5 thermal integration additionally ran the same physical v3 input with
the thermal flag disabled and enabled.  The final H II state and gas internal
energy were bitwise unchanged; the enabled output reported
`dust_ir_transport_semantics=recorded_not_transport_reemitted`, zero thermal
out-of-range hits, and IR energy/power residuals below `1e-5`.

## Remaining promotion gates

The bundle does not approve the Draine mixture as the final source-specific
dust model.  Remaining work includes physical review of the absorption/
depletion normalization, stochastic heating/PAH treatment, IR self-absorption
and recursive transport, dust-gas exchange, source obscuration, native/live
coupling, and matched source/geometry convergence.  These are later bundles;
no synthetic thermal table is a production substitute.
