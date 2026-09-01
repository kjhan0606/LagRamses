# Furlanetto--Stoever secondary-ionization validation

Date: 2026-09-02
Stage status: internal gates PASS; Claude Opus 5 final re-audit PASS

## Physical contract

`snrt_core/secondary.py` uses the Furlanetto & Stoever (2010) Monte-Carlo
fast-electron tables distributed with 21cmFAST. The vendored dataset contains
258 electron energies from 10 to 9937.21 eV at 14 H II fractions from `1e-4`
to `0.999`. Production uses bilinear interpolation in electron energy and H II
fraction, matching the public author implementation. Below 10 eV all energy is
heat; above the table ceiling the terminal fractions are held fixed. The old
100 eV Shull--van Steenberg branch and its discontinuity have been removed.

The table's total ionization energy is divided among H I, He I, and He II using
its tabulated ionization counts and threshold energies. Heating, the three
ionization channels, and excitation are normalized together to exact local
photoelectron-energy closure. Both `multiphysics.py` and
`conservative_primordial.py` carry all three secondary species into their rate
equations and ledgers.

Secondary ionization is opt-in in the library API. The P4/P5 validation runners
enable FS2010 by policy unless their own control explicitly selects OFF. If
H I, He I, or He II is numerically absent, the corresponding tabulated
ionization energy is returned to heat instead of creating a nonexistent
species. Excitation is an explicit escaping-line channel, not gas heat. P4/P5
enforce the resulting multi-group photoelectron ledger and record the
line-escape policy.

Target availability is evaluated from the start-of-step species inventory and
held fixed during the opacity iteration. The FS deposition fractions still use
the iterated H II fraction. This removes an artificial fixed-point
discontinuity in which newly created He II could switch its own secondary
channel on during the same timestep; the 200 eV B2 He II residual falls from
`2.47e-4` to `1.28e-5` without a material change in the converged mean xHII.

The FS2010 composition assumption is explicit: primordial H/He, equal H II and
He II fractions, and negligible He III. The current solver uses the actual H II
fraction as the table coordinate. States far from the tabulated He relation are
therefore an approximation to the published closure; arbitrary-composition
fast-electron transport is not claimed.
Trace, non-primordial helium abundances are likewise outside this closure: the
tracked production input has uniform `nHe/nH=0.0789474`, and P4/P5 fail closed
on all H/He ledgers if that contract is violated.

Interpolation in H II fraction is linear, including the sparse 0.1-to-0.5
interval, to match the pinned public implementation. It is not claimed to be a
new higher-order fit in that gap.

## Data provenance

The 14 tables, the upstream MIT license, a README, and per-file checksums live
in `data/furlanetto_stoever_2010/`. `TABLE_MANIFEST.json` pins upstream
21cmFAST commit `892f98c80cfe985ca6b399ec6b51a3aa95124b11`, the source path, DOI
`10.1111/j.1365-2966.2010.16401.x`, and both upstream and vendored SHA256.
Vendoring adds only a final newline to files whose upstream text ended at EOF;
all numeric rows are unchanged. Runtime B2 and P5 reports bind every table and
the manifest by SHA256.

## Numerical gates

- An independent scalar NumPy interpolation reproduces the JIT production
  result to `1.11e-16` maximum absolute error.
- A pinned 21cmFAST reference at 200 eV and `xHII=0.1` checks interpolated
  `f_ion`, `f_heat`, `f_exc`, and all three ionization counts independently of
  the production implementation.
- The maximum channel change from 99.9 to 100.1 eV over all 14 ionized-fraction
  rows is `3.35457e-4`, below the declared `5e-3` continuity threshold.
- The 9.999/10.001 eV table-floor change is `2.27216e-3`, also below `5e-3`.
- Channel energy sums close to `2.22e-16`; one-cell multiphysics and
  conservative primordial ledgers close to `2.78e-17` relative error.
- The table-floor, table-ceiling, non-negativity, finite-value, and all-three-
  species wiring checks pass.
- Hot neutral gas with no photons remains on the exact zero-electron root;
  seeded photo/collisional chemistry closes its electron density to
  `3.48e-12` relative error. Every root search exposes a bracket-success flag,
  and a resolved, single-root H-only step agrees with an independent host RK4
  integration to `2.38e-10` absolute xHII, or `5.74e-4` of the physical change
  in that step. The exact electron-free multi-root limit separately verifies
  preservation of the stationary neutral branch. Large-step backward-Euler
  accuracy remains an explicit stage-5 timestep-convergence question.

The source-bound canonical result is
`data/furlanetto_stoever_validation.json`, checked by
`tests/secondary_furlanetto_stoever_artifact.py`.

## Independent audit

Claude Opus 5 returned a conditional pass in the initial audit, a conditional
pass after the first remediation round, and a final PASS after the resolved
single-root RK4 fixture and exact-neutral branch test were separated and
independently stress-tested:

- `provenance/claude_opus5_secondary_ionization_audit_2026-09-02.md`;
- `provenance/claude_opus5_secondary_ionization_reaudit_2026-09-02.md`;
- `provenance/claude_opus5_secondary_ionization_final_reaudit_2026-09-02.md`.

The final audit leaves large-step/timestep convergence, 256-cubed memory and
runtime scaling, Ly-alpha trapping/radiation pressure, and coupled H+He front
convergence to their declared later gates. None is promoted by stage 2.

## P5 effect measurement

The matched current-solver 0.1 Myr P5 pair described in
`P5_THERMOCHEMISTRY_VALIDATION.md` passes all fixed-point, H/He ledger, thermal,
and finite-field gates. FS2010 produces a volume-mean change
`delta<xHII>=+2.24653e-8` and `delta<T>=-7.52274e-3 K`; maximum local changes
are `2.43194e-4` and `25.1009 K`. The effect is small for this particular hot,
UV-dominated control. B2's deliberately hard 200 eV primordial fixture remains
the strong secondary-response regression (`delta<xHII>=+0.0190614`).

## Reproduction

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
JAX_PLATFORMS=cpu .venv/bin/python tests/secondary_furlanetto_stoever.py \
  --json-out build/secondary_furlanetto_stoever_validation.json
.venv/bin/python tests/secondary_furlanetto_stoever_artifact.py
.venv/bin/python tests/coupled_photo_collisional_hhe.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
.venv/bin/python tests/b2_multiphysics_artifact.py
```
