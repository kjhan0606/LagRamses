# Independent re-audit request: SNRT stage 2 FS2010 secondary ionization

You are the mandatory independent reviewer for stage 2 of the SNRT
production/publication-readiness program in `/gpfs/kjhan/LRD_JWST`.
The first audit is preserved at
`provenance/claude_opus5_secondary_ionization_audit_2026-09-02.md` and returned
`CONDITIONAL PASS`. Inspect the current dirty worktree directly. Do not edit
anything. Re-run focused tests where useful. Judge physical correctness,
algorithmic wiring, numerical justification, provenance, and whether every
first-audit stage-2 finding is honestly closed or explicitly bounded.

Return exactly one overall verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
`PASS` is allowed only if no production-reachable stage-2 blocker or unresolved
conditional finding remains. Separate genuinely later-gate limitations from
stage-2 defects. Cite files/lines and measured evidence.

## Claimed remediation map

### B1 zero-helium violation and A4 unsafe defaults

- `snrt_core/multiphysics.py` and `snrt_core/thermochemistry.py` now default
  secondary ionization OFF.
- H I/He I/He II channels are enabled only when the start-of-step target
  inventory is non-negligible. Unavailable ionization energy is routed to
  heat, preserving the photoelectron ledger. The availability mask is fixed
  through the opacity iteration to avoid a newly made He II channel switching
  on inside its own fixed point.
- `tests/secondary_furlanetto_stoever.py` contains a zero-He regression that
  requires both secondary-He channels, both He ledgers, and the energy residual
  to be zero.
- `tests/p2_p3_validation.py` now uses primordial `nHe/nH=0.079`, explicitly
  opts in, checks He ledgers, a dimensionless float32 energy ledger, and root
  bracketing. Its xHII remains 0.090666.
- The hard 200 eV B2 He II fixed-point residual changed from 2.47e-4 to
  1.2778e-5 after fixing the within-iteration target switch; species-specific
  residuals are recorded and gated.

### C1 runtime photoelectron ledger and C6 excitation policy

- `ThermochemicalStepResult` in `snrt_core/multiphysics.py` carries exact
  photoelectron energy and its ledger residual.
- `snrt_core/thermochemistry.py`, `tools/p4_run_transport_pilot.py`, and
  `tools/p5_run_thermochemical_pilot.py` accumulate it. P4/P5 fail closed at
  1e-12 in float64 or 1e-5 in float32.
- P4/P5 record
  `excitation_energy_treatment=radiative_line_escape_not_returned_to_gas`.
- The canonical P5 ON multi-group ledger is 1.3120e-17; OFF is exactly zero.

### C2 independent Solver A/B chemistry

- `snrt_core/conservative_hydrogen.py` again uses its independent H-only
  analytic closure, with mean-state electron density and no shared coupled
  H/He root or collisional chemistry.
- Solver A and B share transport/cross sections/case-B H physics but not the
  chemistry iteration. The regenerated B2 mean absolute xHII difference is
  1.52855e-5, accepted below a predeclared 5e-5 bound (about 0.25% of the
  fixture mean). See `tools/validate_multiphysics_b2.py`, canonical JSON, and
  `B2_PRODUCTION_SOLVER_VALIDATION.md`.

### C3 independent interpolation reference and A5 table-floor continuity

- `tests/secondary_furlanetto_stoever.py` retains the independent NumPy path
  and adds a hard-coded reference vector at E=200 eV, xHII=0.1 derived from the
  pinned upstream 21cmFAST tables, checking all fractions and ion counts.
- It explicitly checks 9.999/10/10.001 eV continuity. Maximum floor delta is
  2.27216e-3 below 5e-3; 99.9/100.1 eV delta is 3.35457e-4.
- `data/furlanetto_stoever_2010/TABLE_MANIFEST.json` pins upstream commit,
  DOI, source path, license, and upstream/vendored SHA256 values.

### C4 nearest-root policy

- `snrt_core/implicit.py` returns a per-cell `root_bracket_found` flag. If the
  nearest-root sample scan finds no crossing it solves on the guaranteed full
  physical bracket but reports false; B2/P4/P5 hard-fail any such cell.
- `tests/coupled_photo_collisional_hhe.py` checks the exact neutral hot root,
  seeded electron closure (3.478e-12), every bracket flag, and comparison to an
  independent high-resolution host RK4 H-only integration at a resolved step
  (absolute xHII error 5.541e-5 below 1e-4).
- Large backward-Euler timesteps can select a materially different branch;
  this is disclosed as a stage-5 timestep-convergence requirement, not hidden
  as a stage-2 accuracy claim.

### C5 P5 effect strength and A2 provenance

- `P5_THERMOCHEMISTRY_VALIDATION.md` and
  `SECONDARY_IONIZATION_VALIDATION.md` explain the small effect from the hot,
  already-ionized control (`<xHII>=0.978`, `<T>=5.58e6 K`) and quantify FS
  suppression with xHII.
- Regression bands are now narrow: 1e-8 < delta mean xHII < 5e-8 and
  -0.02 K < delta mean T < -0.002 K.
- The regenerated matched pair gives +2.24653e-8 and -7.52274e-3 K; all three
  secondary channels are nonzero.
- The P5 report binds every `snrt_core/*.py`, the validator, runner, input,
  photon metadata, thermal atlas, manifest, and all 14 tables by SHA256.

### A3 shared thresholds and other declared boundaries

- Threshold energies are imported from `snrt_core/secondary.py` rather than
  duplicated literals.
- FS2010's primordial-composition assumption, use of actual xHII as table
  coordinate, linear upstream-consistent xHII interpolation, high-energy
  clamp, and escaping-line excitation policy are documented. Arbitrary-
  composition fast-electron transport is not claimed.
- The 65-sample nearest-root scan's production memory scaling remains a later
  performance-engineering concern noted by the first audit; determine whether
  it is merely advisory or prevents stage-2 PASS.

## Canonical evidence

- `simulation/snrt/data/furlanetto_stoever_validation.json`
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`
- `simulation/snrt/data/p5_secondary_ionization_validation.json`
- `simulation/snrt/data/helium_case_b_recombination_validation.json`
- `simulation/snrt/SECONDARY_IONIZATION_VALIDATION.md`
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`
- `simulation/snrt/P2_MULTIPHYSICS.md`
- `simulation/snrt/P5_THERMOCHEMISTRY_VALIDATION.md`

The canonical artifact tests bind source hashes and currently pass. The P5
HDF5 pair is intentionally ignored by git but its paths and SHA256 values are
bound in the tracked canonical report.

## Reproduction commands

Run from `simulation/snrt` with the existing `.venv`:

```bash
JAX_PLATFORMS=cpu .venv/bin/python tests/secondary_furlanetto_stoever.py
.venv/bin/python tests/secondary_furlanetto_stoever_artifact.py
JAX_PLATFORMS=cpu .venv/bin/python tests/coupled_photo_collisional_hhe.py
XLA_FLAGS=--xla_force_host_platform_device_count=2 JAX_PLATFORMS=cpu \
  .venv/bin/python tests/p2_p3_validation.py
.venv/bin/python tests/b2_multiphysics_artifact.py
.venv/bin/python tests/p5_secondary_ionization_artifact.py
```

Audit the claims rather than trusting this map. End with a finding-by-finding
closure table for B1, C1-C6, and A2-A5, then the single overall verdict.
