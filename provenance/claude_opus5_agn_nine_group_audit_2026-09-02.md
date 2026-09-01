# Independent stage-4 audit — AGN nine-group source closure

Model: `claude-opus-5`
Mode: read-only
Repository HEAD during audit: `f792791069810c92c3e23c2339cc16370cef15b3`

## Verdict

**CONDITIONAL PASS**

M9 and Mo2 were independently found genuinely closed. The auditor confirmed
that group/ledger wiring and provenance are enforced on the production path,
and that stale field artifacts are demoted or banner-marked rather than
represented as current.

## Independently confirmed

- The former five source columns map bit-identically to new groups 3–7 and the
  new `[2000,10000] eV` group carries `3.529964866772728e51 s^-1`.
- Config, metadata, static HDF5, and transport control all contain nine groups.
- Edge mismatch was tested by mutation and rejected by the production metadata
  loader.
- The exact old leaked cross sections `6.538494382598191e-22`,
  `5.260332312477607e-22`, and `9.150757399259818e-23` are now zero.
- Photon-denominator additivity agrees to `1.11e-16`; H I/He I/He II numerator
  additivity agrees to `0`, `2.22e-16`, and `1.11e-16`. An independent SciPy
  integration agreed with the closure to at most `7e-7`, consistent with
  trapezoid discretization.
- All hard-X diagnostics and cross sections reproduce from the CSV and JSON.
- The failed `0.104502 Myr` probe is preserved and accurately disclosed; it is
  not referenced by a passing gate.
- Required and additional AGN/SED/dust/B2/RSLA tests passed, except the stale He
  provenance gate described below.

## M9/Mo2 closure

| Finding | Status |
| --- | --- |
| M9, five groups ending at 2 keV | **CLOSED**: nine groups through 10 keV in CSV, metadata, static HDF5, and runner control. |
| M9, edge agreement | **CLOSED**: enforced in loader, source rebind, validator, and artifact test. |
| M9, field changes | **DEFERRED HONESTLY** to B3, as permitted by this stage charter. |
| Mo2, threshold endpoint leakage | **CLOSED**: absorber numerators start at their physical thresholds; exact-zero tests pass. |

## Findings requiring remediation

- **F1 — MAJOR:** `helium_case_b_recombination_validation.json` still binds the
  pre-stage-4 `primordial.py` hash, making its fail-closed test red although the
  recombination physics is unchanged. Regenerate it and update its report hash.
- **F2 — MODERATE:** group 2 is serialized as `[5.6,11.2] eV` but the declared
  SED support begins at 10 eV, so it is integrated only on `[10,11.2]`. It is
  incorrectly labelled like a fully supported group, and the limit text does
  not disclose the partial support. Also distinguish the hard group's 4.38%
  share of in-band emitted energy from its approximately 0.79% share of
  bolometric luminosity.
- **F3 — LOW:** hard-group and below-threshold criteria use fixed indices rather
  than deriving masks from configured edges and species thresholds.
- **F4 — LOW:** the source chain depends on gitignored HDF5 files and needs an
  external deposit before publication; current behavior does fail closed.
- **F5 — LOW:** candidate hashes are checked but candidate IDs are not compared
  with ledger IDs, five of ten ledger rows have zero luminosity, and the static
  sidecar's gas-input hash is recorded but not gated.

The audit judged photon and photoelectron ledgers correctly as arithmetic
closure and the species/fixed-point results as solver self-consistency, not as
physical field validation.
