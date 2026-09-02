# Independent SN/AGN feedback physics and implementation audit

Date: 2026-09-02
Model: actual `fable` through `claude -p --model fable`
Mode: read-only
Repository HEAD: `7e6dab63d87707dc4ee1749f242d3a809191cc00`

## Verdict

**BLOCK for production.**

The scaffold is worth keeping, but the compiled feedback runtime is not the
code that passed validation, no physical yield asset exists, and AGN radiation
reaches gas principally through post-processed ledgers. No files were written
or changed by Fable during the audit.

## Executive findings

- The production Makefile builds stellar modules from `patch/lagRamses`, while
  the modules that earned the G1 PASS live in a separate `native/phase0`
  mirror compiled only by the G1 test runner.
- The compiled runtime queries a yield age axis documented in years with an
  age expressed in Gyr; endpoint clamping masks the mismatch.
- Both trees evaluate a forward-shifted cumulative interval. Under variable
  timesteps the increments do not telescope.
- There is no approved wind, AGB, SNII, SNIa, or PISN production grid. Current
  adapters correctly emit zero canonical rows rather than inventing values.
- AGN Bondi/Eddington state and persistence are largely present, but the
  coarse ledger, thermal/jet deposit, and live SNRT source use inconsistent
  luminosity/accreted-mass conventions.
- RT/dust/hydro remains a post-processed ledger chain with limited opt-in live
  photoheating, not complete live feedback.

## Component status

| Component | Status | Evidence highlighted by Fable |
| --- | --- | --- |
| 32-field yield contract, reader, audit | Implemented | `patch/lagRamses/stellar_yield_tables.f90:10-18`, `stellar_yield_audit.f90:52-88` |
| Age units | Broken in compiled tree; fixed only in mirror | See prior-finding closure below |
| Cumulative increments | Partial; interval direction wrong in both trees | See prior-finding closure below |
| IMF/population normalization | Partial | Hard-coded Kroupa and mass windows; no population remnant ledger |
| Out-of-domain policy | Silent clamp in compiled tree; hard error in mirror | `stellar_yield_interpolation.f90:157-158` vs mirror `:171-175` |
| Cell deposition | Implemented, one NGP cell, thermal only | `stellar_ramses_runtime.f90:221-258` |
| Radial momentum | Missing | Contract defines an isotropic source's directed vector as zero |
| Delayed-cooling ownership | Partial but channel-owned | SNII-only in channel mode, total mass loss in legacy mode |
| Star-release restart state | Binary restart only | HDF5 does not persist `tpp`, `mp0`, or `indtab` |
| SNIa/PISN | Disabled and scientifically unapproved | G2 physics contract |
| Physical wind/AGB/SNII grids | Missing | Candidate adapters are review-only and emit zero rows |
| Legacy kinetic-SN energy | Broken | Missing `1e51` factor at `feedback.kjhan3.f90:1319,1987` |
| AGN Bondi/Eddington/retention | Implemented | `sink_particle.kjhan.f90:4012-4035,4456-4497` |
| AGN thermal and jet deposit | Implemented with a different efficiency convention | `sink_particle.kjhan.f90:6902-6929` |
| AGN coarse-state ledger | Implemented before accumulator reset | `sink_particle.kjhan.f90:6172,6303,6312` |
| AGN accumulator restart | Implemented in both checkpoint formats | sink backup/restore paths |
| Ledger deduplication | Missing; duplicates are rejected | `snrt_core/sink_diagnostic.py:148-154` |
| AGN SED/escape/obscuration | Partial and unapproved | Unobscured Sazonov-style pilot, default escape fraction one |
| Nine-group source closure | Implemented; arithmetic closure only | `AGN_NINE_GROUP_VALIDATION.md` |
| Live SNRT AGN source | Partial | Four hard-coded groups, hard-coded ionizing fraction, retained-mass increment |
| Stellar SED | Candidate only | BPASS path clamps most sources and uses escape fraction one |
| Dust absorption/heating | Offline diagnostic implementation | No complete live dust model |
| Scattering/IR re-emission | Missing | Explicitly out of current scope |
| Radiation pressure | Diagnostic only | Not coupled to hydro |
| Thermochemistry atlas | Implemented/B1 PASS; license pending | Metal-only static coupling |
| Live RT to hydro | Partial | Opt-in photoheating to total energy only |
| Live RT to accretion and live stellar RT sources | Missing | Sink-only source loop and no accretion write-back |

## Closure of the four earlier Fable findings

| Earlier finding | Current status | Fable's source evidence |
| --- | --- | --- |
| Years versus Gyr | **OPEN in compiled tree; closed only in mirror.** Not a false alarm. | Table axis is years in `stellar_yield_tables.f90:15,184`; runtime query is Gyr at `stellar_ramses_runtime.f90:168-170,188`; compiled interpolation compares them directly. Mirror converts on load, but `bin/Makefile` builds the patch tree. |
| Cumulative interval direction | **OPEN in both trees.** Not a false alarm. | Runtime uses current end-of-step age as interval start and computes `C(age+dt)-C(age)` in `stellar_source_increment.f90:50-64`. Mirror fixes idempotence, not direction. |
| Missing/duplicated `1e51` | **Missing factor confirmed in legacy Sedov path; channel-mode duplication closed.** | Live Sedov sites omit the factor present in base RAMSES. Phase-0 converts erg once and disables the kinetic call in channel mode. |
| Three versus 11 species/He/NVAR/delayed cooling | **PARTIALLY CLOSED.** | NVAR guard and SNII delayed-cooling ownership close. He initialization still derives from the legacy three-species header; disabled elements and untracked residual are not added to generic metal; the candidate NVAR=18 field map does not match the executed NVAR=17 binary. |

## Severity-ranked findings

| ID | Severity | Finding |
| --- | --- | --- |
| F1 | Critical | No approved physical yield asset for any channel; candidate mass/metallicity gaps include the comparison population. |
| F2 | Critical | Validated mirror and compiled runtime diverge; the mirror has never been linked into the production executable. |
| F3 | Critical | Gyr query against a year axis suppresses release and is hidden by clamping. |
| F4 | Critical | Forward-shifted increments leave gaps/double release under variable timesteps and skip the first interval. |
| F5 | Critical | HDF5 restart drops star birth time, initial mass, and release cursor, so continued release can stop after restart. |
| F6 | High | Legacy Sedov path lacks the `1e51` energy factor. |
| F7 | High | Silent clamping plus a single-mass fixture can assign the 20-solar-mass row to the whole SNII window without a population ledger. |
| F8 | High | Startup is fail-open to the embedded synthetic table; IMF/windows are hard-coded; optional event channels lack production approval gates. |
| F9 | High | Channel mode puts all SN energy thermally into one cell with zero momentum and disables the kinetic scheme. |
| F10 | High | AGN coarse ledger, deposit, MAD jet, and live SNRT source use inconsistent efficiency and accreted-mass conventions. |
| F11 | High | Live SNRT can re-inject a running accreted-mass accumulator when an AGN blast is deferred. |
| F12 | Medium | AGN ledger interval labels/deferred semantics are ambiguous to consumers. |
| F13 | Medium | Deduplication is documentary, converter keys mismatch the reader, and the binary path uses raw rather than effective radiative efficiency. |
| F14 | Medium | He initialization, disabled-element metal loss, untracked residual, and legacy slot reuse remain ambiguous. |
| F15 | Medium | Mirror candidate maps energy to a non-thermal index that can coincide with metal in an NENER=0 build. |
| F16 | Medium | AGN/BPASS SEDs remain unobscured, escape-fraction-one candidates with unresolved normalization/domain sensitivity. |
| F17 | Low | Unused bridges have layout issues; one diagnostic ledger hard-codes sink one; rewind can append conflicting rows. |

## Transitional output-00011 baseline

Fable narrows the allowed interpretation further than the previous handover.
The run used a pre-selector Phase-0 build, nine-row fixture, NVAR 17, no delayed
cooling, no mass loading, no sinks, and no SNRT. The source had the same Gyr
query and fallback behavior. It can support matched hydro/gravity/cooling/star-
formation and layout/decoding regression tests. It cannot support a legacy-
versus-new feedback comparison, enrichment, SN energy, AGN, RT, or dust
comparison. The most accurate label is `effectively_feedback_free_transitional_build`.

## Required acceptance artifacts

1. Compile and test one source of truth: production Makefile objects must be
   the objects exercised by the native contract suite, with a mirror/patch
   divergence gate.
2. Add a compiled-tree age-node test and a variable-timestep telescoping test;
   the first interval must start at age zero and stop/restart must match.
3. Add an HDF5 bitwise round trip for all per-star release state.
4. Approve a cited/checksummed physical full grid with complete coverage,
   age-zero anchor, decay horizon, rotation weighting, wind/terminal ownership,
   license, and population/remnant ledger.
5. Make startup reject a missing external production table and every out-of-
   domain query; make IMF/channel windows explicit configuration.
6. Fix or formally retire the legacy Sedov path and validate one-event energy.
7. Select and validate a physical SN momentum/cooling scheme with an isolated-
   SN resolution study.
8. Unify AGN mass, luminosity, thermal/jet, and live-radiation conventions;
   persist a per-step accreted increment and close energy/mass over restart.
9. Implement duplicate-key coalescing with conflict rejection and fix converter
   field names/effective-efficiency handling.
10. Gate the executed field map, He initialization, disabled-element/untracked
    metal accounting, and total-energy index at startup.

## Priority

- **P0:** F2–F5, F7, F8, F15 in the compiled tree, then F1 physical-asset
  approval. No scientific claim precedes these.
- **P1:** F9 momentum/cooling; F6 fix or retirement; F14 field semantics;
  F10–F13 AGN convention/ledger repairs; then SNIa DTD and PISN gates.
- **P2:** stellar/AGN SED approval and escape/obscuration, BPASS domain,
  normalization sensitivity, dust scattering, IR re-emission, radiation
  pressure, and thermal-atlas license.
- **P3:** live stellar sources, nine-group runtime spectrum, ionization
  write-back, RT-to-accretion, B3 convergence, AMR/MPI determinism, and a
  production rerun.

## Prohibited claims

- That the current compiled lagRamses build implements the validated Phase-0
  contract, or that the G1 PASS applies to it.
- Any stellar mass return, enrichment, SN energy, or physical wind/AGB/SNII/
  SNIa/PISN result from the current build or output-00011.
- That output-00011 is a feedback-on or pure legacy-feedback baseline.
- Any channel-mode SN momentum result.
- Any AGN luminosity/energy/duty-cycle result mixing the incompatible ledger,
  deposit, and live-SNRT conventions.
- Any live stellar/AGN radiation-hydrodynamic, dust-obscuration, IR, radiation-
  pressure, or transported-field claim before the corresponding gates pass.
- Full nine-group RT beyond the source-side arithmetic closure.
