# Independent reproduction of the Fable SN/AGN audit

Date: 2026-09-02
Project root: `/gpfs/kjhan/LRD_JWST`
Repository HEAD: `7e6dab63d87707dc4ee1749f242d3a809191cc00` (worktree contains
uncommitted provenance/B3 work)
Fable source report: [`fable_sn_agn_feedback_audit_2026-09-02.md`](fable_sn_agn_feedback_audit_2026-09-02.md)
Reproduction code: [`reproduce_fable_sn_agn_findings.py`](../simulation/snrt/tools/reproduce_fable_sn_agn_findings.py)
Regression test: [`test_fable_sn_agn_reproduction.py`](../simulation/snrt/tests/test_fable_sn_agn_reproduction.py)

## Result

The independent audit confirms the production **BLOCK** decision. Fifteen
findings are reproduced directly, and two are partially reproduced because
the current code contains some fail-closed protections or require a compiled
dynamic test for the final claim. No Fable finding was dismissed.

This was a read-only audit. It inspected the production Makefile, both source
trees, Fortran source, JSON/CSV metadata, and pure numerical counterexamples.
It did not rebuild or launch a RAMSES executable, so it does not estimate the
observable bias in a full simulation.

Command and result:

```text
python3 simulation/snrt/tests/test_fable_sn_agn_reproduction.py
FABLE_SN_AGN_INDEPENDENT_REPRODUCTION_OK reproduced=15 partial=2
```

The numerical checks independently reproduce two central failures:

- For `C(t)=t^2`, variable steps evaluated as `C(age+dt)-C(age)` give 64,
  while a cumulative ledger must telescope to 36 from age zero.
- A physical age of 1 Gyr passed as `1.0` to an age axis declared in years is
  a coordinate ratio of `1e-9`; this is a wrong-coordinate interpolation
  (not an age-axis clamp).

## Fable finding disposition

| ID | Severity | Independent result | Evidence / reason | Priority |
|---|---|---|---|---|
| F1 | Critical | Reproduced | G2 is `blocked_no_approved_physical_asset`; candidate assets are review-only and the legacy table is `legacy_only`. | P0.5 |
| F2 | Critical | Reproduced | `bin/Makefile` selects `patch/lagRamses`, while the G1 runner builds `native/phase0`; their source objects have different hashes, and the G1 source list does not compile the RAMSES deposition runtime. | P0.1 |
| F3 | Critical | Reproduced | The table documents `age_yr`; the compiled runtime constructs and passes `age_gyr`. | P0.2 |
| F4 | Critical | Reproduced | Both source-increment trees use the forward interval; the numerical counterexample gives 64 instead of 36. | P0.2 |
| F5 | Critical | Reproduced | HDF5 particle backup/restore contains `birth_epoch` but not `tpp`, `mp0`, or `indtab`. | P0.3 |
| F6 | High | Reproduced | lagRamses Sedov expressions omit `1d51`, unlike the base RAMSES expression; channel-mode kinetic feedback is separately disabled. | P1.2 |
| F7 | High | Reproduced | Compiled interpolation silently substitutes endpoint nodes and runtime uses hard-coded mass windows with a one-point fixture. | P0.4 |
| F8 | High | Reproduced | Compiled runtime has an embedded fallback; IMF/windows are source defaults; optional channels are not physical production models. | P0.4 |
| F9 | High | Reproduced | Channel runtime deposits to one target cell and the kinetic path is disabled; radial source-frame momentum is not implemented. | P1.1 |
| F10 | High | Reproduced | Sink blast, coarse ledger, and live SNRT source use different efficiency/accreted-mass paths, including a live-source half-factor. | P1.3 |
| F11 | High | Reproduced | `accounted_mass` resets every coarse step and advances only after `source_ok`, while the sink resets `dMsmbh` only when `ok_blast_agn`; a deferred blast leaves the full accumulator eligible again. | P1.4 |
| F12 | Medium | Reproduced | Rate records expose `aexp/time_code` and instantaneous inflow but not explicit interval endpoints or commit/defer state. | P1.3 |
| F13 | Medium | Partially reproduced | Reader and merger reject some duplicate IDs/collisions; deterministic coalescing and one raw/effective-efficiency schema remain absent. | P1.5 |
| F14 | Medium | Partially reproduced | NVAR=17 is historical output-00011 metadata, while the current Makefile and candidate map both declare NVAR=18; the remaining issue is absent startup validation of He/disabled-element/untracked-metal semantics in the compiled path. | P0.6 |
| F15 | Medium | Reproduced | The compiled path writes raw energy field 5, but the mirror uses `energy_index=inener`; with NENER=0, `inener=imetal=6`, so the mirror field-map validator rejects the production layout. | P0.6 |
| F16 | Medium | Reproduced | BPASS and AGN SEDs retain candidate clamps/unobscured and `fesc=1` assumptions. | P2.1 |
| F17 | Low | Reproduced | Diagnostic writer appends, hard-codes sink ID 1, and has process-local header state, allowing rewind conflicts. | P1.5 |

The full machine-readable result is emitted by the reproduction tool with
`--json`; its source markers are computed from the current working tree rather
than copied from the Fable report.

## Revised implementation priority

The following order is now authoritative for the SN/AGN feedback work. A
lower-priority item cannot promote a higher-priority blocked gate.

### P0 — make the executed stellar runtime trustworthy

1. **P0.1 source-of-truth and build parity (F2, F15).** Make the exact
   `patch/lagRamses` objects used by `bin/Makefile` the objects exercised by
   the native contract suite, or remove the duplicate mirror. Add a source
   identity/hash gate and compile-time report of `NVAR`, `NENER`, `inener`,
   `imetal`, `idelay`, and `ichem`.
2. **P0.2 time and cumulative interval contract (F3, F4).** Convert
   `age_yr` exactly once at the table boundary, resolve the RAMSES `aexp`
   convention, and replace the forward-shifted interval with a declared
   `[previous_age, current_age]` transaction. Add age-node, first-interval,
   variable-dt, repeated-call, and stop/restart tests.
3. **P0.3 complete HDF5 restart state (F5).** Persist and restore `tpp`,
   `mp0`, `indtab`, and any new stellar ledger/progress fields bitwise. Test
   continuation from a checkpoint against an uninterrupted run.
4. **P0.4 fail-closed production input (F7, F8).** Remove the embedded
   synthetic fallback from production, reject out-of-domain queries, and
   move IMF/channel windows and enabled-channel status into versioned,
   approval-checked configuration.
5. **P0.6 field/species/He semantics (F14, F15).** Validate the executed
   field map at startup; prove the NENER=0 index relations and close H, He,
   tracked elements, untracked residual metal, total energy, and
   delayed-cooling ownership without relying on transitional metadata.
6. **P0.5 approve physical yield coverage (F1).** Promote no candidate until
   wind, AGB, and SNII have a complete mass--metallicity--age grid, age-zero
   anchor, release convention, remnant/energy/momentum fields, license, and
   immutable approval sidecar. SNIa/PISN remain separately gated.

P0 exit requires the compiled-tree test matrix to pass. The existing G1 PASS
is retained as a native-mirror result only; it is not a production-runtime
pass until P0.1 is closed.

### P1 — close SN and AGN physical/bookkeeping semantics

1. **P1.1/P1.2 SN deposition and energy.** Select radial momentum, spatial
   deposition, delayed-cooling, and thermal/kinetic ownership; fix or formally
   retire the legacy Sedov path and pass a one-event dimensional `1e51 erg`
   test.
2. **P1.3/P1.4 AGN convention and transaction.** Define one accreted-mass,
   efficiency, luminosity, SED, thermal/jet, photon, and momentum convention;
   write explicit interval/commit/defer fields and make live-source updates
   atomic and exactly once across retry/restart.
3. **P1.5 identity and diagnostic integrity.** Define duplicate-key policy,
   preserve effective efficiency through every converter, replace sink-one
   diagnostics with actual IDs, and make rewind/restart output partition-safe.
4. **P1.6 event channels.** Implement SNIa DTD/binary normalization and the
   explicit PISN/PPISN population gate only after their physical source inputs
   are approved.

### P2 — approve spectra and dusty physics

Approve the stellar and AGN SEDs, low-metallicity/young-age domain, IMF
normalization, obscuration, and escape prescription (F16). Then complete the
dust mixture/depletion model, scattering, IR re-emission, radiation pressure,
and thermal-atlas licensing/coverage tests.

### P3 — live coupling and science qualification

Implement live stellar/AGN sources, nine-group runtime spectrum, ionization
write-back, RT-to-accretion coupling, B3 timestep/spatial convergence,
AMR/MPI ownership and determinism, restart equivalence, and the final
production rerun. No production or publication claim is permitted before
P0--P3 gates are green.

## Current decision

The Fable BLOCK is independently supported. The project should begin at
**P0.1**, not at a new physical-yield import or a production rerun. B3 jobs
already in flight are retained as diagnostic SNRT work and do not alter this
stellar/AGN production gate.
