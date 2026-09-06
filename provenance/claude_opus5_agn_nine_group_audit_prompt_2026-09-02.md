# Independent stage-4 audit: AGN nine-group source closure

Audit the current dirty worktree in `/gpfs/kjhan/LRD_JWST` read-only with
`claude-opus-5`. Do not edit files. This is the mandatory audit after stage 4
of the SNRT production/publication-readiness sequence.

Return exactly one verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`. PASS
requires the original M9 and Mo2 findings to be closed, the group/ledger wiring
and provenance to be hard-gated rather than documentary, and no stale field
artifact to be represented as current. Do not require the next B3 timestep and
spatial-convergence result to pass this source-ledger stage, but do verify that
its necessity and the failed large-step probe are disclosed honestly.

## Original findings

Read the source audit at
`provenance/claude_opus5_rt_audit_2026-09-01.md`, especially M9 and Mo2.

- M9: `config/p0_photon_group_edges_ev.txt` defines nine groups through
  10 keV, while the former production ledger had five groups ending at 2 keV
  and falsely stated that P0 ended there.
- Mo2: closed endpoint sampling leaked non-zero threshold cross sections into
  groups wholly below the H I, He I, and He II thresholds.

## Implementation to inspect

- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/p4_attach_pilot_sources.py`
- `simulation/snrt/tools/validate_agn_nine_group_ledger.py`
- `simulation/snrt/tests/agn_photon_ledger.py`
- `simulation/snrt/tests/p0_sed_closure.py`
- `simulation/snrt/tests/agn_nine_group_artifact.py`
- `simulation/snrt/data/p4_pilot_agn_photon_ledger.csv`
- `simulation/snrt/data/p4_pilot_agn_photon_ledger.json`
- `simulation/snrt/data/agn_nine_group_validation.json`
- `simulation/snrt/data/p4_coeval_static_rt_input_agn9.json`
- `simulation/snrt/AGN_NINE_GROUP_VALIDATION.md`

The common SED quadrature now constructs each group with explicit endpoints,
but integrates each absorber numerator only from
`max(group lower, species threshold)` to
`min(group upper, Verner fit maximum)`. Judge whether this correctly removes
the threshold leak without dropping the final photon-quadrature segment or
double-counting finite measure. Metadata declares `[lower,upper)` except the
final group closed on the right.

The metadata loader now requires top-level `group_edges_ev` and exact equality
of every interval. The canonical validator additionally requires bit-exact
agreement with `config/p0_photon_group_edges_ev.txt`, nine CSV columns, exact
sub-threshold zeros, positive 2–10 keV closure, CSV/metadata totals, source IDs,
and all generation/static-input/runner hashes.

## Canonical measured result

- `[2,10] keV` photon rate: `3.529964866772728e51 s^-1`.
- photon-number ratio to `[0.5,2] keV`: `0.2194249407776824`.
- energy-power ratio to `[0.5,2] keV`: `1.015228072754792`.
- fraction of total emitted energy power: `0.0437845714507376`.
- H I/He I/He II hard-group cross sections:
  `3.1479425e-25`, `9.5803070e-24`, `8.2470181e-24 cm^2`.
- all groups wholly below each species threshold have exact zero cross section.

The canonical source was rebound to an external `(10,9)` coeval static HDF5.
A `0.001 Myr` production P4 runner control with 32 opacity iterations passes
photon ledger `1.42e-16`, H ledger `1.77e-8`, fixed point `2.33e-10`,
photoelectron ledger `8.50e-22`, and zero root failures. A separate full-CFL
`0.1045 Myr` one-step probe failed H ledger at `7.57e-3`; the report discloses
this and assigns field-level timestep/spatial remeasurement to B3.

The old five-group P5 report was moved, not deleted, to
`data/p5_secondary_ionization_validation_legacy5_historical.json`; its test now
verifies historical provenance differs from current AGN9 metadata/core rather
than pretending source closure. B2 and RSLA were fully rerun after the core
change and retain identical physical values with current hashes.

## Required checks

Run from `simulation/snrt`:

```bash
PYTHONPATH=. .venv/bin/python tests/agn_photon_ledger.py
PYTHONPATH=. .venv/bin/python tests/p0_sed_closure.py
PYTHONPATH=. .venv/bin/python tests/agn_nine_group_artifact.py
PYTHONPATH=. .venv/bin/python tests/b2_multiphysics_artifact.py
PYTHONPATH=. .venv/bin/python tests/rsla_refinement_artifact.py
PYTHONPATH=. .venv/bin/python tests/p5_secondary_ionization_artifact.py
PYTHONPATH=. .venv/bin/python tests/p4_dust_runner.py
PYTHONPATH=. .venv/bin/python tests/p5_dust_runner.py
```

Independently recompute the hard-X diagnostics from CSV/JSON, inspect the
actual diff, verify exact edge and source/hash closure including ignored HDF5
artifacts, and distinguish arithmetic ledgers from physical validation. Give a
closure table for M9/Mo2, findings with severity, and one overall verdict.
