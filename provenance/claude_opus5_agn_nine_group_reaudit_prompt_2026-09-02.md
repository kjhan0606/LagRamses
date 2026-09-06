# Focused re-audit: AGN nine-group source closure

Audit the current dirty worktree in `/gpfs/kjhan/LRD_JWST` read-only using
`claude-opus-5`. Do not edit files. This is the mandatory final audit for SNRT
stage 4. Start with exactly one verdict: `PASS`, `CONDITIONAL PASS`, or
`BLOCK`. Stage 4 may close only on `PASS`.

Read the initial audit at
`provenance/claude_opus5_agn_nine_group_audit_2026-09-02.md`, then independently
inspect the actual worktree and test the remediations. Preserve the original
scope: M9 and Mo2 source/group closure belong to this stage; B3 transported
field timestep/spatial convergence belongs to the next stage and must be
honestly disclosed but is not a prerequisite for this source-ledger PASS.

## Initial findings that must now be closed

- F1: the He case-B artifact had a stale `primordial.py` hash.
- F2: `[5.6,11.2] eV` was only supported on `[10,11.2] eV` but labelled as
  fully supported; 4.38% in-band power was not distinguished from bolometric.
- F3: validator indices and threshold slices were hard-coded.
- F4: ignored HDF5 dependencies lacked a stage-specific external asset
  manifest and publication-deposit status.
- F5: candidate IDs/kinds/positions, luminous-row count, gas-input hash, gas
  metadata, and zoom manifest were not all hard-gated.

## Remediation implementation to inspect

- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/p4_attach_pilot_sources.py`
- `simulation/snrt/tools/validate_agn_nine_group_ledger.py`
- `simulation/snrt/tests/agn_photon_ledger.py`
- `simulation/snrt/tests/agn_nine_group_artifact.py`
- `simulation/snrt/data/p4_pilot_agn_photon_ledger.{csv,json}`
- `simulation/snrt/data/p4_coeval_static_rt_input_agn9.json`
- `simulation/snrt/data/agn_nine_group_external_assets.json`
- `simulation/snrt/data/agn_nine_group_validation.json`
- `simulation/snrt/data/helium_case_b_recombination_validation.json`
- `simulation/snrt/HELIUM_RECOMBINATION_VALIDATION.md`
- `simulation/snrt/AGN_NINE_GROUP_VALIDATION.md`

Expected corrected facts:

- groups wholly below 10 eV are zero; `[5.6,11.2] eV` is explicitly partial
  with supported interval `[10,11.2] eV`; groups above it are fully supported;
- the hard group and 0.5--2 keV comparison group are located by exact edge
  matching, and H/He zero-opacity masks derive from the production Verner
  thresholds;
- hard-X power is 4.378457% of the represented 10 eV--10 keV SED power and
  0.794605% of candidate bolometric luminosity;
- there are 10 candidate/ledger rows, exactly 5 luminous; ordered IDs, kinds,
  and positions agree;
- the static sidecar's gas input, gas metadata, zoom manifest, edge table,
  photon files, and rebind tool hashes are verified;
- ignored HDF5 files have exact path, size, SHA256, purpose, and pass/fail
  status in the stage manifest. Publication deposit remains explicitly pending
  final archive; do not mistake local reproducibility for a completed deposit;
- the stale He artifact is regenerated and its physical errors remain
  unchanged.

Run at least these commands from `simulation/snrt`:

```bash
PYTHONPATH=. .venv/bin/python tests/agn_photon_ledger.py
PYTHONPATH=. .venv/bin/python tests/p0_sed_closure.py
PYTHONPATH=. .venv/bin/python tests/agn_nine_group_artifact.py
PYTHONPATH=. .venv/bin/python tests/helium_recombination_artifact.py
PYTHONPATH=. .venv/bin/python tests/b2_multiphysics_artifact.py
PYTHONPATH=. .venv/bin/python tests/rsla_refinement_artifact.py
PYTHONPATH=. .venv/bin/python tests/p5_secondary_ionization_artifact.py
PYTHONPATH=. .venv/bin/python tests/p4_dust_runner.py
PYTHONPATH=. .venv/bin/python tests/p5_dust_runner.py
```

Independently recompute candidate luminosity accounting and both hard-X power
fractions. Check that the validator's dynamic masks are correct, that the
external manifest really fails closed on altered/missing assets, and that no
historical five-group field result is presented as current. Distinguish
arithmetic/self-consistency ledgers from physical field validation. Return a
finding-by-finding F1--F5 table, M9/Mo2 closure, any new findings with severity,
and the single final verdict.
