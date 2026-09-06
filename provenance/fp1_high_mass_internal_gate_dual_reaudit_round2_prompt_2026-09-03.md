# Independent re-audit round 2 — F-P1 40–120 M☉ internal controls

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit, commit, push, download
source data, run full RAMSES, create physical nodes, or enable runtime
deposition. Inspect the live checkout and dirty diff rather than trusting this
prompt or provenance prose.

This round follows the first re-audit in which Fable found three medium defects
that AGY missed. Independently attack the corrected future-promotion path.

Verify at minimum:

1. `audit_fp1_source_node_contract.py` now rejects nonzero direct-collapse
   terminal mass/components/terminal momentum-energy, null failed or direct-
   collapse wind histories, failed/direct-collapse null baryonic remnants,
   malformed rights/identifier types, terminal and wind component non-closure,
   resolver-branch mass-cell overlap, approved branch gaps, and 39/121 M☉
   domain excursions. Check that legitimate branch/source separation is not
   incorrectly treated as overlap.
2. `audit_fp1_physical_package_admission.py` cannot qualify even one gate from
   self-authored evidence and a merely hashed validator file. Current gate
   evidence activation must remain disabled until gate-specific executable
   validators are implemented and code-registered. Confirm there is no config-
   only or absolute-path bypass and that this remaining validator work is
   honestly reported as a physical-package blocker.
3. `convert_yield_rows_to_canonical.py` requires the approved repository node
   contract and matching approval identity, refuses the current zero-node
   review contract, rejects unknown source-node ids, and checks node mass and
   metallicity. Confirm `audit_stellar_yield_asset.py` independently repeats
   repository path, contract audit/approval, mapped-id membership, and
   coordinate checks. Look for substitution or in-process/default-path bypasses
   relevant to the production CLI.
4. Both interval and cumulative Fortran driver entry points refuse channel 3
   above 40 M☉ while the compile-time source-node consumer flag is false, and
   both are now exercised by the native population-ledger test. Confirm the
   explicit 121 M☉ source-node contract negative is in the suite.
5. The roadmap no longer calls contract-only ownership or exactly-once
   declarations a runtime implementation.
6. Re-run bounded F-P1 tests, inspect generated-hash consistency, and search for
   regressions or overclaims. Preserve the distinction between safe internal
   controls and the absent physical source package/runtime consumer.

Return:

- top-level and separate internal/physical/F-P1H-F/production/publication
  verdicts;
- `VERIFIED FIXED`, `PARTIAL`, or `OPEN` for N1--N3 and the cumulative/121/
  wording follow-ups, with exact `file:line` evidence;
- any new finding with severity and a concrete reproduction;
- an explicit answer whether the corrected narrow fail-closed claim is now
  confirmed, while physical admission remains blocked.
