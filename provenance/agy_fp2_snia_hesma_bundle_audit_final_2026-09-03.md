# AGY final bundled audit — F-P2 SNIa/HESMA review package

- Model: `gemini-3.8-flash-high`
- Date: 2026-09-03
- Review scope: DTD/event ledger, HESMA source mirror, profile estimator,
  selection packet, approval sidecar, contract audit, and runner evidence

## Verdict

- Review-only scaffolding: **PASS**
- Production source selection/runtime activation: **BLOCK**

AGY found the checksum chain, fail-closed sidecar, quarantine behavior,
relative audit provenance, self-rooted runner, and canonical promotion-field
requirements coherent.  It did not treat those checks as approval of the
physical SNIa model.

The local audit is authoritative for the profile numbers: `n300c` has relative
discrepancy 6.4104173893 (641.04%), and `n1600c` has 0.05112189345 (5.112%).

## Remaining blockers

The selected model/mixture, DTD normalization and IMF conversion, decay
mapping, returned mass, WD debit/remnant ownership, signed momentum and
deposition policy, realization, thermal coupling, metallicity dependence,
source commit binding, and named approval remain unset.  Runtime activation
therefore remains false.

This record is an independent review result; the final warning/mirror
hardening was subsequently revalidated by the full local F-P2 runner.
