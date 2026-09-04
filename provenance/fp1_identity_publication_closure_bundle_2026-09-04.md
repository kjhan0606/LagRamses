# F-P1 identity and publication closure bundle

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Implementation commit: `25bd05f` (`Close FP1 identity and publication boundaries`)
Plan review: Fable `APPROVE WITH CHANGES`, recorded in
`fable_fp1_identity_publication_closure_plan_audit_2026-09-04.md`

## Outcome

The authorized F-P1 bundle is implemented and locally verified. It closes
the source-package identity, source-node mapping, derived-artifact
publication, and LC18 failed-wind diagnostic boundaries without selecting a
physical yield package or enabling runtime feedback.

The code-owned state remains fail-closed:

- no physical source nodes are admitted (`physical_node_count=0`);
- the unresolved fate seams remain `[0.8, 1.0]` and `[40, 120] M_sun`;
- canonical conversion, runtime deposition, production readiness, and
  publication readiness remain disabled;
- the LC18 derived cross-check remains review-only because source and derived
  redistribution/approval evidence is incomplete;
- LC18 CDS terminal-wind zero endpoints are partitioned as 4 successful
  controls, 3 failed models, and 7 total.

## Implemented controls

1. The validator registry now rejects a passed identity validator without a
   valid package fingerprint. Positive package admission requires all nine
   identity-matched passed gate reports and binds the selected package hash to
   the executable source-identity evidence.
2. A shared canonical serializer and validator owns the source-node mapping
   schema and SHA-256. Admission checks source-node coverage, coordinates,
   approval identities, package hash, asset hash, and source-node-contract
   hash. The converter compares its generated mapping byte-for-byte with the
   admitted mapping before writing any output.
3. Yield conversion has a non-writing proposal mode. It produces evidence for
   review but cannot create an admitted repository asset or contract.
4. A code-owned publication-rights gate requires a hash-locked terms record,
   verified license, attribution, explicit artifact approval, explicit
   derived-artifact redistribution permission, and `review_use_only=false`.
   The LC18 cross-check calls this gate; mutable labels cannot override it.
5. Successful-control and all-model LC18 diagnostics now expose symmetric
   positive/zero endpoint counts and retain the existing fail-closed anomaly
   behavior.

## Verification

The following completed successfully on GPFS:

- Python compilation of all changed tools/tests;
- the full `simulation/snrt/tests/run_g2_preflight.sh` fixture matrix, ending
  in the expected `G2_PREFLIGHT_BLOCKED` state;
- the F-P1 physical admission, source mapping, publication-rights,
  LC18-cross-check, converter, source-node, terminal-deposition, and
  population/fate contract tests;
- direct state invariants for both unresolved intervals, zero physical nodes,
  publication blocking, and the `48/4`, `53/3`, `101/7` LC18 partitions;
- `git diff --check`;
- deterministic regeneration: all 248 files under `simulation/snrt/config`
  and `simulation/snrt/data` retained identical paths and SHA-256 values
  across the preflight rerun.

No simulation, source selection, CDS redistribution, author contact, or
runtime feedback activation was performed. Bundle-end AGY and Claude Opus 5
audits remain pending under the project audit policy.
