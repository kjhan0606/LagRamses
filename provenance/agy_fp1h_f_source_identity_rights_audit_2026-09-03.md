# AGY audit: F-P1 source identity and rights validator

Date: 2026-09-03
Model: `gemini-3.8-flash-high`
Scope: uncommitted source-identity/rights executable-gate bundle relative to
commit `296dd0c`
Mode: read-only implementation and adversarial review

## Verdict

**CONDITIONAL PASS for the narrow registry wiring; FAIL for accepting the
rights gate as established.** Production and publication remain correctly
blocked.

AGY confirmed that the code-owned registry rejects unknown validators,
wrong-gate validators, hash-only evidence, and registry/contract divergence.
It also confirmed that altered source bytes are detected when the declared
manifest and source contract are not changed with them. Those controls are
useful, but the validator does not yet provide an independent trust anchor for
the claimed candidate identity, exact release, file inventory, or rights.

## Findings

1. **High -- candidate identity blockers do not affect the verdict.** A source
   contract candidate-ID or release-root mismatch appends a blocker at
   `simulation/snrt/tools/validate_fp1_source_identity_rights.py:154` and
   `:158`, but `passed` is computed only from requirement booleans at `:318`.
   Blockers are then erased on a pass at `:328`. AGY reproduced a substituted
   candidate identity that returned `passed=true` with an empty blocker list.
2. **Medium -- in-tree symlinks are accepted.** Paths are resolved and checked
   only for escape from the candidate root at `:205-211`. A source entry that
   is a symlink to another file inside the root passes, contradicting the
   claimed immutable-mirror semantics.
3. **Medium -- malformed identity input can raise an uncaught exception.** At
   `:274`, `source.get("article_doi") in manifest_candidate["citation"]`
   raises `TypeError` when the DOI is null. The registry invokes the runner
   directly at `simulation/snrt/tools/fp1_gate_validator_registry.py:74-75`
   without converting unexpected exceptions into a controlled blocked result.
4. **Low -- adversarial coverage is incomplete.** Tests did not exercise
   candidate-ID mismatch, path traversal, duplicate paths, internal symlinks,
   missing files, or malformed/wrongly typed JSON fields.
5. **Low -- the version-record filename is hard-coded.**
   `zenodo_record_19503168.json` at
   `simulation/snrt/tools/validate_fp1_source_identity_rights.py:262` is not
   derived from a pinned release profile.

## Disposition

The generated gate pass is provisional and must not be treated as an accepted
F-P1 physical-package gate. The implementation remains fail-closed at the
higher production boundary: eight other executable gates, four physical
blockers, zero physical nodes, and false canonical conversion, runtime
deposition, production, and publication approvals prevent activation.
