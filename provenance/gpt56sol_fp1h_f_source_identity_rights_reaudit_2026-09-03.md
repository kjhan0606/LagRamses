# Codex gpt-5.6-sol independent re-audit: F-P1 source identity and rights

Date: 2026-09-03
Model: Codex `gpt-5.6-sol`
Scope: the same uncommitted bundle audited independently from AGY
Mode: read-only code inspection and temporary adversarial fixtures

## Verdict

**FAIL.** The five rights requirements are not independently established.
Production and publication are currently blocked, but several admission
invariants are incomplete.

## Independently reproduced findings

1. **High -- candidate-ID and release-root mismatches can pass.** The validator
   records blockers but does not include blocker emptiness in `passed`; this is
   the same concrete bypass identified independently by AGY.
2. **High -- empty inventories can pass.** The implementation compares the
   manifest and contract sets at
   `simulation/snrt/tools/validate_fp1_source_identity_rights.py:166-199`, but
   neither the exact expected set nor non-emptiness is pinned. Empty equal
   inventories produce a syntactically valid composite at `:252-256`.
3. **High -- a self-consistently rewritten package can pass.** Source bytes,
   manifest hashes, source-contract hashes, and their composite can all be
   changed together. The validator recomputes consistency at `:220-256`, but
   has no code-owned or separately locked expected file set, byte counts,
   per-file hashes, and composite fingerprint.
4. **High -- DOI/version and rights are self-asserted.** The manifest, local
   source contract, staged Zenodo JSON, and terms record can be rewritten
   together and pass `:269-305`. A minimal terms object with the expected
   status strings and any non-empty citation is accepted. There is no pinned
   DOI/record/license/attribution trust profile independent of these editable
   inputs.
5. **High -- blocked fate admission does not require publication to be
   false.** `simulation/snrt/tools/audit_fp1_fate_admission.py:270-305` checks
   canonical conversion, runtime deposition, and production consistency but
   does not reject `publication_ready=true` in a blocked sidecar.
6. **Medium -- the claimed immutable mirror is not immutable.** The staged
   external tree is writable and untracked, and internal symlinks are accepted.
   The implementation should use the accurate term `hash-locked local mirror`
   unless operating-system immutability is actually enforced and audited.
7. **Medium -- fate-sidecar artifact paths are not confined.**
   `simulation/snrt/tools/audit_fp1_fate_admission.py:64-70` accepts absolute
   paths and resolves repository-relative paths without checking that they
   remain under the SNRT/repository root. A self-consistent external artifact
   set can therefore be substituted.
8. **Medium -- malformed types and dates are under-validated.** Floating byte
   counts compare equal to integers at
   `validate_fp1_source_identity_rights.py:228-231`; the date regex at
   `:279-280` accepts impossible calendar dates; null DOI raises `TypeError`;
   and an unhashable validator ID can escape registry handling at
   `fp1_gate_validator_registry.py:67`.

## Controls that resisted testing

- Unknown or wrong-gate validators, undeclared hash-only evidence, and
  contract/registry divergence were rejected.
- Absolute and `..` manifest paths, duplicate manifest paths, external
  symlinks, missing declared files, and altered bytes without matching
  metadata changes were rejected.
- The checked-in production, publication, canonical-conversion, and runtime
  deposition fields are currently false.

## Required remediation

- Pin exact candidate/release identity, DOI, record ID and filename, non-empty
  file inventory, integer byte counts, per-file SHA256 values, and composite
  fingerprint in a code-owned or separately locked trust profile.
- Make a pass require both all requirement booleans and an empty blocker list.
- Reject every symlink and non-regular path component; validate strict types
  and actual ISO dates; convert malformed input and unexpected validator
  exceptions into controlled blocked results.
- Bind redistribution to pinned structured rights evidence rather than mutable
  status strings and a merely non-empty citation.
- Require false publication state while blocked and confine all sidecar
  artifact paths to pinned repository-relative locations.
- Add an adversarial mutation matrix covering each bypass above.
