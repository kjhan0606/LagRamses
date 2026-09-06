# Fable F-P1 high-mass internal-gate re-audit, round 2

Date: 2026-09-03
Model: `fable` through `claude -p --model fable`

## Verdict

- Top level: **CONDITIONAL PASS for engineering / BLOCK for physics**.
- Internal controls: **CONDITIONAL PASS**; N1--N3 fixed, with three new latent
  medium findings.
- Physical gap, F-P1H-F, production, and publication: **BLOCK**.

Fable confirmed that the current checkout cannot approve a node, qualify a
package, convert a canonical row, or integrate channel 3 above 40 M☉. It marked
N1--N3, cumulative-driver coverage, 121 M☉ coverage, and wording `VERIFIED
FIXED`. It also regenerated the F-P1 artifacts byte-identically and verified
all recorded hashes.

## New findings

1. **F1, medium:** approved-node status did not require non-null license,
   redistribution/research-use status, approval id, DOI, retrieval date, or
   package fingerprint.
2. **F2, medium:** converter and asset auditor consulted F-P1H-B node approval
   but not F-P1H-E package admission. A self-approved node contract could
   therefore unlock conversion while package admission remained blocked.
3. **F3, medium:** mapping checked only node id, mass, and metallicity. A
   direct-collapse node could map to a channel-3 row with inconsistent age,
   returned mass, remnant, ejecta, and energy.
4. **F4, low:** F-P1H-E evidence artifacts accepted absolute paths outside the
   repository.
5. **F5, low:** `birth_metallicity_value` accepted strings and booleans.

Fable additionally noted that diagnostic energy on a direct-collapse source
record may legitimately remain nonzero as classifier evidence, and that the
half-open cell convention was not enforced at a non-final right edge.
