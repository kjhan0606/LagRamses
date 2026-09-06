# Fable F-P1 high-mass internal-gate re-audit, round 3

Date: 2026-09-03
Model: `fable` through `claude -p --model fable`

Fable marked F1--F5 `VERIFIED FIXED` and confirmed the current internal
fail-closed implementation. It reproduced rejection of null/mistyped rights,
NaN/boolean metallicity, interior right-edge membership, gap/overlap, absolute
and symlink-escaping evidence, declarative gate evidence, mapping/package hash
substitution, wrong channels, pre/post-lifetime payloads, component swaps,
energy/momentum substitutions, and twin-node remapping.

The physical gap, F-P1H-F, production, and publication remain **BLOCK**.

Fable found four lower-severity latent hardening items for the future validator
activation path:

1. Pin F-P1H-E's source-node evidence path to the repository contract and
   compare its hash to the converter/auditor-bound node contract.
2. Bind `selected_package_sha256` to approved nodes' package fingerprints.
3. Restrict approved research-use and redistribution status to an allowed
   vocabulary rather than accepting strings such as `denied`.
4. Type-check conditional binary population axes rather than accepting booleans
   or integers.
