# Targeted final confirmation — F-P1 high-mass internal controls

Read-only audit `/gpfs/kjhan/LRD_JWST`. Do not edit or run full RAMSES. Confirm
the round-3 latent hardening against the live checkout:

1. F-P1H-E accepts only its five exact repository evidence paths and executes
   the exact source-node contract audit; alternate in-repository files,
   absolute paths, `..`, and symlink escapes fail.
2. Converter and asset audit compare F-P1H-E's source-node evidence path/hash to
   the exact source-node contract they bind.
3. A selected package requires valid package and mapping SHA256 values; the
   selected package hash equals every admitted node package fingerprint and the
   converter/sidecar source hash.
4. Approved rights statuses reject `denied`/`forbidden`; conditional binary axes
   reject bool/int.
5. Re-run the bounded F-P1 suite and verify that the current checkout remains
   zero-node, zero-package, conversion-disabled, runtime-disabled, and physically
   BLOCKED.

Return `VERIFIED FIXED/PARTIAL/OPEN` for these four items, any concrete new
bypass, and separate internal versus physical/production/publication verdicts.
