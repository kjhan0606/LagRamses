# Fable final confirmation — F-P1 high-mass internal controls

Date: 2026-09-03
Model: `fable` through `claude -p --model fable`

## Verdict

- Exact evidence paths plus executable source-node audit: **VERIFIED FIXED**.
- Package/node evidence path and hash parity: **VERIFIED FIXED**.
- Package/mapping SHA256 and node/source equality: **VERIFIED FIXED**.
- Rights vocabulary and conditional binary-axis types: **VERIFIED FIXED**.
- New concrete bypass: **none found**.
- Internal controls: **PASS**.
- Physical/production/publication: **BLOCK**, unchanged.

Fable attacked alternate repository paths, absolute and traversal paths,
external symlink targets, extra/missing evidence keys, divergent node-contract
copies and hashes, null/non-hex/different package and mapping hashes, rights
values such as `denied` and `forbidden`, case and whitespace variants, and
boolean/integer/float binary-axis values. All were rejected.

All F-P1 tests regenerated under `/tmp` and passed. The source-node inventory
remained empty, all four candidates remained unqualified, runtime deposition
remained false, and both Fortran high-mass entry points remained unavailable.

Fable noted two non-bypass cleanups for future selection work: type approver/date
fields more strictly, and add a direct regression for the valid state in which
an approved node contract coexists with an intentionally blocked package.
