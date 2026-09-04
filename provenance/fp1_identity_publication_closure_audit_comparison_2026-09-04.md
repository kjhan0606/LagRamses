# F-P1 identity/publication closure audit comparison

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Implementation under audit: `25bd05f`; verification boundary: `5aeb6d3`

## Decisions

| Auditor | Role in this historical review | Decision |
|---|---|---|
| AGY / gemini-3.8-flash-high | historical bundle-end audit | PASS |
| Claude Opus 5 | independent bundle-end audit | CONDITIONAL PASS |
| Grok / grok-4.6-build | replacement for Fable's bundle-start plan review | APPROVE |
| GPT-5.6 Sol | conditional-result adjudication | CONDITIONAL PASS |

AGY was subsequently retired from the active auditor roster. The historical
AGY result remains provenance only. The active workflow is Claude Opus 5 for
each completed implementation step and Grok for the driver's plan before a
new bundle begins. GPT-5.6 Sol remains the conditional/negative-result
adjudication path.

## Consensus

All reviewers agree that the current F-P1 implementation is an integrity
boundary rather than a physical-source approval. No reviewer found a live
package-identity bypass, mapping-hash bypass, publication-rights bypass, or
runtime-feedback activation. The current blocked state is correct:

- zero physical source nodes;
- unresolved fate seams `[0.8,1.0]` and `[40,120] M_sun`;
- false production, publication, canonical-conversion, and runtime-deposition
  flags;
- LC18 terminal-wind anomaly unresolved and publication gate false.

The LC18 partition is independently confirmed as 48 positive/4 parsed-zero
successful controls, 53 positive/3 parsed-zero failed models, and 101/7
overall. The word “parsed” is important: a CDS Table 5 endpoint equal to zero
at `0.01 M_sun` print precision is not an attested physical zero. The four
successful parsed-zero endpoints have positive BR26 Wind sums; the three
failed endpoints retain the unresolved BR26 zero-release anomaly.

## Findings adopted for the next driver plan

1. Add an admitted synthetic converter fixture that reaches the positive
   write-path guard, mutates mapping content/hash independently, and proves
   no output is created on disagreement. The driver reproduced the existing
   mismatch guard in an isolated temporary fixture, so this is evidence
   strengthening rather than a current bypass.
2. Rename/qualify LC18 zero counters as parsed/quantized endpoints and expose
   the quantization limit plus successful/failed release outcome counts.
3. Reorder `run_fp1_population_fate_contract.sh` so generated high-mass
   evidence is validated in the same invocation, or run all dependent audits
   after regeneration; add a freshness/dependency-order assertion.
4. Harden the publication-rights API so the code-owned gate reads and hashes
   the locked terms bytes internally. Keep the current terms and derived
   artifact review-only until explicit rights and approval exist.

These are not yet an implementation authorization. They enter the next
driver-written bundle plan. Before that bundle starts, Grok must audit the
plan and the driver must wait for its decision; after each completed step,
Claude Opus 5 audits the implementation. AGY will not be called.

## Verification context

The driver ran the full G2 preflight, which ended in the expected
`G2_PREFLIGHT_BLOCKED`, all F-P1 focused tests and the population/fate contract,
direct state invariants, deterministic 248-file config/data hash comparison,
and isolated F1/F2 reproduction. No physical source was selected, no job was
launched, and no CDS artifact was redistributed.
