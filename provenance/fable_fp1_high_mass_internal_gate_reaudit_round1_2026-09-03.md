# Fable F-P1 high-mass internal-gate re-audit, round 1

Date: 2026-09-03
Model: `fable` through `claude -p --model fable`
Mode: read-only inspection; temporary artifacts only under `/tmp`

## Verdict

- Top level: **CONDITIONAL PASS for engineering / BLOCK for physics**.
- Internal F-P1H-A--E controls: **CONDITIONAL PASS**.
- Physical gap, F-P1H-F, production, and publication: **BLOCK**.

Fable confirmed that every remediation claimed in the implementation note was
present and that the checkout failed closed. It nevertheless declined the word
“complete” after independently finding three medium future-promotion defects:

1. The node validator accepted nonzero tracked terminal ejecta for direct
   collapse, null wind histories for failed/direct-collapse nodes, failed nodes
   with no remnant, overlapping mass cells, and an integer license id.
2. Nine self-authored pass JSON files plus a hashed one-line “validator” file
   could set a candidate's derived `production_qualified` state to true. Empty
   physical-node inventory still stopped final selection.
3. The converter accepted a source-node id absent from the source-node contract,
   including while the repository contract contained zero nodes and denied
   canonical conversion. The asset audit did not require the repository path or
   membership of mapped ids.

Fable also found low-severity evidence gaps: only the interval driver guard was
tested, the 121 M☉ node-domain negative was manually probed rather than in the
suite, the linked binary predates the guard, one roadmap sentence overstated
contract-only ownership/exactly-once work, and patch/native mirrors retain
unrelated known differences.

All headline W18/N20 numbers reproduced. Fable reported no present activation
bypass; the three medium findings matter when a physical package is introduced.
