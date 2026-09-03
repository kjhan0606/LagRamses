# F-P1 high-mass re-audit round-1 comparison

Date: 2026-09-03
Auditors: AGY (`gemini-3.8-flash-high`) and Fable (`fable`)

Both auditors confirmed current fail-closed behavior and independently retained
the physical/F-P1H-F/production/publication **BLOCK**. AGY accepted the internal
controls as complete. Fable exercised malformed future promotion inputs more
aggressively and found three medium defects that AGY missed.

Local independent probes reproduced all three Fable findings exactly:

- a semantically invalid direct-collapse record passed;
- arbitrary pass JSON and a one-line validator qualified a candidate;
- an unknown source-node id converted against the zero-node contract.

The findings are therefore accepted as implementation defects rather than
optional audit preferences. They were remediated before round 2. The stale
linked binary, absent physical package, absent runtime consumer, and unrelated
mirror differences remain explicit blockers or later integration work.
