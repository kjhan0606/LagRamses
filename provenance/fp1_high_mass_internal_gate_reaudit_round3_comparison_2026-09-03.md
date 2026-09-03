# F-P1 high-mass re-audit round-3 comparison

Date: 2026-09-03
Auditors: AGY (`gemini-3.8-flash-high`) and Fable (`fable`)

Both auditors independently confirmed the current internal fail-closed
implementation and all F1--F5 remediations. AGY found no further bypass. Fable
again inspected the hypothetical future promotion state more deeply and found
four lower-severity identity, rights-vocabulary, and type hardening items.

Those four items are in scope and were implemented before the final targeted
confirmation. They do not alter the scientific status: zero physical nodes,
zero qualified packages, disabled validator execution, absent runtime consumer,
and stale linked-build evidence keep all physical and production gates blocked.
