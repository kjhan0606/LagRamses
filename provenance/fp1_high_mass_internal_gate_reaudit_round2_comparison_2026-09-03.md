# F-P1 high-mass re-audit round-2 comparison

Date: 2026-09-03
Auditors: AGY (`gemini-3.8-flash-high`) and Fable (`fable`)

Both auditors confirmed the round-1 fixes and current fail-closed behavior.
AGY gave the internal controls a full PASS. Fable again performed the stronger
future-state adversarial audit and identified F1--F5, including two ways in
which a future node/package could be declaratively promoted without all
intended semantics being enforced.

The findings were accepted because they are concrete and locally reproducible,
not because an auditor requested them. Remediation binds approved-node rights,
F-P1H-E admission, selected package hash, source-node mapping, and the complete
canonical cumulative payload while prohibiting external evidence paths.
