# AGY audit — F-P1 40–120 M☉ internal gates

Date: 2026-09-03
Auditor: AGY / `gemini-3.8-flash-high`
Mode: read-only inspection of `/gpfs/kjhan/LRD_JWST`
Prompt: `fp1_high_mass_internal_gate_dual_audit_prompt_2026-09-03.md`

## Verdict

- Top level: **CONDITIONAL PASS**.
- F-P1H-A--E engineering controls: **PASS with a runtime-consumer condition**.
- Physical 40--120 M☉ gap: **BLOCK**.
- Production readiness: **BLOCK**.
- Publication readiness: **BLOCK**.

AGY found the contracts and negative gates fail closed: no physical package is
qualified, no physical node exists, compiled review identities are blank, and
runtime deposition remains disabled. Source parity is correctly blocked on
stale production-linked build evidence.

## Principal findings

1. **Critical, missing runtime consumer:** channel 3 has an 8--120 M☉ candidate
   domain and a deposition contract, but no source-node fate/deposition consumer
   equivalent to the SNIa adapter is wired into RAMSES. This is mandatory before
   F-P1H-F/runtime activation.
2. **High, expected fail-closed:** the current dirty source is not represented
   by the recorded production-linked binary. A clean rebuild and new linked
   evidence are required before promotion.
3. **High, physical blocker:** all four candidates have zero passed package
   gates, and the admitted physical-node inventory is empty.
4. **High, physical blocker:** 12 of 18 W18/N20 high-mass outcomes are failed;
   failed nodes have zero source remnant records and N20 has no wind record for
   its five failed nodes.
5. **High, source anomaly:** the 100 M☉ W18/N20 stable wind differs by about
   `1.000001e-4 M☉` in Mg; 12 radioactive-epoch warnings reach 36.89 dex, and
   K-40 appears as a cross-segment duplicate in every common record.

AGY independently reproduced 18 outcomes, 12 failed outcomes, 6 terminal-yield
records, maximum rounded-source relative residual `0.003301399`, four common
wind nodes, 12 radioactive warnings, and the K-40 duplication.

## Mandatory conditions

Before F-P1H-F: acquire a qualified multi-Z/multi-rotation package, resolve the
LC18 wind anomaly, supply age-resolved winds or an approved error-bounded
lumping rule, standardize decay, define injected-energy/deposition semantics,
and populate complete physical nodes.

Before runtime activation: implement the channel-3 fate/deposition adapter,
compile the approved identity into a clean binary, regenerate linked-build
evidence, and obtain a passing source-parity result.

## Candor assessment

AGY judged the local wording honest that the physical gap remains unresolved.
Its engineering PASS applies only to the current fail-closed controls, not to
physical source admission or runtime feedback.
