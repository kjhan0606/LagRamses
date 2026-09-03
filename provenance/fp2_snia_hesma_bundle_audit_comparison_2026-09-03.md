# F-P2 SNIa/HESMA bundle audit comparison

- Date: 2026-09-03
- Independent auditors: AGY (`gemini-3.8-flash-high`) and Fable (`fable`)
- Scope: DTD/event ledger, HESMA 15-model review package, profile estimator comparison, selection packet, approval sidecar, contract audit, and directly relevant tools/tests
- Local verification: `bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh` passed, including native tests, HESMA tests, selection packet test, production mirror compile, and contract audit

## Headline comparison

| Area | AGY | Fable | Assessment |
|---|---|---|---|
| Review-only verdict | PASS for scaffolding | CONDITIONAL PASS for evidence | Same operational conclusion: review-only is acceptable; production is blocked |
| DTD mathematics | Strong, concise confirmation | Strong, with additional cancellation warning | Both agree; Fable is more adversarial on numerical edge cases |
| Hash/wiring chain | Reports complete and consistent | Reports consistent, but flags path dependence and missing commit binding | Fable gives the more useful publication/reproducibility qualification |
| Physical closure | Lists all missing fields | Same, plus normalization/IMF, realization, thermal coupling, metallicity | Fable is more complete for production planning |
| Profile screen | Correctly diagnostic | Correctly diagnostic; flags hardcoded threshold/common bias | Fable provides the stronger statistical qualification |
| n300c/n1600c | Production blockers | n300c should be classified as a source anomaly; n1600c warning needs context | Fable adds an important evidence-status distinction |
| Runtime/fail-closed | Confirms multilayer enforcement | Confirms multilayer enforcement | Agreement is strong |
| Promotion process | Not materially discussed | No defined approval/promotion schema | Fable identifies a real lifecycle gap |
| Ledger/remnant/momentum | Treats SNIa zero remnant as consistent and notes unresolved momentum | Identifies WD remnant-reservoir debit and radial-vs-vector momentum as blockers | Fable exposes deeper coupling to the population ledger; this requires code-level reconciliation before promotion |
| Runner completeness | Accepts negative coverage as present in test files | Notes the shell runner does not invoke those negative tests | Fable is more precise about what the claimed runner proves |

## Common findings to accept into the F-P2 backlog

1. Keep production activation blocked until model/mixture, DTD bounds/normalization, decay mapping, returned mass/remnant, energy, momentum, and population weighting are explicitly approved.
2. Keep n300c and n1600c out of production selection until the profile discrepancies are resolved or formally excluded.
3. **Resolved in the current implementation bundle:** add a stable near-α=−1
   DTD evaluation path and its tests.
4. Make the 5% profile threshold an explicit diagnostic contract field and document the common profile bias.
5. Extend promotion prerequisites to include IMF conversion, realization policy, thermal coupling, and metallicity dependence.
6. Define a versioned promotion/approval schema before enabling any source or warning-bearing model.

## Disagreements requiring engineering resolution

### 1. WD remnant debit

AGY accepted `channel_snia` having zero terminal-remnant ownership. Fable points out that this does not by itself account for a SNIa consuming a previously formed WD reservoir: the current population ledger computes living mass from aggregate returned and remnant fields, while the SNIa event source has no remnant-debit term. This is not an immediate runtime bug because SNIa is disabled, but it is a valid production-design blocker. Resolve it with an explicit ledger contract and conservation test rather than choosing one audit interpretation.

### 2. Momentum representation

AGY confirms that signed momentum is unresolved. Fable makes the needed distinction explicit: a spherical source has zero net vector momentum, while feedback deposition may require a scalar radial momentum budget or a declared derivation from energy and ejecta mass. The production contract must define both the event-frame vector convention and the deposition-layer coupling.

### 3. Hash portability and linking evidence

AGY verifies current-file hash continuity. Fable additionally checks portability and build evidence: absolute paths inside hashed JSON make regeneration checkout-dependent, and the SNIa modules are in the production source list despite a “not linked” wording in the runner/provenance. Treat current hashes as valid for this checkout but not yet publication-grade until these claims are corrected.

## Quality conclusion

Neither audit is sufficient alone for production approval. AGY is better at confirming the intended algorithm and the current fail-closed architecture with a clear top-level verdict. Fable is better at adversarial review of promotion readiness, evidence wording, reproducibility, physical ledger coupling, and test-claim scope. For this project, Fable's additional findings should be retained as conditional blockers/required design clarifications, while AGY's PASS should be read narrowly as “the review-only scaffolding is coherent,” not as a production approval.

## Current disposition

The audit comparison is recorded. No audit finding was silently promoted into runtime code. The F-P2 bundle remains `review_only_selection_pending` / `runtime_activation_allowed=false`; the next implementation work should address the reconciled blocker list in a grouped change, followed by one bundled re-audit.
