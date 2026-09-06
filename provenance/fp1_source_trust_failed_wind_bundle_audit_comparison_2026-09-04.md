# F-P1 source-trust and failed-wind bundle audit comparison

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Audited base: `296dd0c` plus the uncommitted current bundle

## Independent audit outcomes

| Auditor | Model | Verdict | Scope limitation |
|---|---|---|---|
| AGY | `gemini-3.8-flash-high` | PASS | Full read-only inspection and local checks |
| Claude | Opus 5 (`--model opus`) | CONDITIONAL PASS | No shell available; source/config/artifact review and hand checks only |

Grok was not run because the xAI CLI/service is down, consistent with the
user's instruction. A Codex `gpt-5.6-sol` re-audit is not required by the
current rule because neither independent auditor returned FAIL; Opus returned
conditional pass and explicitly allowed the driver to draft the next bundle.

## Consensus

Both auditors independently found that:

- the code-owned source-rights lock profile blocks candidate substitution,
  coherent manifest/contract/bytes rewrites, symlinks, non-regular files, and
  malformed identity data;
- the exact 108/108 join, 52/56 classification, 845 nodes/19 duplicates,
  3--8 phase range, and 96/12 table7 coverage are correct;
- the release-internal wind residual and the cross-source CDS discrepancy are
  kept separate and are not silently reconciled;
- the four physical blockers remain visible;
- physical nodes remain zero and production, publication, canonical conversion,
  and runtime deposition remain disabled;
- the bundle does not constitute physical source approval and is complete enough
  to draft the next implementation bundle.

## Independent reproduction of conditions

The driver reproduced the following claims outside the auditor sessions:

1. **AGY Python portability observation — confirmed.** Host `/usr/bin/python3`
   is Python 3.9.21, and the LC18 test fails at the upstream
   `zip(..., strict=True)` call. The project `.venv` runner is the supported
   Python 3.13 environment and passes. This is an operational portability item,
   not a production-gate bypass.
2. **Opus F-2 — confirmed, high forward risk.**
   `audit_fp1_fate_admission.py` executes and embeds
   `physical_package_report`, but its final production result is only
   `sidecar_ready and fate_report["production_ready"]`; it never requires the
   physical-package report's production/publication readiness or physical-node
   count. With unresolved fate intervals closed, a stale/over-claimed sidecar
   could therefore outvote a blocked F-P1H-E report.
3. **Opus F-3 — confirmed.** The five physical-package evidence SHA256 values
   are declared inside the editable physical-package contract being audited;
   they are not yet code-owned trust-root values.
4. **Opus F-4 — confirmed.** The current contract requires
   `required_birth_metallicity_domain_selected == false` at admission entry,
   while the selected-package branch requires the same field to be true. The
   promotion branch is unreachable under one contract, and its positive path
   is not covered by a test.
5. **Opus F-18 — confirmed directly.** Calling the new
   `_build_phase_histories` with 108 synthetic records having negative phase
   durations returns a diagnostic with 108 age violations and does not raise.
   The current valid-data test still passes, but a standalone drifted-data tool
   invocation exits successfully with a violation report.
6. **Opus F-15 — confirmed.** CDS rights status is read through the mutable
   terms evidence and the cross-check artifact does not label that field
   `authoritative_for_verdict: false`, although the whole result is review-only.

## Disposition

### Accepted into the next implementation bundle

- commit the current audited bundle before the next implementation begins;
- couple fate-admission production/publication to the physical-package verdict;
- move physical-package evidence digests to a code-owned profile;
- make the positive selection path reachable only under a coherent contract and
  add guarded positive-path/negative-guard tests;
- make LC18 monotonicity and terminal-phase violations fail closed at tool level;
- close the missing adversarial cases (root/source-contract symlinks, duplicate
  candidate identity, runner pass-with-blocker, and controlled exception edges);
- mark mutable CDS rights evidence explicitly non-authoritative;
- add signed and relative cross-source residual summaries;
- eliminate or explicitly identify the duplicated phase aggregation;
- compute, rather than literalize, evidence booleans and scope blocker reports.

### Deferred as low-priority hygiene or long-term validation

The following are real observations but do not block the next high-level
feedback bundle: parent-directory symlink symmetry, `md5` FIPS API spelling,
hex-check symmetry, repository-root digest pinning for tracked Fortran mirrors,
the wording of the age test (duration positivity versus independent
monotonicity), duplicate archive-open CRC hygiene, and the source-node rights
status binding when physical nodes are eventually populated. The informational
duplicate-identity check itself needs no change.

### Not accepted as defects

The intentional zero-node and all-false production/publication/deposition state,
the unresolved failed-wind anomaly, the lack of author confirmation, missing
age-resolved winds, missing injected-energy mapping, missing canonical
momentum/deposition, and the absence of CDS redistribution rights are correctly
blocked scientific boundaries, not implementation failures.

## Overall decision

The current bundle is **conditionally accepted as an implementation/evidence
bundle only**. It may proceed to coherent commit/finalization and to the next
bundle plan, but it may not promote a physical source or open runtime feedback.
The next plan is recorded separately and requires Fable's final-purpose and
feasibility approval before implementation.
