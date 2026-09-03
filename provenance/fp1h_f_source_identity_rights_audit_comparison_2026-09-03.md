# F-P1 source identity/rights audit comparison

Date: 2026-09-03

## Combined disposition

The first executable source-rights gate is **not accepted**. AGY rated the
registry wiring conditionally acceptable but found a candidate-substitution
bypass; the independent Codex `gpt-5.6-sol` re-audit returned FAIL and both
reproduced the same core defect. Under the project audit cadence, the shared
finding is mandatory remediation rather than a disputed model opinion.

| Area | AGY `gemini-3.8-flash-high` | Codex `gpt-5.6-sol` | Combined action |
| --- | --- | --- | --- |
| Registry/contract wiring | Narrow conditional pass | Confirmed resistant to basic bypasses | Retain and harden exception boundary |
| Candidate identity | Concrete passing substitution | Independently reproduced | Make blockers verdict-bearing and pin exact identity |
| Package bytes | Detected internal-symlink gap | Found empty and self-consistent rewrite attacks | Pin exact non-empty inventory, bytes, hashes, composite; reject all symlinks |
| DOI/version/rights | Null DOI crash; hard-coded record filename | Editable local records are circular evidence | Add independent code/lock trust profile and strict structured evidence |
| Sidecar boundary | Not highlighted | Publication and external-path invariants missing | Require all blocked approval flags false and pin relative paths |
| Malformed inputs | Requested broader adversarial tests | Reproduced type/date/exception gaps | Add exhaustive mutation matrix and controlled blocked outcomes |

AGY was faster at identifying the principal verdict bug and operational
fail-closed gap. `gpt-5.6-sol` was materially stronger on trust-root analysis,
cross-file mutation attacks, sidecar invariants, and malformed-type coverage.
Their conclusions are complementary rather than conflicting.

The Grok audit is excluded: the xAI service was unavailable after a correct
read-only invocation, and the user directed the project to proceed without
it. No Grok verdict is represented here.

Production safety was not breached. The current source package remains
unqualified, the physical-node inventory is empty, and canonical conversion,
runtime deposition, production, and publication remain disabled. The generated
rights-gate pass is nevertheless invalid as approval evidence until the shared
and additional findings are remediated and re-audited.
