# Claude Opus 5 audit of the independent Fable reproduction

Date: 2026-09-02
Prompt: [`claude_opus5_fable_reproduction_priority_audit_prompt_2026-09-02.md`](claude_opus5_fable_reproduction_priority_audit_prompt_2026-09-02.md)
Scope: read-only review of the Fable report, independent reproduction code/test, source evidence, and revised priority plans.

## Initial verdict

**CONDITIONAL.** The production **BLOCK** and the revised P0--P3 order were
accepted, subject to five factual corrections and two document-consistency
corrections. Opus explicitly stated that the audit did not need to be rerun
after these amendments.

## Required corrections and resolution

1. **F15:** corrected from “compiled runtime uses `ndim+2` and needs a
   compiled matrix” to the actual source evidence: production `patch` writes
   raw energy field 5; the native mirror uses `energy_index=inener`; with
   `NENER=0`, `inener=imetal=6`, so the mirror field map is invalid for that
   layout. This is now statically proven in the reproduction tool/test.
2. **F14:** corrected the historical NVAR comparison. Current
   `bin/Makefile` and the candidate map both declare NVAR=18; NVAR=17 belongs
   to transitional output-00011. F14 is now **partially reproduced** for the
   remaining absent startup validation of He, disabled-element, and untracked
   metal semantics.
3. **F3:** removed the implication that age-axis clamping causes the unit
   error. The report now identifies the direct failure as a `1e-9` wrong
   coordinate passed to interpolation; general endpoint clamping remains a
   separate F7 issue.
4. **F11:** upgraded to **reproduced**. `accounted_mass` resets at each new
   coarse step and advances only after `source_ok`, while the sink clears
   `dMsmbh` only for `ok_blast_agn`; a deferred blast can therefore expose the
   full accumulator again on the next coarse step. A dynamic run is still
   appropriate for measuring the coupled effect, not for establishing the
   source-level retry condition.
5. **Priority ordering:** moved P0.6 field/index semantics ahead of P0.5
   physical-yield promotion.
6. **Plan consistency:** revised the detailed plan's opening G1 wording and
   changed the production-readiness plan's interval prescription to
   `C(current_age)-C(previous_age)`.
7. **Evidence hygiene:** recorded that the G1 runner does not compile the
   RAMSES deposition runtime and made the regression test probe repository
   facts rather than only comparing copied labels.

## Accepted final disposition

After the amendments, the independent reproduction reports 15 findings as
directly reproduced (F1--F12, F15--F17) and two as partially reproduced
(F13--F14); none is dismissed. The accepted order is:

`P0.1 -> P0.2 -> P0.3 -> P0.4 -> P0.6 -> P0.5 -> P1 -> P2 -> P3`.

The existing B3 jobs `330195_2` and `330195_3` remain correctly classified as
diagnostic P3 work in `/gpfs/kjhan/LRD_JWST`; they do not close the stellar/AGN
production gate. No source, runtime, or running job was modified by this
audit.
