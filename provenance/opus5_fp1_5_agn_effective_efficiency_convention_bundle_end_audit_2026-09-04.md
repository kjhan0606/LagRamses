# Claude Opus 5 bundle-end re-audit — F-P1.5-R AGN effective-efficiency convention

Date: 2026-09-04 (KST)  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Mode: read-only; no files edited, jobs launched, runtime flags set, commits,
pushes, or secondary auditor  

## Repair verification

The four repairs from the first `FAIL` are present and current:

1. `patch/lagRamses/snrt_agn_efficiency_smoke.f90:15-16` expects
   `inflow=0.2` for `Bondi=0.2, Eddington=1.0`, matching the shared minimum.
2. The helper SHA256 is
   `34915c8cafc688763e38aa11641e613662ec2bbf420f96f64448cabeaa2bcc01`,
   matching the generated audit report and evidence.
3. The audit regex recognizes the block-form source commit.
4. `retained_seen(isink)=dMsmbh(isink)` is consumed only inside the successful
   transaction block.

## Findings

The shared pure resolver, writer/driver coefficient parity, explicit statuses,
supplied-inflow accounting, one-sided retained check, all-group transaction,
source API, reader/converter, fixture, Makefile prerequisites, and actual
default/SNRT-CUDA build evidence are substantively correct. The fixture is
arithmetically consistent after refresh, and the evidence honestly excludes
live CUDA/hydro, AGN SED, obscuration, escape, jet, dust, and crash-journal
claims.

Two representation details remain to be recorded more explicitly:

- `dMBH_coarse`/`dMEd_coarse` are cumulative rate-based supplied ledgers while
  `accrete_bondi` accumulates the sum of per-fine-step minima.  The minimum of
  cumulative sums can exceed the cumulative sum of minima when Bondi/Eddington
  ordering changes inside a coarse step.  The driver bound also uses the
  instantaneous helper efficiency while accretion applied its per-fine-step
  efficiency, so the retained relation is an upper-bound/approximate
  convention, not an equality.
- At a coarse-step rollover `accounted_inflow` resets while `retained_seen`
  persists.  A retained increment without a successful prior transaction can
  therefore exceed the fresh-step bound and remain skipped until a reset
  rebase; this is fail-closed and must be documented.

The Python `parity=` banner is stronger than its hard-coded algebra oracle:
the native smoke, not that Python test, is the Fortran parity evidence. The
rearm static criterion should also be structurally scoped to the successful
`if (source_ok) ... end if` block, and its no-clamp check is a literal
blacklist rather than a complete semantic proof. Finally, the helper's
zero-Eddington status is deliberately promotable: its ratio/inflow are zero,
and no source is emitted without a positive cumulative supplied increment;
that policy should be stated in the coarse-state documentation.

## Verdict

**CONDITIONAL PASS**

Mandatory bounded follow-up conditions:

1. Add the two supplied-ledger slack sources above to the implementation
   evidence.
2. Record the coarse-step rollover asymmetry and fail-closed consequence.
3. Rename/re-scope the Python `parity=` banner so native smoke is identified as
   the Fortran parity evidence.
4. Scope the rearm regex to the success block and disclose the literal
   blacklist nature of the no-clamp criterion.
5. State the zero-Eddington promotable decision in `AGN_COARSE_STATE.md`.

After these documentation/test-scope amendments, issue one final read-only
bundle audit. No production-code change is required by these conditions.

## Closed items

One pure helper with two production consumers; effective coefficient in Lbol
and all photon groups; no hidden raw clamp; supplied-inflow `idsink` accounting;
commit-ordered marker and post-commit retained cursor; explicit raw/base and
effective ranges/statuses; fresh hash-matched audit; direct Makefile edges;
unity-efficiency API rejection.

## Deferred

AGN SED/obscuration/escape/jet coupling; legacy `accrete_bondi`/`AGN_blast`
parity; live RT-hydro activation; durable restart journal; production-dump
fixture refresh; stellar fate/yield and dust/IR closure.
