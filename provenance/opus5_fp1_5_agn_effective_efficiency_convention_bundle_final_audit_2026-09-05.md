# Claude Opus 5 final bundle-end audit — F-P1.5-R

Date: 2026-09-05 (KST)  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Mode: read-only; no files edited, builds/tests rerun, jobs, runtime flags,
commits, pushes, or secondary auditor

## Verdict

**PASS** — the F-P1.5-R engineering bundle is accepted.

## Verified closure

- `snrt_agn_efficiency.f90:56-165` is one pure RAMSES-independent resolver;
  the writer and driver each call it once, with one MAD transform and no
  duplicated consumer formula.
- The explicit status matrix covers spin-disabled default, spin-enabled
  uninitialized raw efficiency, invalid raw values/rates, zero Eddington,
  invalid `X_floor`, and invalid effective efficiency.  The raw/base
  `(0,1)` and effective `[0,1)` conventions, including the documented
  spin-disabled raw-zero exception, are consistent.
- The writer uses `epsilon_eff` for `Lbol` and emits raw/resolved/effective
  efficiency, status, mode, and contract at the pre-feedback/pre-reset
  boundary.  The driver uses the same effective coefficient for every photon
  group and has no former `0.99` raw clamp.
- Source accounting uses the `idsink`-keyed increment of cumulative
  `min(dMBH_coarse,dMEd_coarse)`.  `dMsmbh` remains a one-sided retained-mass
  check, and the marker/cursor advance only follows a complete all-group
  transaction.
- The source API names its first argument `delta_inflow_mass_code` and rejects
  unity/super-unity efficiency.  The reader, converter, fixture, static audit,
  Makefile graph, and recorded build evidence are coherent.
- F1 is closed: `AgnCoarseState` carries `efficiency_status` and
  `efficiency_contract_ok`; P4 emits both columns and refuses false contracts
  before opening an artifact.  C1 is closed: status-flagged non-promotable
  rows remain readable in the append-only diagnostic ledger but are refused
  by state promotion and P4 conversion.  Tests cover floor-disabled,
  spin-uninitialized, and null raw-efficiency rows.
- Current helper SHA256 is
  `34915c8cafc688763e38aa11641e613662ec2bbf420f96f64448cabeaa2bcc01` and
  matches the machine-readable audit report.  `git diff --check` is clean.

## Non-blocking follow-up before live activation

1. Record the rank-local `accounted_inflow`/`retained_seen` limitation and the
   reset/`nstep_coarse` boundary window; cross-rank sink migration can otherwise
   re-emit a marker-local increment.  MPI/live coupling is already deferred.
2. Rename the driver's `supplied_mass` variable because it briefly receives an
   instantaneous rate before being overwritten by cumulative mass.
3. Strengthen the native smoke from required-bit subset checks to exact status
   equality, and add explicit `raw_nonfinite` and `floor_nonfinite` cases.
4. Consider adding `sink_diagnostic.py` and
   `p4_build_agn_rate_ledger.py` hashes to the machine audit provenance.

## Deferred scope

No AGN SED, obscuration, escape fraction, jet/radiation-pressure, hydro/dust,
legacy `accrete_bondi`/`AGN_blast` parity, live RT-hydro/AMR/MPI production,
durable crash journal, production-dump fixture refresh, stellar fate/yield, or
physical/publication claim follows from this engineering bundle.
