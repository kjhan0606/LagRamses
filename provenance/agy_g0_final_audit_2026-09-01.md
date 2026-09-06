# AGY Gate G0 Final Audit — 2026-09-01

## Scope and provenance

- Project: `/gpfs/kjhan/LRD_JWST`
- Gate: **G0 only** — inventory, identity, migration, provenance, environment,
  repository policy, and production fail-closed enforcement.
- Auditor: AGY (Gemini Antigravity CLI), model `gemini-3.1-pro-high`, CLI
  version `1.1.21`.
- Mode: read-only plan audit; no files were modified, assets copied, jobs
  launched, or commits made by the auditor.
- AGY artifact: `/home/kjhan/.gemini/antigravity-cli/brain/35397f06-302c-43a6-9570-10d65404154f/g0_audit_plan.md`
- Local evidence: `simulation/snrt/data/g0_production_manifest_audit.json`

## Verdict

**BLOCKED.** The G0 fail-closed policy is active and detects the unresolved
production identity and reproducibility prerequisites. G0 therefore permits
continued development and audit work, but does not authorize a production
run.

## AGY findings

1. **CRITICAL — production HDF5 payload is not migrated.**
   `p4_snapshot_hdf5` remains `not_migrated` and has no SHA256 in
   `manifests/lrd_jwst_external_assets.json`. The source payload is the
   approximately 866 GB external file; its identity must be recorded only
   after staging to `/gpfs` and hashing the staged file.
2. **CRITICAL — CUDA smoke executable is missing.**
   `snrt_cuda_smoke_executable` is required by the production overlay but is
   registered as `missing` without a hash.
3. **CRITICAL — repository/source identity is dirty.**
   The production overlay requires a clean tree, while the recorded tree is
   dirty. `lagramses_pinned_source` is also
   `available_external_dirty`; its source/build identity cannot yet serve as
   a production-pinned identity.
4. **HIGH — required physical assets and per-asset approvals are incomplete.**
   The required production IDs for stellar yields, stellar SED, AGN SED, dust,
   and the thermal atlas are not registered as approved production assets.
   Existing registry assets lack the required `license_status`,
   `provenance_status`, and `owner` records. AGY reported 27 such existing
   records; the current local fail-closed report enumerates 30 missing records.
   The count discrepancy does not affect the blocked verdict and should be
   reconciled when metadata is populated.
5. **HIGH — environment is recorded but not locked.**
   `manifests/snrt_cpu_environment_2026-09-01.json` records JAX/JAXLIB
   `0.11.1` on CPU and dependency hashes, but its
   `reproducibility_status` is `recorded_environment_not_yet_locked` rather
   than `locked`.
6. **INFO — the fail-closed controls are correctly wired.**
   `simulation/snrt/tools/audit_production_manifest.py` rejects blocking asset
   states, missing required hashes/registrations, unapproved metadata, an
   unlocked environment, a dirty production tree, and disabled fallback/
   legacy rejection policy. The synthetic regression test confirms pass and
   failure paths.

## Local verification

- `PYTHONPATH=. .venv/bin/python tests/production_manifest.py` — passed:
  `PRODUCTION_MANIFEST_TEST_OK fail_closed=true`.
- The real project audit returned exit code `2`, status `blocked`, and
  `production_gate_pass: false`.
- JSON parsing, Python bytecode compilation, and `git diff --check` passed.

## Minimum closure conditions

1. Stage `p4_snapshot_hdf5` under `/gpfs`, verify the staged payload, and record
   its final SHA256 without overwriting the existing external registry history.
2. Build or provide the CUDA smoke executable, record its binary/build identity
   and SHA256, and run the smoke qualification.
3. Establish a clean, reproducible source/build identity for the production
   code. Existing user changes must be preserved or explicitly committed by
   the project owner; no cleanup was performed automatically.
4. Register the approved physical yield, SED, AGN, dust, and thermal-atlas
   assets and populate complete, evidence-backed per-asset license and
   provenance metadata.
5. Lock the CPU/JAX environment (and separately qualify the CUDA environment
   where required) so the fingerprint is reproducible.
6. Re-run the local G0 audit and commission a fresh G0-only AGY audit. G0 is
   closed only when both the machine-readable audit and the independent audit
   return `PASS`.

No G1 work is promoted by this audit. The next action is G0 blocker closure,
followed by the required G0 re-audit.
