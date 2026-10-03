# SMBH pre-compaction capture ledger

`merge_sink` irreversibly replaces every FOF sink group with one centre-of-mass
sink.  The version-2 capture ledger records the group immediately before that
replacement and commits the batch only after sink-array and sink-particle
compaction finishes on every MPI rank. Logging is rank-0-only and does not change group membership,
sink forces, masses, spins, or the existing merge decision.

## Runtime controls

The `&physics_params` namelist accepts:

```fortran
smbh_capture_ledger = .true.
smbh_capture_ledger_file = 'smbh_capture_ledger_v2.jsonl'
```

Logging is enabled by default.  Give every independent simulation a separate
ledger path; provenance also remains tied to the RAMSES output directory and
its `info`/build metadata.

The writer is invoked only when `smbh = .true.`.  In other sink modes the
ledger controls are ignored.  When enabled for an SMBH run, failure to open,
write, flush, or close the ledger is fatal: RAMSES calls `clean_stop` before
any rank enters the irreversible sink-compaction section.  A partially
written tail from that failed run remains invalid and must not be consumed as
a capture event.

Version 2 uses a two-phase file-level transaction. `batch_begin`, all event
rows, and `batch_prepared` are flushed and closed before compaction. Only after
all ranks finish irreversible compaction does rank 1 append and flush
`batch_commit`. The validator never promotes v2 events from an uncommitted
batch, including a complete early event followed by later I/O failure. A
failure to write the post-compaction commit is fatal and leaves the batch
censored. Legacy standalone version-1 event transactions remain readable.

Each invocation also appends an `attempt_begin`. It carries an execution-unique
`attempt_uid`, a stable `run_uuid`, the configured ledger filename, and either
a null parent (fresh run) or the checkpoint UID selected by the final
MPI-broadcast `nrestart`. Every v2 batch marker carries that attempt UID and a
one-based `committed_batch_seq`; the sequence advances only after the commit
record has been flushed.

Every complete RAMSES output publishes `SMBH_CAPTURE_LINEAGE` immediately
before `COMPLETE`. The sidecar binds the checkpoint UID to its producing
attempt, run UUID, ledger filename, coarse step, output number, and committed
batch high-water mark. Restart with capture logging enabled fails closed unless
both markers exist, their output and coarse-step values match the restored
state, and the referenced parent attempt is present in the configured ledger
with the same run and ledger identities. The runtime gate also streams the
ledger without allocating from the sidecar value, requires contiguous matching
prepare/commit pairs through the checkpoint high-water, and rejects commits
whose coarse step exceeds the checkpoint. Replacing an output removes stale
lineage and completion markers before writing any component.

The runtime scan is a restart-safety gate, not a replacement for catalogue
validation: the Python validator remains authoritative for all event/member/
pair invariants and whole-ledger branch selection.

## Transaction layout

Within a version-2 batch every event is a contiguous JSONL transaction:

1. `event_begin`: integration time, cosmology, code-unit conversions, merge
   radius, FOF group size, classification, group COM, and maximum separation.
2. one `member` row for every original sink: ID, the surviving primary ID,
   a primary flag, mass, position, velocity,
   formation time, accretion/feedback accumulators, BH spin, gas angular
   momentum, and the last available Bondi gas context.
3. one `pair` row for every unordered member pair: minimum-image separation,
   relative velocity, reduced mass, Newtonian two-body energy, angular
   momentum, the current code's legacy `1/r^2` binding proxy, and both binding
   flags.
4. `event_end`: expected member/pair counts and `complete=true`.

Two-member groups are `BINARY`; larger transitive FOF groups are `MULTIPLE`.
No arbitrary binary ordering is inferred for a multiple.

`primary_sink_id` is the global ID retained by `merge_sink`: the most massive
member, with the lowest pre-compaction sink index breaking exact mass ties.
It is stored in both `event_begin` and every `member` row.  Thus each captured
sink row carries the requested `(sink_id, primary_sink_id)` relation without
having to reconstruct the compaction order.

The deterministic event UID contains coarse step, level, minimum/maximum sink
ID, and member count. A restart may append the same committed batch again.
Consumers deduplicate identical batch and event UIDs. If a crash occurs before
`batch_commit`, all v2 events in that batch are censored. A repeated UID with
different content is a provenance conflict, not a valid restart duplicate. A
restart from state after compaction cannot reconstruct a missing commit;
downstream completeness still requires checkpoint/run provenance.

The batch UID additionally binds the pre-compaction sink count, ID sum, and a
position/group-sensitive rotated-XOR identity. This is a compact deterministic
identity, not a collision-proof cryptographic hash. The validator therefore
also compares the complete ordered event digests: if two batches share a UID
but differ in any captured content, it rejects the ledger fail-closed instead
of deduplicating or promoting either conflicting replay.

Validate a ledger with:

```bash
python3 patch/lagRamses/aux/validate_smbh_capture_ledger.py \
  smbh_capture_ledger_v2.jsonl --output-root /path/to/run \
  --checkpoint-uid ATTEMPT-output-NNNNN
```

With `--output-root`, the validator follows the selected leaf attempt back
through COMPLETE checkpoint sidecars. The selected leaf itself and every
ancestor are capped at the corresponding checkpoint's batch high-water and
coarse step, excluding uncheckpointed tails, superseded tails, and siblings.
`--checkpoint-uid` is required when more than one leaf exists or the selected
leaf has more than one COMPLETE checkpoint. No timestamp or append order is
used to choose a physical history. Without `--output-root`, validation checks
transaction structure but cannot certify restart supersession.

For a ledger being written, `--allow-incomplete-tail` reports a final partial
event or uncommitted batch without treating that condition alone as invalid;
its events remain excluded from the reported physical event counts.

## Physical interpretation

`event_begin` is a **numerical-capture event**, not a claim that the SMBHs form
a physical bound binary or coalesce at that time.  The `two_body_bound` field
uses the isolated Newtonian `1/r` pair energy and omits the host potential.  It
is an audit diagnostic for the downstream `kpc_to_pc` state classifier.  The
`legacy_pair_bound` field reproduces the current merge proxy, which divides by
`r^2`; it is recorded only so historical merge decisions can be reconstructed.

The last Bondi context is a local scalar diagnostic, not the stellar/gas/FDM
radial profile needed by the delay model.  Profile extraction and host/galaxy
provenance are separate, non-destructive follow-up products keyed by the sink
IDs and capture time.
