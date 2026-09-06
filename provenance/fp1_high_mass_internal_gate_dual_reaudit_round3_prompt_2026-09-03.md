# Final independent re-audit — F-P1 40–120 M☉ internal controls

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit, commit, push, download
physical data, run full RAMSES, create source nodes, or enable runtime
deposition. Inspect the live checkout and dirty diff.

Rounds 1 and 2 found real latent promotion defects. Attack the final corrected
path rather than trusting the implementation note.

Verify:

1. An `approved_physical_nodes` contract cannot contain null or mistyped rights,
   provenance, population/binary state, approval identity, package fingerprint,
   or metallicity, and source nodes obey half-open mass-cell membership.
2. F-P1H-E evidence cannot use absolute or escaping paths. Non-empty declarative
   gate evidence remains unable to pass until executable validators are
   implemented and code-registered.
3. Converter and asset auditor require both F-P1H-B node approval and F-P1H-E
   admitted-package approval, with matching approval id, selected package/source
   hash, and selected mapping hash. Confirm there is no circular hash dependency
   between package contract and generated mapping.
4. `fp1_source_node_projection.py` binds canonical wind/SNII rows to source-node
   channel, mass, metallicity, age/lifetime, cumulative wind, terminal ejecta,
   remnant, energy kind, and momentum. Reproduce the prior malicious direct-
   collapse row and try wrong channel, pre/post lifetime, component, energy,
   momentum, and mapping substitutions against both converter and asset audit.
5. The current zero-node/zero-package checkout still emits no canonical output,
   keeps runtime deposition false, refuses both interval and cumulative high-
   mass driver paths, and reports physical/F-P1H-F/production/publication BLOCK.
6. Run the bounded F-P1 suite, inspect all hash chains and wording, and report
   any new concrete bypass. Do not demand physical data as an internal-control
   fix; list it separately as the next scientific gate.

Return separate verdicts, `VERIFIED FIXED/PARTIAL/OPEN` for F1--F5, exact
`file:line` evidence, any new finding with a concrete reproduction, and an
explicit statement whether the internal fail-closed implementation is now
confirmed without claiming that the physical high-mass gap is solved.
