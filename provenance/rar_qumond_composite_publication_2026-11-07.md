# Native composite QUMOND publication

Sources f3859bd and 7ca18e0 connect the previously converged native composite
Newtonian and QUMOND candidate solves to full-hierarchy phi, force and source
publication. The opt-in RAR_QUMOND_COMPOSITE_PUBLISH path requires stage-history
activation, serial CPU, periodic fixed two-level topology, static gas, no PIC,
cosmology, restart, subcycling, stars/sinks/cooling/chemistry/RT or other gravity.
It gathers gas plus the declared analytic source directly on leaves, solves
both stages with true relative residual at most 1e-9, constructs conservative
phantom flux, restricts covered coarse fields, and publishes only after both
solves converge. Explicit legacy gravity half kicks are skipped. Optional
first-epoch leaf/subface exports support independent reconstruction.

Build 1275524 compiled successfully. Matched execution 1275530 and default-off/
history analysis 1275531 completed. All sixteen old/new-off gravity arrays are
bitwise identical. Independent publication verifier 1275556 FAILED: the static
gas acquired momentum. Its earlier geometry, physical source, units, both true
Poisson residuals and independent phantom/force reconstruction checks passed,
but that partial success does not establish publication acceptance.

Source inspection identified a second gravity application hidden inside the
legacy cuRamses set_uold wrapper. Static Godunov returns without transport, but
set_uold still calls add_gravity_source_terms. Correction 631678f skips the
whole hydro-update block only in the explicitly frozen-source publication
mode; ordinary evolution is unchanged. This is within the reviewed frozen
source design, not a new physical model. Build 1275598, matched run 1275601,
default-off/history audit 1275602 and independent verifier 1275603 were submitted
on Grammar. All completed: build 2m55s, matched run 19s, history audit 4s and independent verifier 6s. Both audit artifact manifests verify. All sixteen default-off arrays remain bitwise equal. The independent calculation covers 59312 leaves and 180564 subfaces: N residual 9.66689487188552e-10, Q residual 9.651508668088988e-10, phantom relative difference 1.4249643653715343e-15, force relative difference 7.173273467265204e-17, analytic source difference exactly zero. Both saved outputs retain zero gas momentum and thermal energy to relative 1.1293772630057339e-16. Published leaf fields agree within 3.55e-15 absolute, and covered potential restriction within 8.89e-16. This passes the bounded frozen-source publication acceptance. Audit paths: /gpfs/kjhan/RAR_QUMOND/audits/rar_stage_history_1275602 and /gpfs/kjhan/RAR_QUMOND/audits/rar_publication_1275603.
Existing failed outputs are preserved.

The independent verifier reconstructs the native subface operator in NumPy,
checks source identity and MOND units, true N/Q residuals, phantom density and
least-squares force, leaf publication and covered-cell restriction. It also
requires unchanged uniform gas density, zero initial momentum and unchanged
initial thermal energy. This measures the declared discrete frozen-source
operation, not continuum accuracy, particle force consistency, MPI/regrid/
restart, evolving hydro, a cosmological MOND model or a galaxy RAR result.

Next reviewed requirements remain continuous-source/force accuracy at two
resolutions, then actual coupled source deposition and synchronized evolving
publication with lifecycle support. Existing legacy continuous-source studies
must be distinguished from this new native operator. The full RAR programme
and subsequent BIMOND/FDM work remain active.
