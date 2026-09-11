# Separate science plan: galaxy calibration and resolution convergence

Operator-approved separation,2026-09-11. This is not a prerequisite for the
three-step effective-model implementation closeout. Status: planned, no
calibration simulation launched or calibrated parameter claim made.

After a fixed runnable physical configuration is delivered:
1. Predeclare calibration observables, selection functions, observational
   uncertainties and independent holdout observables; distinguish stellar
   mass/size/gas/metal/dust constraints from predictions not used in fitting.
2. Set a bounded parameter prior and computational budget. Use a COLIBRE-type
   simulation design/emulator strategy only where sample size and measured
   predictive error justify it; retain direct simulation verification.
3. Compare identical physical parameters across resolutions (strong
   convergence), then separately assess recalibration (weak convergence).
   Do not assume interpolation of tuned parameters is valid: test at held-out
   resolution with matched phases/ICs where feasible and estimate variance.
4. Report degeneracies, fit/holdout residuals and applicability ranges.

Before launching this science campaign, submit actual observable/data
selection, resolution matrix, run/storage budget and acceptance criteria to
the operator. Do not silently turn it into more implementation gates.
