# Reference test for the Group 434 initial-condition analysis

Grammar job 1578549 applied the prepared geometry and Fourier analysis to the existing Group 205 two-fluid ΛCDM initial conditions. It completed in 22 seconds. The source revision and input checksums are included with the results. The analysis read existing fields and did not generate new initial conditions.

For all four levels, the new audit reproduced the previously recorded gas-density minima, counts of density contrast at or below −1, and refinement-mask counts. The base-grid gas density and gas/CDM velocity statistics reproduced the previous matched-parent analysis to within 10⁻¹². The new output also records the CDM displacement components. [Geometry](geometry_audit.json) · [Field statistics](summary.json) · [Fourier shells](shell_coherence.csv)

The older Group 205 z=49 data contain three level-11 raw array cells with gas contrast δb ≤ −1. One coincides with a cell flagged in the same-level refinement map. The array includes cells that may be buffers or covered at another level. An AMR leaf inventory and the RAMSES reader are needed to determine the physical hydro state. A separate Group 205 z=99 high-wavenumber realization recorded no such raw cells through level 11, with a minimum contrast of −0.447. This comparison does not establish the density range of Group 434.

The Fourier comparison measures the gas density, gas and CDM velocity components, and CDM displacement components on the shared level-8 parent grid. It records shell powers and phase coherence with the parent realization. Higher-level modes and absolute normalization against an independent transfer calculation remain to be checked. The Fourier amplitudes for CDM displacement are not CDM density spectra.
