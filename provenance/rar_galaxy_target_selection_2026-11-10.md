# RAR target selection, 10 November 2026

## Purpose and selection

The pilot sample consists of NGC 2403, NGC 7793, DDO 154, and NGC 2976. NGC 3198 is a supplementary high-inclination test. The sample deliberately spans stellar surface brightness and gas importance. It is a targeted convenience sample, not a representative or statistically complete population. It was not selected using which gravity model fits best. The population residual-slope analysis retains its separately defined SPARC selection.

The local source `data/sparc/SPARC_Lelli2016c.mrt` is the SPARC catalogue. Quality values 1 and 2 mean high and medium quality. The data inventory `reports/leroy_sparc_overlap_1267515/overlap.csv` counts observed rotation-curve points with Vobs > 0 and fractional velocity error <= 0.1 falling within the distance-adjusted endpoints of the Leroy HI radial table. The counts are geometrical coverage only, not a claim of complete gas measurements or an accurate three-dimensional source. The principal population selection uses inclination < 70 degrees. The separate HI inventory does not impose that upper cut. The numbers are not interchangeable with the population-analysis radial counts.

| Galaxy | D (Mpc) | SPARC inclination (deg) | Q | HI interval coverage / qualified | Purpose and limitation |
|---|---:|---:|---:|---:|---|
| NGC 2403 | 3.16 | 63 +/- 3 | 1 | 47/72 | Existing matched force baseline. The adopted FR02 gas reconstruction is not established as the exact SPARC input. |
| NGC 7793 | 3.61 | 47 +/- 9 | 1 | 23/26 | Existing independent spiral test. Inclination uncertainty and molecular gas need propagation. |
| DDO 154 | 4.04 | 64 +/- 3 | 2 | 1/11 | Gas-rich dwarf and low-acceleration contrast. Leroy HI table extends only to 1.3153 kpc after distance adjustment. An extended HI profile and force consistency are prerequisites. |
| NGC 2976 | 3.58 | 61 +/- 10 | 2 | 21/21 | Compact stellar-dominated comparison and gas-source systematics. SPARC gas contribution and THINGS-derived source differ substantially. Resolve provenance, non-circular motions and geometry before accepting a force model. |
| NGC 3198 | 13.80 | 73 +/- 3 | 1 | 22/37 | Supplementary extended spiral. Excluded from the i < 70 degree primary population analysis. Inclination, distance and outer gas coverage require sensitivity checks. |

NGC 2976 has an existing source-provenance audit, not a validated new QUMOND galaxy solution. NGC 7793 and NGC 2403 already have fixed-source compiled QUMOND analyses, with documented numerical and input limitations. DDO 154 and NGC 3198 are proposed targets, not completed simulations. Individual-object literature and data sources include the [SPARC archive](https://astroweb.case.edu/SPARC/), [THINGS rotation curves](https://doi.org/10.1088/0004-6256/136/6/2648), and [Leroy et al. gas profiles](https://doi.org/10.1088/0004-6256/136/6/2782). See `docs/ngc2976_sparc_provenance_2026-11-04.md` and `docs/leroy_sparc_overlap_2026-11-04.md` for the input discrepancies and exact source audits.

Before a new fixed-source solve, obtain an HI surface-density profile with quantified coverage outside the observed radii, check helium convention and distance/inclination, document molecular gas and disc thickness, and compare the baryon-only Newtonian force with the observed mass-model contribution. Use common stellar and gas fields for QUMOND and the CDM mass model. Run resolution, box, boundary, and covariance checks. Do not infer a cosmological zoom match from similarity to any observed galaxy. High-resolution simulations and the SPARC population test remain separate stages.
