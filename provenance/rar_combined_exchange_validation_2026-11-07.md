# RAR combined local exchange validation

Validated source: c27b714933d50e37067f3a8dc324c4a33b062cac. The three modified Fortran files match the archived source of Grammar build 1275118 byte for byte. Build, three-step evolution 1275121, and audit 1275122 completed with exit 0. Physics: baryon-only QUMOND, hydro, original cooling/Haardt-Madau background, ordinary and effective-Ia enrichment, seeded AGN; dust and radiation transport disabled. The fixture uses artificial star formation, 32 cubed uniform cells, one MPI rank and one OpenMP thread.

Measured local accepted exchanges:

- 328 positive stellar-return transactions and 3608 element records. Maximum stellar mass residual: 8.470329472543003e-22 code mass. Maximum element residual: 1.0587911840678754e-21 code mass. The eleven-element auditor requires complete nonaliased mapping and relative small-return accuracy.
- 112 positive AGN receiver events. Requested/accepted energy each total 2.3502550743283316e-7 code energy; deferred energy zero. Maximum event energy residual: 8.271806125530277e-25. Nonzero cap deferral was not exercised by this evolution.
- 98304 effective thermal updates: 97163 heating and 1141 cooling. Net requested -1.0823046254344385e60 erg; net measured -1.0823046254344387e60 erg. Signed residual -1.5344344247033094e44 erg is about 1.4e-16 of the net magnitude. This aggregate does not bound every event's relative error or subsequent correction operators.

Authoritative artifacts are private RAR outputs at `/gpfs/kjhan/RAR_QUMOND/audits/rar_combined_1275122/`; the three audit scripts and three JSON result digests verify against artifacts.sha256. The project-root scripts and launch records preserve the exact audits. Earlier provenance files describing these modules as awaiting validation are superseded only within this fixture's scope by this record.

This evidence establishes these local operator exchanges in the declared fixture. Global mass/element/energy budgets, AGN donor withdrawal and sink generation, later velocity/floor corrections, MPI/AMR synchronization, restart/diagnostic parity, and calibrated galaxy formation remain unverified. It supplies no galaxy RAR comparison, viable relativistic BIMOND cosmology or FDM science result. The first matched Lambda-CDM/QUMOND prescribed-source comparison remains open.
