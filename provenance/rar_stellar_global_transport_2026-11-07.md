# Stellar gross-return transport: global compiled measurement

Source 7984e30 adds opt-in global inventories before thermal_feedback and after
set_uold following hydro reverse exchange at levelmin. The interval includes
hydro transport. In the uniform periodic fixture, density and passive-element
integrals can therefore be compared with accepted stellar material increments;
momentum and energy cannot be interpreted as pure source increments because
forces and hydro work occur in the same interval.

Grammar build 1275316 completed in 2m50s. Three-step one-/two-rank runs
1275317/1275318 completed in 23s/17s, and reconciliation 1275319 completed in
2s with exit 0. Its input and output SHA-256 manifests verify. The parser requires
NVAR=18, 32768 owned leaf cells, complete finite records, valid particle lists,
operator physical mass invariance and positive stellar-to-gas exchange. The
reconciler restricts accepted transactions to uniform level 5 and volume 1/32768.

Both runs have 328 accepted stellar-return transactions. Accepted stellar debit
is 2.5792576764650586e-5 code mass; measured global stellar losses are
2.579257676464986e-5 and 2.579257676465005e-5. Accepted gas mass return is
2.579257676465082e-5; post-transport global changes differ by
1.0390278725293592e-14 and -4.5977321071460214e-15 code mass, consistent with
subtracting totals near unity. Total-metal return is 7.391151587342369e-6;
its residuals are -4.235164736271502e-21 and -3.3881317890172014e-21.

All eleven tracked element fields also pass the declared summation/relative
criterion. Fields 7/8 have large initial backgrounds; their maximum absolute
residual is 2.2940955640140745e-15 code mass. Other tracked-element residuals
are at most 1.9058241313221758e-21. The test uses 256 ulp of the bounding field
scale plus 1e-8 of accepted return, not a fixed absolute allowance that could
hide tiny trace-element losses. Delayed-cooling clocks are excluded from the
element inventory comparison.

This establishes gross material delivery through the measured CPU MPI exchange
and transport path in this synthetic periodic uniform fixture. It does not
certify nuclear net-yield production, stellar/BH chemical reservoirs, arbitrary
AMR interfaces, global cosmic energy or production galaxy calibration. The
SSP/seed/artificial star-formation inputs remain exploratory, with dust/RT off.
The current first scientific comparison remains Lambda-CDM with a CDM halo
versus full-field baryon-only QUMOND using matched baryonic sources.
