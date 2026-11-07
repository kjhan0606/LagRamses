# RAR global physical mass inventory

Opt-in RAR_GLOBAL_INVENTORY at source 89ecd0f measures all conserved gas fields
on owned active leaf cells before/after the coarse AGN operator. It also counts
unique local DM/star/sink-cloud particle entries. Physical mass includes gas,
DM, stars, other particle mass, and the replicated sink catalogue counted once
on rank 1. Cloud masses are reported separately and excluded to avoid counting
both BHs and their quadrature particles. Invalid or duplicated particle-list
entries are exposed for rejection.

Grammar build 1275266 completed in 2m53s. Three-step synthetic QUMOND hydro runs
1275267/1275268 completed in 23s/18s. Global operator audit 1275270 passes both:
all 32768 leaf cells are represented, and physical total mass is unchanged by
each AGN call. Positive global gas-energy increments are measured at steps 2/3
(approximately 1.05012876169175e-7 and 1.30012631263659e-7 code energy).
These measurements alone are not a global energy-conservation claim.

Full coarse-stage mass reconciliation 1275271 uses the one-sink pre-feedback
ledger and explicit spin_bh=false,mad_jet=false fixture. Source
agn_accretion_receipt removes gross gas mass and retains (1-epsilon) times that
mass in the BH, with epsilon=0.1 here. Expected radiated code mass is therefore
accepted retained mass times epsilon/(1-epsilon). The inter-stage measured
physical-mass losses agree with this prediction in both rank decompositions.
At steps 1→2 the expected loss is 2.05504924853669e-8; measured one-/two-rank
losses are 2.055048664573178e-8 and 2.0550494639337558e-8. At steps 2→3 expected
loss is 2.4812636887875422e-8; measured values are 2.4812634169535386e-8 and
2.4812639054516694e-8. Maximum residual over all intervals is 5.84e-15 code
mass. The zero-accretion interval is also consistent with summation rounding.
Both input and artifact hash manifests verify for both audits.

This tests global physical mass, including star formation/return and BH mass
loss, over these three synthetic coarse intervals. It does not establish
stellar nuclear-element production closure, BH elemental inventories, cosmic
energy conservation, general AMR accuracy or production galaxy calibration.
The SSP bundle remains exploratory; dust and RT are disabled. The physical
mass baseline includes the deliberately added seed BH, approximately 1.000001
code mass, rather than assuming that the seed was taken from the initial gas.
Do not use these results to override the older full-run wrapper's missing
stellar/global-energy inventory requirements.
