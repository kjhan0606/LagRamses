# PARSEC low-Z common source: printed-precision reconstruction proposal

Final goal: physically justified simulation-ready RT/stellar-AGN feedback
and dust. Medium groups are preapproved; continue without another operator
approval wait. Q-GOAL FIRST and Q-LEAN SECOND; one Fable plan review, driver
end evaluation. This follows the completed five-Z bundle, not an invented
new finish-line gate. Main target is missing low-Z common high-mass material
and radiation, not full SSP/binaries or an unrelated restart campaign.

## Actual new findings, not a blanket roundoff assertion

Using Decimal80 arithmetic and half the final printed decimal unit of each
Mfin/Mbar field, the negative eleven-element wind/terminal residuals admit
consistent baryonic masses within their PRINTED intervals at every14--600
node except24Msun/Z=1e-6 (infeasible by3.197519368631e-8Msun). No isotope
mass has to be rescaled to show this feasibility. Largest required change
in Mfin/initial mass is3.3772606890714284e-6; Mbar shifts<=1.8433e-6 of
initial mass. This is a consistent representative inside printed intervals,
NOT proof of the author's exact hidden values or a Gaussian noise model.

Separately,150Msun/Z=.001 is labelled PPISN but M_He=64.586Msun is beyond
the existing W17/explicit62--64 bridge. Do not extrapolate/clamp its energy
or omit just that mass node. Neither issue can be erased by normalization.
Proposed explicit source grid therefore contains eleven complete Z branches:
1e-11, .0001, .002, .004, .006, .008, .01, .014, .017, .02, .03.
That is495nodes, within the existing512 cap. The omitted1e-6/.001 are not
author source nodes in this model; queries there use the declared same-age
linear-Z mixture of selected neighbors, not those omitted source records.
State this visibly in metadata/docs; no claim to reproduce their original
fates or the entire thirteen-Z catalogue. Judge whether this selection is
scientifically defensible for an explicit comparison population.

Do NOT normalize all41 tabulated nuclear species as a supposedly complete
inventory: minor unlisted isotopes exist. Trying a full-isotope-sum wind
mass authority, followed by terminal renormalization, changes a weak PPISN
terminal normalization by17.46% atZ=.02. Driver rejects that blanket method.
Preserve published eleven-element gross values; only minimally reconstruct
printed baryonic masses, with explicit interval bounds. Larger physical
source inconsistencies remain outside this numerical remedy.

## Proposed baryonic policy (new explicit model, old defaults untouched)

For each node, M0 is the exact source mass coordinate. W is the sum of the
eleven tracked wind masses; T is the sum of tracked total-minus-wind ejecta.
Choose F=Mfin and R=Mbar inside their original half-last-digit intervals,
requiring M0-F>=W, F-R>=T and nonnegative returns/remnant. Prefer unchanged
F; choose the closest feasible F and retain R when feasible, otherwise the
closest allowed R. FSN/DBH additionally require R=F and exact zero terminal
ejecta; PISN keeps R=0. Reject empty feasible intervals. Gross elements and
fates never change. Returned mass, net columns and wind energy follow the
selected baryonic budgets. Record original/reconstructed values and bounds
for every adjusted node. Do not call this recovery of an exact author model.

The current converter subtracts rounded track masses for wind timing. Tiny
primordial winds may fall below those mass digits. Inspect actual positive
RATE column28 before implementing; its trapezoidal integral in years can
supply an explicit source-based loss SHAPE, normalized to the selected wind
budget. Surface abundances still set relative phase composition; existing
positive endpoint/interval balancing and1e-4 cumulative compression apply.
Do not synthesize uniform release or zero wind if RATE is absent/zero while
the source yields are positive. Actual track M/L/Teff remain untouched.
Judge whether RATE should be the new model's common shape, or only a
declared fallback where mass subtraction is unresolved. Avoid duplicate
pipelines. Fixed and phase-escape speeds are existing named comparisons;
the full selected grid may use fixed if any phase case fails its domain,
with no Gamma cap or invented velocity. All old source models remain intact.

Handle only arithmetic-scale negative residuals introduced by converting
the selected Decimal budgets to binary floats, with explicit compensated
sum/outward-rounding conventions and tests. Never use a tolerance to admit
the known physically unselected1e-6 node. A source-specific model ID binds
this new mass/shape policy to its matched SED; merely relabelling it as the
old published-mass model is forbidden.

## Native connection and lean evidence

Reuse the two existing converters and structural native reader/binder,
shared Z bracket, same per-source mass cells and full .08--600 IMF.
Precheck actual per-node Q5/track timing, Planck feasibility and energy
domain; do not relax those contracts to publish an invalid full grid.
Every material/radiation coordinate, fate, lifetime and model must bind.
No new carrier, RAMSES NML key, MPI format, source atlas framework or
capacity expansion is necessary. Preserve Makefile VPATH. Update both
namelist generators only if an actual new runtime NML key is introduced.

Extend existing native material/radiation fixtures, one actual low-Z
interior query and precise unchanged-node checks. One small MPI2/OMP2
integration/restart with this input model, not one run per Z or IMF.
Report its full effective NML/output budget before launch; after driver
evaluation remove exact raw HDF5 and retain compact evidence and inputs.
No extra external end audit, simulation matrix or per-helper approval.

## Acquisition/precheck update while the single review runs

All proposed added track/photon ZIPs are now staged. Source RATE is positive
at all inspected45nodes/Z, including nine Z=1e-11 nodes14--30Msun whose
printed track MASS differences are exactly zero despite positive winds.
RATE is therefore necessary there. At other nodes the RATE integral/mass
difference can differ appreciably (minimum .8533 atZ=.002); do not replace
the resolved mass-based shape everywhere without acknowledging that change.

Two photon files extend31.4861yr (270Msun/Z=.004) and27.3327yr
(70Msun/Z=.006) beyond the corresponding track deaths, fractional1.33e-5
and7.53e-6. Existing two-Z converter required photon end<track death. These
are measured-time overlaps, not missing spectra: a possible bounded repair
is to truncate their cumulative/rate integral at the actual feedback death,
interpolating positive Q/E between bracketing measured-time rows and never
publishing photons after death. Shorter files retain the existing distinct
full-track Planck tail. Record discarded overlap and do not extrapolate.
This is not a claim that every raw radiation time grid has the same endpoint.

At Z=1e-11, projected surface fractions exceed unity by at most9.821e-12:
major H/He/CNO track fractions have10digits and erase some trace fraction,
while the separately supplied initial heavy fractions are still retained.
Any repair must be confined to the actual printed precision of the major
fraction, preserving trace species and all published endpoint yields. Do
not impose a broad abundance floor or a physical total-Z normalization.
