# Native radioactive material/source connection — proposed bundle

Final goal: physically justified simulation-ready RT, stellar/AGN feedback
and dust. All medium groups are preapproved. Q-GOAL FIRST, Q-LEAN SECOND;
one Fable plan review, driver end evaluation. Do not create another source
survey framework: the existing LC18 adapter,333-isotope decay sensitivity
report and pinned NUBASE2020 are already present.

## Actual missing physics and bounded target

Current LC18 native source aggregates parent elements with no decay. It
therefore leaves Ni56 outside tracked Fe forever, unlike already-decayed
PARSEC/KL16 sources. LC18table8/9 supply actual fresh Al26,Fe60,Co60,Ni56,
Co56 inventories. Initial target: transport these five radioisotope mass
subsets with gas, evolve Al26->Mg26,Fe60->Co60->Ni60,Ni56->Co56->Fe56,
and change the actual gas Mg/Fe elemental fields, NOT just a diagnostic.
Use the same original source masses/Z/fates/channel ownership/IMF and
wind/terminal release times as the selected LC18 source. Do not attach
LC18 isotopes to PARSEC or double-decay already-stable yield tables.
This is a dominant-chain subset, not a claim to evolve all333nuclides.
Judge explicitly whether that bounded target usefully closes a real part
of medium group4 without becoming an instrumentation project.

Nuclear inputs already pinned at
`external/g2_candidates/nuclear_decay/nubase_4.mas20` SHA256
1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859:
half-lives Al26 717ky,Fe60 2.62My,Co60 5.2714y,Ni56 6.075d,Co56 77.236d.
IAEA LiveChart ground-state26Al independently agrees with717ky. These are
laboratory/ground-state rates, not temperature/ionization-dependent EC.
LC18 single Al26 output is treated as the long-lived ground inventory;
its nuclear network separately handles the short-lived isomer. Verify the
source interpretation before conversion; do not invent an isomer split.

Use the existing baryonic-mass convention: transmutation conserves A*m_u,
not atomic rest-mass. No extra rho or metal mass is created by tracer fields.
Ni/Co/Al remain part of the existing untracked-metal reservoir; daughterMg/Fe
move into tracked fields and Fe60 decay moves Fe out to untrackedCo/Ni.
MeV gamma/leptonic energy is OUTSIDE this abundance-only/transparent limit:
no new gas heat, CR source or photons inserted into the <=10keV SNRT groups.
Do not call this a decay-heating or ionization-rate implementation.
Decide in review whether this explicit approximation is adequate for the
abundance task; reject a falsely complete general-energy claim.

## Source timing is essential, not optional instrumentation

Advect five nonnegative mass-density SUBSETS through the existing passive
hydro path. Gas astration removes them in the same fraction as gas; they are
not additional stellar returned mass. This models gas abundances, not decay
emission inside newly formed stars. Use exact positive Bateman-chain updates
in physical elapsed seconds, including daughter buildup, stable limiting
forms for small dt and underflow-safe large dt. Do not force timesteps to
resolve Ni56 days in a galaxy run.

Do NOT simply inject fresh isotopes at timestep end then leave Ni56 intact
for a Myr step. Convolve each source interval with the analytic decay kernel:
piecewise-linear wind cumulative releases give constant emission rates on
their existing knot intervals; terminal ejecta occur at their exact lifetime.
Integrate survival and daughter production to the end of the same step.
Pre-existing gas inventory decays once per local level timestep, separately
from those already-aged new ejecta. No second decay, star-age clock in place
of time-since-ejection, or deletion of gas mass through a delayed source.

## Native wiring and supported combinations

Extend existing source table/history loading with an explicit matched
isotope companion; verify coordinates/age/channel and source identity before
publication. Reuse native IMF cell weights, not a separately normalized
population. Carry radioactive increments through the existing transactional
source into gas. Reserve five passive slots after the complete current
element window and bind layout/rates/source content in existing MPI/HDF5
identity. Preserve disabled defaults and old checkpoints when disabled.
Any new namelist selection must be wired in BOTH mkrun.py and
patch/cuRamses/aux/ramses_nml_generator.py/shared GUI. Preserve VPATH.

Initially admit gas-only, native LC18 channel-resolved stellar feedback.
Reject live dust/CHIMES with isotope evolution: Fe60 incorporated in grains
would transmute lattice composition; assuming unchanged olivine/Fe optics
would silently violate the current composition model. This is not new
approval gating; it is a named unsupported combination. SNRT fixed/HHe and
ordinary hydro/MHD can coexist if passive offsets/energy indices are correct.
Do not claim a new GPU decay kernel unless the actual native path needs it;
five exponentials per cell are not a reason for a dispatch framework.

## Minimum evidence

Reuse native source/bridge tests: analytic parent/daughter conservation,
zero dt, source-time convolution, one vs split steps, actual LC18 node/IMF
material change, source mismatch rejection, and disabled identity. One small
MPI2/OMP2 injection/evolution and exact restart with real source input; test
ordinary initial gas decay without star formation to isolate the analytic
law as part of the same run/fixture work, not another audit campaign.
Retain compact inputs/results, then delete evaluated raw outputs. No
publication claim for all333nuclides, lab EC at all ionizations, MeV feedback,
or isotope-resolved dust. Do not turn those limits into invented source data.
