# AGN source spectrum closure — bounded bundle plan (draft)

Plan initiated 2026-10-31; revised after Astra review 2026-11-02.
Workspace: `/gpfs/kjhan/chimes-cvode-gpu-20261001`,
`kjhan0606/LagRamses`.

## Q-GOAL — ask first

Does adding a literature-based AGN emission shape, comparing its
moment-matched interaction rates to the actual runtime reconstruction, and
coupling it through existing SNRT/gas/dust paths advance the approved
galaxy-scale effective-model goal? Identify scope mismatch before reviewing
implementation details.

## Q-LEAN — ask second

Is the proposed bundle lean enough for that goal? Flag unnecessary diagnostics,
duplicate gates, new test infrastructure, or any work better deferred. The
target is one source-spectrum closure and one bounded coupled confirmation.

## Current state

The successful M5 + stellar/thermal-AGN + dust + CHIMES integrations verify
source routing and coupling under reference controls. They do not validate the
AGN SED. `snrt_group_contract_reference_control_v1.nml` is explicitly marked
`reference_control_not_science_approved`; the nine-group contract is bound to a
pilot ledger, and its represented energy fraction is about 0.1815. The runtime
transports group photon-number and energy moments and reconstructs a positive
within-group spectrum with the shared MaxEnt prior in
[`snrt_band_spectrum.h`](../patch/lagRamses/snrt_band_spectrum.h); CHIMES
interaction rates integrate that reconstructed spectrum in
[`snrt_chimes_spectrum.cpp`](../patch/lagRamses/snrt_chimes_spectrum.cpp).
Thus the relevant comparison is against the actual moment-matched runtime
reconstruction, not against a fixed group coefficient alone.

### Code-path constraints confirmed before implementation

- The source builders can emit source-specific group photon and energy
  moments, but the material interaction path reconstructs a spectrum from the
  transported moments using one shared MaxEnt prior. Candidate SED work must
  compare rates produced by that exact runtime reconstruction against the
  resolved candidate spectrum, with matched group moments. Merely replacing a
  shared coefficient or splitting STAR and AGN labels does not measure or
  remove within-family spectral variation or attenuation-induced hardening.
- `candidate_explicit_sed` currently passes file/value inspection but is not
  runtime-admissible. The only runtime statuses are `approved_production` and
  `reference_control` with `SNRT_ALLOW_REFERENCE_CONTROL=1`; the current
  `partition_reference_v1` namelist path additionally insists on
  `reference_control`. A candidate coupled smoke would require a deliberately
  separate candidate-only admission path, with a fail-closed opt-in and no
  route to production admission.
- The current release-candidate manifest still selects the reference AGN
  group contract and sets `launch_admission_manifest=false`. The latest
  reference-control integration therefore proves wiring only; it is not
  physical-SED or production qualification evidence.

## Astra plan-audit result (read-only)

**Q-GOAL — conditional yes.** A literature-based AGN emission shape followed
through M5, gas chemistry, dust absorption, and the existing energy ledger
supports the effective galaxy-formation goal. A fixed group-coefficient swap
alone would not establish this because the actual runtime uses the moment-based
MaxEnt reconstruction.

**Q-LEAN — conditional.** Keep the initial bundle to one bounded source
contract, one comparison against the runtime reconstruction, and one coupled
confirmation. Reuse the existing runner and conservation evidence. Add
source-resolved mesh state only if the measured shared reconstruction misses a
declared accuracy target and source mixing is shown to cause that miss. Such a
split would not by itself fix within-family variation or spectral hardening.

The audit recommends Su et al. `agnSED` as a shape-only candidate for the first
bounded fixture. Preserve LagRamses's accepted electromagnetic-energy receipt
and mechanical partition; do not derive a second bolometric normalization from
Su and add it to the existing feedback budget. The sink currently carries a
scalar pending-energy reservoir rather than a per-event spectral history, so
the first candidate must declare a fixed template evaluated at consumption
time. Event-weighted spectra would require explicit reservoir, merge, and
restart changes and are not part of this bundle.

For evidence, compare the resolved candidate spectrum with the actual runtime
MaxEnt reconstruction at matched group photon and energy moments for pure
stellar, pure AGN, mixed, and representative attenuated spectra. Summarize the
materially contributing photoionization/heating rates, CHIMES rates, dust
absorbed energy and momentum, and transmitted photon number and energy. Astra
proposed a 10% relative engineering target for materially contributing rates
and transfers, with an absolute contribution floor for near-zero channels and
integration error smaller than that target. The operator explicitly accepted
it on 2026-11-02 after the corrected comparison report had been produced; its
existing results are evaluated retrospectively against the criterion below.
It remains an engineering screen, not a universal physical guarantee.

The coupled candidate run must demonstrate an in-domain black hole and
accepted electromagnetic receipt, distinct stellar source moments reaching
M5 alongside AGN photons, consistent gas/chemistry and dust/IR transfers, and
no duplicated mechanical or absorbed/reprocessed energy. Candidate admission
must fail closed in both AGN and spectral-model loaders and remain separate
from production admission. The audit does not approve production use or clear
external data-use and release-manifest blockers.

## Driver-defined comparison criterion — 2026-11-02

The operator explicitly accepted Astra's proposed **10% engineering
screening target** on 2026-11-02. The corrected comparison report was already
present at 16:27 KST, while this explicit acceptance arrived at about 18:11
KST. Its existing results are therefore being evaluated retrospectively
against the accepted screen; this was not a preregistered blind comparison.
The target is not a physical-accuracy guarantee, calibration result, or
production approval.

For a process with resolved-spectrum reference contributions `R_j`, define
`R_total = sum_j abs(R_j)` across its configured groups/channels. A channel is
material when `max(abs(R_resolved_j), abs(R_runtime_j)) >= 1e-3 * R_total`;
each material channel must have relative discrepancy at most 10%. For a
sub-floor channel, require `abs(R_resolved_j - R_runtime_j) <= 1e-4 *
R_total`, the absolute equivalent of applying the same 10% tolerance at the
contribution floor. Apply this to the declared ionization/heating rates and
group-integrated gas/dust energy and momentum transfers; report signed
quantities with both absolute and budget-normalized residuals. The actual
runtime MaxEnt reconstruction receives exactly the same group photon-number
and energy moments as the resolved candidate SED.

The user-approved numerical integration error limit is strictly below 10%.
For this comparison, the calculation used a stricter working target of at most
1% for material channels and `1e-5 * R_total` for sub-floor channels. The
convergence estimate must come from an independent quadrature refinement, not
from comparison of two representations alone. If the candidate fails the 10%
screen, record the failure and identify which source/rate channels dominate;
do not weaken the threshold after seeing results. Source-resolved mesh state
remains conditional on demonstrating that shared-source reconstruction is the
cause.

### Candidate-data audit (2026-10-31)

The repository does not yet contain a physical AGN SED candidate for this
bundle. `p4_explicit_agn_photon_ledger.json` is explicitly marked
`synthetic_non_physical_wiring_fixture`; its ten-point CSV assigns a constant
energy fraction per eV over 0.01--10,000 eV and therefore places about 80% of
bolometric power in the 2--10 keV group. It is an engineering fixture, not a
plausible adopted AGN spectrum. In `run_snrt_native_spectral_contract.sh`, the
loader's `candidate_explicit_sed` case is made by copying the reference
namelist and changing only its status string. Passing that case validates
status parsing, not a physical candidate or its spectrum.

The other checked-in input, `p4_pilot_agn_photon_ledger.json`, is a
parameterized Sazonov--Ostriker--Sunyaev-style reference control with support
only from 10 eV to 10 keV. It represents 0.18148066885 of bolometric power;
the lower two configured groups are zero and the group crossing 10 eV is only
partially supported. Its own metadata requires Lyman normalization and escape
fraction variation before a production interpretation. The control must not
be promoted by changing its approval string.

This leaves the source spectrum itself as a prerequisite to quantifying the
shared-kernel error: the synthetic flat fixture cannot represent either the
selected high-redshift AGN model or its intended parameter domain. The
candidate derivation must state its black-hole mass/accretion-rate domain and
whether the transported component is intrinsic or escaped, with dust
reprocessing accounted for exactly once. The Obelisk precedent integrates
SED-dependent interaction moments but omits X-rays and assigns 30% of the
bolometric disc luminosity to dust-reprocessed IR that it does not transport;
it therefore cannot supply our full 0.01 eV--10 keV plus live-dust closure
unchanged. The Su et al. accretion-state model depends on black-hole mass and
Eddington-scaled accretion rate, and its own cosmological comparison reports
high-redshift population discrepancies; it is not an automatic LRD-valid
template. These source-domain limits must be carried into the Q-GOAL/Q-LEAN
review.

The authors' [`agnSED` implementation and interpolation data](https://github.com/SuTong1999/agnSED/tree/200dbbf4e606bbf4a7d6d72c7883951940c8f57a)
are available upstream at commit `200dbbf4e606bbf4a7d6d72c7883951940c8f57a`
under the [MIT license](https://github.com/SuTong1999/agnSED/blob/200dbbf4e606bbf4a7d6d72c7883951940c8f57a/LICENSE). The repository README
says its `data/` directory was restored on 2026-03-31. Its API calculates
SEDs from black-hole mass and Eddington-scaled accretion rate, switching
between ADAF/truncated-disk and disk-corona branches at the published critical
rate ([Su et al.](https://arxiv.org/abs/2501.10793)); the paper's model does
not include spin. The package's numerical interpolation axes are branch
specific; they are not the same thing as the black-hole population selected
for the paper's cosmological analysis. The package interpolation data cover
log frequency 8.0--22.9, which contains the current transported 0.01 eV--10
keV interval. However, the package loads pickled SciPy
`RegularGridInterpolator` objects and pins SciPy to <=1.14.1 because of a
documented compatibility issue. If selected, it should be used only by a
pinned, audited offline table generator; the simulation runtime should read
plain, hash-bound numerical arrays instead of importing this Python/pickle
stack. The upstream model is a traceable candidate source, not a
production-approved SED or proof of shared-opacity accuracy. Its high-redshift
LRD applicability still requires review.

### LagRamses-to-model accretion contract audit (2026-10-31)

The pinned `agnSED` API accepts `log10(M_BH/Msun)` and `log10(mdot)` and
returns tabulated `log10(L_nu)`; its mode switch is at `log10(mdot)=-2.9`.
The paper defines `mdot = Mdot_acc/Mdot_Edd` with
`Mdot_Edd=L_Edd/(0.1 c^2)` and assumes zero SMBH spin. Direct inspection of
the pinned interpolation objects gives these numerical axes:

| Component | `log10(M_BH/Msun)` | `log10(mdot)` |
| --- | ---: | ---: |
| ADAF | 0.8--11.5 | -12.2-- -2.9 |
| Truncated disk | 0.8--11.0 | -3.0-- -2.9 |
| Disk-corona | 0.8--11.5 | -3.0--3.0 |
| Standard thin disk | 3.8--11.0 | -7.0--3.2 |

The public `SEDs` helper defaults to a sample cut `mcut=4`, whereas the paper's
cosmological population plots exclude SMBHs below about `10^6 Msun` because of
simulation resolution. Neither selection is a statement that every point in
the wider interpolation arrays has been observationally validated. The
transported 0.01 eV--10 keV band lies within the table's log-frequency support
8.0--22.9. Runtime evaluation must check the selected component's actual axes
and reject nonfinite or out-of-domain inputs; the public wrapper's conversion
of NaNs to `-100` must not be used to hide an invalid interpolation.

The retained job-444355 wiring fixture has `M_BH=2.4527e11 Msun`
(`log10 M_BH=11.39`). That is not outside every numerical interpolation axis,
as an earlier draft incorrectly stated, but it is outside some component
axes. Interpolation support alone does not establish physical validity at
this mass, and this fixture's applicability to the paper's population has
not been established. Therefore it is not a justified candidate fixture. The
bounded candidate comparison uses
`M_BH=1e8 Msun`, `mdot=0.1` (`log10 mdot=-1`), which lies within the
disk-corona interpolation domain and the paper's modeled SMBH/accretion
regime. It must not be extrapolated.

LagRamses does not use that same Eddington-rate convention when spin is active:
`dM_Edd/dt = 4 pi G m_p M_BH/(sigma_T c epsilon_r)`, where `epsilon_r` is
spin-dependent (and defaults to 0.1 without spin). The coarse-state record
already carries the BH mass, Bondi and Eddington rates, capped inflow rate,
base/effective efficiencies, and code bolometric luminosity. Thus the Su-model
input is derivable from audited fields as
`mdot_Su = Mdot_inflow/[L_Edd/(0.1 c^2)]`; if expressed using the code's
uncapped ratio `X=Bondi/Mdot_Edd(code)`, it is
`min(X,1) * 0.1/epsilon_r`. Feeding `X` directly is equivalent only when
`epsilon_r=0.1` and the cap has not changed the inflow. The existing effective
efficiency can also include the MAD jet-mode reduction, while Su et al.'s SED
assumes zero spin and does not encode that reduction.

**Selected first-bundle normalization:** use the Su interpolator for spectral
shape only. Define the spectral-fraction denominator before clipping to the
transport bands: use the model-reported bolometric luminosity only if the
candidate data/code exposes it and its consistency with the spectrum can be
checked. Otherwise use the integral over the complete supported SED grid,
label it `L_grid` rather than `L_bol`, compare it with a paper-reported
bolometric value when available, and carry the unsupported low/high-energy
tails as a normalization uncertainty. For each group,
`f_g = integral_group(L_nu dnu) / L_denominator` and
`E_g = f_g * E_EM,accepted`; then `N_g = E_g / <h nu>_g`. If the denominator
is a supported-grid integral, report the in-band fraction of that grid and do
not call its omitted fraction a bolometric remainder. If it is a checked
bolometric value, report the true out-of-band remainder explicitly. Never
renormalize only the sampled in-band groups to sum to one. Here
`E_EM,accepted` is the energy already accepted into `agn_pending_erg`; do not
count it again as thermal/mechanical feedback. The existing high-accretion
branch allocates 85% of released energy to electromagnetic feedback and 15%
to mechanical feedback; the low branch is jet-only and emits no AGN photons.
The candidate must preserve that partition.

Compute `mdot_Su` from the accepted/capped inflow rate, not the requested
Bondi rate or the code's spin-dependent Eddington ratio. The first fixture is
`M_BH=1e8 Msun`, `mdot_Su=0.1`, zero spin/MAD, and must verify that the actual
`X_floor` selects the electromagnetic branch. Any source-side escape or
obscuration is explicit. For this bounded fixture, set the unresolved source
escape fraction to `f_escape=1`: inject the selected intrinsic SED at the
resolved transport boundary without an unresolved torus/obscurer correction.
This is a declared diagnostic assumption, not a general AGN escape model.
Resolved gas and dust then process the injected radiation once.

### Su-template to RAMSES energy-scale cross-check — 2026-11-02

For the proposed bounded mapping, the complete supported-grid integral is
`L_grid=1.020811918e45 erg/s` at `M_BH=1e8 Msun` and `mdot_Su=0.1`.
With spin disabled, the live accretion path uses `epsilon_r=0.1`; Su's
definition `Mdot_Edd=L_Edd/(0.1 c^2)` therefore makes the pre-partition
radiative power `0.1 L_Edd` at this fixture. The high-state reference
partition sends 85% of that receipt to `agn_pending_erg`, so
`E_EM,accepted/dt=0.085 L_Edd`. Using the constants in
`sink_particle.kjhan.f90` gives `L_Edd=1.246899043e46 erg/s` and
`E_EM,accepted/dt=1.059864186e45 erg/s`, 3.83% above `L_grid`. This is within
the operator-approved 10% engineering comparison target for this exact
fixture. It is an amplitude-consistency cross-check, not a bolometric closure
or a bound on the unmodeled SED tails.

Applying the unrenormalized represented fraction `0.9732287970` to
`E_EM,accepted` injects `1.031490347e45 erg/s`; the remaining
`2.67712030%` of `E_EM,accepted` is recorded as unrepresented. The resulting
represented injection is 1.05% above `L_grid`. Neither result measures
source-side escape. The proposed diagnostic convention remains
`f_escape=1` at the resolved transport boundary, with no unresolved obscurer.
It is recorded in the candidate source-boundary contract. The build generator
checks the sidecar and its moment-ledger hash, compiles its source identity into
the group profile, and the Fortran loader rejects a runtime contract with a
different source ID, ledger hash, commit, status, approval marker, or fraction
semantics. The profile remains runtime-inadmissible, and this arithmetic check
does not measure physical source escape or establish bolometric closure.

### Corrected bounded Su AGN moment-closure result (2026-11-02)

The preliminary calculation is superseded. Astra found that it used 64 nodes
for a dust/CHIMES path that uses 128 and assigned a spurious He I contribution
to the group ending exactly at 24.59 eV. The corrected calculation applies the
runtime's left-limit upper endpoint and repeats both applicable node counts.
He I ionization in group 5 (13.6--24.59 eV) is exactly zero in both the
resolved integral and runtime operator; do not use the earlier nonzero
residual.

The fixture is `log10(M_BH/Msun)=8`, `log10(mdot)=-1`, disk-corona, zero spin,
and no unresolved obscuration. Pinned SciPy interpolators were evaluated
directly, bypassing the public wrapper's NaN-to-`-100` masking. The model's
`log10(L_nu)` was interpolated linearly in `log10(nu)`, matching the upstream
helper. `L_grid = integral(L_nu dnu)` over its complete supported frequency
grid `log10(nu/Hz)=8.0..22.9` is `1.020811921e45 erg/s`; this is a shape
denominator, not a bolometric luminosity. Independent 40,000/80,000-interval
refinement changes `L_grid` by `5.6e-11` relative. The nine transport groups
contain `0.9732287941` of `L_grid`; fractions were not renormalized. The
omitted `0.0267712059` is the unrepresented fraction under this `L_grid`
normalization, not a claimed bolometric tail or measured physical escape. The
later M_N source-receipt implementation records the accepted input,
group-injected and absolute unrepresented energies and checks their arithmetic
closure before clearing the pending receipt. It does not reassign the omitted
energy to mechanical heat or dust.

| Group (eV) | `L_group/L_grid` | Photon-weighted mean (eV) |
| --- | ---: | ---: |
| 0.01--1 | 0.01611252680 | 0.4810330 |
| 1--5.6 | 0.13949771354 | 2.7645220 |
| 5.6--11.2 | 0.16262271554 | 7.9838109 |
| 11.2--13.6 | 0.05741634377 | 12.3348930 |
| 13.6--24.59 | 0.17327787945 | 17.8852513 |
| 24.59--54.42 | 0.14086657979 | 33.8069843 |
| 54.42--500 | 0.18354965189 | 120.8728798 |
| 500--2000 | 0.06073922181 | 871.4394326 |
| 2000--10000 | 0.03914616149 | 3723.1961018 |

The resolved reference integrates `L_nu/(h nu) * sigma(E)` and its
primary-electron excess-energy counterpart. Runtime receives the same group
photon-number and energy moments and uses the exact exponential-tilt MaxEnt
weights, trapezoidal prior and energy nodes from
[`snrt_band_spectrum.h`](../patch/lagRamses/snrt_band_spectrum.h). The 64-node
path is the dust-free direct H/He kernel. The 128-node path is the applicable
D03-dust and CHIMES grid; the H I, He I and He II shell cross sections in the
pinned CHIMES bank match the same Verner kernel at all 128 nodes to at most
`1.5e-15` relative. Thus the 128-node H/He comparison uses the coupled path's
actual nodes and opacity rather than a 64-node proxy.

For each process, `R_total` is its resolved contribution summed over all nine
groups. A group is material when `max(|R_resolved|,|R_runtime|) >= 1e-3
R_total`; a subfloor group must satisfy `|R_resolved-R_runtime| <= 1e-4
R_total`.

| Process | 64-node max material error | 128-node max material error | 128-node max subfloor `|delta|/R_total` | Result |
| --- | ---: | ---: | ---: | --- |
| H I ionization | 19.61% (group 7) | 19.69% (group 7) | `5.24e-9` | Fail |
| H I primary excess energy | 12.93% (group 7) | 12.97% (group 7) | `5.29e-7` | Fail |
| He I ionization | 16.70% (group 7) | 16.75% (group 7) | `2.53e-7` | Fail |
| He I primary excess energy | 7.36% (group 7) | 7.36% (group 7) | `1.18e-5` | Pass |
| He II ionization | 18.31% (group 7) | 18.38% (group 7) | `4.83e-6` | Fail |
| He II primary excess energy | 5.87% (group 7) | 5.95% (group 7) | `1.19e-6` | Pass |

The independent 40,000-to-80,000 integration refinement changes material
resolved rates by at most `1.57e-5` relative and subfloor rates by at most
`2.46e-13 R_total`, below the comparison's stricter 1% working target and
`1e-5 R_total` subfloor target, and below the user-approved 10% integration
error ceiling. This estimates quadrature error for the pinned, interpolated
SED; it does not estimate uncertainty in the upstream SED interpolation.

All paired group contributions, normalized by each process's resolved total,
are in [`su_agn_sed_shape_comparison_2026-11-02.md`](su_agn_sed_shape_comparison_2026-11-02.md).
The 128-node result confirms a within-spectrum closure miss in group 7
(54.42--500 eV), not a demonstrated source-mixing error. Moving from 64 to 128
nodes does not cure it. The 10% screen fails for H I, He I and He II ionization
and H I primary excess energy. The corrected group-5 He I residual is zero;
all subfloor channels pass.

At the time this comparison was completed, no coupled candidate run or syntax
GPU build had been submitted. The nine-group closure fails the now-accepted
screen; the follow-up four-band result supports evaluating an isolated
12-group candidate. This comparison does not approve Su et al. for LRDs,
calibrate absolute luminosity, change the mechanical/EM partition, or clear
other production blockers. See the compact calculation record for methods,
paired group contributions and input hashes.

### Group-7 closure sensitivity (2026-11-02)

The bounded follow-up found that the power-law photon-number prior
`p_j proportional to w_j E_j^alpha` is not a general replacement: although it
reduced the selected unattenuated `logmdot=-1` maximum material error to
0.478%, it left one 10.523% failure in the low-accretion ADAF case and six
large failures after a neutral-H screen. One geometric split of group 7 also
left three ionization-rate failures after that screen.

The bounded diagnostic subdivided `[54.42,500)` eV into four equal-log
sub-bands, with edges `54.42, 94.746113, 164.954539, 287.188561, 500 eV`,
and applied the existing 128-node MaxEnt closure. It passed the 10% material
and `1e-4 R_total` subfloor limits for three sampled accretion states and the
same neutral-H screen. This changes the canonical count from 9 to 12 groups,
which is 33.3% more group-indexed state and loop work by count; actual runtime
and total-memory impacts are unmeasured. Adding an `E^2` moment alone did not
resolve the screened mismatch for the tested quadratic-exponential family.
Detailed evidence is recorded in
[`su_agn_sed_shape_comparison_2026-11-02.md`](su_agn_sed_shape_comparison_2026-11-02.md).

These offline results support four-bin group-7 subdivision as the current
runtime candidate, not as production-approved physics. They do not establish
accuracy for arbitrary columns, full CHIMES response, dust coupling, or
source families outside the sampled cases. Astra independently reproduced
the four-case errors and supports advancing to migration design. Actual
memory/runtime cost remains unmeasured.

### Versioned group-contract migration status (2026-11-02)

The operator explicitly accepted the 10% relative error screen for material
reaction and transfer channels on 2026-11-02; the already-generated comparison
is evaluated retrospectively. Integration error remains separately bounded
below that screen, and near-zero channels use an absolute contribution floor.
The four-way logarithmic split of 54.42--500 eV is being implemented as an
isolated, nonproduction candidate profile. The nine-group profile remains the
default and the only profile currently admissible for runtime use.

The working tree contains a generated profile registry, profile-selected
Fortran/C dimensions, profile-bound PARSEC v3 records, variable-length Fe/HDF5
identity handling, and a 12-group source builder. The 12-group profile is
explicitly marked `runtime_admissible=false`; `mkrun.py` and the RAMSES driver
reject an attempted candidate run.

Syntax Slurm job 447576 completed successfully and built both
`snrt_agn9_v1` and `snrt_agn12_g7log4_v1`, including link checks. This was a
CUDA-enabled RAMSES build (`USE_CUDA=1`) with the CPU CHIMES backend
(`CHIMES_CUDA=0`); it proves compilation/linkage only. The profile-stamp
Makefile edit made after that job was syntax-checked with `make -n`, but was
not part of the successful full build. Neither binary was run as a simulation.

The resolved 12-group Su AGN source ledger and profile-bound CHIMES atomic and
molecular candidate banks have since been generated and passed a C++ runtime
loader probe; this closes only the source-input and loader subgate. The
accepted-screen material response comparison, profile-matched Draine/Fe
tables, dust/material and restart fixtures, and integrated M5+CHIMES coupling
evidence remain outstanding. Do not launch the 12-group candidate as a
simulation until those checks pass.

The following was the original migration design. Current partial completion
and blockers are recorded after the checklist below; compilation does not
mean that all data or physical interfaces have migrated.

At planning time, `snrt_ngroups=9` and `chimes_ng=9` were fixed compile-time
constants, and dust-enabled CPU transport rejected other group counts. The
following was the recommended bounded migration design; refer to the current
status above for what has since been implemented and what remains blocked.

1. Keep `snrt_agn9_v1` as the existing reference build and input contract.
   Add a build variable `SNRT_GROUP_SET` (default `snrt_agn9_v1`) and a
   separately compiled `snrt_agn12_g7log4_v1` candidate profile. It replaces
   `[54.42,500] eV` with exact edges `[54.42, 94.746113497183657,
   164.95453919186340, 287.18856104645209, 500] eV`. Two isolated builds
   avoid imposing the extra group storage on nine-group runs. Use one
   authoritative versioned profile manifest to generate the Fortran and
   C/C++/CUDA constants and feed the table builders. Separate object/build
   directories by profile so stale objects cannot mix. Update `mkrun.py`,
   `ramses_nml_generator.py`, and the GUI selection to carry the profile ID;
   preserve the established Makefile VPATH precedence. A namelist/manifest
   field must identify the selected profile and fail if it disagrees with the
   executable. Absence is accepted only as the current nine-group legacy
   profile; the 12-group candidate must name its profile explicitly. This is
   not dynamic switching inside one binary; that would
   require a broader allocation/interface refactor.
2. Bind the profile ID, 13 exact edges, and SHA-256 in source, dust, CHIMES
   atomic/molecular tables, and run metadata. Regenerate both active PARSEC
   stellar and AGN moments from resolved spectra; never split parent N/E
   moments after the fact. Preserve the PARSEC Q-constrained spectral prior.
   Update `build_parsec_native_sed.py` to emit a new 12-group record version:
   Q slots `1:12`, E slots `13:24`, and the existing total/remainder fields
   `25:27`. Derive offsets from the shared profile and reject a record whose
   explicit profile/hash does not match. The raw legacy PARSEC interval file
   has no edge hash; admit it only with a recognized legacy file version and
   an enclosing stellar contract whose existing edge hash matches
   `snrt_agn9_v1`. A 12-group record requires explicit ID/hash. Candidate
   admission in the stellar loader must allow an explicit `science_candidate`
   source with either the nine-group comparison profile or the 12-group
   candidate profile; it must not relax production approval generally.
   Rebuild D03/Fe dust and CHIMES atomic/molecular tables against the same
   edges and 128-node grids.
3. Make the common profile drive every fixed group count: Fortran
   `snrt_ngroups`, CHIMES Fortran `chimes_ng`, C/C++/CUDA `NG`/edge arrays,
   spectral-bank dimensions, and table builders. Generalize the CHIMES
   spectrum, photo-integration and C-bridge constants. For `NG=12`, the
   bridge identity has `32*(NG+1)=416` doubles plus four metadata values;
   update `chimes_live_identity`, its 324-value output to 420, and every
   writer/reader consumer and offset, including HDF5 backup identity.
   Generalize both `ngroup==9` guards in
   `snrt_openmp_transport.cpp` (dust and H/He column paths). Keep CPU/OpenMP
   and CUDA shapes, group IDs, edge hashes, opacity tables, and left-limit
   endpoint behavior identical.
4. Preserve both restart paths. Stream checkpoints already record group
   count and edge SHA-256 along with source, closure, status, and material
   table identities; M5 HDF5 restart records payload width and contract
   identities. Retain every existing check and add/verify `group_set_id` in
   new HDF5 metadata. An old HDF5 file without this attribute maps to the
   nine-group profile only when its payload width and all existing source,
   closure, status, table, and nine-group edge-hash identities match;
   otherwise reject it. Reject cross-profile restarts in
   both paths. A nine-group checkpoint runs only with `snrt_agn9_v1`; a
   12-group checkpoint runs only with `snrt_agn12_g7log4_v1`. No
   parent-to-subband conversion is allowed.
5. After both profiles build and their spectral/opacity contract checks pass,
   measure two matched pairs: 9 vs 12 groups on CPU-only `grammar/normal`,
   then 9 vs 12 on the same typed GPU partition on `syntax`. For both group
   counts, generate moments from the same resolved candidate AGN and PARSEC
   spectra. The nine-group comparison uses a separate `science_candidate`
   table with the current nine-group profile's ten edge values; it neither
   uses nor replaces the historical reference SED. Within each pair, hold those resolved source
   spectra, IC/run input, mesh, integrated duration, tolerances, CPU/MPI/
   OpenMP allocation, and (for GPU) card type/count fixed.
   Record wall time, peak host RSS, and peak device memory. The 12/9 count
   predicts 33.3% more group-proportional work/state; it is not a total
   memory or wall-time measurement. Preserve effective namelists, logs and
   compact evidence, then remove only evaluated raw outputs.

The checklist does not authorize runtime admission by itself. Keep the
candidate profile blocked until the listed physical data and coupling inputs
are present and match the profile identity. The adopted 10% screen is a
numerical engineering target, not a claim that the candidate SED or its
physical source model is production approved.

### Current migration gate and audit update — 2026-11-02

The explicit 10% acceptance followed the existing offline report, so its
pass/fail labels are retrospective. Astra's read-only audit found the sampled
H/He classifications consistent with the adopted rule. A persistent probe now
reproduces profile-bound HDF5 loading and atomic moments through the production
C++ loader, but it does not reproduce the offline source-vs-MaxEnt rate table
or establish full CHIMES response. Dust absorption/momentum, stellar/AGN
mixing, and transport-generated attenuation are also outside that evidence.

The successful 447576 build does not qualify candidate material tables. Dust
and Fe compatibility paths currently pad old nine-group data into twelve
slots; because group 7 was subdivided, old high-energy group indices no longer
map to the same energy bands. CHIMES molecular phase coupling also rejects a
group count that differs from its nine-group Draine table. Those paths must
remain blocked; changing `runtime_admissible` alone would be unsafe.

The source-to-CHIMES-input portion of the lean gate is complete. The 12-group
Su moments were integrated from the pinned source, both CHIMES banks carry the
same profile ID/edge digest, and the production C++ loader accepted them and
closed reconstructed N/E. Profile-bound D03 and Fe electric-plus-eddy tables
have now been generated against the same exact edges and 128-node grids. The
Fortran build selects uniquely named 12-group includes under the candidate
macro, avoiding source-local nine-group include precedence; both profile
Makefile dry-runs select their intended inputs and preserve the exact `VPATH`
ordering. A direct GNU Fortran candidate module compile/probe passed dimensions,
edge binding and finite-data checks. Fe remains comparison-only (magnetic spin
absorption is missing). Its table is profile-sized, but the current Fe UV
operator rejects nonzero populations in groups 5--12, so the full AGN
candidate spectrum cannot yet use that pathway. A full CUDA-enabled syntax
build was not completed: H100 was blocked by `QOSGrpGRES`, A100 by priority,
and A40 had no free card. All requests were canceled before execution; full
Intel/MPI/CUDA linkage remains unconfirmed.
Hashes and limitations for this source/CHIMES/D03/Fe bundle are recorded in
[`su_agn12_profile_loader_evidence_2026-11-02.md`](su_agn12_profile_loader_evidence_2026-11-02.md)
and the 22-file asset manifest.

At this point in the migration log, the full 10% response/transfer comparison
was still pending. The later bounded atomic and D03 one-crossing results are
recorded in the dated updates below; they do not cover full molecular CHIMES,
Fe absorption/momentum, or coupled transport. The molecular table's coarse
line projection and Habing-shape approximation remain limitations. No
coupled candidate simulation or performance matrix should start until the
remaining profile-compatible inputs and source receipt are implemented. A
coupled check must close `E_EM = E_injected + E_unrepresented` before clearing
the pending source receipt.

### First-pass stellar interaction-moment sensitivity (2026-10-31)

The independent BPASS v2.2.1 radiation-reference source was integrated against
the [Verner et al. (1996) ground-state H/He photoionization fits](https://arxiv.org/abs/astro-ph/9601009)
and the authors' [fit-parameter table](https://www.pa.uky.edu/~verner/dima/photo.dat).
The checked-in source record binds its HDF5 input to SHA256
`b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3`.
For each of its 51 tabulated ages and 13 metallicities, the photon-weighted
group cross section was evaluated as
`sigma_N = integral(q_lambda * sigma(E)) d(lambda) / integral(q_lambda d(lambda))`,
using the native 100,000-point wavelength grid, the BPASS builder's photon
integrand convention, and the nine configured group edges. Groups whose upper
edge is at or below the species threshold have zero support. For this summary,
only rows with `Q_group / Q_total(0.01 eV--10 keV) > 1e-8` are included. The high-redshift
subset `age <= 1000 Myr` contains 403 age-metallicity spectra; 396 contribute
above that relative threshold in 54.42--500 eV. This is a sensitivity estimate
from native-grid trapezoids, not an adaptive-quadrature error bound.

| Species | Group (eV) | BPASS `sigma_N` range (cm^2) | BPASS span | Ratio to current pilot kernel |
|---|---:|---:|---:|---:|
| H I | 13.6--24.59 | `3.070e-18`--`4.100e-18` | 1.34x | 0.86--1.15 |
| H I | 24.59--54.42 | `5.775e-19`--`9.759e-19` | 1.69x | 0.95--1.61 |
| H I | 54.42--500 | `6.827e-20`--`1.222e-19` | 1.79x | 1.50--2.69 |
| He I | 24.59--54.42 | `4.534e-18`--`6.449e-18` | 1.42x | 0.97--1.38 |
| He I | 54.42--500 | `1.055e-18`--`1.681e-18` | 1.59x | 1.45--2.31 |
| He II | 54.42--500 | `9.203e-19`--`1.576e-18` | 1.71x | 1.48--2.53 |

These values establish that even a stellar spectral family can move the
photon-weighted H/He interaction moments appreciably within a broad group.
They do **not** directly qualify the active stellar model: the integrated
profile's runtime receipt identifies `PARSEC_mixed_lowmass_comparison`, while
this BPASS file is explicitly `independent_radiation_reference`. The current
PARSEC runtime knots preserve group photon and energy moments, not a resolved
within-group spectrum or source-specific H/He cross-section moments. Therefore
the BPASS table is a comparison envelope only. Before accepting a shared
source-independent kernel, derive the corresponding moments from the actual
PARSEC Q-constrained source construction and the selected AGN source over their
declared domains. If only Q/E group moments are available, the missing
within-group shape must remain an explicit limitation; it cannot be recovered
uniquely from those two moments.

### Actual active PARSEC-model kernel sensitivity (2026-10-31)

The active source profile was checked directly in both continuous and resumed
logs for job 444355: `parsec_mixed_lowmass_truncated_v1`, with the shared
spectral contract still admitted as `reference_control`. The matching source
package declares `full_ssp=false`, `same_evolution=false`, mass 2--600 Msun,
and 11 metallicities. This is a comparison profile, not a complete SSP or a
production-approved stellar SED.

An offline sweep covered all 979 active source nodes and verified the
runtime package's 283,185-knot extent, then evaluated the pinned PARSEC
cumulative-Q tables at their native uncompressed rows. It followed
the source builder's five-band Q-constrained Planck prior for measured rows
and separately used its actual-track L/Teff Planck tail. Ground-state H I,
He I, and He II interaction moments use the same Verner et al. (1996) fits
as `snrt_core.primordial`. The photon support cut is
`Q_group/Q_total(0.01 eV--10 keV) > 1e-8`. The reported distributions are
unweighted instantaneous knots, not IMF- or time-weighted population means.
It ran as CPU job 1251422 on Grammar `normal` (2 CPUs, 16 GiB requested,
503,444 KiB peak RSS, 4m56s elapsed).
The compact result and all 24 PARSEC input hashes are recorded in
`simulation/snrt/data/parsec_shared_kernel_sensitivity_20261031.json`;
its SHA256 is
`5256c3c5caf298805fb8c4f4b12a620b375dad3fd670df4d59d6a5444460355b`.

| Species | Group (eV) | PARSEC Q-constrained sigma range (cm^2) | Span | Range / current pilot | Best possible max relative error* |
|---|---:|---:|---:|---:|---:|
| H I | 13.6--24.59 | `2.781e-18`--`5.361e-18` | 1.93x | 0.78--1.50 | 31.7% |
| H I | 24.59--54.42 | `4.613e-19`--`1.238e-18` | 2.68x | 0.76--2.04 | 45.7% |
| H I | 54.42--500 | `5.534e-20`--`1.074e-19` | 1.94x | 1.22--2.37 | 32.0% |
| He I | 24.59--54.42 | `3.885e-18`--`7.432e-18` | 1.91x | 0.83--1.59 | 31.3% |
| He I | 54.42--500 | `8.864e-19`--`1.529e-18` | 1.73x | 1.22--2.10 | 26.6% |
| He II | 54.42--500 | `7.571e-19`--`1.403e-18` | 1.85x | 1.22--2.26 | 29.9% |

*For a single fixed value chosen to minimize the maximum relative error over
the measured minimum-to-maximum interval, the lower bound is
`(sigma_max-sigma_min)/(sigma_max+sigma_min)`. This is a range bound on the
sampled model states, not a source-population error estimate.

In the track-Planck-tail sample, spans for H I groups 13.6--24.59,
24.59--54.42, 54.42--500, and 500--2000 eV are 2.05x, 2.65x, 2.76x,
and 1.03x; for He I groups 24.59--54.42, 54.42--500, and 500--2000 eV they
are 1.92x, 2.35x, and 1.03x; for He II groups 54.42--500 and 500--2000 eV
they are 2.59x and 1.03x. The 500--2000 eV group has support above the cut
only in the track-tail sample; the 2--10 keV group has no stellar support
above the cut. Across representative masses 2, 8, 12, 13, 14,
20, 50, 80, 120, 200, 300, and 600 Msun at every metallicity, order-64 and
order-128 cross-section integrations differed by at most `5.9e-13` in either
phase. Thus quadrature error is negligible compared with the within-group
shape range for this assumed source reconstruction.

This does not recover a unique atmosphere spectrum from runtime Q/E moments:
the within-band Planck shape is an explicit source-builder prior, while the
runtime stores only photon and energy moments. The table is conditional on
that prior, and its ratios compare against a pilot reference kernel that is
itself marked `reference_control_not_science_approved`. The operator adopted
the 10% engineering comparison screen on 2026-11-02. The corrected
matched-source comparison is recorded above; it does not approve a new kernel,
source family, or runtime change.

## Closure decision after the Astra plan audit

Retain the existing shared two-moment MaxEnt reconstruction for the first
candidate. Do not replace it with fixed AGN coefficients or introduce a
declared mixed spectrum before comparing the actual runtime reconstruction to
the resolved source spectra at matched photon-number and energy moments. The
measured PARSEC ranges above show that the prior matters; they are conditional
on the stated source-builder reconstruction and do not establish the runtime
error for the candidate AGN source.

Source-resolved mesh channels are a contingent follow-up only if the declared
rate/transfer target is missed and controlled comparisons show that mixing
source families causes the miss. A source split alone does not remove
within-family spectral variation or attenuation hardening. The added primary
state would cost at least 5,184 bytes per M5 slot for order 5 and nine groups,
before transaction and checkpoint overhead, so do not incur that memory and
AMR/restart scope without this evidence.

The bounded coupled smoke does require candidate-only, fail-closed admission
in both the AGN and spectral-model loaders. It must not change the existing
reference-control or production-admission paths.

## Literature comparison and proposed direction

1. **Sazonov, Ostriker & Sunyaev (2004), average-quasar SED.** This is a
   population-average spectrum constrained jointly by composite optical/UV/X-ray
   spectra, the cosmic X-ray background, and cumulative AGN infrared emission.
   It is a useful fixed-template comparison, not an individual high-redshift
   intrinsic disk/corona model. Its hard-X-ray template peaks near 60 keV and
   extends well beyond this code's 10 keV ceiling; its IR component represents
   dust-reprocessed AGN light. Neither the out-of-band hard-X-ray energy nor
   the reprocessed IR component may be silently reassigned to the modeled
   primary groups or to the local dust transaction.
2. **Trebitsch et al. (2021), Obelisk high-redshift AGN treatment.** The
   published prescription builds a Shakura–Sunyaev multicolour disk plus a
   high-energy power law (UV slope -1.5), approximates it piecewise, integrates
   group luminosities, mean photon energies and interaction cross sections,
   then precomputes/interpolates SED-dependent quantities. Its stated aim is
   to capture variation in ionizing luminosity, not to reproduce an observed
   broadband AGN SED; the simulation transports only three ionizing bands
   (13.6–1000 eV). It assumes 30% of bolometric disk luminosity is absorbed
   and re-emitted in IR, but does not model that IR. This is valuable
   high-redshift RT precedent, not a ready-made nine-group intrinsic spectrum
   for our explicit local dust/IR model. The 30% factor must not be copied as
   an additional sink from the radiation budget when absorbed energy is
   already processed by our dust transaction.
3. **Su et al. (2025; arXiv v2 revised 2026), accretion-state-dependent
   disk/corona/slim-disk SED.** It computes spectra from BH mass and
   Eddington-scaled accretion rate, using a modified magnetic-reconnection
   disk-corona model at high rates and ADAF+thin-disk at low rates, including
   super-Eddington photon trapping. Its cosmological population comparison
   reproduces AGN fractions to z≈2.5, but its bolometric/X-ray luminosity
   functions agree mainly to z≈1 and diverge at higher redshift. It is a
   candidate for a state-dependent source model, but that population-level
   comparison does not validate its LRD/high-z applicability; the relevant
   BH-mass/accretion-rate domain and spectral closure still need to be stated.

**Candidate source family for the first bundle:** the intrinsic high-redshift
disk-plus-corona spectrum remains the leading AGN emission model. Preserve its
emitted N/E receipts and compare material rates through the existing shared
MaxEnt closure. Do not add an unresolved torus/IR component: absorbed energy is
handled by the existing gas/dust/IR transaction. Per-sink
mass/accretion-dependent opacity remains out of scope unless the matched-rate
evidence shows it is needed for the declared validity domain.

The primary-paper checks do not establish that the runtime MaxEnt
reconstruction is accurate for both the active stellar source and the selected
AGN spectrum. The matched-moment comparison below is the bounded evidence for
that question; it is not replaced by the Obelisk precedent or by the literature
review alone. Astra's plan verdict and its conditions are recorded above.

## Proposed bounded implementation and evidence

- Add a science-candidate AGN emission contract generated from the selected
  high-redshift disk-plus-corona model; leave `partition_reference_v1` and
  current reference files unchanged for controls. Load it only through the
  candidate-only, fail-closed path.
- Preserve the existing nine-group model as the reference/control. The
  candidate uses a versioned 12-group set that replaces 54.42--500 eV with
  four equal-log sub-bands. Derive source photon/energy moments by documented
  quadrature over the candidate's exact edges, including partially covered
  groups and 2--10 keV. Inspect interpolation validity before upstream
  NaN-to-floor masking; never turn an invalid model point into a plausible
  faint spectrum.
- Normalize with the selected shape-only contract above. Record group edges,
  nonnegative finite moments, H/He thresholds, quadrature convergence, and the
  out-of-band remainder in the existing source-contract evidence. Do not add
  per-cell counters or a parallel test framework.
- Compare the resolved source spectrum to the actual runtime MaxEnt at
  matched candidate-group N/E for the declared representative source cases.
  Apply the operator-adopted 10% engineering screen and keep numerical
  integration error below its declared limits. Do not carry the power-law
  prior into the runtime candidate: its filtered-spectrum failures reject it
  as a general replacement.
- Before the candidate M5+dust+CHIMES integration, close the AGN source
  receipt explicitly: record accepted `E_EM`, represented/injected
  `E_EM*sum(f_g)`, and unrepresented `E_EM*(1-sum(f_g))`, and verify their
  sum equals `E_EM` before clearing the pending receipt. The current driver
  reports fractions only, so do not claim event-level energy closure until
  this is implemented in the existing source receipt.
- Then run the existing M5+dust+CHIMES integration once with the in-domain BH
  fixture and distinct stellar/AGN receipts. Verify gas and dust/IR transfers
  and the no-double-count ledger. Do not capture an additional final HDF5
  dump. After stage evaluation, delete only this run's raw outputs and retain
  the exact cleanup manifest and compact evidence.
- Keep the candidate marked nonproduction until its literature derivation,
  provenance, and applicable data-use terms are recorded in the active release
  manifest. This bundle cannot clear CHIMES, Draine, yield, PARSEC, or other
  independent release blockers.

### Source-receipt status reconciliation — 2026-11-03

The proposed source-receipt closure item above is implemented for the
restricted 25-group M5 path. `snrt_moment_ramses` globally checks accepted
AGN `E_EM` against the group-injected plus group-unrepresented energy to
`1e-10` relative before committing the pending source receipt. The corrected
unit job 448663 passed that path; its logged values close to `1.27e-14`
relative (see the integrated-run record in
[`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md)).
This closes the arithmetic receipt subtask only. The unrepresented source
energy remains explicitly accounted but is not transported or deposited, and
the candidate remains `runtime_admissible=false`; sustained coupled evolution,
broader response validation, and release/data-use qualification remain open.

## References

- Sazonov, Ostriker & Sunyaev (2004), *Quasars: the characteristic spectrum
  and the induced radiative heating*, MNRAS 347, 144–156,
  <https://doi.org/10.1111/j.1365-2966.2004.07184.x>;
  <https://arxiv.org/abs/astro-ph/0305233>.
- Trebitsch et al. (2021), *The Obelisk simulation: galaxies contribute more
  than AGN to HI reionization of protoclusters*, A&A 653, A154,
  <https://doi.org/10.1051/0004-6361/202037698>;
  <https://arxiv.org/abs/2002.04045>.
- Volonteri et al. (2017), *High-redshift galaxies and black holes in the
  eyes of JWST: a population synthesis model from infrared to X-rays*, ApJ
  849, 155, <https://arxiv.org/abs/1704.00753>.
- Su et al. (2025), *Modeling the Spectral Energy Distribution of Active
  Galactic Nuclei: Implications for Cosmological Simulations of Galaxy
  Formation*, arXiv v2 revised 2026-03-20,
  <https://arxiv.org/abs/2501.10793>.

This is a draft plan for review, not model approval or production admission.

## Full atomic CHIMES response screen update — 2026-11-02

The accepted 10% material-channel rule has been applied to all 311
profile-bound atomic reactions for three pinned Su `agnSED` accretion states
at `M_BH=1e8 Msun`. The registered 12-group moments fail the fixed
`log10(mdot)=-1` disk-corona case (145 rate and 107 excess-energy material
channels). A 20-group greedy split tuned to that case passes it but fails the
ADAF and high-rate disk-corona states. A separate greedy split ranked jointly
over all three states finds a 25-group diagnostic with zero material failures
and a passing `1e-4 R_total` subfloor criterion for each state.

That joint result was first measured with child opacity curves resampled from
the existing 12-group bank. The same fixed 25-group edges were then rebuilt
directly from the pinned Verner95/96 and KM93 source fits, without parent-node
resampling. The direct bank also has zero material failures in all three
states. Maximum rate/excess-energy errors are 7.64%/6.96% at
`log10(mdot)=-1`, 5.00%/6.85% at `-4`, and 8.77%/6.42% at `+0.5`. The largest
near-zero response difference is `4.80e-5 R_total`. Independent source-moment
40k-to-80k relative changes in `L_grid` and per-group source moments are below
`5.8e-9`. Separately, 32-to-128 response quadrature changes are below
`3.7e-6` relative for material rates, `3.93e-7` for material excess energy,
and `2.63e-10 R_total` in absolute subfloor response. These pass the 1%
material and `1e-5 R_total` response-integration working limits.

The 25-group result is a greedy diagnostic, not a proven minimum and not a
registered runtime profile. One direct-source atomic bank is retained only as
a diagnostic; no compiled profile, matching D03/Fe or molecular tables, or
runtime migration was produced. The group count alone is 2.08x the 12-group
count; total memory and runtime are unmeasured. Attenuation, transmitted
photon N/E, gas/dust energy and momentum, molecular CHIMES, Fe/D03, and coupled
M5/RAMSES response remain open. Do not run or admit this profile on the
current evidence.

The read-only Astra review accepted the limited atomic-response conclusion,
confirmed the normalization/materiality rule, and required corrections to
quadrature bookkeeping, the subfloor working threshold, and evidence wording.
Those corrections were applied and rechecked. Summary metrics, per-group
maxima, up to 20 largest material failures, exact edges, the direct-source
atomic bank, and Slurm build records are in
[`su_agn12_profile_loader_evidence_2026-11-02.md`](su_agn12_profile_loader_evidence_2026-11-02.md)
and [`snrt_agn25_joint_diagnostic_2026-11-02/`](snrt_agn25_joint_diagnostic_2026-11-02/).

## Fixed-edge gas and dust transfer screen — 2026-11-02

The subsequent bounded D03 one-crossing screen used the same fixed 25 edges,
direct 128-node D03 absorption data, and the actual 128-node runtime MaxEnt
reconstruction. Of 36 homogeneous cases, all 24 cases at the two sampled
disk-corona states (`log10(mdot)=-1,+0.5`) pass the adopted 10% material
criterion; all 12 low-rate ADAF cases (`-4`) fail in dust absorption. The
largest ADAF dust-photon discrepancy is 46.91%; dust absorbed energy and its
absorption-only momentum proxy reach 26.87%. Quadrature and conservation
checks pass, so these failures are not explained by integration noise or a
ledger imbalance. Astra's independent read-only review agreed that the
failure is material and that no additional edge search is justified by this
screen.

The exact proposed fixed source (`M_BH=1e8 Msun`, `log10(mdot)=-1`) passes
both the direct atomic screen and this limited D03 transfer screen. This
supports only that exact sampled template. It does not establish a continuous
accretion-rate interval, and the ADAF branch is outside the passing domain.
For the selected template, 97.3% of its supported `L_grid` falls within the
0.01--10,000 eV group range; the remaining 2.68% stays explicitly
unrepresented. It is not renormalized or labeled as escaped, absorbed, or
reprocessed energy.

This is still diagnostic evidence, not runtime admission. Astra's plan audit
answers **Q-GOAL yes**: the exact passing source template is a useful bounded
step toward RT/dust coupling. **Q-LEAN is yes with a narrow scope**: no more
edge search, broad accretion-state survey, or full molecular/Fe campaign is
justified here. The approved audit recommendation is one restricted
implementation bundle, with production admission deferred:

1. Freeze the experiment to the exact intrinsic, unobscured, zero-spin
   `M_BH=1e8 Msun`, `log10(mdot)=-1` template, the existing 25 edges, and an
   explicit mapping between accepted `E_EM` and the source model's `L_grid`.
   Keep the 2.6771% outside the bands as unrepresented energy; do not
   renormalize it or call it an accuracy bound, escape, absorption, or
   reprocessing.
2. Add explicit experimental 25-group runtime support and matching assets
   only for enabled physics. Atomic CHIMES and D03 absorption are required
   for the bounded fixture. Molecular, Fe, stellar-source, scattering, and
   re-emission paths must either get matching 25-group inputs and declared
   coverage or be explicitly disabled for this fixture. The 12-group profile
   is not a substitute: it fails the fixed-template atomic screen.
3. Complete the event source receipt with accepted, represented/injected,
   and unrepresented energy; preserve the source identity and prevent duplicate
   injection before clearing the pending receipt.
4. Run one bounded integrated M5 + D03 absorption + atomic-CHIMES fixture.
   Check actual injection, transport/material exchange, conservation and
   finite states. Label the result as restricted coupling evidence; it does
   not validate molecular/Fe response, scattering, re-emission, filtered
   spectra, or a disk-corona interval.

At this audit checkpoint the 25-group profile was unregistered and the build
selection supported only the reference 9-group and candidate 12-group
profiles. The registration and restricted build support were added in the
follow-up below; no coupled run is ready to launch. Do not interpolate across
accretion states or include the failing ADAF branch based on these samples.
Astra also noted that the passing
one-crossing result does not constrain spectra after filtering or repeated
interactions, and that absorbed-energy/`c` is only an absorption-momentum
proxy, not an independent scattering-force validation. Detailed methods,
hashes, and the independent review are recorded in
[`su_agn12_profile_loader_evidence_2026-11-02.md`](su_agn12_profile_loader_evidence_2026-11-02.md)
and [`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md).

## Fixed-template molecular response screen — 2026-11-02

For the exact `M_BH=1e8 Msun`, `log10(mdot)=-1` Su template, a line-resolved
optically thin H2/CO response screen compared the pinned SED against the
25-group 128-node MaxEnt reconstruction. All four material responses pass the
operator-approved 10% ceiling: maximum group errors are 1.83% (H2 absorption),
0.174% (H2 dissociation), 1.19% (CO absorption), and 1.19% (CO dissociation).
Integrated total-rate errors are respectively 0.108%, 0.121%, 0.0081%, and
0.400%. Subfloor differences pass `1e-4 R_total`; 16-to-32 point quadrature
changes material rates by at most `7.4e-13` relative. The H2 absorption and
both CO screens exceed the 1% working target in at least one material group,
but remain below the accepted outer ceiling.

This covers only line-resolved H2/CO absorption and dissociation in the
optically thin limit. The separate fixed-source check below evaluates the
scalar Habing-energy mapping for the 30 empirical channels; neither comparison
validates their species-specific AGN cross sections, shielding,
chemistry/CVODE evolution, dust coupling, or source semantics. These results do
not unblock runtime admission. Methods and
hashes are in
[`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md)
and [`molecular_source_response_2026-11-02.json`](snrt_agn25_joint_diagnostic_2026-11-02/molecular_source_response_2026-11-02.json).

The 30 empirical CHIMES channels use a shared Habing-energy-density closure.
In the native 25-group molecular-node operator this is represented by
`sigma_eff(E)=k_CHIMES E/(c U_Habing)` for `6<=E<13.6 eV`, applied to the
current node populations in the same photo-CVODE survival transaction. For
the fixed intrinsic Su `M_BH=1e8 Msun`, `log10(mdot)=-1` source, the direct
FUV energy is `0.207436694446 L_grid`; the production 128-node N/E MaxEnt
reconstruction gives `0.207541061057 L_grid` (0.0503% difference, with direct
integration stable from 40k to 80k log-frequency intervals). This passes the
operator-approved 10% scalar G0 screen for the unfiltered fixed source. It
does not measure species-specific AGN reaction-rate errors or bound
post-transport sub-bin spectral hardening in general. A subsequent bounded
screen of this same source through one, two, and three repeated homogeneous
H/He+D03 absorption crossings passed the accepted 10% FUV-energy-density
criterion in all 12 column/pure-dust-bin controls per crossing; the worst
integrated error was 7.455% after three crossings. This closes only that
fixed-template absorption screen. The legacy grey `G0_parameter` remains a
fixed within-group power-law conversion and must not be conflated with the
node-wise closure. The full method, limitations and result are in the
diagnostic README. Angular M5 transport, scattering, evolved dust, other SED
states, and live chemistry remain unvalidated.

## 25-group asset integration follow-up — 2026-11-02

The operator accepted the 10% relative-error ceiling for materially
contributing reaction-rate and transfer channels and reaffirmed it on
2026-11-03. The applied rule is:
material channels (`>=1e-3 R_total`) have relative error `<=10%`; channels
below that floor have absolute difference `<=1e-4 R_total`; resolved
integration error is `<10%`. The 1% material and `1e-5 R_total` subfloor
limits remain working targets, not the approved outer ceiling.

The fixed 25-group candidate is registered as
`snrt_agn25_joint_diag_20261102_v1`, with `runtime_admissible=false`. Its
source-matched active D03/Fe base tables were regenerated using the exact
photon-weighted means in the fixed Su source ledger; the earlier geometric-
center D03 base is retained only as a diagnostic archive. All 17 active
asset hashes in `simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_asset_manifest.json`
were rechecked. The active package contains the 25-group atomic and coarse
molecular CHIMES banks and matching D03/Fe files.

The production C++ CHIMES loader was exercised with the existing probe against
the active atomic and molecular HDF5 banks. It verified the 25-group profile
identity, HDF5 dimensions/reaction ordering, molecular-bank loading, N/E
reconstruction closure, and finite atomic reaction moments. Result:
`CHIMES_PROFILE_LOADER_OK id=snrt_agn25_joint_diag_20261102_v1 groups=25 reactions=311 shells=682 molecular=loaded N/E=closed moments=finite`.
This is a loader check, not the full material-response criterion, CVODE
evolution, or an integrated simulation. Compilation of the probe used GCC
13.2 with installed HDF5/OpenSSL libraries because the local `h5c++` wrapper
points to a missing Conda compiler path. The bounded asset-generation follow-
up used one CPU on Lageunha; it did not submit a Slurm job or run a simulation.

The M_N source receipt now records accepted input, group-injected, and
unrepresented energy and checks their arithmetic closure across MPI ranks.
That check uses the fractions supplied by the source contract; it does not
resolve whether this intrinsic finite-grid SED is the physically correct
escaped fraction of RAMSES `E_EM`, nor establish bolometric closure. The
generated profile constant remains false, and ordinary production startup
rejects a non-admissible profile. Explicit diagnostic selection has allowed
bounded candidate runs, but does not change production admission.

The fixed `log10(mdot)=-1` candidate retains passing atomic, D03 one-crossing,
and repeated-FUV absorption-screen evidence. The sampled low-rate ADAF D03
result fails the approved 10% screen, so any restricted run must enforce the
exact passing template and exclude that state. The fixed-template
`E_EM`/`L_grid` amplitude cross-check is within 3.83% under the accepted 10%
screen. The `f_escape=1` boundary convention is explicitly declared and
source-identity bound in the candidate runtime contract; it remains an
approximation and not a physical escape measurement. The finite template
support still leaves the SED tails unconstrained.
Fe high-energy response and full species-resolved chemistry/material response
remain blockers. Jobs 448037 and 448079 provide bounded one-rank M5+dust+
native-CHIMES capture/conservation and same-job restart evidence for the exact
fixed source; they do not establish sustained evolution, multi-rank behavior,
or a broader SED domain. Do not enable runtime admission from the registration,
loader, or these restricted runs. The evidence and current
asset hashes are
recorded in [`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md)
and [`su_agn12_profile_loader_evidence_2026-11-02.md`](su_agn12_profile_loader_evidence_2026-11-02.md).

The candidate source-boundary declaration is
`simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_source_boundary_contract_v1.json`.
It binds the exact intrinsic, unobscured, zero-spin M8/mdot=-1 template to the
declared `f_escape=1` resolved-boundary approximation, the accepted `E_EM`
receipt, and the unrenormalized represented/unrepresented group accounting.
The runtime loader enforces the compiled source identity, while the contract
remains `candidate_explicit_sed` and the profile remains
`runtime_admissible=false`. The declaration does not measure physical escape
or establish bolometric closure.

The native contract smoke compiled both the 9-group reference and 25-group
candidate. It confirmed the candidate identity is accepted only as compiled,
rejected mutations to each source identity component, loaded valid 25-group
source moments, rejected altered group means/fractions with unchanged hashes,
and kept candidate runtime/checkpoint admission false. The source JSON reader
rejects duplicate keys. The full native spectral-contract script passed; this
is contract/loader evidence, not a 25-group RAMSES simulation.

### Fe UV group-coverage diagnosis — 2026-11-02

The UV group-indexing gap is now fixed in `dust_iron_photons.f90`. The kernel
selects the contiguous profile groups whose upper edge and representative
photon-weighted mean are both at or below 13.6 eV. Native Fe smoke passes for
the 9-, 12-, and 25-group profiles under GNU bounds checking and Intel `ifx`
bounds checking; the 25-group fixture also verifies photon transfer in each
newly admitted group 5–7. This admits groups 1–7 in the 25-group profile
(including groups 5–7, with means 2.7645, 7.9838, and 12.335 eV, carrying 35.95% of
`L_grid`) while preserving hard rejection of groups 8–25. Those harder groups
carry 59.76% of `L_grid`; the metallic-Fe UV model still has no validated
photoelectron/Auger energy partition for them. Groups 1–4 carry 1.61% of
`L_grid` and have means below 0.77 eV. Thus this resolves an indexing/runtime
guard defect, but does not clear the 25-group Fe blocker or make that profile
runtime-admissible. Do not convert hard-photon absorption to Fe heat and call
it production physics.

### Fe hard-photon scope disposition — 2026-11-03

The operator adopted the 10% outer engineering screen for SED-response
comparisons: material channels contributing at least `1e-3 R_total` must have
relative error `<=10%`; lower-contribution channels use an absolute difference
limit of `1e-4 R_total`; resolved integration error must be `<10%`. The 1%
material and `1e-5 R_total` limits remain working targets only.

The literature check does not provide a table that can be copied into the
current elemental-metallic-Fe receiver. Kimura (2016) defines and tabulates
photoelectric quantum yield, i.e. escaping electrons per absorbed photon for
iron grains at 1, 10, and 100 nm; it does not provide the full absorbed-energy
partition among Fe heating, escaping electrons, and escaping fluorescence over
the candidate's 13.6 eV–10 keV bands ([paper](https://doi.org/10.1093/mnras/stw820),
[full table](https://academic.oup.com/mnras/article/459/3/2751/2595216)). Dwek &
Smith (1996) do tabulate deposited energy from 10 eV to 1 MeV for graphite and
silicate grains of 50 Å–1 μm, including Fe K fluorescence within silicate;
those are not elemental metallic-Fe grains and are a cross-material benchmark,
not a valid Fe replacement ([paper](https://doi.org/10.1086/176933),
[CDS data record](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/459/686)).

The project completion scope is an effective galaxy-formation model; a
microscopic elemental-Fe electron/photon transport calculation is a separate
long-term research item, not a closeout prerequisite. At the operator's
direction, a bounded CPU-only Geant4 Livermore/Penelope Fe-sphere diagnostic
was run on 2026-11-03; its results and limits are recorded in
[`geant4_fe_grain_benchmark_2026-11-03/README.md`](geant4_fe_grain_benchmark_2026-11-03/README.md).
That two-model comparison does not constitute independent experimental
validation or runtime approval. Standalone PENELOPE and NIST XCOM remain
possible references for further long-term work. The 10% material-channel /
transfer screen and `<10%` independent integration-error limit remain the
adopted comparison criteria for any future response-table work.

For the present 25-group diagnostic, metallic-Fe grains are explicitly
disabled (`dust_iron_model='none'`). The existing Fe path continues to admit
only its declared UV domain and hard-reject unsupported groups; no absorbed
hard-photon energy is reassigned to Fe heat.

### Cold-transition coupled confirmation — job 448037

The bounded cold-transition D03 + CHIMES candidate run completed on Syntax
A40 (Slurm exit 0, 36 s) with the exact M8, `log10(mdot)=-1` source identity.
One packet supplied `8.41423054e49 erg`; `8.18897147e49 erg` entered the
represented groups and `2.25259074e48 erg` remained explicitly unrepresented.
M5 order 5 produced positive gas and dust photon/energy capture, a finite
positive CHIMES state-change diagnostic, and a maximum radiation-ledger
relative residual of `1.064e-15`. The zero-source preparatory transaction
had exactly zero energy-ledger terms. No `output_*` dump was written, and the
candidate remains `runtime_admissible=false`.

At this checkpoint the restricted cold-transition live-capture confirmation
was closed, while HDF5 restart equivalence remained open. The run used
absorption-only D03; metallic Fe and PAH were off, and scattering and IR
exchange were not exercised. `CHIMES_RHS_COUNTS` recorded no GPU RHS calls, so
this is not GPU chemistry or performance evidence. Do not enable candidate or
production admission from this one bounded transaction. Its full logs, build
identity, and measurements are recorded in
[`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md).

### Restricted restart-equivalence follow-up — job 448079

Job 448037 left HDF5 restart equivalence open. Job 448079 subsequently ran
continuous and resumed segments through fine step 4 in the same one-rank
Syntax A40 allocation. The candidate-only comparator passed all 626 datasets
at `rtol=1e-10`, `atol=1e-25`, with exact grid/particle identity and no clock
violations. The earlier `rtol=1e-12` comparator failed, so Slurm correctly
retains job 448079 as `FAILED (1:0)`; the post-hoc pass does not rewrite the
job outcome. Ordinary M5 comparisons retain `rtol=1e-12`.

The worst passing tolerance ratio was 0.856. Astra reviewed the retained
evidence and accepted this narrowly scoped restart result. It closes restart
equivalence only for the fixed one-rank AGN25 candidate; `runtime_admissible`
remains false. The raw HDF5 files were checksum-verified and removed after
evaluation. Full details are in
[`snrt_agn25_joint_diagnostic_2026-11-02/hdf5_restart_equivalence_2026-11-03.md`](snrt_agn25_joint_diagnostic_2026-11-02/hdf5_restart_equivalence_2026-11-03.md).

### CHIMES conservation-scope clarification — 2026-11-03

The existing native spectral photo-CVODE path enforces chemistry conservation
at each successful operator call: `snrt_chimes_budget` sums all 157 species
into the 11 tracked elemental inventories and net charge, and `photo_step`
rejects a step if an elemental inventory changes by more than `1e-8` relative
(with a `1e-20` abundance floor) or net charge changes by more than `1e-8`
absolute. The live cold/molecular adapter obtains shielding and pumping
factors before calling the profile-bound molecular spectral receiver, then
advances the dark CHIMES network. The successful cold-transition fixtures
therefore show these operator-level guards passed for their exercised calls;
describing the remaining work as an absence of chemistry/material
conservation checks would overstate the gap.

The unresolved item is broader response validation: independently compare
the 157-species state/rate response and molecular shielding/heating against a
resolved-spectrum reference over the admitted gas states, and close the
chemistry, thermal, grain, and radiation budgets across sustained coupled
hydrodynamic evolution. Existing per-call guards and short capture fixtures
do not establish that broader result or change `runtime_admissible=false`.

### MPI2 restart follow-up — 2026-11-04

The distributed receiver eligibility defect is fixed, and job 449907 shows that
the global sink-stat reduction and sink position update agree across the
continuous/restarted branches at the measured steps. The previous large
`sink_stat_1_3` difference was not reproduced. Strict restart equivalence is
still open because seven level-3 hydro values exceed the unchanged
`rtol=1e-10`, `atol=1e-25` gate; the maximum ratio is 52.81. This restricted
candidate remains `runtime_admissible=false` and is not a production-admission
result.

The MPI2 restart outliers in job 449907 map to CHIMES species slots 10, 12, 15,
and 146 in the HDF5 `hydro/uold_N` passive-field block. Bulk hydro state,
sink statistics, particle state, and SNRT moments pass the unchanged gate.
The remaining issue is bounded chemistry-state drift under repartition; no
restart tolerance has been broadened and candidate admission remains disabled.

### Latest MPI2 restart checkpoint — jobs 449951/449972

The follow-up to canceled job 449759 ran in 449951 and 449972. Both runs
completed their continuous/resumed integration segments and passed live source
accounting, but both remained strict comparator failures. In 449972, 494/625
datasets were bitwise equal; nine values outside tolerance were limited to
five CHIMES species (`HeIII`, `CIII`, `CIV`, `CV`, `CH`), two sink diagnostics,
and two SNRT moment values. The sampled AGN gate inputs and decision matched,
and the selected oct trace was bitwise equal; the previous AGN-boundary
difference was not reproduced there. The failures do not justify relaxing the
restart tolerance or changing the physical deposition path. Candidate runtime
admission remains false. If MPI2 restart closure is pursued further, measure
signed species deltas and their elemental-budget impact at comparator-identified
cells; this is a focused follow-up, not a new source-spectrum survey. Detailed
results are in [`hdf5_restart_equivalence_2026-11-03.md`](snrt_agn25_joint_diagnostic_2026-11-02/hdf5_restart_equivalence_2026-11-03.md).


### Signed chemistry / elemental-budget checkpoint — job 450008

The focused follow-up requested after jobs 449951/449972 is complete. Job
450008 retained the same MPI2 candidate, source, physics and strict comparator.
Both branches reached fine step 4 and passed live source capture, but the
comparator still failed: 494/625 datasets are bitwise identical; 13 hydro
values and four SNRT-state values exceed the unchanged rtol=1e-10, atol=1e-25
limits. The largest SNRT component is persisted HeIII state, not an M5 angular
moment. Signed species changes at the principal outlier cell conserve the
tracked elemental inventories to floating-point precision (maximum relative
H difference 2.06e-16, C difference 1.12e-16). This narrows the interpretation
to small species-state drift while leaving strict restart reproducibility
unclosed. Do not weaken tolerance or repeat broad tracing on this evidence.

The restricted candidate remains runtime_admissible=false. Outstanding
physical/admission limits are unchanged: one fixed M8, log10(mdot)=-1 source
state; the sampled ADAF+D03 response fails the adopted 10% screen; the
f_escape=1 boundary is an explicit approximation; source energy outside the
represented groups is recorded but not transported; hard-photon metallic-Fe
energy partition is unvalidated (Fe remains disabled in this candidate); and
sustained 157-species chemistry/material response, wider source-state coverage,
MPI-scale production behavior, and release/data-use closure are not established.
The exact results, source mappings and hashes are in
[the restart-equivalence record](snrt_agn25_joint_diagnostic_2026-11-02/hdf5_restart_equivalence_2026-11-03.md).

The signed-state and element-budget questions from the earlier runs are
answered. Following Astra's focused review of job 450478, one targeted paired
receiver-input capture is queued as job 450842. It targets the dominant CIII
cell and is a prerequisite to a controlled receiver replay; it does not reopen
a broad or equivalent MPI2 run matrix. After that focused diagnosis, continue
the approved AGN-source closure only where it removes a listed physical
limitation, keeping candidate admission fail-closed until the applicable
response domain and source/release contract are justified.

### Targeted CIII receiver restart replay — job 451068 (2026-11-05)

The queued A40 capture was replaced by a Syntax A100 run after its scheduled
start was delayed. Job 451068 built with the full patch layout and ran both
continuous and restarted MPI2 branches. The target was level 3, child 4, at
`(0.375, 0.375, 0.625)`. Both branches reached the target capture and the
coupled AGN/M5/dust/CHIMES smoke passed; radiation closure was `2.339e-15`.
The output also contains repeated CVODE recoverable-RHS warnings, so this
smoke result does not establish a clean CHIMES integration.

The paired step-2 target inputs are bitwise identical across all captured
species, group photo rates, elements, grains, controls and units. Their CHIMES
outputs are not: 35/157 species differ (maximum relative difference
`8.34e-9`) and one ledger component differs by `6.30e-17`. This locates the
first measured divergence after identical captured inputs, within the
CHIMES update/output path or hidden execution state; its mechanism is not yet
identified. At step 3, species and photon inputs have consequently diverged
at up to `8.34e-9` and `1.12e-8` relative, while the captured elemental inputs
remain bitwise equal.

The integrated comparator still fails at the unchanged restricted tolerance
(`rtol=1e-10`, `atol=1e-25`): 626 datasets compared, 498 bitwise identical,
with 26 values outside tolerance in 11 datasets (22 hydro, two sink cloud
statistics and two SNRT moment datasets). Grid ownership and CPU maps match
exactly; elemental budgets remain conserved to floating-point precision. The
failure profile therefore extends beyond the trace-species-only set reported
for job 450478. Raw HDF5 files were removed after comparison; sizes and
SHA256s remain in the run directory. Keep `runtime_admissible=false`; do not
loosen tolerance or repeat a broad MPI2 matrix.

Astra's read-only technical consultation answers **Q-GOAL: yes** and
**Q-LEAN: yes**, for one target-only replay that separates prior execution
context from RHS backend. The identical captured vector is strong localization,
but does not prove that every effective solver input or call-order dependency
was captured; exact ownership does not prove identical thread order or hybrid
dispatch. Run the same saved input and solver settings once in each cell of a
small 2x2 comparison: warmed/reused versus freshly initialized CHIMES execution
context, crossed with production hybrid versus forced CPU RHS. Restore the
input bytes before each call, retain the actual prior-call history for the warm
case, and compare species, photons, ledgers, return status and existing
target-local dispatch counts. No full evolution or HDF5 dump is needed. This
distinguishes backend sensitivity from context/order dependence without
another broad MPI2 replay; it does not by itself establish which path is more
accurate. Repeated recoverable-RHS warnings are not by themselves proof of an
uncompleted integration or its cause.


### Explicit-Courant MPI2 restart follow-up — job 450478

Job 450478 used the unchanged MPI2 candidate and strict restart tolerance with
an explicitly recorded CHIMES Courant factor of 5.0e-8. Continuous and resumed
branches both reached coarse step 4, restored the same grid ownership, and
passed live AGN/M5/dust/CHIMES source accounting. Three source-bearing
transactions had a maximum relative photon/energy closure residual of
1.34e-15. This confirms bounded coupled execution at the diagnostic step size;
it does not establish production operation at the default Courant factor.

The strict comparator still failed: 626 datasets, 491 bitwise identical, and
45 values in 12 datasets outside rtol=1e-10, atol=1e-25. The differences are
limited to CHIMES trace species and the saved HeIII fraction. H and C elemental
inventories agree to 2.06e-16 and 4.49e-16 relative, respectively, while the
species-level restart gate remains unmet. This reproduces the small chemistry
and HeIII drift seen in 450008; do not loosen tolerance or schedule another
equivalent MPI2 replay. Keep the candidate runtime_admissible=false.

The job's Slurm exit was caused first by a wrapper requiring failure-only
CHIMES input-capture records after a successful solve. That guard is removed;
the exact post-hoc comparator result and raw-output cleanup checksums are
recorded in
[the restart evidence report](snrt_agn25_joint_diagnostic_2026-11-02/hdf5_restart_equivalence_2026-11-03.md#explicit-courant-integrated-restart-result--job-450478--2026-11-05).
No physics or tolerance changed.

## Proposed next bundle: event-weighted, state-dependent AGN SED

**Status: Astra plan review complete (2026-11-05); conditional proceed after
scope narrowing. Implementation not started; supported-state boundary remains
to be resolved.**

### Q-GOAL — ask first

Does carrying the pinned Su et al. source shape from each accretion event into the existing M5, gas and dust source transaction advance the approved galaxy-scale effective-model objective? Does it close a material physical gap left by the current single M8, log10(mdot)=-1 template, without expanding into unresolved spin, obscuration or microscopic source physics? The current 25 transport groups cover 97.32% of `L_grid` at log10(mdot)=-1, 16.01% at -2.95, and 21.80% at -4.00, and only 0.000015% at -12. Is a state-dependent SED scientifically admissible on those edges? Should the implementation instead extend groups and material response, restrict the runtime state domain with a fail-closed coverage criterion, or explicitly declare the outside-band component untransported under a justified approximation? At the upstream ADAF/disk-corona switch, the represented fraction changes from 13.92% at log10(mdot)=-2.9000 to 98.40% at -2.8999, while `L_grid` changes by about 9.3%; is this branch discontinuity part of the accepted source model and how must a runtime table preserve it? Also decide whether the model should retain RAMSES `epsilon_r=0.1` as an explicitly shape-only source contract, or apply Su et al.'s state-dependent efficiency consistently to radiated energy and BH retained mass.

### Q-LEAN — ask second

After answering Q-GOAL, is the smallest physically sufficient response lean enough? Flag unnecessary table dimensions, instrumentation, duplicate gates, or tests. Keep verification to the existing source-response criterion, event/merge/restart energy accounting, and one candidate integrated restart fixture; do not add a broad benchmark matrix.

### Proposed implementation

The current sink path reduces each accepted event's electromagnetic energy into the scalar per-sink agn_pending_erg reservoir. Applying one SED after events with different masses or accretion rates have accumulated would assign the wrong spectrum to earlier energy. The candidate path therefore needs to retain group moments at emission time.

1. Use the already selected, pinned Su agnSED model, with spin fixed to zero and the existing intrinsic/no-obscurer, f_escape=1 boundary declaration. At each accepted AGN accretion event, evaluate the source at that event's black-hole mass and Eddington-scaled accretion rate. Reject out-of-support states; do not extrapolate. The supported mass/accretion domain must be stated explicitly from the model source and the intended simulation, not guessed from the one-cell fixture.
2. Generate an offline, hash-bound numerical lookup from the pinned upstream SED data. It supplies per-state group energy and photon-number fractions plus the explicitly unrepresented fraction. The accepted 25-group edges are not assumed to cover the full SED: implement only after the Q-GOAL decision selects extended groups/material response, a fail-closed restricted state domain, or a declared untransported component. Likewise, freeze the normalization contract before implementation: retaining the accepted RAMSES E_EM receipt is a shape-only approximation, while using state-dependent efficiency changes radiated energy and BH retained mass and must be accounted through the existing mechanical partition. Keep the fixed-template candidate as a comparison.
3. Accumulate per-sink, per-group pending energy and photon number from each event. Keep the total accepted electromagnetic-energy receipt and unrepresented fraction available for closure. Reduce by addition across OpenMP threads, MPI ranks and sink mergers. Consume these arrays only after a successful SNRT transaction; failed transactions preserve the pending receipt.
4. Persist the pending source moments, source-table identity and profile identity in HDF5 sink restart state. Legacy/reference-control paths retain their current representation. A candidate restart with positive scalar pending energy but no matching group moments must fail closed; a zero-pending restart may initialize empty group moments.
5. Keep all state-dependent SED use behind the existing candidate-only opt-in. This bundle must not set runtime_admissible=true or alter reference defaults.

### Bounded verification

- Compare table interpolation against the pinned upstream interpolator at table nodes and independent cell midpoints, including both sides of its ADAF/disk-corona transition. Enforce no extrapolation, nonnegative finite moments, per-event energy closure and the existing accepted material-channel/integration-error limits over the declared domain.
- Reuse the existing native smoke and source-ledger checks for group N/E, unrepresented energy, and no double-counting with mechanical feedback.
- Exercise additive thread/MPI reduction, sink merge, HDF5 write/restore, and exact pending-ledger recovery using the existing restart harness.
- Run one MPI2 candidate M5+dust+CHIMES integrated restart using the existing fixture and output cleanup policy. Do not add a broad trace matrix or a new benchmark framework.
- Retain runtime_admissible=false until state-domain response validation, the MPI2 restart gate, Fe hard-photon treatment, the f_escape boundary approximation, sustained coupled chemistry/material response, and release/data-use closure have each been resolved or explicitly accepted for a declared restricted domain.

### Review request

Please return Q-GOAL first and Q-LEAN second, then assess scientific validity and code feasibility. In particular, address the measured transport-band coverage at log10(mdot)=-1, -2.95 and -4; choose the minimum defensible treatment of outside-band energy; and decide whether the accepted `epsilon_r=0.1` receipt can remain as a declared shape-only approximation or must follow the tabulated state-dependent efficiency. Also check whether event-time state is available and consistent with the accepted E_EM transaction, whether group N/E moments plus the selected material response are sufficient for spectral transport, and whether the HDF5/merge contract is coherent. Treat this as advisory plan review; do not edit files or run jobs.

### Driver source-path review — 2026-11-05 (before auditor reassignment)

This is an implementation-feasibility review, not a plan-audit verdict. At
the time, the previously assigned Claude CLI route reported `loggedIn=false`
on the local host, Syntax, and Lageunha. That authentication state is
historical and no longer blocks the reassigned Astra review.

- At Bondi-event time, `msink(ksink)`, `dMBHoverdt(ksink)`,
  `dMEdoverdt(ksink)`, and the selected efficiency are available. The accepted
  radiative receipt is calculated from the donor-limited `gross_mass` after
  the density-floor cap. Define the SED state convention explicitly: use the
  sink mass at the start of that accretion update and the same Eddington-rate
  ratio/cap convention that selected the accepted flow. Do not infer a new
  accretion state from an individual donor cell's clipped mass. Record how
  `mdot > 1`, the low-rate `f_bondi` correction, and any out-of-table state are
  handled; out-of-support states must fail closed.
- The current code reduces `event_radiated(j)` into a thread-local scalar per
  sink, MPI-reduces it, then stores it in `agn_pending_erg`. Different SEDs
  cannot be reconstructed from that reservoir. Event-weighted group moments
  must be accumulated alongside the accepted event receipt through the
  OpenMP, MPI, sink-merger, and HDF5 paths. Keep the scalar energy receipt as
  an independent closure check.
- The pinned table's group energy fractions are relative to `L_grid`, which
  is explicitly not asserted to be bolometric. Before implementation, state
  the normalization contract that maps the accepted `radiated_erg` to this
  grid-normalized shape, and define the disposition of energy outside both the
  model grid and transport-group edges. The ledger must not label this as a
  measured bolometric correction or escape fraction.
- The source builder currently derives AGN photons using a fixed group mean
  and then reconstructs source energy as `N * mean_energy`. The candidate
  state-dependent path must instead inject table-derived N and E separately;
  convert group E from erg to the transport's eV units and preserve the
  represented/unrepresented energy ledger. The downstream M5 source split
  already accepts the N/E pair.
- The active CHIMES band reconstruction consumes group N/E and reconstructs
  128 frequency nodes per group, so it can retain source-dependent group
  means within the declared reconstruction model. This does not make every
  receiver source-independent: check the selected candidate's PAH, Fe, and
  any non-band/legacy receiver paths, since some eligibility and heating
  checks still read the global fixed `snrt_group_mean_energy_ev`.
- HDF5 currently writes `agn_state_schema=2` and one scalar pending-radiation
  dataset. The proposed group ledger therefore requires an explicit schema
  revision and a restart compatibility rule keyed to the SED-table/profile
  identity; scalar-positive legacy checkpoints cannot be losslessly promoted
  to state-dependent group moments.

These points should be included in the Q-GOAL/Q-LEAN review prompt. In
particular, ask whether the grid-to-accepted-energy normalization can be
justified within the approved effective-model scope before approving a wider
state-dependent SED domain.

### Pinned-model and implementation cross-check — 2026-11-05

The source package and primary paper resolve part of the state-definition
question, while exposing an efficiency-scope decision that still needs review.

- The pinned package defines `logmdot_crit=-2.9`, `logmdot_trunc=-3.0`, and
  `lognu=8.0..22.9` in 0.1-dex increments. Inspecting its serialized
  interpolator axes gives the usable domains: ADAF `log10(M/Msun)=0.8..11.5`,
  `log10(mdot)=-12.2..-2.9`; truncated disk `log10(M/Msun)=0.8..11.0`,
  `log10(mdot)=-3.0..-2.9`; disk-corona `log10(M/Msun)=0.8..11.5`,
  `log10(mdot)=-3.0..3.0`. The runtime Bondi/Eddington cap limits the current
  RAMSES candidate to `mdot<=1`; it cannot exercise the upstream super-Eddington
  table region without a separate accretion-model change.
- Su et al. define `mdot=Mdot_acc/Mdot_Edd` with
  `L_Edd=0.1 Mdot_Edd c^2`, assume zero spin, and calculate model bolometric
  luminosities by integrating the SED. The paper reports efficiency near 0.1
  for `-2.9 < log10(mdot) <= 0`, with lower efficiency in the ADAF and
  super-Eddington regimes. The local artifact therefore remains named
  `L_grid`: it is the converged integral over the pinned finite frequency
  support, not an external bolometric correction or an assertion about flux
  outside that support. Source: [Su et al. (2026), ApJ, DOI
  10.3847/1538-4357/ae41bd](https://doi.org/10.3847/1538-4357/ae41bd),
  [arXiv:2501.10793](https://arxiv.org/abs/2501.10793).
- In the spin-disabled RAMSES path, `epsilon_r=0.1` and
  `dMEdoverdt` is proportional to `M_BH/epsilon_r`; therefore the existing
  `min(dMBHoverdt,dMEdoverdt)/dMEdoverdt` rate convention matches the paper's
  reference Eddington mass rate when MAD efficiency changes are disabled.
  Apply the low-rate `f_bondi` correction before mapping that ratio, and keep
  zero accretion as a zero-source case. The event's donor-density cap still
  scales its accepted energy receipt separately from the sink-level SED state.
- The pinned paper couples SED shape and state-dependent radiative efficiency.
  Current RAMSES receipts use `epsilon_r=0.1` for the candidate's spin-off,
  MAD-off contract. Preserving that receipt while changing only the SED shape
  is a deliberate shape-only approximation; outside the paper's approximately
  0.1-efficiency disk-corona regime it does not reproduce the paper's
  state-dependent bolometric luminosity or BH retained-mass fraction. Astra's
  Q-GOAL/Q-LEAN review must decide whether this project should (a) explicitly
  retain and label the existing energy contract, or (b) include the tabulated
  radiative efficiency in the accretion receipt and mass-retention update.
  The latter is more physically self-consistent but changes AGN mass/energy
  feedback and requires the event ledger to remain conservative through the
  existing mechanical partition.
- The existing SED-moment builder omitted the upstream truncated thin-disk
  component on `-3.0 <= log10(mdot) <= -2.9`. It now sums the ADAF and truncated
  disk luminosities in linear space on that closed interval, matching
  `agnSED.SED(truncation=True)`. A focused parity check at six points
  (`-3.1, -3.0, -2.95, -2.9, -2.899, -1.0`) matched the pinned upstream
  routine to at most `3.6e-15 dex`; this also confirms the existing `-1.0`
  diagnostic profile is unchanged by the correction. The broader dynamic table
  and runtime coupling remain unimplemented pending the design decision and
  plan review.

### Transport-band coverage consequence — 2026-11-05

An integration over the pinned SED's full frequency support quantifies a
separate limitation of reusing the current 25-group edges (`0.01 eV` to
`10 keV`). These are bounded source-domain checks, not an expanded benchmark:

| `log10(mdot)` | SED branch | Energy below 0.01 eV | Energy above 10 keV | Energy inside current groups |
| ---: | --- | ---: | ---: | ---: |
| -1.00 | disk-corona | 0.00000015 | 0.02677105 | 0.97322879 |
| -2.95 | ADAF + truncated disk | 0.00061921 | 0.83926716 | 0.16011363 |
| -4.00 | ADAF | 0.40674672 | 0.37522147 | 0.21803180 |

Thus a dynamic SED lookup on the existing edges can preserve the source ledger,
but it does not provide coupled transport for most of the modeled radiative
energy near the ADAF transition or in the low-rate ADAF branch. Do not describe
those states as fully RT-coupled or production-admitted solely because their
missing fraction is logged. The plan review must decide whether the minimal
valid implementation is to extend the energy groups and matching material
response, declare a restricted runtime domain with fail-closed checks, or
explicitly accept an untransported/escaping component as a model approximation.
The current group coverage result makes that a Q-GOAL decision; Q-LEAN should
then judge the smallest physically sufficient response rather than add a
generic instrumentation matrix.

The same plan review must decide how to treat Su et al.'s state-dependent
radiative efficiency. Keeping the current `epsilon_r=0.1` receipt while using
their low-rate spectral shapes is a shape-only model; using the table's
integrated efficiency would alter the sink retained-mass and radiative-energy
receipts. Do not implement the broad state-dependent path or claim full source
physics until both the band-coverage and energy-normalization choices are
settled.

A bounded 13-state M8 sweep across `-12 <= log10(mdot) <= 0` shows the present
groups' represented fraction ranges from `1.5e-7` at -12 through `0.00138` at
-6, `0.218` at -4, and `0.139` at -2.9, then rises to `0.984` just above the
branch switch and reaches `0.995` at 0. The builder matches the pinned
`agnSED.SED(truncation=True)` pointwise to at most `3.6e-15 dex` at
log10(mdot)=`-2.9001, -2.9000, -2.8999, -2.8990`; the abrupt coverage change
therefore comes from the upstream piecewise SED itself, not this builder.
Preserve the upstream branch as a categorical boundary; do not interpolate or
smooth across it without a separate source-model justification.

`P0_REQUIREMENTS.md` explicitly leaves photons above 10 keV outside the P1
groups pending evidence that they affect the target gas. The current SNRT
source/transport code search found CMB Compton exchange but no AGN-photon
Compton heating operator. Thus the new Su SED's large >10 keV fraction in the
ADAF states is not yet shown to be dynamically negligible or physically
coupled. The Q-GOAL decision should account for this missing response when
choosing between extended transport/Compton response, an explicitly admitted
escape approximation, or a justified runtime state boundary.

The updated builder was exercised through JSON serialization at the transition
state `log10(mdot)=-2.95` using 40,000 and 80,000 quadrature intervals. It
reported 25 groups, represented fraction `0.160113633239`, a group-coverage
closure residual of `-8.38e-10`, and maximum refined quadrature change
`5.35e-9`; all are within the builder's 1% numerical working limit. The
temporary JSON artifact was removed after the check. The three-state coverage
check at `log10(mdot)=-1, -2.95, -4` reproduced the table above.

### Astra plan audit — 2026-11-05

**Q-GOAL: Conditional yes. Q-LEAN: yes after scope narrowing.** The audit
concluded that event-weighted AGN spectral moments can advance the final
SNRT/M5, gas and dust coupling objective, because the existing scalar
`agn_pending_erg` cannot reconstruct different spectra from separate
accretion events. The all-state version of this proposal is not supported by
the available response evidence. Keep the shape-only `epsilon_r=0.1` energy
contract, preserve the unrepresented group fraction without renormalizing or
calling it transported/escaped, and keep `runtime_admissible=false`.

The full integrated evidence still applies only to the exact
`M_BH=1e8 Msun`, `log10(mdot)=-1` source. Fixed-edge offline atomic, neutral-H/He
plus D03 absorption, and H2/CO line-response screens now also pass at two
additional accretion rates at M8 and at four luminosity-weighted mass/rate
points from the P4 pilot source ledger. These are discrete response samples;
they do not cover the 30 empirical CHIMES Habing channels, chemistry/CVODE
evolution, repeated-FUV response, Fe hard-photon response, or integrated
multi-state M5 coupling. Their fixed 25-group coverage leaves 3.64--9.42% of
`L_grid` explicitly unrepresented. The P4 ledger is a single `aexp=0.2085`
pilot snapshot, not an approved production population. It suggests a useful
multi-state design question but cannot define the future production domain.
The `+0.5` M8 sample is outside the current Eddington cap, and the low-rate
ADAF sample fails D03. Before building an interpolated table or changing HDF5
state, derive the intended cosmological run's BH-mass/accretion distribution,
set a fail-closed domain and outside-domain treatment, then validate the
chosen interpolation midpoints and one integrated multi-state use. Keep
`runtime_admissible=false` and preserve the exact-template scalar path until
that contract is justified. The P4 distribution and four-state response
evidence are recorded in
[`snrt_agn25_joint_diagnostic_2026-11-02/README.md`](snrt_agn25_joint_diagnostic_2026-11-02/README.md#p4-pilot-luminosity-weighted-source-state-screen--2026-11-05).

If a multi-state domain is later justified, use sink-level accretion state at
the start of the update and weight group contributions by the accepted
post-mechanical-partition `event_radiated` receipt. Accumulate group energy
and photon number independently (`q_g` in photons/erg), and carry them through
thread/MPI reductions, sink merges, HDF5 restore, and commit-only consumption.
Do not infer state from donor-clipped cell mass. Existing PAH and Fe receivers
that depend on fixed group means must remain disabled unless their source
support is separately connected. The single integrated restart proposed by
the audit is adequate for the eventual bounded implementation; no repeated
MPI2 restart, broad SED grid, edge optimization, Fe microphysics, or general
performance matrix is authorized by this review.

Full verdict: [`astra_agn_sed_bundle_plan_audit_2026-11-05.md`](astra_agn_sed_bundle_plan_audit_2026-11-05.md).

### Astra follow-up review and run-domain check — 2026-11-05

**Q-GOAL: conditional proceed with bounded offline work. Q-LEAN: defer the
runtime/HDF5 transaction until its domain is justified.** The added M8 and P4
samples support continuing the spectral-response assessment across discrete
source states. They do not establish a continuous cosmological runtime domain.
The state-dependent moments address a real limitation of scalar
`agn_pending_erg`, but do not yet justify an interpolated runtime table.

The intended 12.5-cMpc/128³ low-redshift calibration baseline explicitly has
sinks and AGN off ([baseline record](galaxy128_grackle_calibration_baseline_2026-10-08.md));
it therefore supplies no black-hole mass/accretion distribution for this SED
contract. The P4 rate ledger is one `aexp=0.2084978` snapshot from a separate
pilot and is not an approved production population. No AGN-on target-run domain
is defined in the reviewed plan. Keep the M8, `log10(mdot)=-1` case as a fixed
diagnostic template and keep `runtime_admissible=false`.

The minimum next evidence is: (1) a run-specific joint mass/accretion domain,
redshift and rate convention, with positive-energy out-of-domain handling;
(2) an explicit disposition for the 3.64–9.42% `L_grid` fraction outside the
25 transport groups, without inferring acceptability from the 10% material
screen; and (3) boundary plus independent interpolation-midpoint checks in
the chosen domain, including the active empirical CHIMES FUV/Habing response
at one demanding state. Existing offline tools suffice; no persistent table,
new diagnostic framework, broad scan, or runtime/HDF5 change is warranted
before these choices. After a domain is justified and the implementation is
made, use one bounded integrated restart; the strict restart failure remains
open.

### Operator source-domain decision — 2026-11-06

The operator selected: **keep production AGN SED work blocked until a target
AGN-on run is designated**. No run-specific mass/accretion domain is currently
approved. The P4 Run0 snapshot remains diagnostic evidence only, and the
12.5-cMpc/128³ calibration baseline has AGN disabled. Keep
`runtime_admissible=false`; do not build or admit the persistent event-weighted
SED/HDF5 runtime path from those datasets. This decision blocks production
admission only; it does not invalidate the existing fixed-template diagnostic
or the separate, unresolved CHIMES restart-consistency investigation.

Current fail-closed enforcement was rechecked on 2026-11-06. `mkrun.py`
rejects every selected group profile with `runtime_admissible=false`, including
this 25-group profile. The direct RAMSES diagnostic exception requires the
explicit `SNRT_ALLOW_25G_AGN_DIAGNOSTIC=1` opt-in and then validates the
identity-bound source contract plus the restricted M5 order-5, sink/AGN,
CHIMES/D03, Fe-off, dust and bounded-fixture conditions collectively across
MPI ranks. It records a diagnostic permission separately; it does not set
`snrt_spectral_contract_runtime_allowed` or promote the profile. Thus the
fixed-template diagnostic remains available under its narrow fixture while
production use stays blocked.

The existing native admission/contract smoke was executed CPU-only on grammar
`normal`, Slurm job 1271987 (grammar001, 2 CPUs, 6 s), and completed with
`SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK`. It generated the 25-group profile as
`runtime_admissible=false` with source binding enabled; the loader accepted the
matching candidate values for inspection while keeping both runtime and
diagnostic permission closed until the full all-rank fixture gate. It also
rejected legacy grey H/He opacity in the CHIMES+D03 candidate path. The test
does not launch RAMSES or promote the profile. The first attempt, job 1271983,
exposed a missing OpenMP link flag in this existing runner (`GOMP_critical_*`
undefined); adding the existing `openmp_flags` at the affected smoke link
fixed the harness, and the rerun passed. Its temporary build directory was
removed by the runner; both compact Slurm logs and the exact source/input
hashes are retained in `.snrt-admission-gate-20261106/`.

### Syntax cancellation and requeue check — 2026-11-05

Slurm accounting confirms that jobs 451122 (`chimes-ctx-2x2`) and 450842
(`chimes-rx-target`) were cancelled before execution (zero elapsed time, no
allocated node). Their replacement jobs are already in the queue: 451130 is
the A40 four-case fresh/warm × CPU/hybrid target replay, and 451146 is the
A40 integrated restart replay. Both currently show `QOSMaxGRESPerUser`, with
an estimated start of 2026-11-06 09:02 KST. No duplicate submission was made
during this history check. Jobs 450313 and 450676 were cancelled H200/H100
alternatives with corresponding completed runs 450314 (A100) and 450678 (A10),
so they do not need another copy.

Before job 451130 starts, the replay harness now snapshots CHIMES dispatch
counters immediately before and after the captured step-2 call. GPU calls,
CPU fallbacks, batches, errors and uploads are reported as target-local deltas;
peak concurrency, streams and worker mask remain explicitly process-cumulative
because the existing API exposes them as high-water/union counters. The binary
result schema is version 2, and the comparator rejects the old schema. This
keeps warm-up work from being mistaken for target-step dispatch. The submitted
Slurm wrapper reads the driver, make fragment and comparator from the project
tree at execution and records their hashes in its result directory.

The step-2 fixture's recorded `solid_charge` is exactly zero; the replay now
passes that captured value instead of hard-coding zero and checks that the
loaded molecular-table temperature ceiling matches the captured `tmax`. The
451068 effective namelist leaves `dust_relative_motion` false (also the source
default), so the production caller omitted the optional four-phase arrays;
the replay's group-only photo call matches that branch.

Scheduler recheck at 2026-11-06 00:07 KST still found both jobs pending under
`QOSMaxGRESPerUser`; Slurm moved the estimated start from 09:02 to 21:03 KST
on 2026-11-06. Their submit IDs and resource requests are unchanged.

### Final Slurm disposition — 2026-11-06

A fresh accounting query on `syntax` reports job 451130 (`chimes-ctx-a40`)
and job 451146 (`chimes-rx-target`) as terminal `FAILED`, exit `2:0`, on
`syn07`, with elapsed times 4 s and 29 s and reason `QOSMaxGRESPerUser`.
Neither allocated a usable compute interval or produced scientific evidence;
this supersedes the 00:07 pending estimate above. They were not resubmitted:
the same-allocation A/B and subsequent source review supplied no new testable
failure hypothesis, so another replay would repeat the diagnostic matrix
without resolving the remaining uncertainty. Preserve the strict gate and
leave the CHIMES restart cause open.
