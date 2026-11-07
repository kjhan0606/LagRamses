# Matched native AMR disc force comparison

Numerical source f0a5129 connects the established fixed disc sources to the
native composite publication operation. The rho_ana implementation is imported
from explicitly mapped immutable source65eeb24 in
/home/kjhan/BackUP/rar/build/qumond_ngc2403_gpu_47a643d, not inferred from its
folder name. It provides the existing analytic stellar/gas sech-squared discs,
truncated NFW halo and tabulated surface profiles. An additional finite-parameter
check rejects nonfinite analytic disc inputs. Source-1 spherical behavior is
preserved. This is reuse of the previously selected controlled mass model,
within the reviewed spatial-force integration; no new physical model is chosen.

The opt-in native publication now accepts use_mond=false as a matched Newtonian
control, retaining all serial/static/two-level/restart/physics guards. It solves
precisely the same native Poisson operator. Phantom density is exactly zero;
the second potential/history slot mirrors N and does not represent another
physics solve. Log RAR_COMPOSITE_MODEL_V1 records use_mond and gravity_type;
exports have a0=0 and external=0 for Newtonian controls. Default history mode
without publication continues to require QUMOND. Ordinary/default physics paths
are unchanged.

The existing configs/qumond_disc_pilot inputs prescribe shared stellar
(M,R,h)=(0.06,0.055,0.025) and gas(0.02,0.085,0.035) in code units.
Newtonian baryon+halo alone adds truncated NFW(M,rs,rt)=(0.4,0.04,0.35).
Baryon-only Newtonian is a diagnostic for gbar, not the Lambda-CDM comparison.
QUMOND solves the full two-stage field on baryons without a CDM halo.
All cases use identical periodic fixed native topology and a10kpc box. The halo
is a prescribed controlled comparator, not a cosmological realization or fitted
observational halo. These are source prescriptions, not live galaxies.

Grammar build1275755 completed in2m48s. Matched execution1275759 uses levels5/6
and completed in28s; execution1275760 uses levels6/7 and completed in5m30s.
Each case deliberately uses nstepmax=1,noutput=1,aout=1.1,tout=1e100,foutput=1,
fbackup=1e9, with exactly one expected full snapshot. Three cases reserve2GiB
including leaf/face exports (snapshot estimate273318792bytes), and require
>3GiB free before launch. Effective nml paths, source revisions, executable and
input/output hashes are preserved under /gpfs/kjhan/RAR_QUMOND/runs/
rar_native_disc_JOBID/CASE. No login-host calculation or previous output
replacement is involved.

Dependent audit1275827 independently reconstructs geometry, physical density,
both potential residuals, phantom density and forces for all six cases. It
requires matched native cells/faces, bitwise equal baryon density for N/Q, and
matching QUMOND first-stage/N-baryon potentials. Rotation curves use32 common
radii and128 azimuths, symmetric finest z planes averaged to the midplane and
bilinear xy interpolation. Radius minimum is four low-resolution fine widths;
maximum is0.25box. Saved CSV includes physical accelerations, circular speeds,
resolution changes, numerical azimuth variation and full-QUMOND departure from
the algebraic simple-nu radial boost. That algebraic formula is a diagnostic,
not used to generate QUMOND forces or treated as an independent MOND prediction.

The figure presents the high-resolution rotation curves, gbar/gtot relation and
force changes with resolution. Numerical azimuth variation is not population
RAR scatter. Periodic images/background, finite thickness and resolution remain
explicit. No galaxy formation, observed RAR/model-selection, evolving MPI,
restart, or full-project completion is claimed.

Source regression execution1275828 and audits1275829/1275830 all completed.
Both artifact manifests verify; default-off gravity arrays remain bitwise equal and the independent spherical publication check passes.
Matched audit1275827 completed in41s and its full artifact manifest verifies. All six independent publication verifiers pass. At high resolution the QUMOND N/Q residuals are9.82737191350802e-10 /9.88181467168656e-10, phantom relative difference2.444518145032272e-15, force relative difference4.8272154247543485e-17 and analytic-source relative difference2.3086181761171936e-16. Gas momentum remains exactly zero and energy unchanged to1.13e-16 relative. Baryon density is bitwise equal between the N-baryon/QUMOND cases, and their first-stage potentials are exactly identical at both resolutions. The next work uses those results to assess compact disc forces
and then connects actual source evolution; it must not replace the full hydro,
BIMOND/lagCAMB and FDM objectives with this fixed-source comparison.


## Measured disc comparison

The32 common radii cover0.625--2.5kpc. High/low force differences have RMS/max
0.981%/2.491% for Newtonian baryons,1.351%/3.276% for baryons+CDM halo, and
0.888%/2.289% for QUMOND. These are resolution differences, not estimated exact
continuum errors. Sampled baryon source integrals change from
0.07965881582647566 to0.07971474203751605; halo-inclusive integrals from
0.47699264982694833 to0.4789051372804839, all including the uniform1e-6 tracer.
The maximum change lies at the innermost sampled radii. This must remain in
interpretation rather than being erased by selecting favorable radii.

The imposed CDM-halo acceleration is1.779--4.037 times the QUMOND acceleration
over this interval. That exceeds the measured resolution changes but reflects
this chosen, unfitted halo; it is not observational model preference. Full
QUMOND differs from the radial algebraic simple-nu boost by mean-0.01668dex,
RMS0.01810dex and maximum absolute0.03252dex. Finite disc geometry, periodic
prescription and finite resolution remain part of this diagnostic, so it is not
a clean observational discriminator or intrinsic RAR scatter.

Saved compact results: rar_native_disc_result_1275827.json and
rar_native_disc_rotation_rar.csv. Source and numerical verification outputs
remain in /gpfs/kjhan/RAR_QUMOND/audits/rar_native_disc_1275827. Initial figure
layout overlapped the CDM curve and titles; rendering-only job1275859 reads
the unchanged CSV/results and moves the labels into clear areas, without
rerunning scientific calculations. Its artifacts verify and final PNG was
visually inspected. Final PDF: figures/rar_native_disc_rotation_rar.pdf.
The original audit and figure remain preserved. Runners, source oracle,
analysis and independent rendering script are recorded under provenance/tools.
The manuscript checkout already had unrelated local changes and was not edited.
