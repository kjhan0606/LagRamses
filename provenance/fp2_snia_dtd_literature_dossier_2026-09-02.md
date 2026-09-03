# F-P2 SNIa DTD literature dossier — 2026-09-02

Status: **candidate review only; no DTD or SNIa event source selected**.

## Decision-relevant findings

Field-galaxy reconstructions support a continuous DTD close to a power law.
Maoz, Mannucci & Brandt (2012) report a power-law index of about `-1.07`
with an integrated production of about `1.30e-3` SNe Ia per formed solar mass,
under the convention used in their analysis. Their result is a useful field
baseline, but it does not uniquely identify the binary progenitor channel or
remove the need to define the project's IMF and binary normalization.

The cluster analysis of Freundlich & Maoz (2021), using extended star-formation
histories, also finds a slope near `-1.1`, but a higher amplitude for cluster
populations. That amplitude cannot be transplanted to the LRD calculation
without an explicit environment/binary/IMF model.

Dubay, Johnson & Johnson (2024) compare DTD forms inside Milky-Way chemical
evolution models and find that an extended DTD with fewer prompt events than a
fiducial `t^-1` law can fit the abundance constraints better. This is a shape
systematic, not a directly portable event-rate normalization.

Binary population-synthesis calculations show why the DTD cannot be treated as
an ordinary 3--8 M☉ IMF yield channel: common-envelope, mass-transfer,
initial-binary distributions, and progenitor channel choices alter both the
shape and normalization. Some double-degenerate branches approach `t^-1`, but
the model parameters remain part of the physical source definition.

## Candidate disposition

The machine-readable comparison is in
`simulation/snrt/config/fp2_snia_dtd_candidate_matrix_v1.json`.

- **Baseline candidate for sensitivity design:** field observational power
  law, not yet approved.
- **Systematic comparison:** extended DTD with fewer prompt events, not yet
  normalized for this project.
- **Environmental comparison:** cluster power law, not a default LRD source.
- **Physical population alternative:** explicit binary population synthesis,
  requiring project-approved binary distributions and common-envelope physics.

No candidate is written into the production DTD contract. The contract keeps
`minimum_delay_gyr`, `maximum_delay_gyr`, and `events_per_initial_msun` null;
the SNIa yield source, event energy, momentum, and composition remain null as
well. The current `alpha=-1` in the executable test is a mathematical fixture,
not a project physics choice.

## Required approval inputs

Before F-P2 can activate, the project must pin the population basis and IMF,
binary fraction and distributions, metallicity dependence, minimum/maximum
delay, event normalization per initial SSP mass, event-yield source and
license/checksum, energy and momentum semantics, decay/composition convention,
and conversion code checksum. The event ledger and runtime gate already refuse
activation if these inputs are absent.

## Primary literature

- [Maoz, Mannucci & Brandt 2012](https://academic.oup.com/mnras/article/426/4/3282/1017965)
- [Freundlich & Maoz 2021](https://academic.oup.com/mnras/article/502/4/5882/6145511)
- [Dubay, Johnson & Johnson 2024](https://arxiv.org/abs/2404.08059)
- [Ruiter, Belczynski & Fryer 2011](https://academic.oup.com/mnras/article/417/1/408/979905)
