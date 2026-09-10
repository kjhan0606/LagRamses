# Metallic Fe kinetics: pinned physical inputs

Primary references:

- Geometric accretion comparison: [Zhukovska, Henning & Dobbs 2018](https://arxiv.org/abs/1803.01929).
- Fe thermal erosion fit: [Choban et al. 2026, section 3.5.2/Table 3](https://doi.org/10.1093/mnras/stag1020).
- Underlying Fe calculation: [Nozawa, Kozasa & Habe 2006, Figure 2](https://arxiv.org/abs/astro-ph/0605193).

Local source captures in `external/g2_candidates/choban2026_fe_kinetics/`:

- `article.html`, downloaded from `https://oup.silverchair-cdn.com/article-minimal/8698771`,
  SHA256 `ad7c9c828f1478103814f0b0da33b5f24e83a505fcf95016dfb3bac8ed3bd204`.
- `nozawa2006.pdf`, downloaded from `https://arxiv.org/pdf/astro-ph/0605193`,
  SHA256 `01afee4dadacc02e8cfeaefc017b8afb5f89ad133e2391e4ee06b677c8fad998`.

The Fe row gives coefficients (constant first)
`[-156.88,82.110,-18.238,2.0692,-0.11933,0.0027788]` for
`log10(Y/[micron yr^-1 cm^3])` versus `log10(T/K)`.
Native conversion uses `1e-4/31557600` to cm^4/s, `n_H` from the
gas hydrogen budget, and `3*|da/dt|/a` for the fixed-radius mass closure.
No clumping enhancement (`C2=1`). The projectile mixture is the source's
low-metallicity mixture, not a composition-dependent sputtering calculation.

The displayed Nozawa figure spans 1e4--1e10 K. This comparison admits only
1e4--1e9 K for the fit, explicitly neglects erosion below 1e4 K, and rejects
higher temperatures when erosion is active and Fe seeds exist. Do not
extrapolate the polynomial to cold gas. The 1e4 K cutoff is a declared
approximation, not a physical threshold. At 1e6 K the polynomial is
`-6.8845312` (decimal evaluation independent of the Fortran Horner evaluation).

This is not the complete cited galaxy model. It uses fixed 0.01/0.1 micron
carriers already bound to the live Fe optical bank, neutral/unfocused
accretion with explicit sticking, and no Fe coagulation/shattering,
nonthermal sputtering, Coulomb focusing, grain charging, nucleation,
adsorption latent heat, or magnetic opacity. Sensible dust+gas energy and
donor momentum are conserved. No unresolved SN shock efficiency is borrowed.
