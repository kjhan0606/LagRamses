# Matched level-11 initial-field spectra at z = 99

Grammar job 1578622 compared two saved Group 205 ΛCDM initial conditions. Both use the same 25 Mpc/h box, candidate region, redshift, random seeds, and level-11 geometry. The baseline CAMB table ends at 155.88 h/Mpc. The extended table reaches 595.78 h/Mpc and uses 200 samples per logarithmic interval. Sampling and range both differ between the tables.

The selected level-11 rectangular array contains 326 × 286 × 188 cells. A separable Hann window and a window-weighted mean subtraction were applied to each gas density and gas/CDM velocity component before the Fourier transform. Shell powers sum the three velocity components. The real-space and rFFT energies agree under a Parseval check. [Shell measurements](shells.csv)

| Wavenumber band (h/Mpc) | Extended / baseline gas density power | Extended / baseline gas velocity power | Extended / baseline CDM velocity power |
| --- | ---: | ---: | ---: |
| 96–128 | 0.99989 | 0.99986 | 0.99999 |
| 155.88–200 | 0.99055 | 0.98863 | 0.99866 |
| 200–240 | 0.95802 | 0.94527 | 0.99399 |
| 240–257.36 | 0.92580 | 0.90107 | 0.98952 |

The baseline table's upper endpoint falls inside this range. The MUSIC CAMB plugin extrapolates the baseline table beyond that endpoint. The level-11 axial Nyquist wavenumber is 257.36 h/Mpc. At 200–240 h/Mpc the extended-table realization has 4.2% less gas density power and 5.5% less gas velocity power than the baseline realization. The corresponding CDM velocity difference is 0.6%. The measurements are ratios of tapered patch spectra. The window mixes modes and includes padding. The table's pointwise input-power ratio at a shell midpoint is recorded for context and is not a window-convolved prediction.

The figure plots the measured shell ratios, with the baseline table endpoint and axial Nyquist marked. Wavenumbers above the axial Nyquist use only available diagonal modes. Grammar job 1578633 rendered the figure at revision `a60eb1d`. [Fine-grid power comparison PDF](fine_spectra.pdf)

This result describes one matched z = 99 Group 205 realization. It does not isolate table sampling from extrapolation, certify CAMB solver accuracy, measure AMR leaf spectra, or establish the Group 434 initial fields. The planned Group 434 table extends to 1,525 h/Mpc and covers the level-12 diagonal Nyquist wavenumber of 891.5 h/Mpc. Its realized fine-grid spectrum remains to be measured.
