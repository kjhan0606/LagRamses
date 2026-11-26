#!/usr/bin/env python3
"""Plot matched level-11 patch power ratios."""

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


source = Path("/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_g205_z99_windowed_spectra_1578622")
out = Path(os.environ["RAR_PLOT_OUTPUT"])
output_format = os.environ["RAR_PLOT_FORMAT"]
if output_format not in ("pdf", "png"):
    raise ValueError("RAR_PLOT_FORMAT must be pdf or png")
with (source / "shells.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))

plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 9})
fig, ax = plt.subplots(figsize=(6.3, 4.2))
for field, label, color, marker in (
    ("gas_density", "Gas density", "#245f96", "o"),
    ("gas_velocity", "Gas velocity", "#bb632a", "^"),
    ("cdm_velocity", "CDM velocity", "#43815b", "s"),
):
    data = [row for row in rows if row["field"] == field]
    k = np.array([np.sqrt(float(row["k_min_h_mpc"]) * float(row["k_max_h_mpc"]))
                  for row in data])
    ratio = np.array([float(row["extended_over_baseline_power"]) for row in data])
    ax.plot(k, ratio, marker=marker, ms=4, lw=1.25, color=color, label=label)

ax.axhline(1, color="0.55", ls="--", lw=0.8)
ax.axvline(155.8788, color="0.45", ls=":", lw=1)
ax.axvline(257.35927, color="0.45", ls="-.", lw=1)
ax.text(155.8788 * 0.98, 0.59, "Baseline table end",
        rotation=90, ha="right", va="bottom", fontsize=9, color="0.35")
ax.text(257.35927 * 0.98, 0.59, "Axial Nyquist",
        rotation=90, ha="right", va="bottom", fontsize=9, color="0.35")
ax.text(0.04, 0.08, r"Level 11, $z=99$", transform=ax.transAxes,
        ha="left", va="bottom")
ax.set_xscale("log")
ax.set_xlim(16, 450)
ax.set_ylim(0.57, 1.05)
ax.set_xticks([20, 40, 80, 160, 320], labels=["20", "40", "80", "160", "320"])
ax.set_xlabel(r"$k\;(h\,{\rm Mpc}^{-1})$")
ax.set_ylabel(r"$P_{\rm extended}/P_{\rm baseline}$")
ax.tick_params(direction="in", top=True, right=True)
ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.15), frameon=False)
fig.subplots_adjust(left=0.15, right=0.98, bottom=0.13, top=0.98)
fig.savefig(out / f"fine_spectra.{output_format}", dpi=180)
plt.close(fig)
