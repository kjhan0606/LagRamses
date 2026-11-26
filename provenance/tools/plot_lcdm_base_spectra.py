#!/usr/bin/env python3
"""Plot the parent and Group 205 base-grid spectra against CAMB."""

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


source = Path("/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_base_camb_spectra_1578580")
out = Path(os.environ["RAR_PLOT_OUTPUT"])


def read(case, field):
    with (source / case / "shells.csv").open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row["comparison"] == field]
    rows = [row for row in rows if float(row["k_min_h_mpc"]) >= 1.0]
    return {
        name: np.array([float(row[name]) for row in rows])
        for name in ("k_min_h_mpc", "k_max_h_mpc", "measured_over_camb_power",
                     "music_input_over_camb_power", "gaussian_sample_sigma_on_ratio")
    }


plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 9})
fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.7), sharex=True, sharey=True)
panels = (("gas_density", "Gas density"), ("cdm_density", "CDM density"),
          ("gas_velocity", "Gas velocity"), ("cdm_velocity", "CDM velocity"))
for ax, (field, title) in zip(axes.ravel(), panels):
    parent = read("parent", field)
    zoom = read("zoom", field)
    x = np.sqrt(parent["k_min_h_mpc"] * parent["k_max_h_mpc"])
    ax.axhline(1.0, color="0.45", lw=0.8, ls="--")
    ax.plot(x, parent["music_input_over_camb_power"], color="#b34939",
            lw=1.2, ls=":", label="MUSIC input")
    ax.errorbar(x * 0.98, parent["measured_over_camb_power"],
                yerr=parent["gaussian_sample_sigma_on_ratio"], color="#245f96",
                fmt="o", ms=3.5, lw=0.8, capsize=2, label="Parent")
    ax.plot(x * 1.02, zoom["measured_over_camb_power"], color="#a26522",
            marker="^", ms=4, lw=1, label="Group 205 zoom")
    ax.text(0.97, 0.08, title, transform=ax.transAxes, ha="right", va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(0.95, 17)
    ax.set_ylim(0.85, 1.08)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_xticks([1, 2, 4, 8, 16], labels=["1", "2", "4", "8", "16"])
axes[0, 0].legend(loc="upper right", frameon=False)
axes[0, 0].set_ylabel(r"$P_{\rm IC}/P_{\rm CAMB}$")
axes[1, 0].set_ylabel(r"$P_{\rm IC}/P_{\rm CAMB}$")
for ax in axes[1]:
    ax.set_xlabel(r"$k\;(h\,{\rm Mpc}^{-1})$")
fig.subplots_adjust(left=0.12, right=0.985, bottom=0.11, top=0.98, wspace=0, hspace=0)
fig.savefig(out / "base_spectra.pdf")
fig.savefig(out / "base_spectra.png", dpi=180)
plt.close(fig)
