#!/usr/bin/env python3
"""AGN bolometric luminosity function at z~5 across SIDM/AGN runs.

Reads sink_00010.csv (RAMSES sink output) and info_00010.txt for code units.
Columns of sink CSV (output_sink.f90:222):
  id, msink, x, y, z, vx, vy, vz, age, dMBHoverdt   (all code units)

L_bol = eps_r * dMdt * c^2   (eps_r = 0.1)
LF: dN / dlogL / V_comoving [(Mpc/h)^-3 dex^-1]
"""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/gpfs/kjhan/Hydro/Sidm/Agn/Run0"
RUNS = ["run_cdm","run_sidm","run_sidm5","run_sidm10",
        "run_isidm","run_isidm_endo","run_dsidm"]   # idm omitted (only z=10.21)
SNAP = "output_00010"
EPS_R = 0.1
MSUN = 1.989e33   # g
CLIGHT = 2.998e10 # cm/s
MPC = 3.086e24    # cm
YR = 3.156e7      # s


def parse_info(path):
    """Return dict of unit_l, unit_d, unit_t, aexp, H0."""
    out = {}
    with open(path) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                v = v.split()[0].strip()
                if k in ("unit_l","unit_d","unit_t","aexp","H0",
                         "boxlen","omega_m","omega_l","omega_b"):
                    try:
                        out[k] = float(v)
                    except ValueError:
                        pass
    return out


def load_sink(run):
    info = parse_info(os.path.join(BASE, run, SNAP, "info_00010.txt"))
    csv = np.loadtxt(os.path.join(BASE, run, SNAP, "sink_00010.csv"),
                     delimiter=",")
    if csv.ndim == 1:
        csv = csv[None, :]
    msink_code  = csv[:, 1]
    mdot_code   = csv[:, 9]
    unit_m   = info["unit_d"] * info["unit_l"]**3            # g
    unit_mdot= unit_m / info["unit_t"]                       # g/s
    M_BH     = msink_code * unit_m / MSUN                    # Msun
    Mdot     = mdot_code  * unit_mdot                        # g/s
    L_bol    = EPS_R * Mdot * CLIGHT**2                      # erg/s
    L_edd    = 1.26e38 * M_BH                                # erg/s
    h = info["H0"] / 100.0
    aexp = info["aexp"]
    Lbox_com_Mpch = info["unit_l"] / aexp / MPC * h
    Vbox = Lbox_com_Mpch**3
    return dict(M_BH=M_BH, Mdot=Mdot, L_bol=L_bol, L_edd=L_edd,
                Vbox=Vbox, aexp=aexp, z=1/aexp-1, Lbox_Mpch=Lbox_com_Mpch,
                run=run)


def lf(L, V, edges):
    L = L[L > 0]
    H, _ = np.histogram(np.log10(L), bins=edges)
    dlog = np.diff(edges)
    return H / V / dlog, np.sqrt(H) / V / dlog


def main():
    data = {}
    for run in RUNS:
        path = os.path.join(BASE, run, SNAP)
        if not os.path.isdir(path):
            print(f"  [skip] {run}: no {SNAP}")
            continue
        try:
            d = load_sink(run)
        except Exception as e:
            print(f"  [error] {run}: {e}")
            continue
        nact = (d["Mdot"] > 0).sum()
        print(f"  {run:20s}  z={d['z']:.3f}  N_sink={len(d['M_BH']):6d} "
              f"N_active={nact:5d}  Lbox={d['Lbox_Mpch']:.2f} Mpc/h "
              f"L_max={d['L_bol'].max():.2e} erg/s")
        data[run] = d

    edges = np.arange(42.0, 47.6, 0.25)
    centers = 0.5*(edges[:-1]+edges[1:])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = dict(run_cdm="black", run_sidm="C0", run_sidm5="C1",
                  run_sidm10="C3", run_isidm="C2", run_isidm_endo="C4",
                  run_dsidm="C5")
    labels = dict(run_cdm="CDM", run_sidm="SIDM σ=1", run_sidm5="SIDM σ=5",
                  run_sidm10="SIDM σ=10", run_isidm="iSIDM (f=0.3)",
                  run_isidm_endo="iSIDM (f=0)", run_dsidm="dSIDM")

    # uncapped
    ax = axes[0]
    for run, d in data.items():
        phi, ephi = lf(d["L_bol"], d["Vbox"], edges)
        m = phi > 0
        ax.errorbar(centers[m], phi[m], yerr=ephi[m],
                    fmt="o-", color=colors[run], label=labels[run], lw=1.5)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10}(L_{\rm bol}\,[{\rm erg\,s^{-1}}])$")
    ax.set_ylabel(r"$\Phi$ [(Mpc/h)$^{-3}$ dex$^{-1}$]")
    ax.set_title(f"AGN LF at z≈5.31 ($L_{{\\rm bol}}=0.1\\,\\dot M\\,c^2$)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    # Edd-capped
    ax = axes[1]
    for run, d in data.items():
        Lcap = np.minimum(d["L_bol"], d["L_edd"])
        phi, ephi = lf(Lcap, d["Vbox"], edges)
        m = phi > 0
        ax.errorbar(centers[m], phi[m], yerr=ephi[m],
                    fmt="s-", color=colors[run], label=labels[run], lw=1.5)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10}(L_{\rm bol}\,[{\rm erg\,s^{-1}}])$")
    ax.set_ylabel(r"$\Phi$ [(Mpc/h)$^{-3}$ dex$^{-1}$]")
    ax.set_title(r"AGN LF (Eddington-capped)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()
    out_png = os.path.join(BASE, "analysis_z5", "agn_lf_z5.png")
    plt.savefig(out_png, dpi=140)
    print(f"\n  → {out_png}")

    # text dump
    out_txt = os.path.join(BASE, "analysis_z5", "agn_lf_z5.txt")
    with open(out_txt, "w") as f:
        f.write("# logL_bin_center  ")
        for run in data: f.write(f"{labels[run]:>14s}  ")
        f.write("\n")
        phis = {run: lf(d["L_bol"], d["Vbox"], edges)[0] for run, d in data.items()}
        for i, c in enumerate(centers):
            f.write(f"{c:6.2f}  ")
            for run in data: f.write(f"{phis[run][i]:14.4e}  ")
            f.write("\n")
    print(f"  → {out_txt}")


if __name__ == "__main__":
    main()
