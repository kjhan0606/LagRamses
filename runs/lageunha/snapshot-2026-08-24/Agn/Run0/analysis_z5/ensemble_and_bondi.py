#!/usr/bin/env python3
"""Two-pronged analysis for the SIDM × AGN signal at z=5.3:

(1) UNIVERSALITY across top-N halos:
    - Pick CDM top-10 BHs at z=5.3 as anchor halos.
    - Trace each anchor's position through all CDM snapshots.
    - In each sim & each z, find central BH (most massive within R_link)
      → record M_BH and Mdot.
    - Plot ensemble M_BH/M_CDM(z) per sim (individual tracks + median).

(2) BONDI vs WANDERING diagnostic:
    - For the same anchor halos, plot Mdot(z) ensemble.
    - Plot offset distance |Δr| of central BH from CDM anchor as
      a wandering proxy.
    - If M_BH suppression matches Mdot suppression but Δr stays small
      → Bondi-rate suppression (gas thermodynamics / DM-modified gas
        accumulation), NOT wandering.
"""
import os, glob, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/gpfs/kjhan/Hydro/Sidm/Agn/Run0"
RUNS = ["run_cdm","run_sidm","run_sidm5","run_sidm10",
        "run_isidm","run_isidm_endo","run_dsidm"]
MSUN = 1.989e33
MPC = 3.086e24
YR = 3.156e7
N_ANCHOR = 10        # top-10 halos
R_LINK = 1.0         # Mpc/h comoving
COLORS = dict(run_cdm="black", run_sidm="C0", run_sidm5="C1",
              run_sidm10="C3", run_isidm="C2", run_isidm_endo="C4",
              run_dsidm="C5")
LABELS = dict(run_cdm="CDM", run_sidm=r"SIDM $\sigma$=1",
              run_sidm5=r"SIDM $\sigma$=5", run_sidm10=r"SIDM $\sigma$=10",
              run_isidm="iSIDM (f=0.3)", run_isidm_endo="iSIDM (f=0)",
              run_dsidm="dSIDM")


def parse_info(path):
    out = {}
    with open(path) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip(); v = v.split()[0].strip()
                if k in ("unit_l","unit_d","unit_t","aexp","H0","boxlen"):
                    try: out[k] = float(v)
                    except ValueError: pass
    return out


def per_dist(p1, p2, Lbox):
    d = p1 - p2
    d = d - Lbox * np.round(d/Lbox)
    return np.sqrt((d**2).sum(axis=-1))


def collect(run):
    rdir = os.path.join(BASE, run)
    outs = sorted(glob.glob(os.path.join(rdir, "output_*")))
    snaps = []
    for out in outs:
        nstr = out.split("output_")[-1]
        info_p = os.path.join(out, f"info_{nstr}.txt")
        csv_p  = os.path.join(out, f"sink_{nstr}.csv")
        if not (os.path.exists(info_p) and os.path.exists(csv_p)): continue
        try:
            info = parse_info(info_p); csv = np.loadtxt(csv_p, delimiter=",")
        except Exception: continue
        if csv.size == 0: continue
        if csv.ndim == 1: csv = csv[None,:]
        unit_m = info["unit_d"] * info["unit_l"]**3
        unit_mdot = unit_m / info["unit_t"]
        h = info["H0"]/100.0
        Lbox = info["unit_l"]/info["aexp"]/MPC * h
        z = 1.0/info["aexp"] - 1.0
        snaps.append(dict(
            z=z, Lbox=Lbox,
            sid=csv[:,0].astype(int),
            M=csv[:,1] * unit_m / MSUN,
            pos=csv[:, 2:5] * Lbox,
            Mdot=csv[:,9] * unit_mdot / (MSUN/YR)))
    return snaps


def trace_anchor(cdm_snaps, target_sid):
    """Return arrays z, pos for CDM anchor BH (target_sid)."""
    zs, ps = [], []
    for s in cdm_snaps:
        idx = np.where(s["sid"] == target_sid)[0]
        if idx.size == 1:
            zs.append(s["z"]); ps.append(s["pos"][idx[0]])
    if not zs: return None
    order = np.argsort(zs)
    return np.array(zs)[order], np.array(ps)[order]


def find_central(snap, anchor_pos, Lbox, R):
    dists = per_dist(snap["pos"], anchor_pos, Lbox)
    mask = dists < R
    if not mask.any(): return None
    j_in = np.where(mask)[0]
    k = j_in[int(np.argmax(snap["M"][j_in]))]
    return dict(M=float(snap["M"][k]), Mdot=float(snap["Mdot"][k]),
                pos=snap["pos"][k], dist=float(dists[k]),
                sid=int(snap["sid"][k]))


def main():
    print("collecting data ...")
    data = {r: collect(r) for r in RUNS}

    # Pick CDM top-10 anchors at z=5.3 (output_00010)
    cdm_snaps = data["run_cdm"]
    cdm_z = np.array([s["z"] for s in cdm_snaps])
    j_z53 = int(np.argmin(np.abs(cdm_z - 5.31)))
    snap53 = cdm_snaps[j_z53]
    rank = np.argsort(-snap53["M"])[:N_ANCHOR]
    anchor_sids = snap53["sid"][rank].tolist()
    print(f"\nCDM top-{N_ANCHOR} BHs at z={snap53['z']:.2f}:")
    for i, sid in enumerate(anchor_sids):
        k = int(np.where(snap53["sid"]==sid)[0][0])
        print(f"  #{i+1}  sid={sid:6d}  M={snap53['M'][k]:8.3e} Msun  "
              f"pos=({snap53['pos'][k][0]:6.2f},{snap53['pos'][k][1]:6.2f},{snap53['pos'][k][2]:6.2f})")

    # Anchor trajectories (CDM)
    anchors = {}
    for sid in anchor_sids:
        tr = trace_anchor(cdm_snaps, sid)
        if tr is not None: anchors[sid] = tr
    print(f"\n{len(anchors)}/{N_ANCHOR} anchors have CDM trajectories")

    # For each sim, for each anchor, build M_BH(z), Mdot(z), Δr(z)
    # tracks[run][anchor_sid] = {'z':[], 'M':[], 'Mdot':[], 'dr':[]}
    tracks = {r: {} for r in RUNS}
    for run, snaps in data.items():
        for sid, (a_z, a_pos) in anchors.items():
            zs, Ms, Mdots, drs = [], [], [], []
            for s in snaps:
                j = int(np.argmin(np.abs(a_z - s["z"])))
                anc = a_pos[j]
                cen = find_central(s, anc, s["Lbox"], R_LINK)
                if cen is None: continue
                zs.append(s["z"]); Ms.append(cen["M"])
                Mdots.append(cen["Mdot"]); drs.append(cen["dist"])
            if zs:
                tracks[run][sid] = dict(z=np.array(zs), M=np.array(Ms),
                                        Mdot=np.array(Mdots), dr=np.array(drs))

    # Build common z grid for ensembles
    z_grid = np.array([s["z"] for s in cdm_snaps])
    z_grid.sort()

    def at_grid(arr_z, arr_v):
        """Sample arr_v at z_grid via nearest-z within 0.05 tolerance."""
        out = np.full_like(z_grid, np.nan, dtype=float)
        for i, zg in enumerate(z_grid):
            j = int(np.argmin(np.abs(arr_z - zg)))
            if abs(arr_z[j] - zg) < 0.05:
                out[i] = arr_v[j]
        return out

    # Per-run ensemble of ratios M/M_CDM(same anchor, same z)
    ratios = {}     # ratios[run] = (Nanchor x Nz) array
    Ms_grid = {}; Mdots_grid = {}; drs_grid = {}
    for run in RUNS:
        rA = []; mA = []; mdA = []; drA = []
        for sid in anchors:
            if sid not in tracks[run] or sid not in tracks["run_cdm"]:
                continue
            tr_run = tracks[run][sid]; tr_cdm = tracks["run_cdm"][sid]
            mr = at_grid(tr_run["z"], tr_run["M"])
            mc = at_grid(tr_cdm["z"], tr_cdm["M"])
            mdr = at_grid(tr_run["z"], tr_run["Mdot"])
            drr = at_grid(tr_run["z"], tr_run["dr"])
            with np.errstate(divide='ignore', invalid='ignore'):
                rat = mr / mc
            rA.append(rat); mA.append(mr); mdA.append(mdr); drA.append(drr)
        ratios[run] = np.array(rA) if rA else np.zeros((0,len(z_grid)))
        Ms_grid[run] = np.array(mA) if mA else np.zeros((0,len(z_grid)))
        Mdots_grid[run] = np.array(mdA) if mdA else np.zeros((0,len(z_grid)))
        drs_grid[run] = np.array(drA) if drA else np.zeros((0,len(z_grid)))

    # ---------------- plot 1: universality (M_BH/M_CDM ensemble) -----
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))

    ax = axes[0,0]
    for run in RUNS:
        if run == "run_cdm": continue
        R = ratios[run]
        if R.size == 0: continue
        # individual tracks (alpha low)
        for i in range(R.shape[0]):
            ax.plot(z_grid, R[i,:], "-", color=COLORS[run], alpha=0.18, lw=0.8)
        med = np.nanmedian(R, axis=0)
        ax.plot(z_grid, med, "o-", color=COLORS[run],
                label=f"{LABELS[run]}", lw=2.2, ms=6)
    ax.axhline(1, color="black", lw=1, ls="--", alpha=0.5)
    ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$M_{\rm BH} / M_{\rm BH}^{\rm CDM}$ (same halo)")
    ax.set_title(f"Universality: top-{N_ANCHOR} halo ensemble")
    ax.set_ylim(0, 1.6); ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="lower left")

    # ---------------- plot 2: Mdot(z) median ensemble ----------------
    ax = axes[0,1]
    for run in RUNS:
        Md = Mdots_grid[run]
        if Md.size == 0: continue
        med = np.nanmedian(Md, axis=0)
        p16 = np.nanpercentile(Md, 16, axis=0)
        p84 = np.nanpercentile(Md, 84, axis=0)
        ax.fill_between(z_grid, p16, p84, color=COLORS[run], alpha=0.13)
        ax.plot(z_grid, med, "o-", color=COLORS[run], label=LABELS[run],
                lw=2.0, ms=6)
    ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$\dot M_{\rm BH}$ [M$_\odot$/yr]")
    ax.set_title("Central BH accretion rate (median ± 16-84%)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    # ---------------- plot 3: Mdot ratio to CDM -------------------
    ax = axes[1,0]
    cdm_md = Mdots_grid["run_cdm"]
    cdm_med_md = np.nanmedian(cdm_md, axis=0)
    for run in RUNS:
        if run == "run_cdm": continue
        Md = Mdots_grid[run]
        if Md.size == 0: continue
        # per-anchor ratio
        Rmd = []
        for i in range(Md.shape[0]):
            sid = list(anchors.keys())[i] if i < len(anchors) else None
            # find same anchor in cdm
            if sid is None: continue
            with np.errstate(divide='ignore', invalid='ignore'):
                Rmd.append(Md[i,:] / cdm_md[i,:] if cdm_md.shape[0]>i else np.nan)
        if not Rmd: continue
        Rmd = np.array(Rmd)
        med = np.nanmedian(Rmd, axis=0)
        ax.plot(z_grid, med, "o-", color=COLORS[run], label=LABELS[run],
                lw=2.0, ms=6)
    ax.axhline(1, color="black", lw=1, ls="--", alpha=0.5)
    ax.invert_xaxis(); ax.set_yscale("log")
    ax.set_xlabel("z"); ax.set_ylabel(r"$\dot M_{\rm BH} / \dot M_{\rm BH}^{\rm CDM}$")
    ax.set_title("Accretion-rate ratio (same halo)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    ax.set_ylim(0.1, 3)

    # ---------------- plot 4: wandering Δr -----------------------
    ax = axes[1,1]
    for run in RUNS:
        DR = drs_grid[run]
        if DR.size == 0: continue
        med = np.nanmedian(DR, axis=0) * 1000   # → kpc/h comoving
        p84 = np.nanpercentile(DR, 84, axis=0) * 1000
        ax.fill_between(z_grid, np.zeros_like(med), p84, color=COLORS[run], alpha=0.10)
        ax.plot(z_grid, med, "o-", color=COLORS[run], label=LABELS[run],
                lw=2.0, ms=6)
    ax.invert_xaxis(); ax.set_yscale("log")
    ax.set_xlabel("z"); ax.set_ylabel(r"$|\Delta r|$ from CDM anchor [kpc/h, comoving]")
    ax.set_title("BH offset (wandering proxy)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    plt.tight_layout()
    out_png = os.path.join(BASE, "analysis_z5", "ensemble_and_bondi.png")
    plt.savefig(out_png, dpi=140)
    print(f"\n  → {out_png}")

    # ---------------- print summary at z=5.3 ----------------------
    iz53 = int(np.argmin(np.abs(z_grid - 5.31)))
    print("\n=== Summary at z=5.31 (median over top-10 halos) ===")
    print(f"{'Run':18s}  {'M_BH':>10s}  {'M/M_cdm':>8s}  {'Mdot':>10s}  "
          f"{'Mdot/Mdot_cdm':>14s}  {'Δr [kpc/h]':>14s}")
    for run in RUNS:
        Mm = np.nanmedian(Ms_grid[run][:,iz53]) if Ms_grid[run].size else np.nan
        Mr = np.nanmedian(ratios[run][:,iz53]) if ratios[run].size else np.nan
        Md = np.nanmedian(Mdots_grid[run][:,iz53]) if Mdots_grid[run].size else np.nan
        DR = np.nanmedian(drs_grid[run][:,iz53])*1000 if drs_grid[run].size else np.nan
        if run == "run_cdm":
            Mdr = 1.0
        else:
            cdm_md_z = cdm_md[:, iz53]
            run_md_z = Mdots_grid[run][:, iz53]
            with np.errstate(divide='ignore', invalid='ignore'):
                Mdr = np.nanmedian(run_md_z / cdm_md_z)
        print(f"{LABELS[run]:18s}  {Mm:10.3e}  {Mr:8.3f}  {Md:10.3e}  "
              f"{Mdr:14.3f}  {DR:14.3f}")


if __name__ == "__main__":
    main()
