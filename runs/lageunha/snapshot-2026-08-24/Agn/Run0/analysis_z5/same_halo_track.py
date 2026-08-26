#!/usr/bin/env python3
"""Track central BH of the same halo across all simulations.

Anchor halo: position of CDM sid=10 (z=5.3 champion BH).
Strategy:
  1) Build CDM sid=10 position trajectory across all output_NNNNN
  2) For each other sim & each output, take CDM-sid10 position at the
     nearest z as the anchor → pick the most massive BH within R_link
     → that BH = "central BH of the same halo"
  3) Plot M_BH(z) for the same halo across all sims
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
R_LINK = 1.0       # Mpc/h comoving
ANCHOR_SID = 10    # CDM champion at z=5.3
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
    """Return list of dicts: [{z, sid, M, pos, Lbox}] per output."""
    rdir = os.path.join(BASE, run)
    outs = sorted(glob.glob(os.path.join(rdir, "output_*")))
    snaps = []
    for out in outs:
        nstr = out.split("output_")[-1]
        info_p = os.path.join(out, f"info_{nstr}.txt")
        csv_p  = os.path.join(out, f"sink_{nstr}.csv")
        if not (os.path.exists(info_p) and os.path.exists(csv_p)): continue
        try:
            info = parse_info(info_p)
            csv = np.loadtxt(csv_p, delimiter=",")
        except Exception: continue
        if csv.size == 0: continue
        if csv.ndim == 1: csv = csv[None,:]
        unit_m = info["unit_d"] * info["unit_l"]**3
        h = info["H0"]/100.0
        Lbox = info["unit_l"]/info["aexp"]/MPC * h
        z = 1.0/info["aexp"] - 1.0
        sid = csv[:,0].astype(int)
        M = csv[:,1] * unit_m / MSUN
        pos = csv[:, 2:5] * Lbox
        snaps.append(dict(z=z, sid=sid, M=M, pos=pos, Lbox=Lbox))
    return snaps


def main():
    print("collecting snapshots ...")
    data = {r: collect(r) for r in RUNS}

    # CDM sid=10 trajectory (z, pos)
    cdm_track = []   # list of (z, pos)
    for s in data["run_cdm"]:
        idx = np.where(s["sid"] == ANCHOR_SID)[0]
        if idx.size == 1:
            cdm_track.append((s["z"], s["pos"][idx[0]], s["M"][idx[0]]))
    cdm_track.sort(key=lambda t: -t[0])  # high z first
    print(f"\nCDM sid={ANCHOR_SID} trajectory ({len(cdm_track)} snapshots):")
    for z, p, M in cdm_track:
        print(f"  z={z:.3f}  pos=({p[0]:6.2f},{p[1]:6.2f},{p[2]:6.2f})  M={M:.3e} Msun")
    cdm_z = np.array([t[0] for t in cdm_track])
    cdm_pos = np.array([t[1] for t in cdm_track])

    # for each sim, for each output_z, find anchor position by
    # nearest CDM z, then most massive BH within R_LINK
    tracks = {}    # {run: [(z, M_BH, dist_to_anchor, sid)]}
    for run, snaps in data.items():
        tr = []
        for s in snaps:
            j = int(np.argmin(np.abs(cdm_z - s["z"])))
            anchor = cdm_pos[j]
            dists = per_dist(s["pos"], anchor, s["Lbox"])
            mask = dists < R_LINK
            if not mask.any():
                tr.append((s["z"], np.nan, np.nan, -1))
                continue
            jcand = np.where(mask)[0]
            kbest = jcand[np.argmax(s["M"][jcand])]
            tr.append((s["z"], float(s["M"][kbest]),
                       float(dists[kbest]), int(s["sid"][kbest])))
        tracks[run] = tr
        zs = [t[0] for t in tr]
        ms = [t[1] for t in tr]
        sids = [t[3] for t in tr]
        print(f"\n  {run:18s}  z={min(zs):.2f}-{max(zs):.2f}  "
              f"central M_BH at z~5.3: {[f'{m:.2e}' for z,m,_,_ in tr if 5.2<z<5.4]}")
        for t in tr:
            print(f"    z={t[0]:5.2f}  M={t[1]:8.2e}  Δ={t[2]:.4f} Mpc/h  sid={t[3]}")

    # ------- plot -------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for run in RUNS:
        tr = np.array([(t[0], t[1]) for t in tracks[run] if not np.isnan(t[1])])
        if len(tr) == 0: continue
        ax.plot(tr[:,0], tr[:,1], "o-", color=COLORS[run],
                label=LABELS[run], lw=2.0, ms=6)
    ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$M_{\rm BH}^{\rm central}$ [M$_\odot$]")
    ax.set_title("Central BH of CDM's most massive halo (z=5.3 anchor)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)

    ax = axes[1]
    # ratio to CDM
    cdm_arr = np.array([t[1] for t in tracks["run_cdm"]])
    cdm_zs  = np.array([t[0] for t in tracks["run_cdm"]])
    for run in RUNS:
        if run == "run_cdm": continue
        tr = np.array([(t[0], t[1]) for t in tracks[run] if not np.isnan(t[1])])
        if len(tr) == 0: continue
        # match to nearest CDM z
        ratios = []
        for zv, mv in tr:
            j = int(np.argmin(np.abs(cdm_zs - zv)))
            if not np.isnan(cdm_arr[j]) and cdm_arr[j] > 0:
                ratios.append((zv, mv/cdm_arr[j]))
        ratios = np.array(ratios)
        ax.plot(ratios[:,0], ratios[:,1], "o-", color=COLORS[run],
                label=LABELS[run], lw=2.0, ms=6)
    ax.axhline(1.0, color="black", lw=1, ls="--", alpha=0.6)
    ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$M_{\rm BH} / M_{\rm BH}^{\rm CDM}$")
    ax.set_title("Central BH mass relative to CDM (same halo)")
    ax.set_ylim(0.0, 1.4)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10)

    plt.tight_layout()
    out_png = os.path.join(BASE, "analysis_z5", "same_halo_track.png")
    plt.savefig(out_png, dpi=140)
    print(f"\n  → {out_png}")


if __name__ == "__main__":
    main()
