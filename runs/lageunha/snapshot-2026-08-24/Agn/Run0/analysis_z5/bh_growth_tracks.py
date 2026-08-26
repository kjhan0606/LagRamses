#!/usr/bin/env python3
"""Trace central BH growth M_BH(z) per simulation.

For each run:
  - read every output_NNNNN/sink_NNNNN.csv
  - build {sink_id: [(z, M_BH_Msun, Mdot_Msun_yr), ...]}

Plots:
  1) Champion BH (max M_BH at each z) per sim
  2) Mean M_BH of top-100 BHs at each z
  3) Top-5 most massive BHs at z~5.3 (output_00010) traced BACK to seeding
"""
import os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/gpfs/kjhan/Hydro/Sidm/Agn/Run0"
RUNS = ["run_cdm","run_sidm","run_sidm5","run_sidm10",
        "run_isidm","run_isidm_endo","run_dsidm"]
MSUN = 1.989e33
COLORS = dict(run_cdm="black", run_sidm="C0", run_sidm5="C1",
              run_sidm10="C3", run_isidm="C2", run_isidm_endo="C4",
              run_dsidm="C5")
LABELS = dict(run_cdm="CDM", run_sidm=r"SIDM $\sigma$=1",
              run_sidm5=r"SIDM $\sigma$=5", run_sidm10=r"SIDM $\sigma$=10",
              run_isidm="iSIDM (f=0.3)", run_isidm_endo="iSIDM (f=0)",
              run_dsidm="dSIDM")
SEL_SNAP = "output_00010"   # selection epoch z~5.3


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


def collect_history(run):
    rdir = os.path.join(BASE, run)
    outs = sorted(glob.glob(os.path.join(rdir, "output_*")))
    snap_z = []
    sink_at = []   # list of dict {id: (M_BH_Msun, Mdot_Msun_yr)}
    for out in outs:
        nstr = out.split("output_")[-1]
        info_path = os.path.join(out, f"info_{nstr}.txt")
        csv_path  = os.path.join(out, f"sink_{nstr}.csv")
        if not (os.path.exists(info_path) and os.path.exists(csv_path)):
            continue
        try:
            info = parse_info(info_path)
            csv = np.loadtxt(csv_path, delimiter=",")
        except Exception:
            continue
        if csv.size == 0:
            continue
        if csv.ndim == 1: csv = csv[None, :]
        unit_m = info["unit_d"] * info["unit_l"]**3
        unit_mdot = unit_m / info["unit_t"]
        ids   = csv[:, 0].astype(int)
        M     = csv[:, 1] * unit_m / MSUN
        Mdot  = csv[:, 9] * unit_mdot / (MSUN / 3.156e7)  # Msun/yr
        z = 1.0/info["aexp"] - 1.0
        snap_z.append(z)
        sink_at.append({int(i): (M[k], Mdot[k]) for k, i in enumerate(ids)})
    return np.array(snap_z), sink_at, outs


def champion(snap_z, sink_at):
    return np.array([max(d.values(), key=lambda x: x[0])[0] for d in sink_at])


def topk_mean(snap_z, sink_at, k=100):
    out = []
    for d in sink_at:
        masses = np.array([v[0] for v in d.values()])
        masses.sort()
        sel = masses[-k:] if len(masses) >= k else masses
        out.append(sel.mean())
    return np.array(out)


def trace_back(sink_at, snap_z, target_ids):
    """Return dict {id: (z_arr, M_arr)}  -- earlier outputs only included if id present."""
    tracks = {tid: ([], []) for tid in target_ids}
    for z, d in zip(snap_z, sink_at):
        for tid in target_ids:
            if tid in d:
                tracks[tid][0].append(z)
                tracks[tid][1].append(d[tid][0])
    return {tid: (np.array(zs), np.array(ms)) for tid, (zs, ms) in tracks.items()}


def main():
    histories = {}
    for run in RUNS:
        print(f"  reading {run} ...")
        snap_z, sink_at, outs = collect_history(run)
        histories[run] = dict(z=snap_z, sink_at=sink_at,
                              champ=champion(snap_z, sink_at),
                              top100=topk_mean(snap_z, sink_at, 100))
        # top-5 selection at z~5.3 (output_00010 -> index where filename matches)
        sel_idx = None
        for i, o in enumerate(outs):
            if o.endswith(SEL_SNAP):
                sel_idx = i; break
        if sel_idx is None:
            histories[run]["top5"] = None
            continue
        d_sel = sink_at[sel_idx]
        items = sorted(d_sel.items(), key=lambda kv: -kv[1][0])
        top5_ids = [k for k, _ in items[:5]]
        histories[run]["top5"] = trace_back(sink_at, snap_z, top5_ids)
        histories[run]["sel_z"] = snap_z[sel_idx]
        print(f"    z range: {snap_z.min():.2f} - {snap_z.max():.2f}, "
              f"champ at z={snap_z[sel_idx]:.2f}: {histories[run]['champ'][sel_idx]:.2e} Msun")

    # ============ plot ============
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    ax = axes[0]
    for run, h in histories.items():
        ax.plot(h["z"], h["champ"], "o-", color=COLORS[run],
                label=LABELS[run], lw=2)
    ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$M_{\rm BH}^{\rm max}$ [M$_\odot$]")
    ax.set_title("Champion BH (most massive at each z)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    ax = axes[1]
    for run, h in histories.items():
        ax.plot(h["z"], h["top100"], "s-", color=COLORS[run],
                label=LABELS[run], lw=2)
    ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$\langle M_{\rm BH} \rangle_{\rm top\,100}$ [M$_\odot$]")
    ax.set_title("Mean of top-100 BHs")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    ax = axes[2]
    for run, h in histories.items():
        if h.get("top5") is None: continue
        for j, (tid, (zs, ms)) in enumerate(h["top5"].items()):
            lab = LABELS[run] if j == 0 else None
            ax.plot(zs, ms, "-", color=COLORS[run], alpha=0.7, lw=1.3, label=lab)
    ax.axvline(5.3, color="grey", ls=":", alpha=0.6)
    ax.text(5.3, 5e4, "selection\nz=5.3", color="grey", fontsize=8, ha="center")
    ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("z"); ax.set_ylabel(r"$M_{\rm BH}$ [M$_\odot$]")
    ax.set_title("Top-5 BHs at z=5.3 traced back")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower left")

    plt.tight_layout()
    out_png = os.path.join(BASE, "analysis_z5", "bh_growth_tracks.png")
    plt.savefig(out_png, dpi=140)
    print(f"  → {out_png}")

    # ============ side info: total BH mass density ===========
    out_txt = os.path.join(BASE, "analysis_z5", "bh_growth_summary.txt")
    with open(out_txt, "w") as f:
        f.write("# z   ")
        for run in RUNS:
            f.write(f"  champ[{LABELS[run]:>16s}]  top100mean[{LABELS[run]:>16s}]")
        f.write("\n")
        # use cdm z as reference
        zref = histories["run_cdm"]["z"]
        for i, z in enumerate(zref):
            f.write(f"{z:6.3f}")
            for run in RUNS:
                h = histories[run]
                # match z by index if same length, else nearest
                if len(h["z"]) > i and abs(h["z"][i]-z) < 0.05:
                    j = i
                else:
                    j = int(np.argmin(np.abs(h["z"]-z)))
                f.write(f"  {h['champ'][j]:18.4e}  {h['top100'][j]:22.4e}")
            f.write("\n")
    print(f"  → {out_txt}")


if __name__ == "__main__":
    main()
