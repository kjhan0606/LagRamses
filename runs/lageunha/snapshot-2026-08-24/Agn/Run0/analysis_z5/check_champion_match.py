#!/usr/bin/env python3
"""Verify whether the 'champion BH' at z=5.3 is the same halo across simulations.

For each run, at output_00010 (z~5.3):
  1. Find champion BH (most massive)
  2. Print its sink_id, mass, comoving position
  3. Cross-check: take CDM champion's position as anchor, find nearest BH
     in every other sim — is it also their champion?
"""
import os, numpy as np
BASE = "/gpfs/kjhan/Hydro/Sidm/Agn/Run0"
RUNS = ["run_cdm","run_sidm","run_sidm5","run_sidm10",
        "run_isidm","run_isidm_endo","run_dsidm"]
SNAP = "output_00010"
MSUN = 1.989e33
MPC = 3.086e24


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


def load(run):
    info = parse_info(os.path.join(BASE, run, SNAP, "info_00010.txt"))
    csv  = np.loadtxt(os.path.join(BASE, run, SNAP, "sink_00010.csv"),
                      delimiter=",")
    if csv.ndim == 1: csv = csv[None,:]
    unit_m = info["unit_d"] * info["unit_l"]**3
    h = info["H0"]/100.0
    Lbox_Mpch = info["unit_l"]/info["aexp"]/MPC * h    # comoving Mpc/h
    M = csv[:,1] * unit_m / MSUN
    pos = csv[:, 2:5]                                  # code units (0..1)
    pos_Mpch = pos * Lbox_Mpch                         # comoving Mpc/h
    sid = csv[:,0].astype(int)
    return dict(sid=sid, M=M, pos=pos_Mpch, Lbox=Lbox_Mpch)


def per_dist(p1, p2, Lbox):
    d = p1 - p2
    d -= Lbox * np.round(d/Lbox)
    return np.sqrt((d**2).sum(axis=-1))


def main():
    data = {r: load(r) for r in RUNS}

    # champion of each sim
    print("=== Champion BH at z=5.3 (output_00010) ===")
    champs = {}
    for r, d in data.items():
        i = int(np.argmax(d["M"]))
        champs[r] = dict(idx=i, sid=d["sid"][i], M=d["M"][i], pos=d["pos"][i])
        print(f"  {r:18s}  sid={d['sid'][i]:6d}  M={d['M'][i]:8.3e} Msun  "
              f"pos=({d['pos'][i][0]:6.2f},{d['pos'][i][1]:6.2f},{d['pos'][i][2]:6.2f}) Mpc/h")

    # distance to CDM champion
    cdm = champs["run_cdm"]
    Lbox = data["run_cdm"]["Lbox"]
    print(f"\n=== Champion BH separations (relative to CDM champion sid={cdm['sid']}) ===")
    print(f"  CDM champ pos = ({cdm['pos'][0]:6.2f},{cdm['pos'][1]:6.2f},{cdm['pos'][2]:6.2f}) Mpc/h\n")
    for r, c in champs.items():
        dist = per_dist(c["pos"], cdm["pos"], Lbox)
        same = "SAME ID" if c["sid"] == cdm["sid"] else "       "
        print(f"  {r:18s}  Δ={float(dist):6.3f} Mpc/h   sid={c['sid']:6d}  {same}")

    # Find each sim's BH nearest to CDM champion position
    print(f"\n=== Nearest BH to CDM champion's position in each sim ===")
    for r, d in data.items():
        dists = per_dist(d["pos"], cdm["pos"], Lbox)
        j = int(np.argmin(dists))
        is_champ = "YES" if j == champs[r]["idx"] else "NO  "
        # mass rank of this sink
        rank = int(np.sum(d["M"] > d["M"][j])) + 1
        print(f"  {r:18s}  nearest sid={d['sid'][j]:6d}  Δ={float(dists[j]):.4f} Mpc/h  "
              f"M={d['M'][j]:8.3e}  rank={rank:3d}/{len(d['M']):5d}  is_champion={is_champ}")

    # Take top-5 in each sim and report their pairwise overlap with CDM top-5
    print(f"\n=== Top-5 BHs: position overlap with CDM top-5 ===")
    cdm_top5 = np.argsort(-data["run_cdm"]["M"])[:5]
    cdm_top5_pos = data["run_cdm"]["pos"][cdm_top5]
    cdm_top5_sid = data["run_cdm"]["sid"][cdm_top5]
    print(f"  CDM top-5 sid: {cdm_top5_sid.tolist()}")
    print(f"  CDM top-5 positions:")
    for i, (s, p) in enumerate(zip(cdm_top5_sid, cdm_top5_pos)):
        print(f"    #{i+1} sid={s:6d}  ({p[0]:6.2f},{p[1]:6.2f},{p[2]:6.2f}) Mpc/h")
    print()
    for r, d in data.items():
        top5 = np.argsort(-d["M"])[:5]
        sid5 = d["sid"][top5]
        pos5 = d["pos"][top5]
        # for each, find nearest CDM top-5
        report = []
        for k in range(5):
            dists = per_dist(pos5[k], cdm_top5_pos, Lbox)
            jmin = int(np.argmin(dists))
            report.append((sid5[k], cdm_top5_sid[jmin], float(dists[jmin])))
        n_close = sum(1 for _,_,dd in report if dd < 1.0)
        print(f"  {r:18s}  matches within 1 Mpc/h to CDM top-5: {n_close}/5")
        for k,(my, cdms, dd) in enumerate(report):
            tag = "★" if dd < 0.5 else ("·" if dd < 2 else "✗")
            print(f"      #{k+1} my_sid={my:6d}  → CDM_sid={cdms:6d}  Δ={dd:6.3f} Mpc/h  {tag}")


if __name__ == "__main__":
    main()
