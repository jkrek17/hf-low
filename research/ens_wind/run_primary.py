"""Registered analysis from count tables (see PREREGISTRATION.md). usage: run_primary.py COUNTS.npz OUT_PREFIX"""
import sys, json, numpy as np
import stats as S

def load(path):
    z = np.load(path); K = z["keys"]; n = z["n"]; y = z["y"]
    return K, n, y

def season_of(ymd):
    y, m = int(ymd[:4]), int(ymd[4:6]); return y if m >= 7 else y - 1

def cell(K, n, y, truth, basin, thr, lead):
    """-> N, {season: n[k]}, {season: y[k]}"""
    sel = [i for i, k in enumerate(K) if k[1] == truth and k[2] == basin and int(k[3]) == thr and int(k[5]) == lead]
    N = int(K[sel[0]][6]); nb, yb = {}, {}
    for i in sel:
        s = season_of(K[i][4]); nb[s] = nb.get(s, 0) + n[i][:N + 1]; yb[s] = yb.get(s, 0) + y[i][:N + 1]
    return N, nb, yb

if __name__ == "__main__":
    K, n, y = load(sys.argv[1]); out = sys.argv[2]
    res = {}; rows = []
    for basin in ("natl", "npac"):
        for thr, truth in ((48, "ascat"), (64, "opc")):
            for lead in (24, 48, 72, 96):
                N, nb, yb = cell(K, n, y, truth, basin, thr, lead)
                r = S.summarize(nb, yb, N); r.update(basin=basin, thr=thr, truth=truth, lead=lead); res[f"{basin}_{thr}_{truth}_{lead}"] = r
                rows.append(r)
    ok = [r for r in rows if r["status"] == "ok"]
    rej = S.bh([r["p_overconf"] for r in ok]) if ok else []
    for r, j in zip(ok, rej): r["bh_reject"] = bool(j)
    # secondary: 48 kt vs OPC
    sec = {}
    for basin in ("natl", "npac"):
        for lead in (24, 48, 72, 96):
            N, nb, yb = cell(K, n, y, "opc", basin, 48, lead)
            r = S.summarize(nb, yb, N); r.update(basin=basin, thr=48, truth="opc", lead=lead); sec[f"{basin}_48_opc_{lead}"] = r
    json.dump(dict(primary=res, secondary_48_opc=sec), open(out + ".json", "w"), default=float, indent=1)
    def line(r):
        if r["status"] != "ok": return f'{r["basin"]} {r["thr"]}kt {r["truth"]} {r["lead"]:>2}h: {r["status"]} (cases {r["cases"]})'
        rate = [f'{b[1]/b[0]*100:.1f}%' if b[0] else "-" for b in zip(r["bins_n"], r["bins_y"])]
        return (f'{r["basin"]} {r["thr"]}kt {r["truth"]} {r["lead"]:>2}h: cases {r["cases"]:,} events {r["events"]:,} | rate by bin {rate} | slope {r["slope"]:.2f} '
                f'({r["slope_ci"][0]:.2f} to {r["slope_ci"][1]:.2f}) p={r["p_overconf"]:.3f} BH={r.get("bh_reject")} | BSS raw {r["bss_raw"]:+.3f} cal {r["bss_cal"]:+.3f} gain {r["cal_gain"]:+.3f} ({r["cal_gain_ci"][0]:+.3f} to {r["cal_gain_ci"][1]:+.3f})')
    with open(out + ".txt", "w") as f:
        f.write("PRIMARY (48 kt vs ASCAT, 64 kt vs OPC HF neighbourhood)\n"); [f.write(line(r) + "\n") for r in rows]
        f.write("\nSECONDARY (48 kt vs OPC storm-or-stronger neighbourhood, no test)\n"); [f.write(line(r) + "\n") for r in sec.values()]
    print(open(out + ".txt").read())
