"""Secondary checks S1-S7 of PREREGISTRATION.md (own BH family; none changes a primary result). ERA5 PROXY, pipeline A.

usage: secondary.py OUTDIR [NPERM] [NBOOT]
"""
import json
import os
import sys

import numpy as np
import pandas as pd

import chanlib as C
import channels as CH


def slope_test(D, name, idx, nperm, nboot, rng):
    """Single-outcome slope: estimate, clustered SE, season-block permutation p, bootstrap interval."""
    rows = np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in idx])
    cl = np.repeat(np.arange(len(idx)), C.NWEEK)
    b = D.slopes(rows, D.x, [name])[name]
    ns = len(idx)
    xs = D.x[rows].reshape(ns, C.NWEEK)
    null = np.empty(nperm)
    for i in range(nperm):
        xf = D.x.copy(); xf[rows] = xs[rng.permutation(ns)].reshape(-1)
        null[i] = D.slopes(rows, xf, [name])[name]
    bt = np.empty(nboot)
    for i in range(nboot):
        pick = rng.integers(0, ns, ns)
        r = np.concatenate([np.arange(idx[j] * C.NWEEK, (idx[j] + 1) * C.NWEEK) for j in pick])
        bt[i] = D.slopes(r, D.x, [name])[name]
    p = (1 + int((np.abs(null) >= abs(b) - 1e-12).sum())) / (1 + nperm)
    se = C.joint_se([(D.X(name, rows), D.Y[name][rows])], cl)
    return dict(est=b, se=se, p=p, lo=np.percentile(bt, 2.5), hi=np.percentile(bt, 97.5))


def rows_for(tag, basin, res, extra=None):
    out = []
    for k in ("T1", "T2", "T3"):
        r = res[k]
        out.append(dict(check=tag, basin=basin, test=k, est=np.exp(r["est"]), lo=np.exp(r["lo"]), hi=np.exp(r["hi"]), p=r["p"]))
    f = res["f"]
    out.append(dict(check=tag, basin=basin, test="f", est=f["est"], lo=f["lo"], hi=f["hi"], p=np.nan))
    return out


def main():
    outdir = sys.argv[1]
    nperm = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    nboot = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(C.SEED + 1)
    Tfull = CH.prev_all_tracks()
    M = C.month_dummies()
    allidx = list(range(len(C.SEASONS)))
    half1, half2 = list(range(0, 11)), list(range(11, 22))
    rows = []
    for b in C.BASINS:
        x, y_arch = C.index(b)
        Y, LP, win = CH.build(b, "hf", None, Tfull)
        D = CH.Design(M, x, Y, LP)
        # S1 fixed-depth HF
        Yd, LPd, _ = CH.build(b, "hfd", None, Tfull)
        res, _ = CH.analyse(CH.Design(M, x, Yd, LPd), allidx, CH.OUTS_CORE, nperm, nboot, rng)
        rows += rows_for("S1 fixed-depth HF", b, res)
        # S2a local split, descriptive RR(N)
        ent = win.entrant == 1
        gen = win[(~ent) & (win.mature == 0)]; mat = win[(~ent) & (win.mature == 1)]
        Yg = dict(N=C.weekly_counts(win, (~ent) & (win.mature == 0)), M=C.weekly_counts(win, (~ent) & (win.mature == 1)))
        for nm, lab in (("N", "S2a local, observed genesis (>=1000 hPa)"), ("M", "S2a local, mature at detection (<1000 hPa)")):
            t_all = Tfull[Tfull.basin == b]
            m = (t_all.entrant == 0) & ((t_all.mature == 0) if nm == "N" else (t_all.mature == 1))
            lp = np.log1p(C.prev_week_counts(t_all.assign(_m=m), m, b))
            Dg = CH.Design(M, x, {nm: Yg[nm]}, {nm: lp})
            r = slope_test(Dg, nm, allidx, nperm, nboot, rng)
            rows.append(dict(check=lab, basin=b, test="T1", est=np.exp(r["est"]), lo=np.exp(r["lo"]), hi=np.exp(r["hi"]), p=r["p"]))
        # S2b no previous-week term
        res, _ = CH.analyse(CH.Design(M, x, Y, LP, use_prev=False), allidx, CH.OUTS_CORE, nperm, nboot, rng)
        rows += rows_for("S2b no previous-week term", b, res)
        # S3 leave out peak latitude > 60N
        Y3, LP3, _ = CH.build(b, "hf", lambda t: t.peak_lat <= 60, Tfull)
        res, _ = CH.analyse(CH.Design(M, x, Y3, LP3), allidx, CH.OUTS_CORE, nperm, nboot, rng)
        rows += rows_for("S3 peak latitude <= 60N", b, res)
        # S4 halves
        for lab, ids in (("S4a seasons 2004-14", half1), ("S4b seasons 2015-25", half2)):
            res, _ = CH.analyse(D, ids, CH.OUTS_CORE, nperm, nboot, rng, loso=False)
            rows += rows_for(lab, b, res)
        # S5 archive reference (total only)
        w = pd.read_csv(C.WEEKLY)
        w = w[w.season >= 2004].reset_index(drop=True)
        Ya = dict(A=w[f"y_{b}"].values.astype(float))
        Da = CH.Design(M, x, Ya, dict(A=np.log1p(w[f"prev_{b}"].values)))
        r = slope_test(Da, "A", allidx, nperm, nboot, rng)
        rows.append(dict(check="S5 archive weekly HF lows (total only)", basin=b, test="T3", est=np.exp(r["est"]), lo=np.exp(r["lo"]), hi=np.exp(r["hi"]), p=r["p"]))
        # S6 season trend
        E = ((np.repeat(C.SEASONS, C.NWEEK) - 2015) / 10.0)[:, None]
        res, _ = CH.analyse(CH.Design(M, x, Y, LP, E=E), allidx, CH.OUTS_CORE, nperm, nboot, rng)
        rows += rows_for("S6 add season trend", b, res)
        print(f"{b} S1-S6 done", flush=True)
    # S7: 1979-80 to 2000-01, frozen primary-split pattern, depth HF, within-era only
    fz = json.load(open(os.path.join(C.ERA, "hemispheric", "results", "frozen_primary.json")))
    w = pd.read_csv(C.WEEKLY)
    seas7 = list(range(1979, 2001))
    C.SEASONS = seas7
    import importlib
    T7 = C.load_tracks(first=1979, last=2000)
    T7["entrant"] = 0; T7["mature"] = 0; T7["missing_fix"] = 0
    M7 = C.month_dummies()
    idx7 = list(range(len(seas7)))
    pcs_cols = [f"{v}{i}" for v in "zus" for i in range(1, 11)]
    s7 = []
    for b in C.BASINS:
        beta = np.array(fz["basins"][b]["models"]["P"]["beta"])
        ww = w[(w.season >= 1979) & (w.season <= 2000)].sort_values(["season", "week"]).reset_index(drop=True)
        z = ww[pcs_cols].values @ beta[8:]
        x7 = (z - z.mean()) / z.std()
        Y, LP, win = CH.build(b, "hfd", None, T7)
        D7 = CH.Design(M7, x7, Y, LP)
        # reproduction check of PR 41 S5 on the depth total (month + previous week + index)
        rep = slope_test(D7, "H", idx7, nperm, 200, rng)
        s7.append(dict(basin=b, repro_rr=float(np.exp(rep["est"])), target=1.245 if b == "atl" else 1.158))
        res, _ = CH.analyse(D7, idx7, CH.OUTS_CORE, nperm, nboot, rng)
        rows += rows_for("S7 1979-2000 depth HF, frozen pattern", b, res)
        print(f"{b} S7 done; reproduction RR/SD {np.exp(rep['est']):.3f} (PR 41 S5 {s7[-1]['target']})", flush=True)
    S = pd.DataFrame(rows)
    pm = S.p.notna()
    S.loc[pm, "q_secondary"] = C.bh(S.p[pm].values)
    S.to_csv(os.path.join(outdir, "secondary.csv"), index=False)
    json.dump(s7, open(os.path.join(outdir, "s7_reproduction.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
