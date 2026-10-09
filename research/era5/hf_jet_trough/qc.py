"""Outcome-free detection QC Q1, Q2 (see PREREGISTRATION.md). Loads track, time, basin, season, lat, lon, heading only.

usage: ERA5_WORK=DIR qc.py        writes results/qc.txt, results/features.csv.gz (features + meta, no outcome)
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from features import J2, T2, COLS
WORK = os.environ.get("ERA5_WORK", "")
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")


def load():
    X = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in sorted(glob.glob(f"{WORK}/feat/*.csv"))], ignore_index=True)
    M = pd.read_csv(FIX, usecols=["track", "time", "basin", "season", "lat", "lon", "heading"], dtype={"time": str})
    return M.merge(X, on=["track", "time"], how="inner", validate="one_to_one")


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty(len(p))
    r = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


if __name__ == "__main__":
    D = load()
    D.to_csv(os.path.join(HERE, "results", "features.csv.gz"), index=False, float_format="%.5g")
    mon = pd.to_datetime(D.time, format="%Y%m%d%H").dt.month
    S = D[mon.isin([10, 11, 12, 1, 2, 3, 4])]
    out = [f"fixes with features {len(D)}; Oct-Apr model sample {len(S)}; tracks {S.track.nunique()}; seasons {S.season.nunique()}"]
    ok = True
    for b in ("atl", "pac", "both"):
        s = S if b == "both" else S[S.basin == b]
        j = s[s.vmax2500 >= 60]
        q1j = (j.nojet == 0).mean()
        q1t_all = ((s.notrough == 0) & (s.nohead == 0)).mean()
        q1t_head = (s[s.nohead == 0].notrough == 0).mean()
        out.append(f"Q1 {b}: streak found {q1j:.3f} of {len(j)} fixes with vortex-removed 250 hPa max >= 60 kt (need 0.80); "
                   f"trough found {q1t_all:.3f} of all {len(s)} fixes (need 0.70), {q1t_head:.3f} of the {int((s.nohead == 0).sum())} fixes with a heading; "
                   f"no heading {(s.nohead == 1).mean():.3f}")
        if b != "both":
            ok &= (q1j >= 0.80) and (q1t_all >= 0.70)
    # Q2
    rng = np.random.default_rng(20261009)
    rows, p = [], []
    for b in ("atl", "pac"):
        s = S[(S.basin == b) & (S.nojet == 0)]
        s = s[s.RE + s.LX + s.LE + s.RX > 0]
        A = (s.RE + s.LX) > 0
        us = np.sort(s.season.unique())
        agg = {u: (s.div300mean1000[(s.season == u) & A].sum(), ((s.season == u) & A).sum(),
                   s.div300mean1000[(s.season == u) & ~A].sum(), ((s.season == u) & ~A).sum()) for u in us}
        arr = np.array([agg[u] for u in us], float)

        def stat(a):
            return a[:, 0].sum() / a[:, 1].sum() - a[:, 2].sum() / a[:, 3].sum()
        d0 = stat(arr)
        bs = np.array([stat(arr[rng.integers(0, len(us), len(us))]) for _ in range(2000)])
        pv = float((bs <= 0).mean())
        rows.append((b, d0, np.percentile(bs, 5), np.percentile(bs, 95), int(arr[:, 1].sum()), int(arr[:, 3].sum())))
        p.append(max(pv, 1 / 2001))
    q = bh(p)
    for (b, d0, lo, hi, nA, nB), qq in zip(rows, q):
        out.append(f"Q2 {b}: mean 300 hPa divergence within 1000 km, (RE+LX) minus (LE+RX) = {d0:+.3f} e-5/s "
                   f"[90% {lo:+.3f}, {hi:+.3f}], n {nA} vs {nB}, one-sided q {qq:.4f}")
        ok &= (d0 > 0) and (qq < 0.05)
    quad = S[S.nojet == 0]
    out.append("quadrant shares of fixes with a streak (RE LE RX LX none): " + " ".join(
        f"{b}:" + "/".join(f"{x:.3f}" for x in [quad[quad.basin == b][c].mean() for c in ("RE", "LE", "RX", "LX")] +
                           [1 - quad[quad.basin == b][["RE", "LE", "RX", "LX"]].sum(1).mean()]) for b in ("atl", "pac")))
    out.append("Q1 and Q2 " + ("PASS: features may be used for the skill test (Q3 is a by-eye review)" if ok else "FAIL: features are not used for the skill test"))
    open(os.path.join(HERE, "results", "qc.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))
