"""Outcome-free detection QC Q1, Q2 for the omega features (see PREREGISTRATION.md).

usage: ERA5_WORK=DIR qc_w.py        writes results/qc.txt, results/features_w.csv.gz (features + meta, no outcome)
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "hf_jet_trough"))
import numpy as np, pandas as pd
from qc import bh
WORK = os.environ.get("ERA5_WORK", "")
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
JT = os.path.join(HERE, "..", "hf_jet_trough", "results", "features.csv.gz")

if __name__ == "__main__":
    X = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in sorted(glob.glob(f"{WORK}/featw/*.csv"))], ignore_index=True)
    M = pd.read_csv(FIX, usecols=["track", "time", "basin", "season", "lat", "lon", "heading"], dtype={"time": str})
    D = M.merge(X, on=["track", "time"], how="inner", validate="one_to_one")
    D.to_csv(os.path.join(HERE, "results", "features_w.csv.gz"), index=False, float_format="%.5g")
    mon = pd.to_datetime(D.time, format="%Y%m%d%H").dt.month
    S = D[mon.isin([10, 11, 12, 1, 2, 3, 4])]
    out = [f"fixes with features {len(D)}; Oct-Apr model sample {len(S)}; tracks {S.track.nunique()}; seasons {S.season.nunique()}"]
    ok = True
    for b in ("atl", "pac", "both"):
        s = S if b == "both" else S[S.basin == b]
        q1a = (np.fmax(s.a500_1000, s.a500_x750) >= 0.2).mean()
        sh = {v: (s["nasc" + v] >= 2).mean() for v in ("", "_m", "_o")}
        out.append(f"Q1 {b}: ascent maximum >= 0.2 Pa/s within 2500 km at {q1a:.3f} of {len(s)} fixes (need 0.80); fixes with >= 2 centres: "
                   f"primary {sh['']:.3f} (need 0.10 to 0.70), 0.6 Pa/s {sh['_m']:.3f}, first-written rule {sh['_o']:.3f}")
        if b != "both":
            ok &= (q1a >= 0.80) and (0.10 <= sh[""] <= 0.70)
    # Q2: lift higher in RE and LX than in LE and RX (PR 149 quadrants; no outcome)
    Q = pd.read_csv(JT, dtype={"time": str}, usecols=["track", "time", "nojet", "RE", "LE", "RX", "LX"])
    S = S.merge(Q, on=["track", "time"], how="left")
    rng = np.random.default_rng(20261011)
    rows, p = [], []
    for b in ("atl", "pac"):
        s = S[(S.basin == b) & (S.nojet == 0)]
        s = s[s.RE + s.LX + s.LE + s.RX > 0]
        A = (s.RE + s.LX) > 0
        us = np.sort(s.season.unique())
        arr = np.array([(s.a500_1000[(s.season == u) & A].sum(), ((s.season == u) & A).sum(),
                         s.a500_1000[(s.season == u) & ~A].sum(), ((s.season == u) & ~A).sum()) for u in us], float)
        stat = lambda a: a[:, 0].sum() / a[:, 1].sum() - a[:, 2].sum() / a[:, 3].sum()
        d0 = stat(arr)
        bs = np.array([stat(arr[rng.integers(0, len(us), len(us))]) for _ in range(2000)])
        rows.append((b, d0, np.percentile(bs, 5), np.percentile(bs, 95), int(arr[:, 1].sum()), int(arr[:, 3].sum())))
        p.append(max(float((bs <= 0).mean()), 1 / 2001))
    for (b, d0, lo, hi, nA, nB), qq in zip(rows, bh(p)):
        out.append(f"Q2 {b}: mean 500 hPa ascent within 1000 km, (RE+LX) minus (LE+RX) = {d0:+.4f} Pa/s [90% {lo:+.4f}, {hi:+.4f}], "
                   f"n {nA} vs {nB}, one-sided q {qq:.4f}")
        ok &= (d0 > 0) and (qq < 0.05)
    out.append("Q1 and Q2 " + ("PASS: features may be used for the skill test (Q3 is a by-eye review)" if ok else "FAIL: features are not used for the skill test"))
    open(os.path.join(HERE, "results", "qc.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))
