"""Descriptive companion (no model): ascent at HF onset vs the storm-force-only peak (anchors from ../hf_vs_storm/results/storms.csv),
matched on basin x month (stratified difference weighted by HF counts), season-block bootstrap, BH over the 4 measures per basin.
Only anchors that fall on an extracted fix (00/12 UTC, Oct-Apr, seasons 2004-05..2021-22) count. usage: python3 desc_w.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "hf_jet_trough"))
import numpy as np, pandas as pd
from qc import bh
S = pd.read_csv(os.path.join(HERE, "..", "hf_vs_storm", "results", "storms.csv"), dtype={"onset_time": str, "peak_time": str})
X = pd.read_csv(os.path.join(HERE, "results", "features_w.csv.gz"), dtype={"time": str})
hf = S[S.grp == "HF"].dropna(subset=["onset_time"]).assign(time=lambda d: d.onset_time.str[:10]).merge(X, on=["track", "time"], suffixes=("", "_f"))
sf = S[S.grp == "SF"].assign(time=lambda d: d.peak_time.astype(str).str[:10]).merge(X, on=["track", "time"], suffixes=("", "_f"))
meas = {"a500_1000": "500 hPa ascent max within 1000 km (Pa/s)", "a500_x750": "ascent max 750-2500 km (Pa/s)",
        "aarea500": "ascent area >= 0.3 Pa/s within 3000 km (1e6 km2)", "couple_w": "share with a second ascent centre (primary rule)"}
rng = np.random.default_rng(20261011)
out, rows, pv = [], [], []
for b in ("atl", "pac"):
    h, s = hf[hf.basin == b], sf[sf.basin == b]
    seasons = np.sort(np.union1d(h.season.unique(), s.season.unique()))
    out.append(f"{b}: HF onset anchors on an extracted fix {len(h)}; storm-force peak anchors {len(s)}")
    for m, label in meas.items():
        def stat(hs, ss):
            tot = 0.0; n = 0
            for mo, g in hs.groupby("mon"):
                sg = ss[ss.mon == mo]
                if len(sg):
                    tot += g[m].mean() - sg[m].mean() if False else (g[m].mean() - sg[m].mean()) * len(g); n += len(g)
            return tot / n if n else np.nan
        d0 = stat(h, s)
        bs = []
        for _ in range(1000):
            pick = rng.choice(seasons, len(seasons))
            bs.append(stat(pd.concat([h[h.season == u] for u in pick]), pd.concat([s[s.season == u] for u in pick])))
        bs = np.array(bs); bs = bs[np.isfinite(bs)]
        p = float(2 * min((bs <= 0).mean(), (bs >= 0).mean())) or 1 / 500
        rows.append((b, label, d0, np.percentile(bs, 5), np.percentile(bs, 95), h[m].mean(), s[m].mean())); pv.append(max(p, 1 / 500))
q = bh(pv)
for (b, label, d0, lo, hi, mh, ms), qq in zip(rows, q):
    out.append(f"  {b} {label}: HF onset {mh:.3f}, storm-force peak {ms:.3f}, stratified difference {d0:+.3f} [90% {lo:+.3f}, {hi:+.3f}], q {qq:.3f}")
out.append("Descriptive only (HF is defined by the gust index, so HF storms differ from storm-force ones in strength by construction); not a test.")
open(os.path.join(HERE, "results", "descriptive.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
