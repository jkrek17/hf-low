"""Outcome-free test of collocate.py + analyse.py: a synthetic scatterometer built from the ERA5 field itself.
wind = gust / GF (+ noise), saturating near SAT kt, with random swath gaps and times jittered within the window.
Runs on real ERA5 storm fixes (HF-index storms plus a random sub-HF sample). Nothing here is evidence about ASCAT.
usage: test_synthetic.py [n_sub] [gust_factor]"""
import os, sys, numpy as np, pandas as pd
import era5side as E, collocate as C, analyse as A

n_sub = int(sys.argv[1]) if len(sys.argv) > 1 else 700
GF = float(sys.argv[2]) if len(sys.argv) > 2 else 1.3
SAT = 50.0
SYNTH_OUT = os.environ.get("SYNTH_OUT", "synth_footprints.csv")      # scratch output, not committed
rng = np.random.default_rng(1)

S = C.storm_sample()
idx = E.cache_index()
S = S[S.h1900.isin(idx)]
S = pd.concat([S[S.gust800_kt >= 71.7], S[S.gust800_kt < 71.7].sample(n_sub, random_state=0)])
S["fr"] = [idx[h] for h in S.h1900]
S = S.sort_values("fr")
ocean = E.lsm_ocean()
rows = []
for _, r in S.iterrows():
    g = E.field(*r.fr)
    th = rng.uniform(0, np.pi)
    off = rng.uniform(-600, 600)
    R, Cc = E.disc_cells(r.peak_lat, r.peak_lon, 900)
    la, lo = E.LAT[R] - 0.125, E.LON[Cc] + 0.125     # scatterometer-grid cell centres (x.125)
    x = (lo - r.peak_lon + 180) % 360 - 180
    xk = x * 111.2 * np.cos(np.radians(r.peak_lat)); yk = (la - r.peak_lat) * 111.2
    s = xk * np.cos(th) + yk * np.sin(th) - off
    m = (np.abs(s) < 500) & ocean[R, Cc]
    if m.sum() < 20:
        continue
    gv = g[R[m], Cc[m]]
    w = gv / GF + rng.normal(0, 3, m.sum())
    w = np.minimum(w, SAT + 0.1 * np.maximum(w - SAT, 0))
    t = pd.Timestamp("1900-01-01") + pd.to_timedelta(r.h1900 + rng.uniform(-1.2, 1.2), unit="h")
    rows.append(pd.DataFrame(dict(sensor="synth", time=t.isoformat(), lat=la[m], lon=lo[m], wind_kt=np.clip(w, 0, None), ok=1)))
cells = pd.concat(rows)
F, P = C.collocate(cells, S=S, ocean=ocean, idx=idx)
F.to_csv(SYNTH_OUT, index=False)
print("footprints", len(F), "accepted", int((F.cover >= C.MIN_COVER).sum()), "mean cover", round(F.cover.mean(), 2))
print("cell pairs", len(P), "pointwise corr gust vs wind", round(np.corrcoef(P.g, P.wind_kt)[0, 1], 3))
A.NB = 300
R, _ = A.main(SYNTH_OUT)
q = R.test.str.startswith("P3")
se = (R.loc[q, "hi95"].iloc[0] - R.loc[q, "lo95"].iloc[0]) / 3.92
print("P3 se (kt)", round(se, 2), " MDE at 80% power, alpha .05 ~", round(2.8 * se, 2), "kt")
