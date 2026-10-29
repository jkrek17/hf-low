"""Planted-effect test, no network: synthetic cases where GEFS and ENS each see half of the signal, so the blend must beat GEFS;
and a case where ENS is pure noise, so it must not. Run: python3 -I test_blend.py"""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import blend as Bl

def make(rng, ens_informative, nd=37, ncell=3000):
    dates = [f"2025{m:02d}{d:02d}" for m, d in [(10 + (i * 5 + 14) // 30 % 3 if False else 10, 15)] * 0] ; dates = []
    import datetime as dt
    d0 = dt.date(2025, 10, 15)
    dates = [(d0 + dt.timedelta(days=5 * i)).strftime("%Y%m%d") for i in range(nd)]
    n = np.zeros((nd, 32, 52), np.int64); y = np.zeros_like(n)
    for i in range(nd):
        z = rng.normal(size=ncell); s1 = z + rng.normal(size=ncell); s2 = (z + rng.normal(size=ncell)) if ens_informative else rng.normal(size=ncell)
        ev = rng.random(ncell) < 1 / (1 + np.exp(-(2.2 * z - 4.5)))
        kg = np.clip(np.round(31 / (1 + np.exp(-(s1 - 2.0)))), 0, 31).astype(int); ke = np.clip(np.round(51 / (1 + np.exp(-(s2 - 2.0)))), 0, 51).astype(int)
        np.add.at(n[i], (kg, ke), 1); np.add.at(y[i], (kg, ke), ev.astype(int))
    return dates, n, y

rng = np.random.default_rng(1)
for inf in (True, False):
    dates, n, y = make(rng, inf); k = np.arange(32); pG = np.clip(np.cumsum(np.r_[0, np.full(31, 0.03)]), 0, 1)
    # calibrate pG from the pooled marginals (stand-in for the committed table)
    from stats import iso_fit
    pG = iso_fit(n.sum((0, 2)), y.sum((0, 2)))
    weeks, res, clim, ws, bins, ev, cs = Bl.run_cell(dates, n, y, pG)
    den = sum(clim.values()); g = 1 - sum(res["gefs"].values()) / den; b = 1 - sum(res["blend"].values()) / den
    print(f"ENS informative={inf}: BSS gefs {g:+.3f} blend {b:+.3f} diff {b-g:+.3f}")
    assert (b - g > 0.01) if inf else (b - g < 0.01)
print("ok")
