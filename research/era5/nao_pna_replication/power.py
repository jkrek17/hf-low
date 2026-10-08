"""Design power (season-block permutation of the index series; outcome marginals only). usage: power.py OUT_TXT [NPERM]"""
import sys
import numpy as np
import common as C

NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 400
T, I = C.load_tracks(), C.indices()
rng = np.random.default_rng(C.core.SEED)
OUT = []


def say(s=""):
    print(s); OUT.append(s)


for name, basin, s0, s1 in (("P1", "atl", 1979, 2000), ("S5", "atl", 2001, 2013), ("S6", "pac", 1979, 2000)):
    D, E = C.table(T, I, basin, "NAO", "PNA", s0, s1)
    seasons = sorted(D.season.unique()); ns = len(seasons)
    nd = D.day.max() + 1
    Aa = np.zeros((ns, nd)); Ba = np.zeros((ns, nd))
    for i, s in enumerate(seasons):
        d = D[D.season == s]; Aa[i, d.day.values] = d.A.values; Ba[i, d.day.values] = d.B.values
    sidx = {s: i for i, s in enumerate(seasons)}
    Es, Ed = E.season.map(sidx).values, E.day.values
    cov = C.core.covariates(E, False).values
    g = []
    for r in range(NPERM):
        perm = rng.permutation(ns)
        Ae, Be = Aa[perm[Es], Ed], Ba[perm[Es], Ed]
        b, *_ = np.linalg.lstsq(np.column_stack([Ae, Be, Ae * Be, cov]), E[["lon", "lat"]].values, rcond=None)
        g.append(b[2])
    sd = np.array(g).std(0)
    say("%s %s winters %d-%d: %d tracks, %d winters; null SE of gamma (lon, lat) = (%.2f, %.2f) deg; detectable at 80%% power = (%.2f, %.2f); well-powered null (1.96 SE <= 0.55): lon %s"
        % (name, basin, s0, s1, len(E), ns, sd[0], sd[1], 2.8 * sd[0], 2.8 * sd[1], "possible" if 1.96 * sd[0] <= 0.55 else "not possible"))
    say("   outcome SD (lon, lat) = (%.1f, %.1f); corr(NAO, PNA) on analysis days = %.2f" % (E.lon.std(), E.lat.std(), np.corrcoef(D.A, D.B)[0, 1]))
import os
os.makedirs(os.path.dirname(sys.argv[1]), exist_ok=True)
open(sys.argv[1], "w").write("\n".join(OUT) + "\n")
