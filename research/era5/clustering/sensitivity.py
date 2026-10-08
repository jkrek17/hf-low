"""Pre-registered sensitivities for the archive tier (no decision weight; see PREREGISTRATION.md, M1).

variants: base; harm (two-harmonic cycle instead of month dummies); shift3 (block grid shifted by 3 days);
octmar (Oct-Mar window, M1 only); clean (archive ID-quality flags and the second member of each
collision pair removed).  usage: sensitivity.py OUT_DIR [--nsim N]
"""
import argparse
import json
import os

import numpy as np

import cluster as C


def one(ev, s0, s1, rng, nsim, nperm, harm=False, off=0, nd=C.NDAY, m23=True):
    S = s1 - s0 + 1
    tr = np.arange(S) - (S - 1) / 2
    y = C.daily(ev, s0, s1)[:, :nd]
    X = C.make_X("N1", S, tr, None, None, nd, harm)
    D = C.dispersion(y, X, (7, 30), rng, nsim, off)
    out = dict(E_W7=float(D["E"][0]), p_W7=float(D["p"][0]), E_W30=float(D["E"][1]), p_W30=float(D["p"][1]), n=len(ev))
    if m23:
        r = C.m2(ev, D["mu"], s0, rng, nsim)
        out.update(G48_ratio=float(r["ratio_G48"]), p_G48=float(r["p"]["G48"]))
        K = C.m3(ev, rng, nperm)
        out.update(knox_ratio=float(K["ratio"][1, 1]), p_knox=float(K["p"][1, 1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--nsim", type=int, default=1000)
    ap.add_argument("--nperm", type=int, default=2000)
    a = ap.parse_args()
    rng = np.random.default_rng(C.SEED + 5)
    rows, cfg = C.tier_data("arch")
    s0, s1 = cfg["s0"], cfg["s1"]
    res = {}
    for basin in ("atl", "pac"):
        ev, _ = C.window(rows, basin, s0, s1)
        clean = [e for e in ev if e["idok"] and not e["second"]]
        res[basin] = {
            "base": one(ev, s0, s1, rng, a.nsim, a.nperm),
            "harm": one(ev, s0, s1, rng, a.nsim, a.nperm, harm=True),
            "shift3": one(ev, s0, s1, rng, a.nsim, a.nperm, off=3),
            "octmar": one(ev, s0, s1, rng, a.nsim, a.nperm, nd=182, m23=False),
            "clean": one(clean, s0, s1, rng, a.nsim, a.nperm),
        }
        print(basin, json.dumps(res[basin]), flush=True)
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, "arch_sensitivity.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
