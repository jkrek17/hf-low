"""2025-26 only: GEFS vs ECMWF ENS vs AIFS ENS on identical cases (leads 48 and 96), raw member fractions, no calibration.
BSS against the season base rate; paired differences with a block bootstrap over ISO weeks of run date (2,000 resamples, seed 20261028).
usage: compare.py gefs25.npz ens.npz aifs.npz OUT"""
import sys, datetime as dt, numpy as np
import stats as S
from run_primary import load

def per_date(K, n, y, truth, basin, thr, lead):
    out = {}
    for i, k in enumerate(K):
        if k[1] == truth and k[2] == basin and int(k[3]) == thr and int(k[5]) == lead:
            N = int(k[6]); p = np.arange(N + 1) / N
            out[k[4]] = dict(N=N, n=int(n[i][:N + 1].sum()), y=int(y[i][:N + 1].sum()), b=S.brier(n[i][:N + 1], y[i][:N + 1], p), nk=n[i][:N + 1], yk=y[i][:N + 1])
    return out

if __name__ == "__main__":
    M = {"GEFS": load(sys.argv[1]), "ENS": load(sys.argv[2]), "AIFS ENS": load(sys.argv[3])}
    rng = np.random.default_rng(S.SEED); lines = []; res = {}
    for basin in ("natl", "npac"):
        for thr, truth in ((48, "ascat"), (64, "opc")):
            for lead in (48, 96):
                D = {m: per_date(*M[m], truth, basin, thr, lead) for m in M}
                dates = sorted(D["GEFS"]); wk = {d: dt.datetime.strptime(d, "%Y%m%d").isocalendar()[1] for d in dates}
                weeks = sorted(set(wk.values())); ev = sum(D["GEFS"][d]["y"] for d in dates); cases = sum(D["GEFS"][d]["n"] for d in dates)
                tag = f"{basin} {thr}kt {truth} {lead}h"
                if ev < S.MIN_EVENTS: lines.append(f"{tag}: too few events ({ev}) to compare"); continue
                base = ev / cases
                clim = {d: D["GEFS"][d]["y"] * (1 - base) ** 2 + (D["GEFS"][d]["n"] - D["GEFS"][d]["y"]) * base ** 2 for d in dates}
                W = {w: [d for d in dates if wk[d] == w] for w in weeks}
                agg = lambda m, ws: sum(D[m][d]["b"] for w in ws for d in W[w]); cl = lambda ws: sum(clim[d] for w in ws for d in W[w])
                bss = {m: 1 - agg(m, weeks) / cl(weeks) for m in M}
                boots = {m: [] for m in M}; diffs = {("ENS", "GEFS"): [], ("AIFS ENS", "GEFS"): [], ("AIFS ENS", "ENS"): []}
                for _ in range(S.NBOOT):
                    ws = list(rng.choice(weeks, len(weeks))); c = cl(ws)
                    if c <= 0: continue
                    b = {m: 1 - agg(m, ws) / c for m in M}
                    for m in M: boots[m].append(b[m])
                    for a, bb in diffs: diffs[(a, bb)].append(b[a] - b[bb])
                ci = lambda v: np.percentile(v, [2.5, 97.5])
                txt = f"{tag}: runs {len(dates)}, cases {cases:,}, events {ev:,} | " + "; ".join(f"{m} BSS {bss[m]:+.3f} ({ci(boots[m])[0]:+.3f} to {ci(boots[m])[1]:+.3f})" for m in M)
                txt += " | " + "; ".join(f"{a} minus {b} {bss[a]-bss[b]:+.3f} ({ci(v)[0]:+.3f} to {ci(v)[1]:+.3f})" for (a, b), v in diffs.items())
                # reliability by bin per model
                for m in M:
                    nk = sum(D[m][d]["nk"] for d in dates); yk = sum(D[m][d]["yk"] for d in dates); nb, yb = S.reliability(nk, yk, D[m]["20251015"]["N"] if "20251015" in D[m] else D[m][dates[0]]["N"])
                    txt += f"\n    {m} rate by bin: " + ", ".join(f"{int(a)}:{(b/a*100 if a else float('nan')):.1f}%" for a, b in zip(nb, yb))
                lines.append(txt)
    open(sys.argv[4], "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
