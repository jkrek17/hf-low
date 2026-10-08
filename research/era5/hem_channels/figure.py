"""Figure: channel effects per SD of the pattern index, and the absolute decomposition. ERA5 proxy, pipeline A.
usage: figure.py RESULTS_DIR"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]
T = pd.read_csv(os.path.join(R, "tests.csv")); D = pd.read_csv(os.path.join(R, "decomposition.csv"))
fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), gridspec_kw=dict(width_ratios=[1.5, 1]))
lab = {"T1": "more cyclones", "T2": "higher HF share", "T3": "HF lows (total)", "T4": "entrants", "T5": "local storms",
       "T6": "share, entrants", "T7": "share, local"}
col = {"atl": "#1b6ca8", "pac": "#c2571a"}
for r, b in enumerate(["atl", "pac"]):
    a = ax[r, 0]
    t = T[(T.basin == b) & T.test.isin(list(lab))].reset_index(drop=True)
    y = np.arange(len(t))[::-1]
    a.errorbar(t.est, y, xerr=[t.est - t.lo, t.hi - t.est], fmt="o", color=col[b], capsize=3)
    a.axvline(1, color="k", lw=0.8)
    a.set_yticks(y); a.set_yticklabels([lab[k] for k in t.test])
    a.set_xlabel("rate ratio per SD of the pattern index (95% season-block bootstrap)")
    a.set_title({"atl": "Atlantic", "pac": "Pacific"}[b] + " (pipeline A proxy, 2004-05 to 2025-26)", loc="left", fontsize=10)
    a = ax[r, 1]
    d = D[D.basin == b].set_index("term")
    names = ["count_e_frac", "count_l_frac", "share_e_frac", "share_l_frac", "inter_e_frac", "inter_l_frac"]
    nl = ["count\nentrant", "count\nlocal", "share\nentrant", "share\nlocal", "inter.\nentrant", "inter.\nlocal"]
    v = d.loc[names, "est"].values
    lo, hi = d.loc[names, "lo"].values, d.loc[names, "hi"].values
    a.bar(range(6), v, color=col[b], yerr=[v - lo, hi - v], capsize=3)
    a.set_xticks(range(6)); a.set_xticklabels(nl, fontsize=8)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("fraction of the change in HF lows")
    a.set_title("what makes up +1 SD", loc="left", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(R, "channels.png"), dpi=110)
