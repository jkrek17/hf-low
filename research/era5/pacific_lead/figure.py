"""Lag profile figure (S5). Partial correlation of the Pacific index at lag k with the Atlantic index, controls month + Atlantic
index of the previous week. k > 0: the Pacific leads."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
P = pd.read_csv(os.path.join(R, "lag_profile.csv"))
N = json.load(open(os.path.join(R, "posthoc_naive_xcorr.json")))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for a, (lab, d) in zip(ax, P.groupby("data", sort=False)):
    a.fill_between(d.lag, d.null_lo, d.null_hi, color="0.85", label="permutation null (95%)")
    a.errorbar(d.lag, d.partial_r, yerr=[d.partial_r - d.boot_lo, d.boot_hi - d.partial_r], fmt="o-", color="C0",
               capsize=3, label="partial r (95% season bootstrap)")
    a.plot(sorted(map(int, N[lab])), [N[lab][str(k)] for k in sorted(map(int, N[lab]))], "s--", color="C3", ms=4,
           label="month-only correlation (post hoc)")
    a.axhline(0, color="k", lw=0.5); a.axvline(0, color="k", lw=0.5, ls=":")
    a.set_title(lab, fontsize=9); a.set_xlabel("lag k weeks (k > 0: Pacific leads the Atlantic)")
ax[0].set_ylabel("correlation with the Atlantic index"); ax[0].legend(fontsize=7, loc="lower right")
fig.suptitle("RA-5: lead-lag of the Pacific and Atlantic HF pattern indices (ERA5 proxy)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(R, "lag_profile.png"), dpi=130)
