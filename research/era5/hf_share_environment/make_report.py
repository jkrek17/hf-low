"""Summary tables, verdicts under the pre-registered rule, and the figure.

    python3 -I make_report.py <results dir>
"""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]
T = pd.read_csv(f"{R}/tests.csv")
P = pd.read_csv(f"{R}/monthly_profile.csv")
PRIMARY = ["jet250", "eady", "sstgrad", "sst_t500", "flux", "tcwv", "div300", "vadv500"]
BASIN = {"atl": "Atlantic", "pac": "Pacific"}


def g(variant, basin, quantity, ing):
    r = T[(T.variant == variant) & (T.basin == basin) & (T.quantity == quantity) & (T.ingredient.fillna("") == ing)]
    return r.iloc[0] if len(r) else None


def verdict(variant, basin, x):
    d, b, f = g(variant, basin, "delta_rise", x), g(variant, basin, "beta", x), g(variant, basin, "F_rise", x)
    if d is None or b is None or f is None:
        return ""
    if d.q < 0.05 and b.q < 0.05 and d.estimate * b.estimate > 0 and f.lo > 0:
        return "lines up"
    if f.hi < 0.15 or ((d.q >= 0.05 or b.q >= 0.05) and b.lo > -0.10 and b.hi < 0.10):
        return "does not line up"
    return "cannot tell"


def fmt(r, nd=2):
    return f"{r.estimate:+.{nd}f} ({r.lo:+.{nd}f} to {r.hi:+.{nd}f})"


out = []
for variant, title in (("primary", "PRIMARY: first deepening fix"), ("S1_entry", "S1: entry fix (first 00/12 fix)"),
                       ("S2_maxdeep", "S2: fix of greatest 24 h deepening (uses the future of the track)"),
                       ("S3_fixlevel", "S3: fix level, HF onset within 24 h"),
                       ("S6a_commoncc", "S6a: common complete-case sample"), ("S6b_deepening_only", "S6b: deepening cyclones only")):
    for basin in ("atl", "pac"):
        out.append(f"\n== {title} | {BASIN[basin]} ==")
        s = g(variant, basin, "share_ON", "")
        if s is not None:
            sh = {k: g(variant, basin, f"share_{k}", "").estimate for k in ("ON", "DJF", "MA")}
            out.append(f"HF share ON {sh['ON']*100:.1f}%  DJF {sh['DJF']*100:.1f}%  MA {sh['MA']*100:.1f}%   n = {int(g(variant, basin, 'n', '').estimate)}")
        c0 = g(variant, basin, "c0_rise", PRIMARY[0])
        out.append("ingredient | Delta DJF-ON (SD) | q | beta per SD (log-odds) | q | F_rise | unique F | verdict")
        for x in PRIMARY:
            d, b, f, u = (g(variant, basin, q, x) for q in ("delta_rise", "beta", "F_rise", "Funique_rise"))
            out.append(f"{x:9s} | {fmt(d)} | {d.q:.3f} | {fmt(b, 3)} | {b.q:.3f} | {fmt(f)} | "
                       f"{fmt(u) if u is not None else '-'} | {verdict(variant, basin, x)}")
        j = g(variant, basin, "F_rise", "JOINT")
        if j is not None:
            out.append(f"JOINT F_rise {fmt(j)}   F_fall {fmt(g(variant, basin, 'F_fall', 'JOINT'))}   n = {int(g(variant, basin, 'n', 'JOINT').estimate)}")
        if variant == "primary":
            out.append("F_fall (DJF vs MA): " + "; ".join(f"{x} {fmt(g(variant, basin, 'F_fall', x))}" for x in PRIMARY))
            out.append("Delta DJF-MA (SD), q: " + "; ".join(f"{x} {g(variant, basin, 'delta_fall', x).estimate:+.2f} (q {g(variant, basin, 'delta_fall', x).q:.3f})" for x in PRIMARY))
            out.append("Interaction (X x DJF) on log-odds: " + "; ".join(f"{x} {fmt(g(variant, basin, 'inter', x), 3)} q {g(variant, basin, 'inter', x).q:.2f}" for x in PRIMARY))
            out.append("Two stages (log-odds DJF-ON): share " + fmt(g(variant, basin, "share_rise", ""), 3) + "; P(deepening fix) " +
                       fmt(g(variant, basin, "deep_rise", ""), 3) + "; P(HF | deepening) " + fmt(g(variant, basin, "hfgd_rise", ""), 3))
            out.append("  P(deepening fix) ON/DJF/MA: " + "/".join(f"{g(variant, basin, f'pdeep_{k}', '').estimate:.3f}" for k in ("ON", "DJF", "MA")) +
                       "; P(HF | deepening): " + "/".join(f"{g(variant, basin, f'hfgd_{k}', '').estimate:.3f}" for k in ("ON", "DJF", "MA")))
for basin in ("atl", "pac"):
    out.append(f"\n== S5 secondary ingredients | {BASIN[basin]} ==")
    for x in ["lat", "B", "VTL", "VTU", "sst"]:
        d, b, f = (g("S5_secondary", basin, q, x) for q in ("delta_rise", "beta", "F_rise"))
        out.append(f"{x:5s} | Delta {fmt(d)} q {d.q:.3f} | beta {fmt(b, 3)} q {b.q:.3f} | F {fmt(f)}")
    j = g("S5_secondary", basin, "F_rise", "JOINT")
    out.append(f"JOINT (secondary 5) F_rise {fmt(j)}")
open(f"{R}/summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))

# figure: monthly mean z of each ingredient and share; explained fractions
order = [10, 11, 12, 1, 2, 3, 4]
fig, ax = plt.subplots(2, 3, figsize=(15, 8))
for i, basin in enumerate(("atl", "pac")):
    a = ax[i, 0]
    sh = P[(P.basin == basin) & (P.ingredient == "jet250")].set_index("month").loc[order]
    a.bar(range(7), sh.hf_share * 100, color="#bbbbbb")
    a.set_xticks(range(7)); a.set_xticklabels(["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"])
    a.set_ylabel("HF share of cyclones, %"); a.set_title(f"{BASIN[basin]}: share (pipeline A, proxy)")
    a = ax[i, 1]
    for x in PRIMARY:
        z = P[(P.basin == basin) & (P.ingredient == x)].set_index("month").loc[order]
        a.plot(range(7), z.mean_z, marker="o", ms=3, label=x)
    a.axhline(0, color="k", lw=0.5); a.set_xticks(range(7)); a.set_xticklabels(["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"])
    a.set_ylabel("mean at first deepening fix, SD"); a.set_title(f"{BASIN[basin]}: environment by genesis month")
    if i == 0:
        a.legend(fontsize=7, ncol=2)
    a = ax[i, 2]
    xs = PRIMARY + ["JOINT"]
    f = [g("primary", basin, "F_rise", x) for x in xs]
    a.errorbar([r.estimate for r in f], range(len(xs)), xerr=[[r.estimate - r.lo for r in f], [r.hi - r.estimate for r in f]], fmt="o", color="#225")
    a.axvline(0, color="k", lw=0.5); a.axvline(1, color="k", lw=0.5, ls=":")
    a.set_yticks(range(len(xs))); a.set_yticklabels(xs); a.invert_yaxis()
    a.set_xlabel("fraction of the Dec-Feb vs Oct-Nov log-odds rise explained"); a.set_title(f"{BASIN[basin]}: explained fraction, 95% season-block interval")
plt.tight_layout(); plt.savefig(f"{R}/ingredients.png", dpi=110)
