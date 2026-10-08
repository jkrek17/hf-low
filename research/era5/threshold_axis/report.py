"""Assemble RA-23 tables and apply the decision rules of PREREGISTRATION.md. usage: report.py RESULTS_DIR"""
import json, os, sys
import numpy as np, pandas as pd
R = sys.argv[1]
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q
BN = {"atl": "Atlantic", "pac": "Pacific"}
TS = {"t68": 68.0, "t717": 71.7, "t75": 75.0}
rows = []
def add(study, basin, quantity, t, est, lo, hi, p, se):
    rows.append(dict(study=study, basin=basin, quantity=quantity, T=TS[t], t=t, est=est, lo=lo, hi=hi, p=p, se_log=se))
for t in TS:
    d = pd.read_csv(os.path.join(R, f"pr14_{t}.csv"))
    for _, r in d.iterrows():
        add("PR14", r.basin, "RR(HF_T)", t, np.exp(r.b_hf), np.exp(r.boot_lo_hf), np.exp(r.boot_hi_hf), r.pperm_hf, r.se_hf)
        add("PR14", r.basin, "RR(share_T)", t, np.exp(r.b_share), np.exp(r.boot_lo_share), np.exp(r.boot_hi_share), r.pperm_share, r.se_share)
        add("PR14", r.basin, "f", t, r.f_share, r.f_boot_lo, r.f_boot_hi, np.nan, np.nan)
        add("PR14", r.basin, "n_hf", t, r.n_hf, np.nan, np.nan, np.nan, np.nan)
    ts = pd.read_csv(os.path.join(R, f"channels_{t}", "tests.csv"))
    dec = pd.read_csv(os.path.join(R, f"channels_{t}", "decomposition.csv"))
    for b in ("atl", "pac"):
        for tn, qn in (("T3", "RR(HF_T)"), ("T2", "RR(share_T)")):
            x = ts[(ts.basin == b) & (ts.test == tn)].iloc[0]
            add("PR64", b, qn, t, x.est, x.lo, x.hi, x.p, x.se_log)
        x = dec[(dec.basin == b) & (dec.term == "share_frac")].iloc[0]
        add("PR64", b, "share channel part of change", t, x.est, x.lo, x.hi, np.nan, np.nan)
    tt = pd.read_csv(os.path.join(R, f"pr63_{t}", "tests.csv")).set_index("test")
    for tn, qn in (("P1", "P1 share without terrain-type fixes"), ("P3", "P3 mediation"), ("F4_GH_share_orig", "F4 Greenland high on share")):
        x = tt.loc[tn]
        add("PR63", "atl", qn, t, x.RR, x.RR_lo, x.RR_hi, x.p, x.se_boot)
H = pd.DataFrame(rows)
famA = ((H.study == "PR14") & H.quantity.isin(["RR(HF_T)", "RR(share_T)"]) | (H.study == "PR64") & H.quantity.isin(["RR(HF_T)", "RR(share_T)"]) |
        (H.study == "PR63") & H.quantity.isin(["P3 mediation", "F4 Greenland high on share"]) |
        (H.study == "PR63") & (H.quantity == "P1 share without terrain-type fixes") & (H.t == "t75")) & H.t.isin(["t68", "t75"])
assert famA.sum() == 21, famA.sum()
H["q_A"] = np.nan
H.loc[famA, "q_A"] = bh(H.p[famA].values)
ref = H[H.t == "t717"].set_index(["study", "basin", "quantity"])
H["est_717"] = [ref.est.get((s, b, q), np.nan) for s, b, q in zip(H.study, H.basin, H.quantity)]
H["mde_rr"] = np.exp(2.8 * H.se_log)
def verdict(r):
    if pd.isna(r.q_A):
        return ""
    same = np.sign(np.log(r.est)) == np.sign(np.log(r.est_717))
    if same and r.q_A < 0.05:
        return "holds"
    if r.q_A >= 0.05 and r.mde_rr <= np.exp(abs(np.log(r.est_717))):
        return "collapse (well powered)"
    return "can't tell" if same else "sign flip"
H["verdict"] = H.apply(verdict, axis=1)
H.to_csv(os.path.join(R, "headline.csv"), index=False, float_format="%.5g")
P = []
for design, fn in (("PR14 daily", "thr_daily.csv"), ("PR64 weekly", "thr_weekly.csv")):
    c = pd.read_csv(os.path.join(R, fn))
    for b in ("atl", "pac"):
        for st in ("t75-t68", "t75-t717", "t717-t68"):
            x = c[(c.basin == b) & (c.stat == st)].iloc[0]
            P.append(dict(design=design, basin=b, diff=st, est=x.est, lo=x.lo, hi=x.hi, p=x.p_perm, mde_rr=x.mde_rr))
P = pd.DataFrame(P)
P["q_C"] = bh(P.p.values)
P["verdict"] = np.where(P.lo > 1, "grows", np.where(P.hi < 1, "shrinks", "flat"))
P.to_csv(os.path.join(R, "paired.csv"), index=False, float_format="%.5g")
say("RA-23: do the share results move with the gust threshold? ERA5 proxy, pipeline A. Per +1 SD of the index; Oct-Apr 2004-05..2025-26.")
say("Family A: 21 headline tests at 68 and 75 kt, BH. Family C: 12 distinct paired differences, BH. 71.7 kt is the control (original).")
for (study, quantity) in [("PR14", "RR(HF_T)"), ("PR14", "RR(share_T)"), ("PR14", "f"), ("PR64", "RR(HF_T)"), ("PR64", "RR(share_T)"), ("PR64", "share channel part of change"),
                          ("PR63", "P1 share without terrain-type fixes"), ("PR63", "P3 mediation"), ("PR63", "F4 Greenland high on share")]:
    say(f"== {study} {quantity}")
    for b in (("atl",) if study == "PR63" else ("atl", "pac")):
        for t in TS:
            r = H[(H.study == study) & (H.quantity == quantity) & (H.basin == b) & (H.t == t)].iloc[0]
            qs = f" q_A {r.q_A:.4f}" if not pd.isna(r.q_A) else ""
            ps = f" p {r.p:.4f}" if not pd.isna(r.p) else ""
            ms = f" detectable {r.mde_rr:.3f}" if not pd.isna(r.mde_rr) and t != "t717" else ""
            say(f"  {BN[b]:9s} {TS[t]:>5}: {r.est:.3f} [{r.lo:.3f}, {r.hi:.3f}]{ps}{qs}{ms}  {r.verdict}")
say("== Tracks flagged HF (daily-design window)")
for b in ("atl", "pac"):
    say("  ", BN[b], {TS[t]: int(H[(H.study == 'PR14') & (H.quantity == 'n_hf') & (H.basin == b) & (H.t == t)].est.iloc[0]) for t in TS})
say("== Rank: Atlantic share RR above Pacific at each T")
for study in ("PR14", "PR64"):
    for t in TS:
        a = H[(H.study == study) & (H.quantity == "RR(share_T)") & (H.basin == "atl") & (H.t == t)].est.iloc[0]
        p = H[(H.study == study) & (H.quantity == "RR(share_T)") & (H.basin == "pac") & (H.t == t)].est.iloc[0]
        say(f"  {study} {TS[t]}: Atlantic {a:.3f} vs Pacific {p:.3f} -> {'kept' if a > p else 'NOT kept'}")
say("== Paired differences RR(HF_a)/RR(HF_b) (= share differences), family C")
for _, r in P.iterrows():
    say(f"  {r.design:12s} {BN[r.basin]:9s} {r['diff']}: {r.est:.3f} [{r.lo:.3f}, {r.hi:.3f}] p {r.p:.4f} q_C {r.q_C:.3f} detectable {r.mde_rr:.3f} -> {r.verdict}")
open(os.path.join(R, "summary.txt"), "w").write("\n".join(OUT) + "\n")
