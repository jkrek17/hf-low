"""Assemble the RA-22 tables from the committed results and apply the decision rules of PREREGISTRATION.md.
usage: report.py RESULTS_DIR      writes RESULTS_DIR/headline.csv, contrasts.csv, paired.csv, secondary.csv, summary.txt
"""
import json, os, sys
import numpy as np, pandas as pd

R = sys.argv[1]
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


BN = {"atl": "Atlantic", "pac": "Pacific"}
rows = []        # one row per (study, basin, quantity, k), k = 1, 2, 3


def add(study, basin, quantity, k, est, lo, hi, p, se_log, log_scale=True, pipeline="A proxy"):
    rows.append(dict(study=study, basin=basin, quantity=quantity, k=k, est=est, lo=lo, hi=hi, p=p, se_log=se_log, pipeline=pipeline))


# ---- PR 14 (daily genesis counts, NAO Atlantic / PNA Pacific)
for k in (1, 2, 3):
    d = pd.read_csv(os.path.join(R, f"pr14_k{k}.csv"))
    for _, r in d.iterrows():
        b = r.basin
        add("PR14", b, "RR(HF_k)", k, np.exp(r.b_hf), np.exp(r.boot_lo_hf), np.exp(r.boot_hi_hf), r.pperm_hf, r.se_hf)
        add("PR14", b, "RR(share_k)", k, np.exp(r.b_share), np.exp(r.boot_lo_share), np.exp(r.boot_hi_share), r.pperm_share, r.se_share)
        add("PR14", b, "f", k, r.f_share, r.f_boot_lo, r.f_boot_hi, np.nan, np.nan)
        add("PR14", b, "RR(all)", k, np.exp(r.b_all), np.exp(r.boot_lo_all), np.exp(r.boot_hi_all), r.pperm_all, r.se_all)
# ---- PR 64 (weekly pattern index)
for k in (1, 2, 3):
    t = pd.read_csv(os.path.join(R, f"channels_k{k}", "tests.csv"))
    dec = pd.read_csv(os.path.join(R, f"channels_k{k}", "decomposition.csv"))
    for b in ("atl", "pac"):
        for tn, qn in (("T3", "RR(HF_k)"), ("T2", "RR(share_k)"), ("T1", "RR(all)")):
            x = t[(t.basin == b) & (t.test == tn)].iloc[0]
            add("PR64", b, qn, k, x.est, x.lo, x.hi, x.p, x.se_log)
        for term, qn in (("f_share_of_log_RR_HF", "f"), ("share_frac", "share channel part of change"), ("count_frac", "count channel part of change")):
            x = dec[(dec.basin == b) & (dec.term == term)].iloc[0]
            add("PR64", b, qn, k, x.est, x.lo, x.hi, np.nan, np.nan)
# ---- PR 41 (archive weekly, same index)
ca = pd.read_csv(os.path.join(R, "contrast_archive_withm.csv"))
for k in (1, 2, 3):
    for b in ("atl", "pac"):
        x = ca[(ca.basin == b) & (ca.stat == f"h{k}")].iloc[0]
        add("PR41", b, "RR(HF_k) archive", k, x.est, x.lo, x.hi, x.p_perm, x.se_log, pipeline="archive")
# ---- PR 63 (Atlantic)
for k in (1, 2, 3):
    t = pd.read_csv(os.path.join(R, f"pr63_k{k}", "tests.csv")).set_index("test")
    h = json.load(open(os.path.join(R, f"pr63_k{k}", "headline.json")))
    for tn, qn in (("P1", "P1 share without terrain-type fixes (RR)"), ("P2", "P2 retained: share without terrain-type fixes / share (descriptive)"),
                   ("P3", "P3 mediation (share given Greenland high / share)"), ("F4_GH_share_orig", "F4 Greenland high on share, per SD")):
        x = t.loc[tn]
        add("PR63", "atl", qn, k, x.RR, x.RR_lo, x.RR_hi, x.p, x.se_boot)
    add("PR63", "atl", "NAO share given Greenland high (RR)", k, h["RR_share_NAO_given_GH"], h["RR_share_NAO_given_GH_ci"][0], h["RR_share_NAO_given_GH_ci"][1], np.nan, np.nan)
    add("PR63", "atl", "NAO share, original (RR)", k, h["orig_RR_share"], h["orig_RR_share_ci"][0], h["orig_RR_share_ci"][1], np.nan, np.nan)
H = pd.DataFrame(rows)

# ---- family A: 26 headline tests at k = 2, 3
FAM_A = [("PR14", "RR(HF_k)"), ("PR14", "RR(share_k)"), ("PR64", "RR(HF_k)"), ("PR64", "RR(share_k)"), ("PR41", "RR(HF_k) archive"),
         ("PR63", "P1 share without terrain-type fixes (RR)"), ("PR63", "P3 mediation (share given Greenland high / share)"),
         ("PR63", "F4 Greenland high on share, per SD")]
isA = H.apply(lambda r: (r.study, r.quantity) in FAM_A and r.k in (2, 3), axis=1)
assert isA.sum() == 26, isA.sum()
H["q_A"] = np.nan
H.loc[isA, "q_A"] = bh(H.p[isA].values)
# k = 1 reference for sign and size
ref = H[H.k == 1].set_index(["study", "basin", "quantity"])
H["est_k1"] = [ref.est.get((s, b, q), np.nan) for s, b, q in zip(H.study, H.basin, H.quantity)]
H["se_k1"] = [ref.se_log.get((s, b, q), np.nan) for s, b, q in zip(H.study, H.basin, H.quantity)]
H["same_sign_as_k1"] = np.sign(np.log(H.est)) == np.sign(np.log(H.est_k1))
H["mde_rr"] = np.exp(2.8 * H.se_log)
# paired differences from the contrast tables
pair = {}
for design, fn, key in (("PR14", "contrast_daily_withm.csv", None), ("PR64", "contrast_weekly_withm.csv", None), ("PR41", "contrast_archive_withm.csv", None)):
    c = pd.read_csv(os.path.join(R, fn))
    for b in ("atl", "pac"):
        for k in (2, 3):
            x = c[(c.basin == b) & (c.stat == f"h{k}-h1")].iloc[0]
            pair[(design, b, k)] = (x.est, x.lo, x.hi, x.p_perm)


def classify(r):
    if r.k == 1 or pd.isna(r.q_A):
        return ""
    if r.same_sign_as_k1 and r.q_A < 0.05:
        v = "holds"
    elif r.q_A >= 0.05 and r.mde_rr <= np.exp(abs(np.log(r.est_k1))):
        v = "lost, well powered"
    else:
        v = "can't tell"
    key = (r.study, r.basin, r.k)
    if key in pair and r.quantity in ("RR(HF_k)", "RR(share_k)", "RR(HF_k) archive"):
        e, lo, hi, p = pair[key]
        v += "; vs k=1: " + ("larger" if lo > 1 else "smaller" if hi < 1 else "same within noise")
    elif r.study == "PR63":
        v += "; vs k=1: paired difference not computed"
    return v


H["verdict"] = H.apply(classify, axis=1)
H.to_csv(os.path.join(R, "headline.csv"), index=False, float_format="%.5g")

# ---- contrasts (family B) and paired differences (family C)
cont = []
for design, fn in (("PR14 daily", "contrast_daily_withm.csv"), ("PR64 weekly", "contrast_weekly_withm.csv"), ("PR41 archive", "contrast_archive_withm.csv")):
    c = pd.read_csv(os.path.join(R, fn))
    for b in ("atl", "pac"):
        x = c[(c.basin == b) & (c.stat == "h3-b")].iloc[0]
        cont.append(dict(design=design, basin=b, est=x.est, lo=x.lo, hi=x.hi, p=x.p_perm, se_log=x.se_log, mde_rr=x.mde_rr))
        for st in ("b", "m", "h3"):
            y = c[(c.basin == b) & (c.stat == st)].iloc[0]
            cont[-1][f"RR_{st}"] = y.est
C = pd.DataFrame(cont)
C["q_B"] = bh(C.p.values)
C["verdict"] = np.where((C.est > 1) & (C.q_B < 0.05), "favours sustained", np.where((C.est < 1) & (C.q_B < 0.05), "favours brief", "no clear preference"))
C.to_csv(os.path.join(R, "contrasts.csv"), index=False, float_format="%.5g")
P = []
for design, fn in (("PR14 daily", "contrast_daily_withm.csv"), ("PR64 weekly", "contrast_weekly_withm.csv"), ("PR41 archive", "contrast_archive_withm.csv")):
    c = pd.read_csv(os.path.join(R, fn))
    for b in ("atl", "pac"):
        for st in ("h2-h1", "h3-h1"):
            x = c[(c.basin == b) & (c.stat == st)].iloc[0]
            P.append(dict(design=design, basin=b, diff=st, est=x.est, lo=x.lo, hi=x.hi, p=x.p_perm))
P = pd.DataFrame(P)
P["q_C"] = bh(P.p.values)
P.to_csv(os.path.join(R, "paired.csv"), index=False, float_format="%.5g")

# ---- secondary: consecutive run and count-matched gust cut (PR 14 and PR 64)
S = []
for tag in ("k", "run", "gcut"):
    for k in (2, 3):
        d = pd.read_csv(os.path.join(R, f"pr14_{tag}{k}.csv"))
        t = pd.read_csv(os.path.join(R, f"channels_{tag}{k}", "tests.csv"))
        for _, r in d.iterrows():
            b = r.basin
            t3 = t[(t.basin == b) & (t.test == "T3")].iloc[0]; t2 = t[(t.basin == b) & (t.test == "T2")].iloc[0]
            S.append(dict(label={"k": "total fixes", "run": "consecutive run", "gcut": "count-matched gust cut"}[tag], k=k, basin=b, n_hf=int(r.n_hf),
                          pr14_RR_HF=np.exp(r.b_hf), pr14_RR_HF_lo=np.exp(r.boot_lo_hf), pr14_RR_HF_hi=np.exp(r.boot_hi_hf), pr14_RR_share=np.exp(r.b_share),
                          pr14_p_share=r.pperm_share, pr64_RR_HF=t3.est, pr64_RR_HF_lo=t3.lo, pr64_RR_HF_hi=t3.hi, pr64_RR_share=t2.est, pr64_p_share=t2.p))
S = pd.DataFrame(S)
S.to_csv(os.path.join(R, "secondary.csv"), index=False, float_format="%.5g")

# ---- text summary
say("RA-22: do the share results hold for sustained HF lows? ERA5 proxy (pipeline A) except PR41 (archive). Per +1 SD of the index.")
say(f"Family A: {int(isA.sum())} headline tests (k = 2, 3), BH. Family B: 6 sustained-minus-brief contrasts, BH. Family C: 12 paired differences, BH.")
say("")
for (study, quantity) in FAM_A + [("PR14", "f"), ("PR64", "f"), ("PR64", "share channel part of change"), ("PR64", "RR(all)"), ("PR14", "RR(all)"),
                                  ("PR63", "P2 retained: share without terrain-type fixes / share (descriptive)"), ("PR63", "NAO share, original (RR)"), ("PR63", "NAO share given Greenland high (RR)")]:
    say(f"== {study} {quantity}")
    for b in (("atl",) if study == "PR63" else ("atl", "pac")):
        for k in (1, 2, 3):
            r = H[(H.study == study) & (H.quantity == quantity) & (H.basin == b) & (H.k == k)].iloc[0]
            qs = f" q_A {r.q_A:.4f}" if not pd.isna(r.q_A) else ""
            ps = f" p {r.p:.4f}" if not pd.isna(r.p) else ""
            ms = f" detectable {r.mde_rr:.3f}" if not pd.isna(r.mde_rr) and k > 1 else ""
            say(f"  {BN[b]:9s} k={k}: {r.est:.3f} [{r.lo:.3f}, {r.hi:.3f}]{ps}{qs}{ms}  {r.verdict}")
say("")
say("== Sustained (>=3 fixes) minus brief (exactly 1 fix): ratio of RR per SD, family B")
for _, r in C.iterrows():
    say(f"  {r.design:13s} {BN[r.basin]:9s} RR brief {r.RR_b:.3f}, middle {r.RR_m:.3f}, sustained {r.RR_h3:.3f}; sustained/brief {r.est:.3f} [{r.lo:.3f}, {r.hi:.3f}] p {r.p:.4f} q_B {r.q_B:.3f} detectable {r.mde_rr:.3f}  -> {r.verdict}")
nsus = int((C.verdict == "favours sustained").sum()); nbri = int((C.verdict == "favours brief").sum()); npos = int((C.est > 1).sum())
overall = ("favours sustained" if nsus >= 4 and nbri == 0 else "favours brief" if nbri >= 4 and nsus == 0 else
           ("no clear preference; leans sustained" if npos >= 4 else "no clear preference; leans brief" if npos <= 2 else "no clear preference"))
say(f"  cells favouring sustained {nsus}, brief {nbri}, point estimates above 1: {npos} of 6 -> overall: {overall}")
say("")
say("== Paired differences RR(HF_k)/RR(HF_1), family C")
for _, r in P.iterrows():
    say(f"  {r.design:13s} {BN[r.basin]:9s} {r['diff']}: {r.est:.3f} [{r.lo:.3f}, {r.hi:.3f}] p {r.p:.4f} q_C {r.q_C:.3f}")
say("")
say("== Secondary (S1 consecutive run, S2 count-matched gust cut; own family, none changes a verdict)")
for _, r in S.iterrows():
    say(f"  {r.label:23s} k={r.k} {BN[r.basin]:9s} n {r.n_hf}: PR14 RR(HF) {r.pr14_RR_HF:.3f} [{r.pr14_RR_HF_lo:.3f}, {r.pr14_RR_HF_hi:.3f}] share {r.pr14_RR_share:.3f} (p {r.pr14_p_share:.4f}); "
        f"PR64 RR(HF) {r.pr64_RR_HF:.3f} [{r.pr64_RR_HF_lo:.3f}, {r.pr64_RR_HF_hi:.3f}] share {r.pr64_RR_share:.3f} (p {r.pr64_p_share:.4f})")
open(os.path.join(R, "summary.txt"), "w").write("\n".join(OUT) + "\n")
