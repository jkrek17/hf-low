"""Skill tests for the W / C omega groups (see PREREGISTRATION.md). Reuses ../hf_jet_trough/skill.py and ../intensity_extra/evaluate.py.

usage: skill_w.py OUTDIR      (needs results/features_w.csv.gz from qc_w.py; results/qc.txt must say PASS)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hf_jet_trough")); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import skill as S0
import evaluate as E
from features import J2, T2

BASE = S0.BASE
WM = ["a500_1000", "a700_1000", "a500_500", "a700_500", "aarea500", "adx", "ady", "afwd500", "nohead_w"]
WX = ["a500_x750"]
CC = ["nasc", "asc2", "asc2dist", "asc2cos", "asc2sin", "couple_w"]
CONT = {"a500_1000", "a700_1000", "a500_500", "a700_500", "aarea500", "adx", "ady", "afwd500", "a500_x750", "asc2", "asc2dist", "asc2cos", "asc2sin"}


def build():
    D, q75 = S0.build()
    X = pd.read_csv(os.path.join(HERE, "results", "features_w.csv.gz"), dtype={"time": str})
    X["month"] = pd.to_datetime(X.time, format="%Y%m%d%H").dt.month
    allc = [c for c in X.columns if c not in ("track", "time", "basin", "season", "lat", "lon", "heading", "month")]
    for c in allc:                                   # basin x month standardisation over every extracted fix (predictors only)
        if c.split("_o")[0] in CONT or c.replace("_m", "").replace("_o", "") in CONT:
            g = X.groupby(["basin", "month"])[c]
            X[c] = (X[c] - g.transform("mean")) / (g.transform("std") + 1e-9)
    X = X[["track", "time"] + allc]
    D = D.merge(X, on=["track", "time"], how="left", validate="one_to_one")
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    for L in (12, 24, 48):
        key = (t - pd.Timedelta(hours=L)).dt.strftime("%Y%m%d%H")
        Z = X.rename(columns={c: f"{c}_l{L}" for c in allc}).rename(columns={"time": "ktime"})
        Zl = pd.DataFrame({"track": D.track, "ktime": key}).merge(Z, on=["track", "ktime"], how="left")
        for c in allc:
            D[f"{c}_l{L}"] = Zl[f"{c}_l{L}"].values
    D["da500_1000_24"] = D.a500_1000 - D.a500_1000_l24
    D["dasc2_24"] = D.asc2 - D.asc2_l24
    D["dasc2_o_24"] = D.asc2_o - D.asc2_o_l24
    D["dasc2_m_24"] = D.asc2_m - D.asc2_m_l24
    return D, q75


def grp(cols, lags, ch=()):
    return S0.group(cols, lags, "", ch if 24 in lags else ())


def cgrp(suf, lags):
    cols = [c + suf for c in CC]
    return S0.group(cols, lags, "", (f"dasc2{suf}_24",) if 24 in lags else ())


LAGS = {"0": (), "0,24": (24,), "all": (12, 24, 48)}


def models_w(lagset, base_extra=()):
    lags = LAGS[lagset]
    W_, C_ = grp(WM, lags, ("da500_1000_24",)), cgrp("", lags)
    return {"W": W_, "C": C_, "W+C": list(dict.fromkeys(W_ + C_))}


def run(D, base, models, target, mask=None):
    m = np.ones(len(D), bool) if mask is None else mask
    S, P, C = E.loso(D, {"base": base, **models}, target, m)
    return S, P, C, S[target].astype(int).values


def fam(S, P, C, y, tname, names, mask=None):
    T, _ = E.gains(S, P, C, y, tname, names, mask=mask)
    T["q_family"] = E.bh(T.p.values)
    return T


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    assert "PASS" in open(os.path.join(HERE, "results", "qc.txt")).read(), "QC did not pass"
    D, q75 = build()
    print(len(D), "fixes", D.track.nunique(), "tracks", D.season.nunique(), "seasons; hf24", int(D.hf24.sum()),
          "NaN share", {c: round(float(D[c].isna().mean()), 3) for c in ("a500_1000", "adx", "asc2dist", "a500_1000_l24")}, flush=True)
    JT = list(dict.fromkeys(S0.group(J2, (24,), "", ("dvmaxp_24",)) + S0.group(T2, (24,), "", ("dtdepth_24", "dtamp_24"))))
    BASE2 = BASE + JT
    mods = {}
    for ls in LAGS:
        for k, v in models_w(ls).items():
            mods[f"{k} [lags {ls}]"] = BASE + v
    mods["couple_w"] = BASE + ["couple_w"]
    mods["Wx [lags 0,24]"] = BASE + grp(WX, (24,))
    mods["C first-written rule [lags 0,24]"] = BASE + cgrp("_o", (24,))
    mods["C 0.6 Pa/s [lags 0,24]"] = BASE + cgrp("_m", (24,))
    S, P, C = E.loso(D, {"base": BASE, **mods}, "hf24", np.ones(len(D), bool))
    y = S["hf24"].astype(int).values
    np.savez_compressed(f"{out}/hf24_probs.npz", track=S.track.values, time=S.time.values.astype(str), season=S.season.values, y=y,
                        **{k.replace(" ", "_"): v for k, v in P.items()})
    L24 = [f"{k} [lags 0,24]" for k in ("W", "C", "W+C")]
    rows = [fam(S, P, C, y, "hf24 [primary]", L24),
            fam(S, P, C, y, "hf24 lags 0", [f"{k} [lags 0]" for k in ("W", "C", "W+C")]),
            fam(S, P, C, y, "hf24 lags 0,12,24,48", [f"{k} [lags all]" for k in ("W", "C", "W+C")]),
            fam(S, P, C, y, "hf24 onset fixes", L24, mask=~S.hf_now.values.astype(bool)),
            fam(S, P, C, y, "hf24 couple_w alone / x750 ring / other centre rules", ["couple_w", "Wx [lags 0,24]", "C first-written rule [lags 0,24]", "C 0.6 Pa/s [lags 0,24]"])]
    for b in ("atl", "pac"):
        rows.append(fam(S, P, C, y, f"hf24 [primary] {b} (descriptive)", L24, mask=(S.basin == b).values))
    m24 = {f"{k} [lags 0,24]": BASE + v for k, v in models_w("0,24").items()}
    for tg, tname, msk in (("rd", "rapid deepening", (D.rd >= 0).values), ("hf48i", "hf48", (D.hf48i >= 0).values)):
        S2, P2, C2, y2 = run(D, BASE, m24, tg, msk)
        rows.append(fam(S2, P2, C2, y2, tname, list(m24)))
    # incremental to the jet and trough features (base = BASE + J2 + T2 at lags 0,24)
    m24i = {f"{k} [lags 0,24] over J2+T2": BASE2 + [c for c in v if c not in BASE2] for k, v in models_w("0,24").items()}
    for tg, tname, msk in (("hf24", "hf24 over J2+T2", None), ("rd", "rapid deepening over J2+T2", (D.rd >= 0).values)):
        S3, P3, C3, y3 = run(D, BASE2, m24i, tg, msk)
        rows.append(fam(S3, P3, C3, y3, tname, list(m24i)))
    T = pd.concat(rows, ignore_index=True)
    T["verdict"] = [S0.verdict(r) for r in T.itertuples()]
    main = T[~T.target.str.contains("descriptive")]
    T["q_all"] = np.nan
    T.loc[main.index, "q_all"] = E.bh(main.p.values)
    T.to_csv(f"{out}/gains.csv", index=False, float_format="%.5g")
    passing = T[(T.verdict == "ADDS SKILL") & (T.target == "hf24 [primary]")]
    pd.DataFrame([dict(model=r.model, **dict(zip(("d_hss", "lo90", "hi90", "hss_base", "hss"), E.hss_gain(S, P, C, y, r.model))))
                  for r in passing.itertuples()]).to_csv(f"{out}/hss.csv", index=False, float_format="%.4g")
    with open(f"{out}/skill.txt", "w") as f:
        f.write(f"Sample: {len(D)} Oct-Apr fixes, {D.track.nunique()} tracks, 18 seasons 2004-05..2021-22; hf24 events {int(D.hf24.sum())}\n"
                f"Base (PR 12 'full', refit on this sample) BSS {T.bss_base.iloc[0]:+.4f}. Gain = BSS(base+group) - BSS(base), 90% CI over seasons.\n")
        for tg, g in T.groupby("target", sort=False):
            f.write(f"\n== {tg}  (base BSS {g.bss_base.iloc[0]:+.4f})\n")
            for r in g.itertuples():
                f.write(f"  {r.model:48s} gain {r.gain:+.4f} [{r.lo90:+.4f}, {r.hi90:+.4f}] BSS {r.bss:+.4f} p {r.p:.4f} q(fam) {r.q_family:.3f} "
                        f"q(all) {r.q_all:.3f} seasons better {r.seasons_better}/{r.n_seasons}  -> {r.verdict}\n")
        n = len(main)
        f.write(f"\nTests passing 'adds skill' (non-descriptive): {int((T.loc[main.index].verdict == 'ADDS SKILL').sum())} of {n}; "
                f"q(all) < 0.05: {int((T.loc[main.index].q_all < 0.05).sum())} of {n}\n")
    print(open(f"{out}/skill.txt").read())
