"""Skill tests for the J2 / T2 groups (see PREREGISTRATION.md). Reuses ../intensity_extra/evaluate.py (LOSO, gains, sign-flip, BH).

usage: skill.py OUTDIR      (needs results/features.csv.gz from qc.py, and qc.txt must say PASS)
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "intensity_extra"))
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import numpy as np, pandas as pd
import evaluate as E
import model as M
from features import J2, T2, EXTRA

INT = os.path.join(HERE, "..", "intensity", "results")
CONT = ["vmaxp", "Lhalf", "s", "n", "dvds", "div300max", "div300dist", "tdepth", "tdist", "tbear_cos", "tbear_sin", "tamp", "ttilt", "tphase"]
VARIANTS = ["r1000", "r0", "x500", "x1000"]
BASE = M.SETS["full"]


def vcols(cols, suf=""):
    return [c + ("_" + suf if suf and c not in ("nohead",) else "") for c in cols]


def build():
    D = M.load(os.path.join(INT, "fixes_2004.csv.gz"), os.path.join(INT, "env_2004.csv.gz"))
    D = D[D.season <= 2021].reset_index(drop=True)
    X = pd.read_csv(os.path.join(HERE, "results", "features.csv.gz"), dtype={"time": str})
    X = X.drop(columns=["basin", "season", "lat", "lon", "heading"])
    X = X.merge(D[["track", "time", "basin", "month"]], on=["track", "time"])
    allc = [c for c in X.columns if c not in ("track", "time", "basin", "month")]
    # basin x month standardisation over every extracted fixes (predictors only)
    for c in allc:
        base = c.split("_")[0] if c.split("_")[0] in CONT else c
        base = c.rsplit("_", 1)[0] if (c.rsplit("_", 1)[-1] in VARIANTS) else c
        if base in CONT:
            g = X.groupby(["basin", "month"])[c]
            X[c] = (X[c] - g.transform("mean")) / (g.transform("std") + 1e-9)
    X = X.drop(columns=["basin", "month"])
    D = D.merge(X, on=["track", "time"], how="left", validate="one_to_one")
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    have = D[["track", "time"]].assign(_h=1)
    for L in (12, 24, 48):
        key = (t - pd.Timedelta(hours=L)).dt.strftime("%Y%m%d%H")
        Z = X.rename(columns={c: f"{c}_l{L}" for c in allc}).rename(columns={"time": "ktime"})
        K = pd.DataFrame({"track": D.track, "ktime": key})
        Zl = K.merge(Z, on=["track", "ktime"], how="left")
        for c in allc:
            D[f"{c}_l{L}"] = Zl[f"{c}_l{L}"].values
        D[f"miss{L}"] = Zl["vmax2500_l%d" % L].isna().astype(float).values
    D["dvmaxp_24"] = D.vmaxp - D.vmaxp_l24
    D["dtdepth_24"] = D.tdepth - D.tdepth_l24
    D["dtamp_24"] = D.tamp - D.tamp_l24
    # coupled indicator, thresholds outcome-free; tdepth/tdist were standardised above, so recompute from the raw file
    R = pd.read_csv(os.path.join(HERE, "results", "features.csv.gz"), dtype={"time": str}, usecols=["track", "time", "tdepth", "tdist", "RE", "LX"])
    q75 = R.tdepth.quantile(0.75)
    R["couple"] = (((R.RE + R.LX) > 0) & (R.tdist <= 2) & (R.tdepth > q75)).astype(float)
    D = D.merge(R[["track", "time", "couple"]], on=["track", "time"], how="left")
    D["hf24"] = D.hf24.astype(int)
    D["rd"] = np.where(D.cls >= 0, (D.cls == E.RD).astype(int), -1)
    D["hf48i"] = D.hf48.fillna(-1).astype(int)
    D = D[D.month.isin([10, 11, 12, 1, 2, 3, 4])].reset_index(drop=True)
    return D, q75


def group(cols0, lags, suf="", changes=()):
    out = list(vcols(cols0, suf))
    for L in lags:
        out += [f"{c}_l{L}" for c in vcols(cols0, suf)]
        out.append(f"miss{L}")
    out += list(changes)
    return list(dict.fromkeys(out))


def models_for(lagset, suf_j="", suf_t=""):
    lags = {"0": (), "0,24": (24,), "all": (12, 24, 48)}[lagset]
    chj = ("dvmaxp_24",) if 24 in lags else ()
    cht = ("dtdepth_24", "dtamp_24") if 24 in lags else ()
    J = group(J2, lags, suf_j, chj)
    T = group(T2, lags, suf_t, cht)
    return {"J2": BASE + J, "T2": BASE + T, "J2+T2": BASE + list(dict.fromkeys(J + T))}


def run(D, models, target, mask=None, tname=None, njobs=4):
    m = np.ones(len(D), bool) if mask is None else mask
    mods = {"base": BASE, **models}
    S, P, C = E.loso(D, mods, target, m)
    y = S[target].astype(int).values
    return S, P, C, y


def family(S, P, C, y, tname, names, mask=None):
    T, bss = E.gains(S, P, C, y, tname, names, mask=mask)
    T["q_family"] = E.bh(T.p.values)
    return T


def verdict(r):
    if r.gain >= E.THRESH and r.lo90 > 0 and r.q_family < 0.05:
        return "ADDS SKILL"
    if r.hi90 < E.THRESH:
        return "well-powered null (gain >= 0.005 excluded)"
    return "can't tell"


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    assert "PASS" in open(os.path.join(HERE, "results", "qc.txt")).read(), "QC did not pass"
    D, q75 = build()
    D.to_pickle(os.path.join(os.environ.get("SCRATCH", "/tmp"), "jt_D.pkl"))
    print(len(D), "fixes", D.track.nunique(), "tracks", D.season.nunique(), "seasons; hf24", int(D.hf24.sum()),
          "tdepth q75", round(q75, 1), flush=True)
    miss = {c: float(D[c].isna().mean()) for c in ("vmaxp", "tdepth", "ttilt", "tphase", "dvds", "vmaxp_l24", "tdepth_l24")}
    print("NaN share", miss, flush=True)
    allm, fams = {}, []
    sets = {"0": models_for("0"), "0,24": models_for("0,24"), "all": models_for("all")}
    mods = {}
    for ls, ms in sets.items():
        for k, v in ms.items():
            mods[f"{k} [lags {ls}]"] = v
    mods["couple"] = BASE + ["couple"]
    S, P, C, y = run(D, mods, "hf24")
    np.savez_compressed(f"{out}/hf24_probs.npz", track=S.track.values, time=S.time.values.astype(str), season=S.season.values, y=y,
                        **{k.replace(" ", "_"): v for k, v in P.items()})
    rows = []
    prim = family(S, P, C, y, "hf24 [primary]", [f"{k} [lags 0,24]" for k in ("J2", "T2", "J2+T2")])
    rows.append(prim)
    rows.append(family(S, P, C, y, "hf24 lags 0", [f"{k} [lags 0]" for k in ("J2", "T2", "J2+T2")]))
    rows.append(family(S, P, C, y, "hf24 lags 0,12,24,48", [f"{k} [lags all]" for k in ("J2", "T2", "J2+T2")]))
    on = ~S.hf_now.values.astype(bool)
    rows.append(family(S, P, C, y, "hf24 onset fixes", [f"{k} [lags 0,24]" for k in ("J2", "T2", "J2+T2")], mask=on))
    rows.append(family(S, P, C, y, "hf24 couple", ["couple"]))
    for b in ("atl", "pac"):
        rows.append(family(S, P, C, y, f"hf24 [primary] {b} (descriptive)", [f"{k} [lags 0,24]" for k in ("J2", "T2", "J2+T2")], mask=(S.basin == b).values))
    ms = {k: v for k, v in sets["0,24"].items()}
    for tg, tname, msk in (("rd", "rapid deepening", (D.rd >= 0).values), ("hf48i", "hf48", (D.hf48i >= 0).values)):
        S2, P2, C2, y2 = run(D, {f"{k} [lags 0,24]": v for k, v in ms.items()}, tg, msk)
        rows.append(family(S2, P2, C2, y2, tname, [f"{k} [lags 0,24]" for k in ms]))
        if tg == "rd":
            S2.to_pickle(os.path.join(os.environ.get("SCRATCH", "/tmp"), "jt_S_rd.pkl"))
    T = pd.concat(rows, ignore_index=True)
    T["verdict"] = [verdict(r) for r in T.itertuples()]
    T["q_all"] = np.nan
    main = T[~T.target.str.contains("descriptive")]
    T.loc[main.index, "q_all"] = E.bh(main.p.values)
    T.to_csv(f"{out}/gains.csv", index=False, float_format="%.5g")
    # yes/no scores for any passing group
    passing = T[(T.verdict == "ADDS SKILL") & (T.target == "hf24 [primary]")]
    hs = []
    for r in passing.itertuples():
        d, lo, hi, hb, hk = E.hss_gain(S, P, C, y, r.model)
        hs.append(dict(model=r.model, hss_base=hb, hss=hk, d_hss=d, lo90=lo, hi90=hi))
    pd.DataFrame(hs).to_csv(f"{out}/hss.csv", index=False, float_format="%.4g")
    with open(f"{out}/skill.txt", "w") as f:
        f.write(f"Sample: {len(D)} Oct-Apr fixes, {D.track.nunique()} tracks, 18 seasons 2004-05..2021-22; hf24 events {int(D.hf24.sum())}\n"
                f"Base (PR 12 'full', refit on this sample) BSS {T.bss_base.iloc[0]:+.4f}. Gain = BSS(base+group) - BSS(base), 90% CI over seasons.\n")
        for tg, g in T.groupby("target", sort=False):
            f.write(f"\n== {tg}  (base BSS {g.bss_base.iloc[0]:+.4f})\n")
            for r in g.itertuples():
                f.write(f"  {r.model:22s} gain {r.gain:+.4f} [{r.lo90:+.4f}, {r.hi90:+.4f}] BSS {r.bss:+.4f} p {r.p:.4f} q(fam) {r.q_family:.3f} "
                        f"q(all) {r.q_all:.3f} seasons better {r.seasons_better}/{r.n_seasons}  -> {r.verdict}\n")
        main = T[~T.target.str.contains("descriptive")]
        nmain = len(main)
        f.write(f"\nTests passing 'adds skill' (non-descriptive): {int((main.verdict == 'ADDS SKILL').sum())} of {nmain}; "
                f"q(all) < 0.05: {int((main.q_all < 0.05).sum())} of {nmain}\n")
    print(open(f"{out}/skill.txt").read())
