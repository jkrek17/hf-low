"""Score the 24 h gale / storm / HF outlook (models S, A, B of results/outlook_models.json, frozen) on GFS forecast fixes.
Plan: PREREGISTRATION.md and its Deviations. ERA5 pipeline A is a proxy and is the truth here; it is not observation.
Truth: gcls of the nearest pipeline A fix at the same valid time within 500 km (bomb_env/score.load(): peak g800 at +12 and +24 h on the
same track binned at 38/54/71.7 kt); no match = dropped and counted. Score: RPSS vs basin x month climatology of the training seasons.
Season-block bootstrap (5 GFS seasons 2021-22..2025-26), 2,000 draws, seed 20261013, 90 % intervals, BH over G1-G5.
usage: score_outlook.py [OUTDIR] [ENVTABLE]       writes OUTDIR/outlook.txt and OUTDIR/outlook.json
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gfs_test"))
sys.path.insert(0, os.path.join(HERE, "..", "bomb_env"))
import a_score as AS                    # gfs_features, match_era5, hav
import score as BS                      # load, rps, bh, NAMES

RES = os.path.join(HERE, "results")
OUT = sys.argv[1] if len(sys.argv) > 1 else RES
ENVT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "gfs_outlook_env.csv.gz")
SEEDV = 20261013
NB = 2000
SEASONS = [2021, 2022, 2023, 2024, 2025]
LEADS = [0, 12, 24, 36, 48]
LINES = []


def out(s=""):
    print(s)
    LINES.append(s)


# ------------------------------------------------------------------ models
def softmax_predict(m, X):
    X = np.asarray(X, float)
    med, lo, hi, mu, sd = (np.array(m["prep"][k]) for k in ("med", "lo", "hi", "mu", "sd"))
    Z = (np.clip(np.where(np.isnan(X), med, X), lo, hi) - mu) / sd
    L = Z @ np.array(m["coef"]).T + np.array(m["intercept"])
    return L


def sm(L):
    L = L - L.max(1, keepdims=True)
    e = np.exp(L)
    return e / e.sum(1, keepdims=True)


def run_model(models, key, X, use_dp):
    """X: DataFrame with all columns; rows with use_dp take model key, others key_nodp12. Returns logits (n,4)."""
    L = np.zeros((len(X), 4))
    for flag, name in ((True, key), (False, key + "_nodp12")):
        mm = use_dp == flag
        if mm.any():
            m = models["models"][name]
            L[mm] = softmax_predict(m, X.loc[mm, m["predictors"]].values)
    return L


# ------------------------------------------------------------------ bootstrap
def draws(k):
    return np.random.default_rng(SEEDV).integers(0, k, (NB, k))


def sums(v, seas, useas):
    return np.array([v[seas == s].sum() for s in useas])


def stat(num_losses, den_loss, seas, useas, B, sign=1):
    """skill-type statistic sum_i sign_i*loss_i / sum(den): list of (sign, loss) combined. Returns point, lo, hi, boot."""
    S = sum(sg * sums(l, seas, useas) for sg, l in num_losses)
    D = sums(den_loss, seas, useas)
    pt = S.sum() / D.sum()
    bt = S[B].sum(1) / D[B].sum(1)
    return pt, np.percentile(bt, 5), np.percentile(bt, 95), bt


def rpss(l_m, l_c, seas, useas, B):
    return stat([(-1, l_m), (1, l_c)], l_c, seas, useas, B)       # 1 - sum m / sum c  = (sum c - sum m)/sum c


def gain(l_new, l_old, l_c, seas, useas, B):
    return stat([(-1, l_new), (1, l_old)], l_c, seas, useas, B)


def p_le(bt, thr=0.0):
    return float((1 + np.sum(bt <= thr)) / (1 + len(bt)))


def fmt(t):
    return f"{t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]"


def rel_table(p, y, nb=10):
    q = np.unique(np.quantile(p, np.linspace(0, 1, nb + 1)))
    b = np.clip(np.searchsorted(q[1:-1], p, side="right"), 0, len(q) - 2)
    rows = [dict(bin=int(i), n=int((b == i).sum()), forecast=float(p[b == i].mean()), observed=float(y[b == i].mean()))
            for i in range(len(q) - 1) if (b == i).sum()]
    r = pd.DataFrame(rows)
    slope = float(np.polyfit(r.forecast, r.observed, 1, w=np.sqrt(r.n))[0]) if len(r) > 1 else float("nan")
    return rows, slope


# ------------------------------------------------------------------ main
def main():
    models = json.load(open(os.path.join(RES, "outlook_models.json")))
    clim = {(r["basin"], r["month"]): np.array(r["probs_smoothed"]) for r in models["climatology_gcls"]}
    F, n_all = AS.gfs_features()                  # atl/pac, valid 00/12 UTC, leads 0..48 step 12, dp12 and speed along the GFS track
    n_track_rows = n_all
    F = F[F.lead.isin(LEADS)].copy().reset_index(drop=True)
    D = BS.load()
    E = D[D.season.isin(SEASONS)].copy().reset_index(drop=True)
    E["time_i"] = E.time.astype(int)
    mi, md = AS.match_era5(F, E)
    F["era_i"] = mi
    F["matched"] = mi >= 0
    F["y"] = np.where(F.matched, E.gcls.reindex(np.where(mi >= 0, mi, 0)).values, -1).astype(int)
    out("GFS outlook test (S, A, B fed GFS forecast fixes). Truth = ERA5 pipeline A, a PROXY (not observation).")
    out(f"GFS track rows (all, 2021-22..2025-26): {n_track_rows}; atl/pac fixes at valid 00/12 UTC, leads 0/12/24/36/48: {len(F)}")
    out("Fixes by lead: " + ", ".join(f"{l}: {int((F.lead == l).sum())}" for l in LEADS))
    out("Fixes by season: " + ", ".join(f"{s}: {int((F.season == s).sum())}" for s in SEASONS))
    out(f"ERA5 fixes in seasons 2021-2025 with gcls: {len(E)}")
    nm = F[~F.matched]
    out("No pipeline A fix within 500 km at the same valid time (dropped): " + ", ".join(
        f"lead {l}: {int(((F.lead == l) & ~F.matched).sum())}" for l in LEADS) + f"; total {len(nm)} of {len(F)}")
    Fm = F[F.matched].copy().reset_index(drop=True)
    Fm["basin_era"] = E.basin.reindex(Fm.era_i).values
    out(f"Matched: {len(Fm)}; matched ERA5 fix in the other basin: {int((Fm.basin_era != Fm.basin).sum())}; "
        f"match distance median {np.nanmedian(md[mi >= 0]):.0f} km")
    out("Matched truth class counts (below gale, gale, storm, HF) by lead: " + "; ".join(
        f"{l}: {np.bincount(Fm.y[Fm.lead == l], minlength=4).tolist()}" for l in LEADS))
    out("Matched ERA5 fixes whose track ends within 24 h (left24, class 0 by convention): " + str(int(E.left24.reindex(Fm.era_i).sum())))
    # environment table
    T = pd.read_csv(ENVT, dtype={"cycle": str})
    out(f"Environment table rows: {len(T)} (leads {sorted(T.lead.unique().tolist())})")
    Fm["cycle"] = Fm.cycle.astype(str)
    Fm = Fm.merge(T[["cycle", "lead", "track", "jet250", "vadv500", "eady"]], on=["cycle", "lead", "track"], how="left")
    Fm["env_ok"] = Fm[["jet250", "vadv500", "eady"]].notna().all(axis=1)
    for l in (0, 24):
        g = Fm[Fm.lead == l]
        out(f"Lead {l}: matched {len(g)}; with all three environment terms {int(g.env_ok.sum())}; no env row {int(g.jet250.isna().sum())}; "
            f"eady NaN {int(g.eady.isna().sum())}")
    # models
    Fm["use_dp"] = ((Fm.lead >= 12) & Fm.dp12.notna()).values
    out("Share using the no-dp12 fallback, lead >= 12: " + ", ".join(
        f"{l}: {float((~Fm.use_dp[Fm.lead == l]).mean()):.3f}" for l in (12, 24, 36, 48)) + "; speed missing: " + ", ".join(
        f"{l}: {float(Fm.speed[Fm.lead == l].isna().mean()):.3f}" for l in LEADS))
    Fm["doy_c"] = Fm.doy_c; Fm["doy_s"] = Fm.doy_s
    C = np.array([clim[(b, m)] for b, m in zip(Fm.basin, Fm.month)])
    P = {"clim": C}
    for k in ("S", "A", "B"):
        P[k] = sm(run_model(models, k, Fm, Fm.use_dp.values))
        if k != "S":                                              # A, B undefined without env terms: NaN rows
            P[k][~Fm.env_ok.values] = np.nan
    # ERA5-fed S on the same matched fixes
    Xe = pd.DataFrame(dict(msl=E.msl.reindex(Fm.era_i).values, dp12=E.dp12.reindex(Fm.era_i).values, lat=E.lat.reindex(Fm.era_i).values,
                           speed=E.speed.reindex(Fm.era_i).values, g800=E.g800.reindex(Fm.era_i).values))
    young_e = E.young.reindex(Fm.era_i).values == 1
    Xe.loc[young_e, "dp12"] = np.nan
    for c in ("pac", "doy_c", "doy_s"):
        Xe[c] = Fm[c].values
    use_e = Fm.use_dp.values & ~young_e
    P["S_era"] = sm(run_model(models, "S", Xe, use_e))
    out(f"ERA5-fed S: dp12 model used for {int(use_e[Fm.lead.values >= 12].sum())} of {int((Fm.lead >= 12).sum())} lead >= 12 fixes")
    y = Fm.y.values
    seas = Fm.season.values
    L = {k: BS.rps(P[k], y, 4) for k in P}
    B5 = draws(5)
    res = {}
    # ---- family
    def sel(lead, mask=None):
        m = (Fm.lead == lead).values
        return m if mask is None else m & mask
    fam = []
    m24 = sel(24)
    g1 = rpss(L["S"][m24], L["clim"][m24], seas[m24], SEASONS, B5)
    fam.append(("G1", "S fed GFS, lead 24, RPSS vs clim", g1, p_le(g1[3]), int(m24.sum())))
    d = stat([(1, L["S_era"][m24]), (-1, L["S"][m24])], L["clim"][m24], seas[m24], SEASONS, B5)   # GFS-fed minus ERA5-fed RPSS
    p2 = max(p_le(d[3], -0.05), float((1 + np.sum(d[3] >= 0.05)) / (1 + len(d[3]))))   # TOST: P(diff <= -0.05), P(diff >= +0.05)
    fam.append(("G2", "S lead 24, GFS-fed minus ERA5-fed RPSS (equivalence +-0.05)", d, p2, int(m24.sum())))
    mE = sel(24, Fm.env_ok.values)
    g3 = gain(L["A"][mE], L["S"][mE], L["clim"][mE], seas[mE], SEASONS, B5)
    g4 = gain(L["B"][mE], L["S"][mE], L["clim"][mE], seas[mE], SEASONS, B5)
    fam.append(("G3", "A minus S, GFS-fed, lead 24, RPSS gain", g3, p_le(g3[3]), int(mE.sum())))
    fam.append(("G4", "B minus S, GFS-fed, lead 24, RPSS gain", g4, p_le(g4[3]), int(mE.sum())))
    m0, m48 = sel(0), sel(48)
    r0 = rpss(L["S"][m0], L["clim"][m0], seas[m0], SEASONS, B5)
    r24 = g1
    r48 = rpss(L["S"][m48], L["clim"][m48], seas[m48], SEASONS, B5)
    dd = r0[3] - r48[3]
    g5 = (r0[0] - r48[0], np.percentile(dd, 5), np.percentile(dd, 95), dd)
    fam.append(("G5", "S RPSS lead 0 minus lead 48 (paired)", g5, p_le(dd), int(m0.sum()) + int(m48.sum())))
    q = BS.bh([f[3] for f in fam])
    out("\n== Registered family G1-G5 (RPSS; 90 % season-block interval, 5 seasons, 2,000 draws; BH q within the 5)")
    rows = []
    for (tid, lab, t, p, n), qq in zip(fam, q):
        if tid == "G2":
            v = "equivalent" if (t[1] > -0.05 and t[2] < 0.05) else ("degrades" if t[2] < -0.05 else "cannot tell")
        elif tid == "G5":
            ordered = r0[0] >= r24[0] >= r48[0]
            v = ("ordered" if (ordered and t[1] > 0 and qq < 0.05) else ("flat" if t[1] <= 0 <= t[2] else "reversed"))
        else:
            v = "yes" if (qq < 0.05 and t[1] > 0) else "cannot tell"
        out(f"  {tid} {lab}: {fmt(t)}  n {n}  p {p:.4f}  q {qq:.4f}  -> {v}")
        rows.append(dict(test=tid, label=lab, point=float(t[0]), lo90=float(t[1]), hi90=float(t[2]), n=n, p=p, q=float(qq), verdict=v))
    res["family"] = rows
    out(f"  S RPSS by lead (G5 detail): 0 {fmt(r0)}; 24 {fmt(r24)}; 48 {fmt(r48)}")
    res["G5_leads"] = {"0": [float(x) for x in r0[:3]], "24": [float(x) for x in r24[:3]], "48": [float(x) for x in r48[:3]]}
    s_era = rpss(L["S_era"][m24], L["clim"][m24], seas[m24], SEASONS, B5)
    out(f"  Lead 24: S ERA5-fed RPSS {fmt(s_era)}; S GFS-fed RPSS {fmt(g1)} (same {int(m24.sum())} matched fixes)")
    res["S_era_lead24"] = [float(x) for x in s_era[:3]]
    # ---- secondary
    out("\n== Secondary (labelled; not in the family; no multiplicity control)")
    res["secondary"] = {}
    sec = res["secondary"]
    sec["S_by_lead"] = {}
    for l in LEADS:
        m = sel(l)
        t = rpss(L["S"][m], L["clim"][m], seas[m], SEASONS, B5)
        sec["S_by_lead"][str(l)] = [float(x) for x in t[:3]] + [int(m.sum())]
        te = rpss(L["S_era"][m], L["clim"][m], seas[m], SEASONS, B5)
        out(f"  S RPSS lead {l:2d}: GFS-fed {fmt(t)}  ERA5-fed {fmt(te)}  n {int(m.sum())}")
    sec["S_by_basin_lead24"], sec["AB_by_basin_lead24"] = {}, {}
    for b in ("atl", "pac"):
        mb = m24 & (Fm.basin == b).values
        t = rpss(L["S"][mb], L["clim"][mb], seas[mb], SEASONS, B5)
        sec["S_by_basin_lead24"][b] = [float(x) for x in t[:3]] + [int(mb.sum())]
        out(f"  S RPSS lead 24 {b}: {fmt(t)}  n {int(mb.sum())}")
        mbe = mE & (Fm.basin == b).values
        for k in ("A", "B"):
            t = gain(L[k][mbe], L["S"][mbe], L["clim"][mbe], seas[mbe], SEASONS, B5)
            sec["AB_by_basin_lead24"][f"{k}_{b}"] = [float(x) for x in t[:3]] + [int(mbe.sum())]
            out(f"  {k} minus S lead 24 {b}: {fmt(t)}  n {int(mbe.sum())}")
    m0e = sel(0, Fm.env_ok.values)
    sec["AB_lead0"] = {}
    for k in ("A", "B"):
        t = gain(L[k][m0e], L["S"][m0e], L["clim"][m0e], seas[m0e], SEASONS, B5)
        sec["AB_lead0"][k] = [float(x) for x in t[:3]] + [int(m0e.sum())]
        out(f"  {k} minus S lead 0: {fmt(t)}  n {int(m0e.sum())}")
    # recalibrated outlook: class offsets fitted on 2021 (2021-22) lead 24, applied 2022..2025
    mfit = m24 & (seas == 2021)
    mapp = m24 & (seas >= 2022)
    LG = run_model(models, "S", Fm, Fm.use_dp.values)
    yf = y[mfit]

    def nll(b):
        bb = np.r_[0.0, b]
        Z = LG[mfit] + bb
        Z = Z - Z.max(1, keepdims=True)
        return -np.sum(Z[np.arange(len(yf)), yf] - np.log(np.exp(Z).sum(1)))
    r = minimize(nll, np.zeros(3), method="BFGS")
    bb = np.r_[0.0, r.x]
    Prc = sm(LG + bb)
    Lrc = BS.rps(Prc, y, 4)
    B4 = draws(4)
    S4 = [2022, 2023, 2024, 2025]
    t = gain(Lrc[mapp], L["S"][mapp], L["clim"][mapp], seas[mapp], S4, B4)
    t_un = rpss(L["S"][mapp], L["clim"][mapp], seas[mapp], S4, B4)
    t_rc = rpss(Lrc[mapp], L["clim"][mapp], seas[mapp], S4, B4)
    out(f"  Recalibrated S (offsets fitted on 2021-22 lead 24, n {int(mfit.sum())}: {np.round(bb, 3).tolist()}), applied 2022-2025 lead 24 n {int(mapp.sum())}: "
        f"RPSS uncal {fmt(t_un)}, recal {fmt(t_rc)}, gain {fmt(t)} (4 seasons)")
    sec["recalibration_lead24"] = dict(offsets=bb.tolist(), n_fit=int(mfit.sum()), n_apply=int(mapp.sum()),
                                       rpss_uncal=[float(x) for x in t_un[:3]], rpss_recal=[float(x) for x in t_rc[:3]], gain=[float(x) for x in t[:3]])
    # reliability of HF probability, lead 24, S
    rows_, slope = rel_table(P["S"][m24][:, 3], (y[m24] == 3).astype(float))
    out(f"  Reliability of the GFS-fed S hurricane-force probability, lead 24 (deciles of forecast), slope {slope:.2f} (1 = perfect), base rate {float((y[m24] == 3).mean()):.4f}, mean forecast {float(P['S'][m24][:, 3].mean()):.4f}:")
    for rw in rows_:
        out(f"    bin {rw['bin']}: n {rw['n']:5d}  forecast {rw['forecast']:.4f}  observed {rw['observed']:.4f}")
    sec["reliability_HF_lead24"] = dict(slope=slope, bins=rows_)
    # class frequencies: forecast vs observed
    out("  Lead 24 mean forecast class probability vs observed frequency (S GFS-fed): " + "; ".join(
        f"{BS.NAMES[c]} {P['S'][m24][:, c].mean():.3f} vs {float((y[m24] == c).mean()):.3f}" for c in range(4)))
    res["counts"] = dict(track_rows=int(n_track_rows), fixes_scored_domain=int(len(F)), matched=int(len(Fm)), unmatched=int(len(nm)),
                         by_lead={str(l): dict(domain=int((F.lead == l).sum()), matched=int((Fm.lead == l).sum()),
                                               env_ok=int(((Fm.lead == l) & Fm.env_ok).sum())) for l in LEADS},
                         seed=SEEDV, nboot=NB)
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, "outlook.txt"), "w").write("\n".join(LINES) + "\n")
    json.dump(res, open(os.path.join(OUT, "outlook.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
