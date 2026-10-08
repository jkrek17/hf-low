"""Stage A scoring (A1.1-A1.3, A2.1-A2.4) of PREREGISTRATION.md exactly as written plus Deviations 8 and 9.
Pipeline A (ERA5) is a proxy. GFS seasons 2021-2025 are inside the gust-threshold calibration block (2021-22..2025-26).
Frozen model: results/gfs_model.json (never refitted). Truth A1: OPC archive, category HF only. Truth A2 primary: pipeline A hf24
of the nearest ERA5 fix at the same valid time within 500 km (else 0); secondary: archive HF fix within 600 km in (t, t+24 h].
Season-block bootstrap, 2,000 resamples, seed 20261008, 95 % percentile intervals. One-sided p floor 1/(1+2000).
usage: a_score.py [--nsim N]   (writes results/stageA.txt and results/stageA.json)
"""
import os, sys, json, time
import numpy as np, pandas as pd
from scipy.stats import rankdata
from scipy.special import expit, logit
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import model as M   # intensity/model.py: predict()

SEED = 20261008
NBOOT = 2000
SEASONS = [2021, 2022, 2023, 2024, 2025]
LEADS_A1 = [0, 24, 48]
LEADS_A2 = [0, 12, 24, 36, 48]
BASINS = ["atl", "pac"]
HF_KT = 71.7


# ---------------------------------------------------------------- generic scoring helpers
def brier(p, y):
    return (p - y) ** 2


def auc(s, y):
    y = np.asarray(y).astype(bool)
    n1 = y.sum(); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def sums_by_season(v, seas, useas):
    return np.array([np.sum(v[seas == s]) for s in useas], dtype=float)


def draws(k, nboot=NBOOT, seed=SEED):
    return np.random.default_rng(seed).integers(0, k, (nboot, k))


def pct(x, lo=2.5, hi=97.5):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return [float(np.percentile(x, lo)), float(np.percentile(x, hi))] if len(x) else [np.nan, np.nan]


def p_le(boot, thr=0.0):
    """one-sided bootstrap p for H0: stat <= thr (floor 1/(1+B))."""
    b = np.asarray(boot, float); b = b[~np.isnan(b)]
    return float((1 + np.sum(b <= thr)) / (1 + len(b)))


def bss_from_sums(sm, sr, B):
    return 1 - sm[B].sum(1) / sr[B].sum(1)


def bss_stat(pm, pr, y, seas, useas, B):
    """BSS of pm vs reference pr, point estimate, 95 % interval, one-sided p (BSS<=0), bootstrap vector."""
    sm = sums_by_season(brier(pm, y), seas, useas); sr = sums_by_season(brier(pr, y), seas, useas)
    obs = 1 - sm.sum() / sr.sum()
    bs = bss_from_sums(sm, sr, B)
    return dict(bss=float(obs), ci=pct(bs), p=p_le(bs, 0.0)), bs


def auc_boot(s, y, seas, useas, B):
    idx = [np.where(seas == u)[0] for u in useas]
    out = np.empty(len(B))
    for i, row in enumerate(B):
        ii = np.concatenate([idx[j] for j in row])
        out[i] = auc(s[ii], y[ii])
    return out


def contingency(fc, ob):
    fc = np.asarray(fc, bool); ob = np.asarray(ob, bool)
    return np.array([np.sum(fc & ob), np.sum(fc & ~ob), np.sum(~fc & ob), np.sum(~fc & ~ob)], float)


def cat_scores(abcd):
    a, b, c, d = abcd
    n = a + b + c + d
    ex = ((a + b) * (a + c) + (c + d) * (b + d)) / n if n else np.nan
    den = n - ex
    return dict(pod=float(a / (a + c)) if a + c else np.nan, far=float(b / (a + b)) if a + b else np.nan,
                bias=float((a + b) / (a + c)) if a + c else np.nan,
                hss=float((a + d - ex) / den) if n and den > 0 else np.nan, hits=int(a), fa=int(b), miss=int(c), cn=int(d))


def hss_boot(fc, ob, seas, useas, B):
    cs = np.array([contingency(fc[seas == s], ob[seas == s]) for s in useas])   # (k,4)
    out = np.empty(len(B))
    for i, row in enumerate(B):
        a, b, c, d = cs[row].sum(0)
        n = a + b + c + d
        ex = ((a + b) * (a + c) + (c + d) * (b + d)) / n
        out[i] = (a + d - ex) / (n - ex) if n - ex > 0 else np.nan
    return out


def best_hss_threshold(p, y):
    """probability threshold (forecast yes if p >= thr) maximising HSS; candidates are the observed values."""
    o = np.argsort(-p, kind="stable")
    ps, ys = p[o], y[o].astype(float)
    ca, cb = np.cumsum(ys), np.cumsum(1 - ys)
    last = np.r_[ps[1:] != ps[:-1], True]            # last index of each tied value
    a, b = ca[last], cb[last]
    P, N = ys.sum(), len(ys); c = P - a; d = N - P - b
    ex = ((a + b) * (a + c) + (c + d) * (b + d)) / N
    h = (a + d - ex) / (N - ex)
    if not np.isfinite(h).any():
        return np.nan
    k = np.nanargmax(h)
    return float(ps[last][k])


# ---------------------------------------------------------------- archive
BOX = {"atl": lambda la, lo: (la >= 30) & (la <= 67) & ((lo >= 262) | (lo <= 10)),
       "pac": lambda la, lo: (la >= 27) & (la <= 67) & (lo >= 135) & (lo <= 240)}


def load_archive():
    root = os.path.join(HERE, "..", "..", "..", "data", "hf_lows")
    parts, cats, bad = [], {}, []
    for tag in ("Atl", "Pac"):
        d = pd.read_csv(os.path.join(root, f"HF_Data_-_{tag}.csv"), dtype={"date": str})
        cats[tag] = {str(k): int(v) for k, v in d.Category.fillna("(nan)").value_counts().items()}
        ok = d.date.str.fullmatch(r"\d{10}", na=False) & pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce").notna()
        for _, r in d[~ok].iterrows():
            bad.append(dict(file=tag, **{k: (None if pd.isna(v) else v) for k, v in r.items()}))
        d = d[ok].copy()
        d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H")
        d["file"] = tag
        parts.append(d)
    A = pd.concat(parts, ignore_index=True)
    A["lon"] = A.Longitude % 360
    A["lat"] = A.Latitude
    return A, cats, bad


def hf_fix_times(A, basin):
    h = A[(A.Category == "HF")]
    m = BOX[basin](h.lat.values, h.lon.values)
    return np.sort(h.t.values[m].astype("datetime64[h]").astype(np.int64))      # hours since epoch


def label_within(times_h, valid):
    v = pd.to_datetime(valid.astype(str), format="%Y%m%d%H").values.astype("datetime64[h]").astype(np.int64)
    lo = np.searchsorted(times_h, v - 3, side="left"); hi = np.searchsorted(times_h, v + 3, side="right")
    return (hi > lo).astype(int)


# ---------------------------------------------------------------- Stage A1
def doy_terms(valid):
    doy = pd.to_datetime(pd.Series(valid).astype(str), format="%Y%m%d%H").dt.dayofyear.values
    return np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)


def _lr(X, y):
    return LogisticRegression(C=np.inf, max_iter=2000).fit(X, y)


def loso_probs(g, c, s, y, seas):
    """leave-one-season-out probabilities: G model (G standardised on training fold + cos, sin) and month-only model."""
    pg = np.zeros(len(y)); pm = np.zeros(len(y))
    for u in np.unique(seas):
        te = seas == u; tr = ~te
        mu, sd = g[tr].mean(), g[tr].std() + 1e-9
        Xg = np.c_[(g - mu) / sd, c, s]; Xm = np.c_[c, s]
        pg[te] = _lr(Xg[tr], y[tr]).predict_proba(Xg[te])[:, 1]
        pm[te] = _lr(Xm[tr], y[tr]).predict_proba(Xm[te])[:, 1]
    return pg, pm


def planted_power(g, c, s, y, seas, B, nsim, seed, target_bss=0.05):
    """Labels from the month-only rates plus a G effect with expected BSS target_bss; probability the A1.1 criterion passes."""
    Xm = np.c_[c, s]
    p0 = _lr(Xm, y).predict_proba(Xm)[:, 1]
    z = (g - g.mean()) / (g.std() + 1e-9)
    l0 = logit(np.clip(p0, 1e-6, 1 - 1e-6))

    def ebss(beta):
        p1 = expit(l0 + beta * z)
        return 1 - np.sum(p1 * (1 - p1)) / np.sum(p1 * (1 - p1) + (p0 - p1) ** 2)
    lo, hi = 0.0, 6.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if ebss(mid) < target_bss else (lo, mid)
    beta = (lo + hi) / 2
    p1 = expit(l0 + beta * z)
    rng = np.random.default_rng(seed)
    useas = np.unique(seas)
    passes = 0
    for _ in range(nsim):
        ys = (rng.random(len(p1)) < p1).astype(int)
        if ys.sum() < 2 or ys.sum() > len(ys) - 2:
            continue
        pg, pm = loso_probs(g, c, s, ys, seas)
        sg = sums_by_season(brier(pg, ys), seas, useas); sm = sums_by_season(brier(pm, ys), seas, useas)
        obs = 1 - sg.sum() / sm.sum()
        lb = np.percentile(bss_from_sums(sg, sm, B), 2.5)
        passes += int(obs > 0 and lb > 0)
    return dict(beta=float(beta), expected_bss=float(ebss(beta)), nsim=nsim, power=float(passes / nsim))


def stage_A1(nsim, log):
    G = pd.read_csv(os.path.join(RES, "gfs_basin_gust.csv.gz"))
    A, cats, bad = load_archive()
    out = {"archive_categories": cats, "archive_malformed_rows": bad}
    cyc_season = pd.to_datetime(G.cycle.astype(str), format="%Y%m%d%H")
    cs = np.where(cyc_season.dt.month >= 6, cyc_season.dt.year, cyc_season.dt.year - 1)
    out["season_col_equals_cycle_season_share"] = float(np.mean(cs == G.season.values))
    B5, B4 = draws(5), draws(4)
    res = {}
    # sensitivity: malformed Pac row (20241101018) read as 2024110118
    A_fix = A.copy()
    for b in bad:
        if b["date"] == "20241101018" and b["file"] == "Pac":
            row = dict(ID=b["ID"], date="2024110118", Latitude=float(b["Latitude"]), Longitude=float(b["Longitude"]),
                       Category=b["Category"], Pressure=b["Pressure"], t=pd.Timestamp("2024-11-01 18:00"), file="Pac")
            row["lon"] = row["Longitude"] % 360; row["lat"] = row["Latitude"]
            A_fix = pd.concat([A_fix, pd.DataFrame([row])], ignore_index=True)
    pw = {}
    for basin in BASINS:
        th = hf_fix_times(A, basin); th2 = hf_fix_times(A_fix, basin)
        d0 = G[(G.basin == basin) & G.lead.isin(LEADS_A1)].copy()
        d0["y"] = label_within(th, d0.valid_time.values)
        d0["y_fix"] = label_within(th2, d0.valid_time.values)
        out.setdefault("label_changed_by_malformed_row_fix", {})[basin] = int((d0.y != d0.y_fix).sum())
        d0["season"] = d0.season.astype(int)
        per = {}
        aucs_boot = {}
        for lead in LEADS_A1:
            d = d0[d0.lead == lead].sort_values("cycle").reset_index(drop=True)
            c, s = doy_terms(d.valid_time.values)
            g, y, seas = d.gmax_kt.values, d.y.values, d.season.values
            pg, pm = loso_probs(g, c, s, y, seas)
            st, bs = bss_stat(pg, pm, y, seas, SEASONS, B5)
            ab = auc_boot(pg, y, seas, SEASONS, B5)
            aucs_boot[lead] = ab
            per[lead] = dict(n=int(len(d)), label_freq=float(y.mean()), n_pos=int(y.sum()),
                             label_freq_by_season={int(u): float(y[seas == u].mean()) for u in SEASONS},
                             bss=st["bss"], bss_ci=st["ci"], bss_p_le0=st["p"],
                             auc=float(auc(pg, y)), auc_ci=pct(ab), auc_rawG=float(auc(g, y)),
                             auc_nan_boot=int(np.isnan(ab).sum()))
            if lead == 24:
                per[lead]["power_planted_bss0.05"] = planted_power(g, c, s, y, seas, B5, nsim, SEED + (1 if basin == "pac" else 0))
            per[lead]["_arrays"] = (d, g, y, seas)
        # A1.2
        a0, a24, a48 = per[0]["auc"], per[24]["auc"], per[48]["auc"]
        diff = aucs_boot[0] - aucs_boot[48]
        d_ci = pct(diff); d_p = p_le(diff, 0.0)
        ordered = bool(a0 >= a24 >= a48)
        a12 = dict(auc0=a0, auc24=a24, auc48=a48, ordered_point=ordered, diff_0_48=float(a0 - a48),
                   diff_ci=d_ci, p_le0=d_p, passes=bool(ordered and d_ci[0] > 0),
                   diff_0_24_ci=pct(aucs_boot[0] - aucs_boot[24]), diff_24_48_ci=pct(aucs_boot[24] - aucs_boot[48]))
        # A1.3
        d0l, g0, y0, s0 = per[0]["_arrays"]
        m21 = s0 == 2021
        freq21 = float(y0[m21].mean())
        gs = np.sort(g0[m21])[::-1]
        k = int(round(freq21 * m21.sum()))
        thr = float(gs[k]) if k < len(gs) else -np.inf
        d24, g24, y24, s24 = per[24]["_arrays"]
        te = s24 != 2021
        fc = g24[te] > thr; ob = y24[te].astype(bool)
        cat = cat_scores(contingency(fc, ob))
        hb = hss_boot(fc, ob, s24[te], [2022, 2023, 2024, 2025], B4)
        a13 = dict(label_freq_2021_lead0=freq21, n_cycles_2021_lead0=int(m21.sum()), k_above=k, threshold_kt=thr,
                   share_above_2021_lead0=float(np.mean(g0[m21] > thr)), n=int(te.sum()), n_pos=int(ob.sum()), **cat,
                   hss_ci=pct(hb), p_hss_lt_0p30=float((1 + np.sum(hb < 0.30)) / (1 + np.sum(~np.isnan(hb)))),
                   passes=bool(cat["hss"] >= 0.30 and pct(hb)[0] >= 0.30))
        for lead in LEADS_A1:
            per[lead].pop("_arrays")
        res[basin] = dict(per_lead=per, A1_1=dict(bss=per[24]["bss"], ci=per[24]["bss_ci"], p=per[24]["bss_p_le0"],
                                                  passes=bool(per[24]["bss"] > 0 and per[24]["bss_ci"][0] > 0)),
                          A1_2=a12, A1_3=a13)
    out["basins"] = res
    return out


# ---------------------------------------------------------------- Stage A2
def hav(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def gfs_features():
    T = pd.concat([pd.read_csv(os.path.join(RES, f"gfs_tracks_{s}.csv.gz")).assign(season=s) for s in SEASONS], ignore_index=True)
    T["lead"] = T.lead.astype(int)
    key = ["cycle", "track"]
    P12 = T[key + ["lead", "msl"]].assign(lead=lambda x: x.lead + 12).rename(columns={"msl": "msl_m12"})
    P6 = T[key + ["lead", "lat", "lon"]].assign(lead=lambda x: x.lead + 6).rename(columns={"lat": "lat_m6", "lon": "lon_m6"})
    T = T.merge(P12, on=key + ["lead"], how="left").merge(P6, on=key + ["lead"], how="left")
    T["dp12"] = T.msl - T.msl_m12
    dlat = T.lat - T.lat_m6
    dlon = ((T.lon - T.lon_m6 + 180) % 360 - 180) * np.cos(np.radians(T.lat))
    T["speed"] = np.hypot(dlat, dlon) * 60 / 6            # intensity/fixes.py: deg -> nmi per h = kt, dh = 6
    n_all = len(T)
    F = T[T.basin.notna() & T.basin.isin(BASINS) & (T.valid_time % 100).isin([0, 12]) & T.lead.isin(LEADS_A2)].copy()
    F["pac"] = (F.basin == "pac").astype(float)
    F["doy_c"], F["doy_s"] = doy_terms(F.valid_time.values)
    F["month"] = (F.valid_time // 10000 % 100).astype(int)
    return F.reset_index(drop=True), n_all


def match_era5(F, E, rad=500.0):
    """nearest ERA5 fix at the same valid time within rad km; returns index into E (-1 none) and distance."""
    E = E.reset_index(drop=True)
    eg = {t: g for t, g in E.groupby("time_i")}
    mi = np.full(len(F), -1); md = np.full(len(F), np.nan)
    for t, g in F.groupby("valid_time"):
        e = eg.get(t)
        if e is None:
            continue
        D = hav(g.lat.values[:, None], g.lon.values[:, None], e.lat.values[None, :], e.lon.values[None, :])
        j = D.argmin(1); dm = D[np.arange(len(g)), j]
        ok = dm <= rad
        pos = F.index.get_indexer(g.index)
        mi[pos[ok]] = e.index.values[j[ok]]
        md[pos] = dm
    return mi, md


def secondary_truth(F, A, rad=600.0):
    H = A[A.Category == "HF"].copy()
    H["th"] = H.t.values.astype("datetime64[h]").astype(np.int64)
    H = H.sort_values("th")
    th = H.th.values
    v = pd.to_datetime(F.valid_time.astype(str), format="%Y%m%d%H").values.astype("datetime64[h]").astype(np.int64)
    out = np.zeros(len(F), int)
    lo = np.searchsorted(th, v + 1, side="left"); hi = np.searchsorted(th, v + 24, side="right")
    la, lo_ = H.lat.values, H.lon.values
    for i in range(len(F)):
        if hi[i] > lo[i]:
            d = hav(F.lat.values[i], F.lon.values[i], la[lo[i]:hi[i]], lo_[lo[i]:hi[i]])
            out[i] = int((d <= rad).any())
    return out


def apply_model(F, X, models):
    """X: DataFrame with msl, dp12, lat, speed, g800, pac, doy_c, doy_s. dp12 model if lead>=12 and F.dp12 present."""
    use_dp = (F.lead.values >= 12) & F.dp12.notna().values
    p = np.zeros(len(F))
    if use_dp.any():
        p[use_dp] = M.predict(models["model_dp12"], X[use_dp])
    if (~use_dp).any():
        p[~use_dp] = M.predict(models["model_nodp12"], X[~use_dp])
    return p, use_dp


def reliability_table(p, y, nb=10):
    q = np.quantile(p, np.linspace(0, 1, nb + 1))
    b = np.clip(np.searchsorted(q[1:-1], p, side="right"), 0, nb - 1)
    return [dict(bin=i, n=int((b == i).sum()), mean_fc=float(p[b == i].mean()), obs=float(y[b == i].mean()),
                 p_lo=float(p[b == i].min()), p_hi=float(p[b == i].max())) for i in range(nb) if (b == i).sum()]


def stage_A2(log):
    models = json.load(open(os.path.join(RES, "gfs_model.json")))
    clim = {(r["basin"], r["month"]): r["rate_smoothed"] for r in models["climatology_hf24"]}
    F, n_all = gfs_features()
    E = pd.read_csv(os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz"), dtype={"time": str})
    E = E[E.season.isin(SEASONS)].copy(); E["time_i"] = E.time.astype(int)
    E = E.reset_index(drop=True)
    A, _, _ = load_archive()
    out = dict(n_gfs_track_rows_total=int(n_all), n_scored_fixes=int(len(F)))
    X = F[["msl", "dp12", "lat", "speed", "g800", "pac", "doy_c", "doy_s"]].copy()
    F["p"], use_dp = apply_model(F, X, models)
    F["used_dp12"] = use_dp
    F["pclim"] = [clim[(b, m)] for b, m in zip(F.basin, F.month)]
    mi, md = match_era5(F, E)
    F["era_i"] = mi; F["era_dist"] = md
    F["matched"] = mi >= 0
    F["y"] = np.where(F.matched, E.hf24.astype(int).reindex(np.where(mi >= 0, mi, 0)).values, 0)
    Xe = pd.DataFrame(dict(msl=np.nan, dp12=np.nan, lat=np.nan, speed=np.nan, g800=np.nan), index=F.index)
    mm = F.matched.values
    for c in ["msl", "dp12", "lat", "speed", "g800"]:
        Xe.loc[mm, c] = E[c].reindex(mi[mm]).values
    Xe["pac"] = F.pac.values; Xe["doy_c"] = F.doy_c.values; Xe["doy_s"] = F.doy_s.values
    Xe_m = Xe[mm]
    pe = np.full(len(F), np.nan)
    Fm = F[mm]
    use_dp_e = (Fm.lead.values >= 12) & Fm.dp12.notna().values       # same model as the GFS fix
    pe_m = np.zeros(mm.sum())
    if use_dp_e.any():
        pe_m[use_dp_e] = M.predict(models["model_dp12"], Xe_m[use_dp_e])
    if (~use_dp_e).any():
        pe_m[~use_dp_e] = M.predict(models["model_nodp12"], Xe_m[~use_dp_e])
    pe[mm] = pe_m
    F["pe"] = pe
    F["y2"] = secondary_truth(F, A)
    Fl12 = F[F.lead >= 12]
    out["counts"] = dict(
        fixes_by_basin_lead={f"{b}_{l}": int(((F.basin == b) & (F.lead == l)).sum()) for b in BASINS for l in LEADS_A2},
        no_era5_match_within_500km=int((~F.matched).sum()), matched=int(F.matched.sum()),
        no_match_by_basin_lead={f"{b}_{l}": int(((F.basin == b) & (F.lead == l) & ~F.matched).sum()) for b in BASINS for l in LEADS_A2},
        share_lead_ge12_nodp12_fallback=float((~Fl12.used_dp12).mean()), n_lead_ge12=int(len(Fl12)),
        share_lead_ge12_nodp12_by_lead={int(l): float((~g.used_dp12).mean()) for l, g in Fl12.groupby("lead")},
        speed_missing_share_by_lead={int(l): float(g.speed.isna().mean()) for l, g in F.groupby("lead")},
        era5_dp12_missing_when_dp12_model=int((Xe_m.dp12.isna() & use_dp_e).sum()),
        matched_basin_differs=int((E.basin.reindex(mi[mm]).values != Fm.basin.values).sum()))
    B5, B4 = draws(5), draws(4)
    S4 = [2022, 2023, 2024, 2025]
    res = {}
    for basin in BASINS:
        res[basin] = {}
        for lead in LEADS_A2:
            f = F[(F.basin == basin) & (F.lead == lead)].reset_index(drop=True)
            y, p, pc, seas = f.y.values, f.p.values, f.pclim.values, f.season.values
            st, bs = bss_stat(p, pc, y, seas, SEASONS, B5)
            ab = auc_boot(p, y, seas, SEASONS, B5)
            r = dict(n=int(len(f)), base_rate=float(y.mean()), n_pos=int(y.sum()), clim_mean=float(pc.mean()), pred_mean=float(p.mean()),
                     bss=st["bss"], bss_ci=st["ci"], bss_p_le0=st["p"], auc=float(auc(p, y)), auc_ci=pct(ab),
                     bss_by_season={int(u): float(1 - brier(p[seas == u], y[seas == u]).sum() / brier(pc[seas == u], y[seas == u]).sum())
                                    if y[seas == u].sum() else None for u in SEASONS})
            # threshold chosen on 2021 only
            m21 = seas == 2021
            thr = best_hss_threshold(p[m21], y[m21]) if y[m21].sum() > 0 else np.nan
            if np.isfinite(thr):
                te = ~m21
                fc, ob = p[te] >= thr, y[te].astype(bool)
                cat = cat_scores(contingency(fc, ob))
                hb = hss_boot(fc, ob, seas[te], S4, B4)
                r["cat"] = dict(threshold=thr, n=int(te.sum()), n_pos=int(ob.sum()), **cat, hss_ci=pct(hb),
                                hss_2021_insample=cat_scores(contingency(p[m21] >= thr, y[m21].astype(bool)))["hss"])
                # POD/FAR CIs
                r["cat"]["pod_ci"] = None
            # A2.2 matched subset
            mk = f.matched.values
            if mk.sum() > 0 and f.pe.notna().sum() > 0:
                ym, pm_, pem, pcm, sm_ = y[mk], p[mk], f.pe.values[mk], pc[mk], seas[mk]
                sg = sums_by_season(brier(pm_, ym), sm_, SEASONS); se = sums_by_season(brier(pem, ym), sm_, SEASONS)
                sr = sums_by_season(brier(pcm, ym), sm_, SEASONS)
                bg, be = 1 - sg.sum() / sr.sum(), 1 - se.sum() / sr.sum()
                dbs = (1 - sg[B5].sum(1) / sr[B5].sum(1)) - (1 - se[B5].sum(1) / sr[B5].sum(1))
                p_eq = float(max((1 + np.sum(dbs <= -0.05)) / (1 + len(dbs)), (1 + np.sum(dbs >= 0.05)) / (1 + len(dbs))))
                r["A2_2"] = dict(n=int(mk.sum()), base_rate=float(ym.mean()), bss_gfs_fed=float(bg), bss_era_fed=float(be),
                                 diff=float(bg - be), diff_ci=pct(dbs), abs_diff_le_0p05=bool(abs(bg - be) <= 0.05),
                                 ci_within_pm0p05=bool(pct(dbs)[0] >= -0.05 and pct(dbs)[1] <= 0.05), p_equiv_tost=p_eq,
                                 gfs_fed_ci=pct(bss_from_sums(sg, sr, B5)), era_fed_ci=pct(bss_from_sums(se, sr, B5)),
                                 g800_bias_gfs_minus_era=float((f.g800.values[mk] - E.g800.reindex(f.era_i.values[mk]).values).mean()),
                                 g800_rmse=float(np.sqrt(np.mean((f.g800.values[mk] - E.g800.reindex(f.era_i.values[mk]).values) ** 2))),
                                 msl_bias_gfs_minus_era=float((f.msl.values[mk] - E.msl.reindex(f.era_i.values[mk]).values).mean()))
            # secondary truth
            y2 = f.y2.values
            st2, _ = bss_stat(p, pc, y2, seas, SEASONS, B5)
            r["secondary_archive_truth"] = dict(base_rate=float(y2.mean()), n_pos=int(y2.sum()), bss_vs_hf24_clim=st2["bss"], bss_ci=st2["ci"],
                                                auc=float(auc(p, y2)))
            thr2 = best_hss_threshold(p[m21], y2[m21]) if y2[m21].sum() > 0 else np.nan
            if np.isfinite(thr2):
                te = ~m21
                c2 = cat_scores(contingency(p[te] >= thr2, y2[te].astype(bool)))
                r["secondary_archive_truth"].update(threshold=thr2, **{k: c2[k] for k in ("pod", "far", "hss", "bias")},
                                                    hss_ci=pct(hss_boot(p[te] >= thr2, y2[te].astype(bool), seas[te], S4, B4)))
            # A2.4 recalibration
            if y[m21].sum() >= 5:
                lg = lambda q: logit(np.clip(q, 1e-6, 1 - 1e-6))
                rc = _lr(lg(p[m21])[:, None], y[m21])
                te = ~m21
                pr = rc.predict_proba(lg(p[te])[:, None])[:, 1]
                b_before, bb = bss_stat(p[te], pc[te], y[te], seas[te], S4, B4)
                b_after, ba = bss_stat(pr, pc[te], y[te], seas[te], S4, B4)
                r["A2_4"] = dict(intercept=float(rc.intercept_[0]), slope=float(rc.coef_[0][0]), n_fit=int(m21.sum()), n_pos_fit=int(y[m21].sum()),
                                 n_apply=int(te.sum()), n_pos_apply=int(y[te].sum()),
                                 bss_before=b_before["bss"], bss_before_ci=b_before["ci"], bss_after=b_after["bss"], bss_after_ci=b_after["ci"],
                                 gain=float(b_after["bss"] - b_before["bss"]), gain_ci=pct(ba - bb))
            else:
                r["A2_4"] = dict(skipped="fewer than 5 positives in 2021", n_pos_fit=int(y[m21].sum()))
            if lead == 24:
                r["reliability_deciles"] = reliability_table(p, y)
            res[basin][lead] = r
        bl = [res[basin][l]["bss"] for l in LEADS_A2]
        res[basin]["A2_3"] = dict(bss_by_lead=bl, ci_by_lead=[res[basin][l]["bss_ci"] for l in LEADS_A2],
                                  non_increasing=bool(all(bl[i] >= bl[i + 1] for i in range(4))),
                                  steps=[float(bl[i + 1] - bl[i]) for i in range(4)])
        r24 = res[basin][24]
        res[basin]["A2_1"] = dict(bss=r24["bss"], ci=r24["bss_ci"], p=r24["bss_p_le0"], passes=bool(r24["bss"] > 0 and r24["bss_ci"][0] > 0))
        a22 = r24.get("A2_2", {})
        res[basin]["A2_2_lead24"] = dict(a22, passes=bool(a22.get("abs_diff_le_0p05", False)), p=a22.get("p_equiv_tost"))
    out["basins"] = res
    # diagnostics pooled
    mm = F.matched.values
    gd = (F.g800.values[mm] - E.g800.reindex(F.era_i.values[mm]).values)
    out["g800_gfs_minus_era_matched_all"] = dict(mean=float(gd.mean()), n=int(mm.sum()),
                                                 by_lead={int(l): float(gd[(F.lead.values[mm] == l)].mean()) for l in LEADS_A2},
                                                 by_basin={b: float(gd[(F.basin.values[mm] == b)].mean()) for b in BASINS})
    out["match_distance_km_median"] = float(np.nanmedian(F.era_dist[F.matched]))
    return out


# ---------------------------------------------------------------- BH
def bh(pvals):
    names = list(pvals); p = np.array([pvals[k] for k in names]); m = len(p)
    o = np.argsort(p); q = np.empty(m)
    run = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]
        run = min(run, p[i] * m / rank)
        q[i] = run
    return {k: float(q[i]) for i, k in enumerate(names)}


def fmt_ci(c):
    return "[%+.3f, %+.3f]" % (c[0], c[1]) if c and c[0] is not None else "[n/a]"


def main(nsim=300):
    t0 = time.time()
    lines = []
    w = lines.append
    a1 = stage_A1(nsim, w)
    a2 = stage_A2(w)
    P = {}
    crit = {}
    for b in BASINS:
        x = a1["basins"][b]
        P[f"A1.1 {b}"] = x["A1_1"]["p"]; crit[f"A1.1 {b}"] = x["A1_1"]["passes"]
        P[f"A1.2 {b}"] = x["A1_2"]["p_le0"]; crit[f"A1.2 {b}"] = x["A1_2"]["passes"]
        P[f"A1.3 {b}"] = x["A1_3"]["p_hss_lt_0p30"]; crit[f"A1.3 {b}"] = x["A1_3"]["passes"]
        y = a2["basins"][b]
        P[f"A2.1 {b}"] = y["A2_1"]["p"]; crit[f"A2.1 {b}"] = y["A2_1"]["passes"]
        P[f"A2.2 {b}"] = y["A2_2_lead24"]["p"]; crit[f"A2.2 {b}"] = y["A2_2_lead24"]["passes"]
    q = bh(P)
    summary = {k: dict(p=P[k], q_bh_10=q[k], passes=crit[k]) for k in P}

    w("GFS Stage A scoring. Pipeline A (ERA5) is a proxy; GFS seasons 2021-2025 (Oct-Apr) lie inside the gust-threshold calibration block 2021-22..2025-26.")
    w("Bootstrap: 2,000 season-block resamples, seed %d, 95%% percentile intervals; p one-sided, floor 1/2001." % SEED)
    w("")
    w("== Registered criteria (p: one-sided bootstrap; q: Benjamini-Hochberg over the 10 owned here; parent combines with Stage B's 4)")
    for k in P:
        w("  %-8s pass=%-5s p=%.4f q=%.4f" % (k, crit[k], P[k], q[k]))
    w("  passing: %d of 10" % sum(crit.values()))
    w("")
    w("== Archive categories: " + json.dumps(a1["archive_categories"]))
    w("   malformed rows (dropped): " + json.dumps(a1["archive_malformed_rows"]))
    w("   labels changed if the malformed Pac row is read as 2024110118: " + json.dumps(a1["label_changed_by_malformed_row_fix"]))
    w("")
    for b in BASINS:
        x = a1["basins"][b]
        w(f"== A1 {b}")
        for lead in LEADS_A1:
            r = x["per_lead"][lead]
            w("  lead %2d n=%d freq=%.3f (pos %d) BSS %+.4f %s p=%.4f | AUC %.3f %s (raw G AUC %.3f)" % (
                lead, r["n"], r["label_freq"], r["n_pos"], r["bss"], fmt_ci(r["bss_ci"]), r["bss_p_le0"], r["auc"], fmt_ci(r["auc_ci"]).replace("+", ""), r["auc_rawG"]))
        pw = x["per_lead"][24]["power_planted_bss0.05"]
        w("  power (A1.1 criterion passes, planted expected BSS %.3f, beta %.3f, %d sims): %.3f" % (pw["expected_bss"], pw["beta"], pw["nsim"], pw["power"]))
        w("  A1.1 passes=%s" % x["A1_1"]["passes"])
        t = x["A1_2"]
        w("  A1.2 AUC0 %.3f AUC24 %.3f AUC48 %.3f ordered=%s diff0-48 %+.3f %s p=%.4f passes=%s" % (
            t["auc0"], t["auc24"], t["auc48"], t["ordered_point"], t["diff_0_48"], fmt_ci(t["diff_ci"]), t["p_le0"], t["passes"]))
        t = x["A1_3"]
        w("  A1.3 2021 lead0 freq %.3f (n=%d) threshold %.2f kt (share above %.3f); 2022-25 lead 24 n=%d pos=%d POD %.3f FAR %.3f bias %.3f HSS %.3f %s p(HSS<0.30)=%.4f passes=%s" % (
            t["label_freq_2021_lead0"], t["n_cycles_2021_lead0"], t["threshold_kt"], t["share_above_2021_lead0"], t["n"], t["n_pos"], t["pod"], t["far"], t["bias"], t["hss"], fmt_ci(t["hss_ci"]), t["p_hss_lt_0p30"], t["passes"]))
        w("")
    c = a2["counts"]
    w("== A2 counts: scored fixes %d (of %d GFS track rows); no ERA5 match within 500 km %d; matched %d; median match dist %.0f km" % (
        a2["n_scored_fixes"], a2["n_gfs_track_rows_total"], c["no_era5_match_within_500km"], c["matched"], a2["match_distance_km_median"]))
    w("   lead>=12 fixes %d, share using nodp12 fallback %.4f; by lead %s" % (c["n_lead_ge12"], c["share_lead_ge12_nodp12_fallback"], json.dumps(c["share_lead_ge12_nodp12_by_lead"])))
    w("   speed missing share by lead %s; ERA5 dp12 missing where dp12 model used %d; matched basin differs %d" % (json.dumps(c["speed_missing_share_by_lead"]), c["era5_dp12_missing_when_dp12_model"], c["matched_basin_differs"]))
    w("   g800 GFS minus ERA5 at matched fixes: %s" % json.dumps(a2["g800_gfs_minus_era_matched_all"]))
    w("   no-match by basin_lead " + json.dumps(c["no_match_by_basin_lead"]))
    w("")
    for b in BASINS:
        y = a2["basins"][b]
        w(f"== A2 {b} (BSS vs frozen basin-month climatology; primary truth)")
        for lead in LEADS_A2:
            r = y[lead]
            w("  lead %2d n=%d base=%.4f clim=%.4f meanp=%.4f BSS %+.4f %s p=%.4f AUC %.3f %s" % (
                lead, r["n"], r["base_rate"], r["clim_mean"], r["pred_mean"], r["bss"], fmt_ci(r["bss_ci"]), r["bss_p_le0"], r["auc"], fmt_ci(r["auc_ci"]).replace("+", "")))
            if "cat" in r:
                k = r["cat"]
                w("          thr(2021)=%.4f 2022-25 n=%d pos=%d POD %.3f FAR %.3f HSS %.3f %s (2021 in-sample HSS %.3f)" % (
                    k["threshold"], k["n"], k["n_pos"], k["pod"], k["far"], k["hss"], fmt_ci(k["hss_ci"]), k["hss_2021_insample"]))
            if "A2_2" in r:
                k = r["A2_2"]
                w("          A2.2 matched n=%d GFS-fed BSS %+.4f ERA-fed %+.4f diff %+.4f %s |d|<=0.05:%s CI within:%s p_equiv=%.4f | g800 bias %+.2f kt rmse %.2f msl bias %+.2f" % (
                    k["n"], k["bss_gfs_fed"], k["bss_era_fed"], k["diff"], fmt_ci(k["diff_ci"]), k["abs_diff_le_0p05"], k["ci_within_pm0p05"], k["p_equiv_tost"], k["g800_bias_gfs_minus_era"], k["g800_rmse"], k["msl_bias_gfs_minus_era"]))
            k = r["A2_4"]
            if "skipped" in k:
                w("          A2.4 skipped: " + k["skipped"])
            else:
                w("          A2.4 recal a=%+.3f b=%.3f (2021 n=%d pos=%d) 2022-25 BSS before %+.4f %s after %+.4f %s gain %+.4f %s" % (
                    k["intercept"], k["slope"], k["n_fit"], k["n_pos_fit"], k["bss_before"], fmt_ci(k["bss_before_ci"]), k["bss_after"], fmt_ci(k["bss_after_ci"]), k["gain"], fmt_ci(k["gain_ci"])))
            s2 = r["secondary_archive_truth"]
            w("          secondary archive truth base=%.4f (pos %d) BSS vs hf24 clim %+.4f %s AUC %.3f%s" % (
                s2["base_rate"], s2["n_pos"], s2["bss_vs_hf24_clim"], fmt_ci(s2["bss_ci"]), s2["auc"],
                (" | thr %.4f POD %.3f FAR %.3f HSS %.3f %s" % (s2["threshold"], s2["pod"], s2["far"], s2["hss"], fmt_ci(s2["hss_ci"]))) if "threshold" in s2 else ""))
        t = y["A2_3"]
        w("  A2.3 BSS by lead %s non-increasing=%s" % (" ".join("%+.4f" % v for v in t["bss_by_lead"]), t["non_increasing"]))
        w("  A2.1 passes=%s ; A2.2 (lead 24) passes=%s" % (y["A2_1"]["passes"], y["A2_2_lead24"]["passes"]))
        w("  reliability lead 24 (deciles of forecast probability): bin n mean_fc obs")
        for rr in y[24]["reliability_deciles"]:
            w("    %d n=%d fc=%.4f obs=%.4f (p %.4f-%.4f)" % (rr["bin"], rr["n"], rr["mean_fc"], rr["obs"], rr["p_lo"], rr["p_hi"]))
        w("")
    w("elapsed %.0f s" % (time.time() - t0))
    txt = "\n".join(lines) + "\n"
    open(os.path.join(RES, "stageA.txt"), "w").write(txt)

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o
    json.dump(clean(dict(summary=summary, A1=a1, A2=a2)), open(os.path.join(RES, "stageA.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    ns = 300
    if "--nsim" in sys.argv:
        ns = int(sys.argv[sys.argv.index("--nsim") + 1])
    main(ns)
