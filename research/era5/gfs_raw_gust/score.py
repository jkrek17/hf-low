"""Raw GFS gust against OPC HF designations (PREREGISTRATION.md, commit fcfc9eb). Pipeline A detector on GFS = PROXY fixes.
Writes results/raw_gust.json and results/raw_gust.txt.   usage: score.py"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GT = os.path.join(HERE, "..", "gfs_test")
sys.path.insert(0, GT); sys.path.insert(0, os.path.join(GT, "..", "intensity"))
import a_score as S

RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)
SEED = 20261013; NB = 2000
SEAS = S.SEASONS; BAS = S.BASINS
B5 = S.draws(5, NB, SEED)
KT_NOM = 64.0; KT_ERA = 71.7


def label_conc(F, th_pos, rad=600.0, win=3):
    """Y0: archive HF fix within +-win h of valid time and within rad km. th_pos: (hours, lat, lon) sorted by hours."""
    th, la, lo = th_pos
    v = pd.to_datetime(F.valid_time.astype(str), format="%Y%m%d%H").values.astype("datetime64[h]").astype(np.int64)
    a = np.searchsorted(th, v - win, side="left"); b = np.searchsorted(th, v + win, side="right")
    out = np.zeros(len(F), int); flat, flon = F.lat.values, F.lon.values
    for i in np.where(b > a)[0]:
        out[i] = int((S.hav(flat[i], flon[i], la[a[i]:b[i]], lo[a[i]:b[i]]) <= rad).any())
    return out


def hf_pos(A):
    H = A[A.Category == "HF"].copy(); H["th"] = H.t.values.astype("datetime64[h]").astype(np.int64)
    H = H.sort_values("th")
    return H.th.values, H.lat.values, H.lon.values


def pooled(fc, ob, seas):
    return S.cat_scores(S.contingency(fc, ob))


def stat_boot(fc, ob, seas, B=B5, useas=SEAS):
    """HSS point, interval; POD/FAR/bias point; bootstrap vector of HSS."""
    cs = S.cat_scores(S.contingency(fc, ob)); hb = S.hss_boot(fc, ob, seas, useas, B)
    return dict(cs, hss_ci=S.pct(hb)), hb


def match_T(g, y, k=None):
    """threshold T with #(g > T) == #(y==1) (count match); k = number of positives."""
    k = int(y.sum()) if k is None else k
    gs = np.sort(g)[::-1]
    return float(gs[k]) if 0 <= k < len(gs) else float("-inf")


def heldout_T(F, ycol, lead_ref=0):
    """per (basin, held-out season): T from the other four seasons at lead_ref, on label ycol."""
    T = {}
    for b in BAS:
        for s in SEAS:
            d = F[(F.basin == b) & (F.lead == lead_ref) & (F.season != s)]
            T[(b, s)] = match_T(d.G.values, d[ycol].values)
    return T


def apply_fc(F, T, thr_col=None, nominal=None):
    fc = np.zeros(len(F), bool)
    for (b, s), t in T.items():
        m = ((F.basin == b) & (F.season == s)).values
        fc[m] = F.G.values[m] > t
    return fc


def main():
    A, _, _ = S.load_archive()
    pos = hf_pos(A)
    F, _ = S.gfs_features()
    F["season"] = F.season.astype(int)
    F["G"] = F.g800
    F["Y0"] = label_conc(F, pos)
    F["Y24"] = S.secondary_truth(F, A)
    models = json.load(open(os.path.join(S.RES, "gfs_model.json")))
    X = F[["msl", "dp12", "lat", "speed", "g800", "pac", "doy_c", "doy_s"]].copy()
    F["p"], _ = S.apply_model(F, X, models)
    out = dict(n_fixes=int(len(F)), n_Y0=int(F.Y0.sum()), n_Y24=int(F.Y24.sum()))
    out["counts"] = {b: {int(l): dict(n=int(((F.basin == b) & (F.lead == l)).sum()), y0=int(F[(F.basin == b) & (F.lead == l)].Y0.sum()),
                                      y24=int(F[(F.basin == b) & (F.lead == l)].Y24.sum())) for l in S.LEADS_A2} for b in BAS}

    # held-out count-matched thresholds, Y0, lead-0 reference
    T0 = heldout_T(F, "Y0")
    out["T_heldout_Y0"] = {b: {int(s): T0[(b, s)] for s in SEAS} for b in BAS}
    F["fc_nom"] = F.G > KT_NOM
    F["fc_match"] = apply_fc(F, T0)
    F["fc_era"] = F.G > KT_ERA
    tests = {}
    prof = {}
    for b in BAS:
        for l in S.LEADS_A2:
            d = F[(F.basin == b) & (F.lead == l)]
            y = d.Y0.values.astype(bool); s = d.season.values
            row = {}
            for nm, col in (("nominal64", "fc_nom"), ("matched", "fc_match"), ("era717", "fc_era")):
                st, hb = stat_boot(d[col].values, y, s)
                row[nm] = st
                if l in (0, 24) and nm in ("nominal64", "matched"):
                    tests[f"{'R1' if nm == 'nominal64' else 'R2'}_{b}_lead{l}"] = dict(
                        hss=st["hss"], ci=st["hss_ci"], pod=st["pod"], far=st["far"], bias=st["bias"],
                        usable=bool(st["hss"] >= 0.30 and st["hss_ci"][0] >= 0.30),
                        p_hss_lt_0p30=float((1 + np.sum(hb < 0.30)) / (1 + np.sum(~np.isnan(hb)))))
            # AUC of raw G and model p on both labels
            row["auc_G_Y0"] = float(S.auc(d.G.values, d.Y0.values)); row["auc_G_Y24"] = float(S.auc(d.G.values, d.Y24.values))
            row["auc_p_Y24"] = float(S.auc(d.p.values, d.Y24.values)); row["auc_p_Y0"] = float(S.auc(d.p.values, d.Y0.values))
            prof[f"{b}_{l}"] = row
    out["profile"] = prof

    # lead-specific matched threshold (sensitivity S1)
    sens = {}
    for b in BAS:
        for l in S.LEADS_A2:
            Tl = heldout_T(F, "Y0", lead_ref=l)
            d = F[(F.basin == b) & (F.lead == l)]
            fc = np.zeros(len(d), bool)
            for s in SEAS:
                m = (d.season == s).values; fc[m] = d.G.values[m] > Tl[(b, s)]
            st, _ = stat_boot(fc, d.Y0.values.astype(bool), d.season.values)
            sens[f"{b}_{l}"] = dict(hss=st["hss"], ci=st["hss_ci"], pod=st["pod"], far=st["far"], bias=st["bias"])
    out["lead_specific_matched"] = sens

    # R3: ERA5 proxy 71.7 vs GFS lead-0 matched, Y0, same seasons/months/boxes
    E = pd.read_csv(os.path.join(GT, "..", "intensity", "results", "fixes_2004.csv.gz"), dtype={"time": str})
    E = E[E.season.isin(SEAS) & E.basin.isin(BAS)].copy()
    E["valid_time"] = E.time.astype(int)
    E = E[(E.valid_time % 100).isin([0, 12]) & (E.valid_time // 10000 % 100).isin([10, 11, 12, 1, 2, 3, 4])].reset_index(drop=True)
    E["Y0"] = label_conc(E, pos)
    out["era5_n_fixes"] = int(len(E)); out["era5_n_Y0"] = int(E.Y0.sum())
    r3 = {}; r3b = {}
    for b in BAS:
        e = E[E.basin == b]; g0 = F[(F.basin == b) & (F.lead == 0)]
        ste, hbe = stat_boot((e.g800 > KT_ERA).values, e.Y0.values.astype(bool), e.season.values)
        stg, hbg = stat_boot(g0.fc_match.values, g0.Y0.values.astype(bool), g0.season.values)
        sg7, _ = stat_boot(g0.fc_era.values, g0.Y0.values.astype(bool), g0.season.values)
        diff = hbg - hbe
        r3b[b] = diff
        ci = S.pct(diff)
        r3[b] = dict(era5=ste, gfs_lead0_matched=stg, gfs_lead0_717=sg7, diff_hss=float(stg["hss"] - ste["hss"]), diff_ci=ci,
                     same=bool(ci[0] >= -0.10 and ci[1] <= 0.10),
                     p_diff_outside=float(min(1, 2 * min((1 + np.sum(diff <= 0)) / (1 + len(diff)), (1 + np.sum(diff >= 0)) / (1 + len(diff))))))
    out["R3"] = r3

    # R4: lead 24, Y24, AUC model minus AUC raw
    r4 = {}
    for b in BAS:
        d = F[(F.basin == b) & (F.lead == 24)].reset_index(drop=True)
        s = d.season.values
        ap = S.auc_boot(d.p.values, d.Y24.values, s, SEAS, B5); ag = S.auc_boot(d.G.values, d.Y24.values, s, SEAS, B5)
        diff = ap - ag
        r4[b] = dict(auc_model=float(S.auc(d.p.values, d.Y24)), auc_raw=float(S.auc(d.G.values, d.Y24)),
                     diff=float(S.auc(d.p.values, d.Y24) - S.auc(d.G.values, d.Y24)), diff_ci=S.pct(diff),
                     raw_not_worse=bool(S.pct(diff)[1] <= 0.05), n=int(len(d)), n_pos=int(d.Y24.sum()))
    out["R4"] = r4
    out["tests"] = tests

    # BH over the 12 primary tests: one-sided p for R1/R2 (HSS<0.30), R3 p (two-sided diff), R4 p (diff>0.05 one-sided)
    ps = {}
    for k, v in tests.items():
        ps[k] = v["p_hss_lt_0p30"]
    for b in BAS:   # R3: TOST against +-0.10 (Deviation 1); the two-sided p for diff != 0 stays in r3[b]["p_diff_outside"]
        bb = r3b[b]
        ps[f"R3_{b}"] = float(max((1 + np.sum(bb <= -0.10)) / (1 + len(bb)), (1 + np.sum(bb >= 0.10)) / (1 + len(bb))))
    for b in BAS:
        d = F[(F.basin == b) & (F.lead == 24)].reset_index(drop=True)
        ap = S.auc_boot(d.p.values, d.Y24.values, d.season.values, SEAS, B5); ag = S.auc_boot(d.G.values, d.Y24.values, d.season.values, SEAS, B5)
        diff = ap - ag
        ps[f"R4_{b}"] = float((1 + np.sum(diff >= 0.05)) / (1 + len(diff)))   # non-inferiority, H0: model - raw >= 0.05 (Deviation 1)
    keys = list(ps); pv = np.array([ps[k] for k in keys]); o = np.argsort(pv); m = len(pv)
    q = np.empty(m); run = 1.0
    for rank, i in zip(range(m, 0, -1), o[::-1]):
        run = min(run, pv[i] * m / rank); q[i] = run
    out["bh"] = {k: dict(p=float(ps[k]), q=float(q[i])) for i, k in enumerate(keys)}

    # S3 event level
    out["S3"] = event_level(F, A, T0)
    json.dump(out, open(os.path.join(RES, "raw_gust.json"), "w"), indent=1, default=float)
    F[["cycle", "lead", "valid_time", "track", "basin", "season", "lat", "lon", "msl", "G", "p", "Y0", "Y24"]].to_csv(
        os.path.join(RES, "raw_gust_fixes.csv.gz"), index=False)
    print("done", json.dumps(out["tests"], default=float)[:400])


def event_level(F, A, T0):
    H = A[A.Category == "HF"].copy(); H = H.sort_values("t")
    first = H.groupby("ID").first().reset_index()
    first["m"] = first.t.dt.month
    first["season"] = np.where(first.m >= 6, first.t.dt.year, first.t.dt.year - 1)
    ev = first[first.season.isin(SEAS) & first.m.isin([10, 11, 12, 1, 2, 3, 4])].copy()
    rows = []
    for _, e in ev.iterrows():
        t = e.t; v = t.floor("12h") if t.hour % 12 < 6 else t.floor("12h") + pd.Timedelta(hours=12)
        vi = int(v.strftime("%Y%m%d%H"))
        c = F[F.valid_time == vi]
        if len(c) == 0:
            rows.append(dict(ID=e.ID, covered=False)); continue
        d = S.hav(e.lat, e.lon, c.lat.values, c.lon.values); near = c[d <= 600.0]
        b = near.basin.iloc[0] if len(near) else None
        if len(near) == 0:
            rows.append(dict(ID=e.ID, covered=True, tracked=False)); continue
        bas = "atl" if e.ID and A[A.ID == e.ID].file.iloc[0] == "Atl" else "pac"
        Ts = T0[(bas, int(e.season))]
        r0 = near[near.lead == 0]
        hit0 = bool(len(r0) and (r0.G > Ts).any())
        ex = near[near.G > Ts]
        rows.append(dict(ID=e.ID, covered=True, tracked=True, basin=bas, hit_lead0=hit0, any_exceed=bool(len(ex)),
                         max_lead_exceed=int(ex.lead.max()) if len(ex) else -1))
    R = pd.DataFrame(rows)
    out = dict(n_events=int(len(R)), covered=int(R.covered.sum()), tracked=int(R.get("tracked", pd.Series(dtype=bool)).fillna(False).sum()))
    if out["tracked"]:
        Tr = R[R.tracked == True]
        out.update(hit_lead0=int(Tr.hit_lead0.sum()), any_exceed=int(Tr.any_exceed.sum()),
                   max_lead_exceed_median=float(Tr[Tr.max_lead_exceed >= 0].max_lead_exceed.median()),
                   lead_dist={int(k): int(v) for k, v in Tr.max_lead_exceed.value_counts().sort_index().items()})
    return out


if __name__ == "__main__":
    main()
