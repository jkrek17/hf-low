"""Registered tests of PREREGISTRATION_GFS.md. ERA5 pipeline A is a PROXY (training record and truth).
usage: [GFS_WORK=DIR] [STORM_PATCH_DIR=DIR] score_gfs.py   (defaults: gfs/work/ and /mnt/project-files/storm-patches)      -> results/gfs_skill.txt, gfs_tests.csv
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); RE = os.path.join(HERE, "..", "..")
WORK = os.environ.get("GFS_WORK") or os.path.join(HERE, "work")
for p in ("gfs_outlook", "gfs_test", "intensity"):
    sys.path.insert(0, os.path.join(RE, p))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
_argv = sys.argv; sys.argv = [_argv[0]]
import score_outlook_full as SF
sys.argv = _argv
from frozen import Frozen
SO, AS, BS = SF.SO, SF.AS, SF.BS
SEAS = SF.SEASONS; RES = os.path.join(HERE, "results")
LINES = []
def out(t=""):
    print(t, flush=True); LINES.append(t)


def load_patches():
    P = {}
    for f in glob.glob(os.path.join(WORK, "gfs_patches", "*.npz")):
        z = np.load(f)
        for k, x in zip(z["keys"], z["x"]):
            P[str(k)] = x
    return P


def main():
    S = pd.read_csv(os.path.join(RES, "gfs_fixes_lead24.csv.gz"), dtype={"cycle": str})
    Pm = load_patches()
    S["pkey"] = S.cycle + "|" + S.track.astype(str)
    have = S.pkey.isin(Pm)
    out(f"GFS lead 24 scored sample of PR 185: {len(S)}; with a patch {int(have.sum())}; cycles missing: {S.cycle[~have].nunique()}")
    S = S[have].reset_index(drop=True)
    X = np.stack([Pm[k] for k in S.pkey]).astype(np.float32)
    nanfrac = np.isnan(X[:, :4]).mean()
    out(f"NaN fraction in anomaly channels: {nanfrac:.5f}")
    S["month"] = (S.valid_time // 10000 % 100).astype(int); S["hour"] = (S.valid_time % 100).astype(int)
    y = S.y.astype(int).values; seas = S.season.values
    out("Truth classes (below gale, gale, storm, HF): " + str(np.bincount(y, minlength=4).tolist()) + "; by season: " + ", ".join(f"{s}: {int((seas == s).sum())}" for s in SEAS))
    frz = {"full": Frozen(os.path.join(RES, "frozen_pca_full.npz")), "out": Frozen(os.path.join(RES, "frozen_pca_out.npz"))}
    models = json.load(open(os.path.join(RES, "frozen_models.json")))
    mfull = json.load(open(os.path.join(RE, "gfs_outlook", "results", "outlook_models_full.json")))
    models["models"]["FULL"] = mfull["models"]["FULL"]; models["models"]["FULL_nodp12"] = mfull["models"]["FULL_nodp12"]
    clim = {(r["basin"], r["month"]): np.array(r["probs_smoothed"]) for r in mfull["climatology_gcls"]}
    C = np.array([clim[(b, m)] for b, m in zip(S.basin, S.month)])
    use_dp = S.dp12.notna().values
    def add_pcs(df, Xp, tags=("full", "out")):
        df = df.copy()
        for t in tags:
            pc = frz[t].pcs(Xp)
            for i in range(pc.shape[1]):
                df[f"pc{i+1}_{t}"] = pc[:, i]
        return df
    G = add_pcs(S, X)
    P = {"clim": C}
    for k in ("FULL", "FULL7", "FULLPC7", "FULLPC7_OUT"):
        P[k] = SO.sm(SO.run_model(models, k, G, use_dp))
    # placebo: permute patches among fixes of the same basin, month, valid hour
    rng = np.random.default_rng(20261014); perm = np.arange(len(S))
    for _, ix in S.groupby(["basin", "month", "hour"]).indices.items():
        perm[ix] = rng.permutation(ix)
    Gp = add_pcs(S, X[perm], tags=("full",))
    P["PLACEBO"] = SO.sm(SO.run_model(models, "FULLPC7", Gp, use_dp))
    L = {k: BS.rps(P[k], y, 4) for k in P}
    B5 = SO.draws(5); cl = L["clim"]
    def R(k): return SO.rpss(L[k], cl, seas, SEAS, B5)
    def GN(a, b, m=None):
        m = np.ones(len(S), bool) if m is None else m
        return SO.gain(L[a][m], L[b][m], cl[m], seas[m], SEAS, B5)
    out("\n== RPSS vs basin-month climatology, lead 24, GFS-fed (point [90% season-block])")
    for k in ("FULL", "FULL7", "FULLPC7", "FULLPC7_OUT", "PLACEBO"):
        t = R(k); out(f"  {k:12s} {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]")
    rows = []
    t1 = GN("FULLPC7", "FULL7"); s1 = GN("FULLPC7_OUT", "FULL7")
    verdict = "survives" if (t1[0] >= 0.010 and t1[1] > 0) else ("no useful gain" if t1[2] < 0.005 else "cannot tell")
    q = BS.bh([SO.p_le(t1[3]), SO.p_le(s1[3])])
    out(f"\n== T1 (primary): FULLPC7 minus FULL7, GFS-fed: {t1[0]:+.4f} [{t1[1]:+.4f}, {t1[2]:+.4f}]  p {SO.p_le(t1[3]):.4f} q {q[0]:.4f}  verdict (rule): {verdict}")
    out(f"== S1 (outer pattern, inner 300 km removed): FULLPC7_OUT minus FULL7: {s1[0]:+.4f} [{s1[1]:+.4f}, {s1[2]:+.4f}]  p {SO.p_le(s1[3]):.4f} q {q[1]:.4f}")
    rows += [dict(test="T1", gain=t1[0], lo90=t1[1], hi90=t1[2], q=q[0], verdict=verdict), dict(test="S1", gain=s1[0], lo90=s1[1], hi90=s1[2], q=q[1])]
    pl = GN("PLACEBO", "FULL7")
    out(f"\n== Descriptive. Placebo (permuted GFS patches) minus FULL7: {pl[0]:+.4f} [{pl[1]:+.4f}, {pl[2]:+.4f}]")
    out("  T1 gain by season: " + ", ".join(f"{s}-{str(s+1)[2:]} {(L['FULL7'][seas==s].sum()-L['FULLPC7'][seas==s].sum())/cl[seas==s].sum():+.4f}" for s in SEAS))
    for lab, m in (("Atlantic", (S.basin == "atl").values), ("Pacific", (S.basin == "pac").values)):
        t = GN("FULLPC7", "FULL7", m); out(f"  T1 gain {lab:9s} n {m.sum():6d}: {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]")
    cum = lambda Pm_: np.cumsum(Pm_, 1)[:, :3]; Yc = np.cumsum(np.eye(4)[y], 1)[:, :3]
    sh = [(((cum(P["FULL7"])[:, j] - Yc[:, j]) ** 2).sum() - ((cum(P["FULLPC7"])[:, j] - Yc[:, j]) ** 2).sum()) / cl.sum() for j in range(3)]
    out(f"  share of the T1 gain by threshold (>=gale, >=storm, >=HF): {sh[0]:+.4f}, {sh[1]:+.4f}, {sh[2]:+.4f}")
    yh = (y == 3).astype(float)
    for k in ("FULL7", "FULLPC7"):
        bs_ = 1 - ((P[k][:, 3] - yh) ** 2).sum() / ((C[:, 3] - yh) ** 2).sum(); out(f"  HF Brier skill {k}: {bs_:+.3f}")
    # neighbour tracks
    T = AS.gfs_features()[0]; T = T[T.lead == 24]
    byv = {str(v): g[["track", "lat", "lon"]].values for v, g in T.groupby("cycle")}
    nb = np.zeros(len(S), bool)
    for i, r in enumerate(S.itertuples()):
        a = byv.get(str(r.cycle))
        if a is None: continue
        m = a[:, 0] != r.track
        if m.any(): nb[i] = bool((AS.hav(r.lat, r.lon, a[m, 1].astype(float), a[m, 2].astype(float)) < 800).any())
    k_ = ~nb
    t = GN("FULLPC7", "FULL7", k_); out(f"  T1 gain on fixes with no other GFS track within 800 km (kept {int(k_.sum())}, flagged {nb.mean():.3f}): {t[0]:+.4f} [{t[1]:+.4f}, {t[2]:+.4f}]")
    # ERA5-fed comparison, overlap season 2021-22
    D = BS.load(); E = D[D.season.isin(SEAS)].copy().reset_index(drop=True); E["time_i"] = E.time.astype(int)
    z = np.load(os.path.join(os.environ.get("STORM_PATCH_DIR", "/mnt/project-files/storm-patches"), "patches_2021_22.npz"), allow_pickle=True)
    EP = {str(k): x for k, x in zip(z["keys"], z["x"])}
    ek = (S.era_track.astype(int).astype(str) + "|" + S.era_time.astype(int).astype(str)).values
    ov = (seas == 2021) & np.array([k in EP for k in ek])
    out(f"\n== ERA5-fed comparison, overlap season 2021-22 (one season; illustration only): GFS fixes with an ERA5 patch for the matched ERA5 fix: {int(ov.sum())} of {int((seas == 2021).sum())}")
    if ov.sum() > 200:
        So = S[ov].reset_index(drop=True); Xe = np.stack([EP[k] for k in ek[ov]]).astype(np.float32)
        Ee = E.reindex(So.era_i.values.astype(int)).reset_index(drop=True)
        young = (Ee.young == 1).values
        XE = Ee[["msl", "dp12", "lat", "speed", "g800"] + SF.ALLENV + ["nosst"]].copy(); XE.loc[young, "dp12"] = np.nan
        for c in ("pac", "doy_c", "doy_s"): XE[c] = So[c].values
        XE = add_pcs(XE, Xe); use_e = So.dp12.notna().values & ~young
        yo = So.y.astype(int).values; Co = C[ov]; seo = np.full(len(So), 2021)
        Pe = {k: SO.sm(SO.run_model(models, k, XE, use_e)) for k in ("FULL7", "FULLPC7")}
        Pg = {k: P[k][ov] for k in ("FULL7", "FULLPC7")}
        cl_o = BS.rps(Co, yo, 4)
        ge = (BS.rps(Pe["FULL7"], yo, 4).sum() - BS.rps(Pe["FULLPC7"], yo, 4).sum()) / cl_o.sum()
        gg = (BS.rps(Pg["FULL7"], yo, 4).sum() - BS.rps(Pg["FULLPC7"], yo, 4).sum()) / cl_o.sum()
        out(f"  gain of the pattern terms, ERA5-fed {ge:+.4f}, GFS-fed {gg:+.4f} on the same {len(So)} fixes (ratio {gg / ge:.2f}); no interval (one season)")
        rp = lambda Pm_: 1 - BS.rps(Pm_, yo, 4).sum() / cl_o.sum()
        out("  RPSS on the same fixes: FULL7 ERA5-fed %+.4f, GFS-fed %+.4f; FULLPC7 ERA5-fed %+.4f, GFS-fed %+.4f" % (rp(Pe["FULL7"]), rp(Pg["FULL7"]), rp(Pe["FULLPC7"]), rp(Pg["FULLPC7"])))
        rows.append(dict(test="overlap_2021_22", gain=gg, lo90=np.nan, hi90=np.nan, note=f"ERA5-fed {ge:.4f}"))
        out("  Bias diagnostic (GFS minus ERA5 patch, matched fixes of 2021-22), per channel: mean difference, SD of ERA5 patch values, centre pixel mean difference")
        for i, nme in enumerate(["a_z500 (m)", "a_msl (hPa)", "a_u250 (m/s)", "a_v250 (m/s)"]):
            g_, e_ = X[ov][:, i].astype(float), Xe[:, i].astype(float)
            out(f"    {nme:13s} mean diff {np.nanmean(g_ - e_):+.2f}, ERA5 SD {np.nanstd(e_):.2f}, centre diff {np.nanmean(g_[:, 6, 6] - e_[:, 6, 6]):+.2f}, corr {np.corrcoef(np.nan_to_num(g_).ravel(), np.nan_to_num(e_).ravel())[0, 1]:.3f}")
    pd.DataFrame(rows).to_csv(os.path.join(RES, "gfs_tests.csv"), index=False)
    open(os.path.join(RES, "gfs_skill.txt"), "w").write("\n".join(LINES) + "\n")
    np.savez_compressed(os.path.join(WORK, "gfs_oof.npz"), keys=S.pkey.values, y=y, **{k: P[k] for k in P})


if __name__ == "__main__":
    main()
