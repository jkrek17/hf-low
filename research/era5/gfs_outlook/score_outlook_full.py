"""Score the FULL / FULL_NOSF 24 h outlook models (results/outlook_models_full.json, frozen) and the gain retention of A, B, FULL when fed GFS versus ERA5
on the same fixes. Plan: PREREGISTRATION_FULL.md and its Deviations. ERA5 pipeline A is a proxy and is the truth here; it is not observation.
Reuses score_outlook.py (part 1) for matching, truth, dp12/fallback rule, RPSS and the season-block bootstrap (5 GFS seasons 2021-22..2025-26,
2,000 draws, seed 20261013, 90 % intervals, BH over the 6 registered tests).
ERA5-fed side: the matched pipeline A fix supplies msl, dp12 (young -> fallback), lat, speed, g800 and every environment term from env_2004.csv.gz;
season terms (pac, doy) are the GFS fix's. No ERA5 field is pulled.
usage: score_outlook_full.py [OUTDIR] [ENVTABLE_FULL] [ENVTABLE_PART1]      writes OUTDIR/outlook_full.txt and outlook_full.json
"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_argv = sys.argv
sys.argv = [_argv[0]]                       # score_outlook reads sys.argv at import; its own OUT/ENVT are unused here
import score_outlook as SO
sys.argv = _argv
AS, BS = SO.AS, SO.BS
RES = SO.RES
OUT = sys.argv[1] if len(sys.argv) > 1 else RES
ENVF = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "gfs_outlook_env_full.csv.gz")
ENV1 = sys.argv[3] if len(sys.argv) > 3 else os.path.join(RES, "gfs_outlook_env.csv.gz")
SEASONS, LEADS, SEEDV = SO.SEASONS, [0, 24], SO.SEEDV
LINES = []
HART = ["B", "VTL", "VTU"]
ENVC = ["jet250", "div300", "vadv500", "eady", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]
ALLENV = HART + ENVC
REQ = ["B", "VTL", "VTU", "div300", "tcwv", "jet250", "vadv500", "eady"]      # sst, sstgrad, sst_t500 are legitimately NaN (nosst); flux also at lead 24


def out(s=""):
    print(s)
    LINES.append(s)


def ratio_stat(num_g, num_e, den_g, den_e, clim, seas, B):
    """gain retention: (gain over S, GFS-fed) / (gain over S, ERA5-fed); both gains are sum-of-loss differences over the same clim denominator,
    so the ratio is (sum S_g - sum X_g) / (sum S_era - sum X_era) per season-block draw (same draw for numerator and denominator)."""
    n = SO.sums(den_g - num_g, seas, SEASONS)
    d = SO.sums(den_e - num_e, seas, SEASONS)
    pt = n.sum() / d.sum()
    dn, dd = n[B].sum(1), d[B].sum(1)
    bt = dn / dd
    return pt, np.percentile(bt, 5), np.percentile(bt, 95), bt, float((dd <= 0).sum())


def verdict_ratio(t):
    return "survives" if t[1] >= 0.5 else ("shrinks" if t[2] < 0.5 else "cannot tell")


def main():
    models = json.load(open(os.path.join(RES, "outlook_models.json")))
    mfull = json.load(open(os.path.join(RES, "outlook_models_full.json")))
    allm = {"models": {**models["models"], **mfull["models"]}}
    clim = {(r["basin"], r["month"]): np.array(r["probs_smoothed"]) for r in models["climatology_gcls"]}
    assert models["climatology_gcls"] == mfull["climatology_gcls"]
    F, n_all = AS.gfs_features()
    F = F[F.lead.isin([0, 12, 24, 36, 48])].copy().reset_index(drop=True)
    D = BS.load()
    E = D[D.season.isin(SEASONS)].copy().reset_index(drop=True)
    E["time_i"] = E.time.astype(int)
    mi, md = AS.match_era5(F, E)
    F["era_i"] = mi
    F["matched"] = mi >= 0
    F["y"] = np.where(F.matched, E.gcls.reindex(np.where(mi >= 0, mi, 0)).values, -1).astype(int)
    out("GFS FULL-model outlook test (part 2). Truth = ERA5 pipeline A, a PROXY (not observation); ERA5-fed side = pipeline A proxy tables.")
    out(f"GFS track rows (all, 2021-22..2025-26): {n_all}; atl/pac fixes at valid 00/12 UTC, leads 0/12/24/36/48: {len(F)}")
    out(f"No pipeline A fix within 500 km at the same valid time (dropped): total {int((~F.matched).sum())} of {len(F)}")
    Fm = F[F.matched].copy().reset_index(drop=True)
    out(f"Matched: {len(Fm)}; lead 0: {int((Fm.lead == 0).sum())}; lead 24: {int((Fm.lead == 24).sum())}")
    Fm["cycle"] = Fm.cycle.astype(str)
    T1 = pd.read_csv(ENV1, dtype={"cycle": str})
    T2 = pd.read_csv(ENVF, dtype={"cycle": str})
    out(f"Part 1 environment table rows {len(T1)}; part 2 table rows {len(T2)} (leads {sorted(T2.lead.unique().tolist())})")
    Fm = Fm.merge(T1[["cycle", "lead", "track", "jet250", "vadv500", "eady"]], on=["cycle", "lead", "track"], how="left")
    Fm = Fm.merge(T2[["cycle", "lead", "track"] + HART + ["div300", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]], on=["cycle", "lead", "track"], how="left")
    Fm["nosst"] = Fm.sst.isna().astype(float)
    flux_ok = (Fm.lead == 0) | Fm.flux.notna()
    Fm["gfs_ok"] = Fm[REQ].notna().all(axis=1) & flux_ok
    # ERA5-fed environment from the matched ERA5 fix
    ei = Fm.era_i.values
    ecols = ALLENV + ["nosst"]
    Ee = pd.DataFrame({c: E[c].reindex(ei).values for c in ecols + ["msl", "dp12", "lat", "speed", "g800"]})
    Fm["era_ok"] = Ee[REQ + ["flux"]].notna().all(axis=1).values
    Fm["ok"] = Fm.gfs_ok & Fm.era_ok
    out("\n== Counts at each filter step (matched fixes)")
    out(f"{'lead':>5} {'matched':>8} {'part1 env':>10} {'part2 row':>10} {'GFS terms ok':>13} {'ERA5 terms ok':>14} {'both = scored sample':>21}")
    for l in LEADS:
        g = Fm[Fm.lead == l]
        out(f"{l:5d} {len(g):8d} {int(g.jet250.notna().sum()):10d} {int((g.B.notna() | g.tcwv.notna()).sum()):10d} "
            f"{int(g.gfs_ok.sum()):13d} {int(g.era_ok.sum()):14d} {int(g.ok.sum()):21d}")
    out("Scored sample by season, lead 24: " + ", ".join(f"{s}: {int(((Fm.lead == 24) & Fm.ok & (Fm.season == s)).sum())}" for s in SEASONS))
    out("Matched lead 24 fixes without a part 2 row (missing file or no fix pulled), by season: " + ", ".join(
        f"{s}: {int(((Fm.lead == 24) & Fm.B.isna() & Fm.tcwv.isna() & (Fm.season == s)).sum())}" for s in SEASONS))
    m24 = ((Fm.lead == 24) & Fm.ok).values
    y24 = Fm.y[m24]
    out("Truth classes, scored sample lead 24 (below gale, gale, storm, HF): " + str(np.bincount(y24, minlength=4).tolist()))
    out(f"nosst (no ocean cell within 500 km): GFS {float(Fm.nosst[m24].mean()):.3f}, ERA5 {float(Ee.nosst[m24].mean()):.3f}")
    out("Descriptive means on the scored sample, lead 24 (GFS-fed fix vs matched ERA5 fix; NOT a validation): " + "; ".join(
        f"{c} {float(np.nanmean(Fm[c][m24])):.3f}/{float(np.nanmean(Ee[c][m24])):.3f}" for c in ALLENV))
    # ---- models
    Fm["use_dp"] = ((Fm.lead >= 12) & Fm.dp12.notna()).values
    young_e = E.young.reindex(ei).values == 1
    use_e = Fm.use_dp.values & ~young_e
    XE = Ee.copy()
    XE.loc[young_e, "dp12"] = np.nan
    for c in ("pac", "doy_c", "doy_s"):
        XE[c] = Fm[c].values
    XG = Fm.copy()          # GFS lead 0 has no flux: left NaN -> Prep imputes the ERA5 training median (documented, secondary only)
    C = np.array([clim[(b, m)] for b, m in zip(Fm.basin, Fm.month)])
    P = {"clim": C}
    for k in ("S", "A", "B", "FULL", "FULL_NOSF"):
        P[k + "_g"] = SO.sm(SO.run_model(allm, k, XG, Fm.use_dp.values))
        P[k + "_e"] = SO.sm(SO.run_model(allm, k, XE, use_e))
    y = Fm.y.values
    seas = Fm.season.values
    L = {k: BS.rps(P[k], y, 4) for k in P}
    B5 = SO.draws(5)
    clim_l = L["clim"]
    out(f"ERA5-fed dp12 model used for {int(use_e[m24].sum())} of {int(m24.sum())} scored fixes; GFS-fed {int(Fm.use_dp.values[m24].sum())}")

    def RP(k, m):
        return SO.rpss(L[k][m], clim_l[m], seas[m], SEASONS, B5)

    def GN(a, b, m):
        return SO.gain(L[a][m], L[b][m], clim_l[m], seas[m], SEASONS, B5)

    def RT(x, m):                                    # retention of X over S
        return ratio_stat(L[x + "_g"][m], L[x + "_e"][m], L["S_g"][m], L["S_e"][m], clim_l[m], seas[m], B5)

    # ---- registered family
    fam = []
    f1 = RP("FULL_g", m24)
    fam.append(("F1", "FULL fed GFS, lead 24, RPSS vs clim", f1, SO.p_le(f1[3]), int(m24.sum())))
    f2 = GN("FULL_g", "S_g", m24)
    fam.append(("F2", "FULL minus S, GFS-fed, lead 24, RPSS gain", f2, SO.p_le(f2[3]), int(m24.sum())))
    f3 = SO.stat([(1, L["FULL_e"][m24]), (-1, L["FULL_g"][m24])], clim_l[m24], seas[m24], SEASONS, B5)
    p3 = max(SO.p_le(f3[3], -0.05), float((1 + np.sum(f3[3] >= 0.05)) / (1 + len(f3[3]))))
    fam.append(("F3", "FULL lead 24, GFS-fed minus ERA5-fed RPSS (equivalence +-0.05)", f3, p3, int(m24.sum())))
    rt = {}
    for tid, x in (("F4A", "A"), ("F4B", "B"), ("F4F", "FULL")):
        r = RT(x, m24)
        rt[x] = r
        fam.append((tid, f"gain retention {x}: (gain over S, GFS-fed) / (gain over S, ERA5-fed), lead 24", r, SO.p_le(r[3], 0.5), int(m24.sum())))
    q = BS.bh([f[3] for f in fam])
    out("\n== Registered family F1-F4F (RPSS or ratio; 90 % season-block interval, 5 seasons, 2,000 draws; BH q within the 6)")
    out("   (p for F4: one-sided bootstrap p for 'ratio <= 0.5'; the F4 verdict uses the registered interval rule, not q)")
    rows = []
    for (tid, lab, t, p, n), qq in zip(fam, q):
        if tid == "F3":
            v = "equivalent" if (t[1] > -0.05 and t[2] < 0.05) else ("degrades" if t[2] < -0.05 else "cannot tell")
        elif tid.startswith("F4"):
            v = verdict_ratio(t)
        else:
            v = "yes" if (qq < 0.05 and t[1] > 0) else "cannot tell"
        extra = f" (draws with non-positive ERA5-fed gain: {int(t[4])})" if tid.startswith("F4") else ""
        fm = (lambda z: f"{z[0]:+.4f} [{z[1]:+.4f}, {z[2]:+.4f}]") if not tid.startswith("F4") else (lambda z: f"{z[0]:.3f} [{z[1]:.3f}, {z[2]:.3f}]")
        out(f"  {tid} {lab}: {fm(t)}  n {n}  p {p:.4f}  q {qq:.4f}  -> {v}{extra}")
        rows.append(dict(test=tid, label=lab, point=float(t[0]), lo90=float(t[1]), hi90=float(t[2]), n=n, p=p, q=float(qq), verdict=v))
    res = dict(family=rows)
    # ---- gains themselves
    out("\n== Secondary (labelled; not in the family; no multiplicity control). Lead 24, scored sample, same fixes for every number")
    sec = res["secondary"] = {}
    out("  RPSS vs basin-month climatology, GFS-fed | ERA5-fed | GFS minus ERA5:")
    sec["rpss"] = {}
    for k in ("S", "A", "B", "FULL", "FULL_NOSF"):
        g, e = RP(k + "_g", m24), RP(k + "_e", m24)
        d = SO.stat([(1, L[k + "_e"][m24]), (-1, L[k + "_g"][m24])], clim_l[m24], seas[m24], SEASONS, B5)
        out(f"    {k:10s} {SO.fmt(g)} | {SO.fmt(e)} | {SO.fmt(d)}")
        sec["rpss"][k] = dict(gfs=[float(x) for x in g[:3]], era5=[float(x) for x in e[:3]], diff=[float(x) for x in d[:3]])
    out("  Gain over S (RPSS difference), GFS-fed | ERA5-fed | retention ratio (verdict by the F4 rule):")
    sec["gain_over_S"] = {}
    for k in ("A", "B", "FULL", "FULL_NOSF"):
        g, e = GN(k + "_g", "S_g", m24), GN(k + "_e", "S_e", m24)
        r = RT(k, m24)
        out(f"    {k:10s} {SO.fmt(g)} | {SO.fmt(e)} | ratio {r[0]:.3f} [{r[1]:.3f}, {r[2]:.3f}] -> {verdict_ratio(r)}")
        sec["gain_over_S"][k] = dict(gfs=[float(x) for x in g[:3]], era5=[float(x) for x in e[:3]], ratio=[float(x) for x in r[:3]], verdict=verdict_ratio(r))
    out("  FULL_NOSF: does dropping sst, sstgrad, sst_t500, flux, nosst cost anything?  (FULL minus FULL_NOSF RPSS)")
    sec["full_minus_nosf"] = {}
    for tag in ("g", "e"):
        t = GN("FULL_" + tag, "FULL_NOSF_" + tag, m24)
        out(f"    {'GFS-fed' if tag == 'g' else 'ERA5-fed'}: {SO.fmt(t)}")
        sec["full_minus_nosf"][tag] = [float(x) for x in t[:3]]
    # ---- basins
    sec["by_basin"] = {}
    for b in ("atl", "pac"):
        mb = m24 & (Fm.basin == b).values
        g1 = RP("FULL_g", mb)
        g2 = GN("FULL_g", "S_g", mb)
        e2 = GN("FULL_e", "S_e", mb)
        r = RT("FULL", mb)
        out(f"  {b} (n {int(mb.sum())}): FULL GFS-fed RPSS {SO.fmt(g1)}; FULL minus S GFS-fed {SO.fmt(g2)}, ERA5-fed {SO.fmt(e2)}, retention {r[0]:.3f} [{r[1]:.3f}, {r[2]:.3f}] -> {verdict_ratio(r)}")
        sec["by_basin"][b] = dict(n=int(mb.sum()), rpss_full_gfs=[float(x) for x in g1[:3]], gain_gfs=[float(x) for x in g2[:3]],
                                  gain_era5=[float(x) for x in e2[:3]], ratio=[float(x) for x in r[:3]])
    # ---- lead 0 (flux missing on GFS: imputed with the ERA5 training median by Prep)
    m0 = ((Fm.lead == 0) & Fm.ok).values
    out(f"\n  Lead 0 (n {int(m0.sum())}; GFS-fed FULL has flux = ERA5 training median, no GFS flux exists at lead 0; ERA5-fed FULL uses the ERA5 flux of the matched fix; no-dp12 fallback throughout):")
    sec["lead0"] = {}
    for k in ("S", "A", "B", "FULL", "FULL_NOSF"):
        g, e = RP(k + "_g", m0), RP(k + "_e", m0)
        out(f"    {k:10s} RPSS GFS-fed {SO.fmt(g)} | ERA5-fed {SO.fmt(e)}")
        sec["lead0"][k] = dict(gfs=[float(x) for x in g[:3]], era5=[float(x) for x in e[:3]])
    for k in ("A", "B", "FULL", "FULL_NOSF"):
        g, e = GN(k + "_g", "S_g", m0), GN(k + "_e", "S_e", m0)
        r = RT(k, m0)
        out(f"    {k:10s} gain over S GFS-fed {SO.fmt(g)} | ERA5-fed {SO.fmt(e)} | ratio {r[0]:.3f} [{r[1]:.3f}, {r[2]:.3f}] -> {verdict_ratio(r)}")
        sec["lead0"][k + "_gain"] = dict(gfs=[float(x) for x in g[:3]], era5=[float(x) for x in e[:3]], ratio=[float(x) for x in r[:3]])
    # ---- part 1 reference: A/B gains on all part-1 env_ok lead 24 fixes
    m1 = ((Fm.lead == 24) & Fm[["jet250", "vadv500", "eady"]].notna().all(axis=1)).values
    out(f"\n  Reference (part 1 sample rule, lead 24 fixes with jet250/vadv500/eady, n {int(m1.sum())}): A minus S GFS-fed {SO.fmt(GN('A_g', 'S_g', m1))}; B minus S {SO.fmt(GN('B_g', 'S_g', m1))}")
    # ---- reliability of FULL HF probability
    for tag, lab in (("g", "GFS-fed"), ("e", "ERA5-fed")):
        rows_, slope = SO.rel_table(P["FULL_" + tag][m24][:, 3], (y[m24] == 3).astype(float))
        out(f"  Reliability of the {lab} FULL hurricane-force probability, lead 24 (deciles of forecast), slope {slope:.2f}, base rate {float((y[m24] == 3).mean()):.4f}, mean forecast {float(P['FULL_' + tag][m24][:, 3].mean()):.4f}:")
        for rw in rows_:
            out(f"    bin {rw['bin']}: n {rw['n']:5d}  forecast {rw['forecast']:.4f}  observed {rw['observed']:.4f}")
        sec["reliability_HF_FULL_" + tag] = dict(slope=slope, bins=rows_)
    out("  Lead 24 mean forecast class probability vs observed frequency (FULL GFS-fed): " + "; ".join(
        f"{BS.NAMES[c]} {P['FULL_g'][m24][:, c].mean():.3f} vs {float((y[m24] == c).mean()):.3f}" for c in range(4)))
    res["counts"] = dict(track_rows=int(n_all), fixes_domain=int(len(F)), matched=int(len(Fm)),
                         lead24=dict(matched=int((Fm.lead == 24).sum()), gfs_ok=int(((Fm.lead == 24) & Fm.gfs_ok).sum()),
                                     era_ok=int(((Fm.lead == 24) & Fm.era_ok).sum()), scored=int(m24.sum())),
                         lead0=dict(matched=int((Fm.lead == 0).sum()), scored=int(m0.sum())), seed=SEEDV, nboot=SO.NB)
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, "outlook_full.txt"), "w").write("\n".join(LINES) + "\n")
    json.dump(res, open(os.path.join(OUT, "outlook_full.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
