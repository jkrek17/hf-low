"""H1 (skill against the archive) and H2 (where gust and depth disagree).  ERA5 PROXY, pipeline A.
Plan: PREREGISTRATION.md (committed before this was run).

usage: diagnose.py ALL_TRACKS_CSV_GZ ARCHIVE_JSON OUT_DIR [NBOOT] [MATCH_KM]
"""
import json, os, sys, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "freq_split"))
import match, split

warnings.filterwarnings("ignore")
THR = 71.7
FIT = list(range(2021, 2026)); TEST = list(range(2004, 2021))
SEED = 20261009


def season_of(start):
    s = start.astype(str)
    y, m = s.str[:4].astype(int), s.str[4:6].astype(int)
    return np.where(m >= 6, y, y - 1)


def load(tracks, archive):
    T = pd.read_csv(tracks, dtype={"start": str, "end": str, "peak_time": str})
    T["sea"] = season_of(T.start)
    T["pk_month"] = T.peak_time.str[4:6].astype(int)
    T["real"] = ((T.minp <= 1000) & (T.n_fix >= 8)).astype(int)
    ev = match.load_archive(archive)
    if os.environ.get("HF_ONLY"):          # sensitivity: archive events whose peak category is HF
        ev = [e for e in ev if e["peak"] == "HF"]
    M = match.match_all(T, ev, 2004)
    return T, ev, M


def auc(score, y):
    r = rankdata(score)
    n1, n0 = y.sum(), (~y).sum()
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else np.nan


def contingency(flag, obs, unmatched):
    a = int((flag & obs).sum()); b = int((flag & ~obs).sum())
    c = int((~flag & obs).sum()) + unmatched; d = int((~flag & ~obs).sum())
    n = a + b + c + d
    e = ((a + b) * (a + c) + (c + d) * (b + d)) / n
    return dict(a=a, b=b, c=c, d=d, pod=a / max(a + c, 1), far=b / max(a + b, 1), csi=a / max(a + b + c, 1),
                hss=(a + d - e) / (n - e), bias=(a + b) / max(a + c, 1), n_fc=a + b, n_obs=a + c)


def fit_thr(x, n_obs, sign):
    """threshold on x (sign +1: x >= thr flags, -1: x <= thr flags) whose flagged count is closest to n_obs."""
    v = np.sort(x)[::-1] if sign > 0 else np.sort(x)
    k = int(min(max(n_obs, 1), len(v)))
    return float(v[k - 1])


class Skill:
    """Everything needed to score one basin, with season-block resampling of the test seasons."""

    def __init__(self, T, M, basin, G_pub=THR):
        self.T = T[(T.basin == basin) & T.sea.between(2004, 2025)].reset_index(drop=True)
        m = M[(M.basin == basin) & (M.cls == "low")]
        self.M = m
        pos = set(m.track.dropna().astype(int))
        self.T["obs"] = self.T.track.isin(pos)
        self.T["pos_cls"] = self.T.track.isin(pos)
        # unmatched archive events by season
        self.unm = m[m.track.isna()].groupby("season").size().to_dict()
        self.nobs_fit = int(((m.season.isin(FIT))).sum())
        fit = self.T[self.T.sea.isin(FIT)]
        self.G_star = fit_thr(fit.gust800_kt.values, self.nobs_fit, +1)
        self.D_star = fit_thr(fit.minp.values, self.nobs_fit, -1)
        self.G_pub = G_pub

    def flags(self, d):
        return {"G71.7": d.gust800_kt >= self.G_pub, "G*": d.gust800_kt >= self.G_star, "D*": d.minp <= self.D_star}

    def stats(self, seasons, real_only=False):
        d = self.T[self.T.sea.isin(seasons)]
        if real_only:
            d = d[d.real == 1]
        out = {}
        unm = sum(self.unm.get(s, 0) for s in seasons)
        for k, f in self.flags(d).items():
            out[k] = contingency(f.values, d.obs.values, unm)
        out["auc_G"] = auc(d.gust800_kt.values, d.obs.values)
        out["auc_D"] = auc(-d.minp.values, d.obs.values)
        return out

    def boot(self, nboot, rng, real_only=False):
        d = self.T[self.T.sea.isin(TEST)]
        if real_only:
            d = d[d.real == 1]
        groups = {s: g for s, g in d.groupby("sea")}
        res = []
        for _ in range(nboot):
            pick = rng.choice(TEST, len(TEST), replace=True)
            dd = pd.concat([groups[s] for s in pick if s in groups])
            unm = sum(self.unm.get(s, 0) for s in pick)
            aG, aD = auc(dd.gust800_kt.values, dd.obs.values), auc(-dd.minp.values, dd.obs.values)
            cG = contingency((dd.gust800_kt >= self.G_star).values, dd.obs.values, unm)
            cD = contingency((dd.minp <= self.D_star).values, dd.obs.values, unm)
            res.append((aG - aD, cG["hss"] - cD["hss"], aG, aD, cG["hss"], cD["hss"]))
        return np.array(res)


def pval(b, obs):
    """two-sided bootstrap p for a difference: 2 * min(P(b<=0), P(b>=0)), centred on sign of the estimate."""
    p = 2 * min((b <= 0).mean(), (b >= 0).mean())
    return float(min(max(p, 1.0 / (len(b) + 1)), 1.0))


def main():
    tracks, archive, out = sys.argv[1:4]
    nboot = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    mkm = float(sys.argv[5]) if len(sys.argv) > 5 else 800
    match.MATCH_KM = mkm
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    T, ev, M = load(tracks, archive)
    M.to_csv(os.path.join(out, "matches.csv"), index=False)
    rows, L = [], []
    say = lambda s="": (print(s), L.append(s))
    say(f"Archive-to-track matching: {mkm:.0f} km; archive events 2004-2025: {len(M)}; matched {M.track.notna().mean():.3f}")
    say(M.groupby(["basin", "cls"]).apply(lambda g: f"  {g.name}: n={len(g)}, matched={g.track.notna().mean():.2f}").to_string())

    # ------------------------------------------------ H1
    sk = {}
    for basin in ("pac", "atl"):
        s = Skill(T, M, basin)
        sk[basin] = s
        say(f"\n=== {basin}: fit seasons {FIT[0]}-{FIT[-1]}, {s.nobs_fit} archive 'low' events; G*={s.G_star:.1f} kt, D*={s.D_star:.1f} hPa (published G {THR})")
        for tag, ro in (("all tracks", False), ("real cyclones", True)):
            st = s.stats(TEST, ro)
            say(f"  test {TEST[0]}-{TEST[-1]}, {tag}:  AUC gust {st['auc_G']:.3f}  AUC depth {st['auc_D']:.3f}")
            for k in ("G71.7", "G*", "D*"):
                c = st[k]
                say(f"    {k:6s} n_obs={c['n_obs']} n_fc={c['n_fc']} bias={c['bias']:.2f} POD={c['pod']:.2f} FAR={c['far']:.2f} CSI={c['csi']:.2f} HSS={c['hss']:.3f}")
            bt = s.boot(nboot, rng, ro)
            for j, nm in ((0, "dAUC"), (1, "dHSS")):
                if ro and nm == "dHSS":
                    continue
                est = (st["auc_G"] - st["auc_D"]) if nm == "dAUC" else (st["G*"]["hss"] - st["D*"]["hss"])
                lo, hi = np.percentile(bt[:, j], [2.5, 97.5])
                rows.append(dict(family="H1", basin=basin, test=f"T1{'b' if ro else ('a' if nm == 'dAUC' else 'c')}", what=f"{nm} G-D ({tag})",
                                 est=est, lo=lo, hi=hi, p=pval(bt[:, j], est)))
                say(f"    {nm} (G - D) = {est:+.3f}  95% [{lo:+.3f}, {hi:+.3f}]  p={pval(bt[:, j], est):.4f}")
        # leave-one-season-out on the test seasons for dAUC, HSS
        lo_a, lo_h = [], []
        for hold in TEST:
            st = s.stats([x for x in TEST if x != hold])
            lo_a.append(st["auc_G"] - st["auc_D"]); lo_h.append(st["G*"]["hss"] - st["D*"]["hss"])
        say(f"  leave-one-test-season-out: dAUC {min(lo_a):+.3f}..{max(lo_a):+.3f}, dHSS {min(lo_h):+.3f}..{max(lo_h):+.3f}")
        # fit-season (in-sample for thresholds) for reference
        st = s.stats(FIT)
        say(f"  fit seasons (thresholds fitted here): AUC gust {st['auc_G']:.3f} depth {st['auc_D']:.3f}; HSS G* {st['G*']['hss']:.3f} D* {st['D*']['hss']:.3f}")
        # tip jets / centreless: recall among matched
        Mo = M[(M.basin == basin) & (M.cls != "low")]
        mt = Mo.dropna(subset=["track"]).merge(T[["track", "gust800_kt", "minp"]], on="track")
        if len(Mo):
            say(f"  non-'low' archive events (tip jet / centreless): {len(Mo)}, matched {len(mt)}; of matched, gust>=71.7: {int((mt.gust800_kt >= THR).sum())}, minp<=D*: {int((mt.minp <= s.D_star).sum())}")

    # ------------------------------------------------ H2 regression
    say("\n=== H2a: gust at fixed depth (OLS, real cyclones, seasons 2004-2025, season-clustered SE)")
    for basin in ("pac", "atl"):
        d = T[(T.basin == basin) & T.sea.between(2004, 2025) & (T.real == 1)].copy()
        west = (d.peak_lon < 180) if basin == "pac" else (d.peak_lon < 315)
        X = pd.DataFrame({"minp": d.minp - 975, "minp2": (d.minp - 975) ** 2, "lat": d.peak_lat - 50,
                          "west": west.astype(float), "cold": d.pk_month.isin([11, 12, 1, 2]).astype(float),
                          "logn": np.log(d.n_fix) - np.log(d.n_fix).mean(), "trend": (d.sea - 2015) / 10, "const": 1.0})
        r = sm.OLS(d.gust800_kt.values, X).fit(cov_type="cluster", cov_kwds={"groups": d.sea.values})
        say(f" {basin}: n={len(d)}, R2={r.rsquared:.3f}")
        for nm, tst, unit in (("lat", "T2a", "kt per degree latitude"), ("west", "T2b", "kt, west sector vs east"),
                              ("cold", "T2c", "kt, Nov-Feb vs other months"), ("logn", "T2d", "kt per unit log(n_fix)")):
            b, se, p = r.params[nm], r.bse[nm], r.pvalues[nm]
            say(f"    {nm:5s} {b:+.3f} kt (se {se:.3f}, p {p:.2g})  [{unit}]")
            rows.append(dict(family="H2", basin=basin, test=tst, what=f"gust at fixed depth: {nm} ({unit})", est=b, lo=b - 1.96 * se, hi=b + 1.96 * se, p=p))
        say(f"    minp slope at 975 hPa: {r.params['minp']:+.3f} kt per hPa")

    # ------------------------------------------------ H2 groups
    say("\n=== H2b: G-only vs D-only (published thresholds, genesis Oct-Apr cut as in freq_split; full-year tracks 2004-2025)")
    T["hf"] = (T.gust800_kt >= THR).astype(int)
    T["gen"] = pd.to_datetime(T.start.str[:8], format="%Y%m%d")
    matched = set(M.track.dropna().astype(int))
    grp_rows = []
    for basin in ("pac", "atl"):
        cut, n_hf, n_cut = split.depth_cut(T, basin)
        d = T[(T.basin == basin) & T.sea.between(2004, 2025)].copy()
        west = (d.peak_lon < 180) if basin == "pac" else (d.peak_lon < 315)
        d["west"] = west.astype(float); d["cold"] = d.pk_month.isin([11, 12, 1, 2]).astype(float)
        d["G"] = d.gust800_kt >= THR; d["D"] = d.minp <= cut
        d["arch"] = d.track.isin(matched).astype(float)
        d["grp"] = np.select([d.G & d.D, d.G & ~d.D, ~d.G & d.D], ["both", "G-only", "D-only"], "neither")
        say(f" {basin}: published depth cut {cut:.1f} hPa; counts: " + ", ".join(f"{k}={v}" for k, v in d.grp.value_counts().items()))
        for gname, g in d.groupby("grp"):
            if gname == "neither":
                continue
            say(f"    {gname:7s} n={len(g):4d}  mean lat {g.peak_lat.mean():5.1f}  west {g.west.mean():.2f}  cold {g.cold.mean():.2f}  median n_fix {g.n_fix.median():.0f}  archive-matched {g.arch.mean():.2f}  mean gust {g.gust800_kt.mean():.1f} kt  mean minp {g.minp.mean():.1f}")
            grp_rows.append(dict(basin=basin, grp=gname, n=len(g), lat=g.peak_lat.mean(), west=g.west.mean(), cold=g.cold.mean(),
                                 nfix=g.n_fix.median(), arch=g.arch.mean(), gust=g.gust800_kt.mean(), minp=g.minp.mean()))
        seas = np.arange(2004, 2026)
        by = {s: d[d.sea == s] for s in seas}

        def diffs(pick):
            x = pd.concat([by[s] for s in pick])
            g, dd = x[x.grp == "G-only"], x[x.grp == "D-only"]
            return np.array([g.peak_lat.mean() - dd.peak_lat.mean(), g.west.mean() - dd.west.mean(), g.cold.mean() - dd.cold.mean(),
                             g.n_fix.median() - dd.n_fix.median(), g.arch.mean() - dd.arch.mean()])
        est = diffs(seas)
        bs = np.array([diffs(rng.choice(seas, len(seas), replace=True)) for _ in range(nboot)])
        for j, (nm, tst) in enumerate((("peak latitude (deg)", "T2e"), ("west-sector share", "T2f"), ("cold-season share", "T2g"),
                                       ("median n_fix", "T2h"), ("archive-matched share", "T2i"))):
            lo, hi = np.nanpercentile(bs[:, j], [2.5, 97.5])
            p = pval(bs[~np.isnan(bs[:, j]), j], est[j])
            say(f"    G-only minus D-only, {nm}: {est[j]:+.3f}  [{lo:+.3f}, {hi:+.3f}]  p={p:.4f}")
            rows.append(dict(family="H2", basin=basin, test=tst, what=f"G-only minus D-only: {nm}", est=est[j], lo=lo, hi=hi, p=p))
        d[["track", "basin", "sea", "grp", "gust800_kt", "minp", "peak_lat", "peak_lon", "pk_month", "n_fix", "arch"]].to_csv(
            os.path.join(out, f"groups_{basin}.csv"), index=False)
    pd.DataFrame(rows).to_csv(os.path.join(out, "diagnose_tests.csv"), index=False, float_format="%.5g")
    pd.DataFrame(grp_rows).to_csv(os.path.join(out, "diagnose_groups.csv"), index=False, float_format="%.4g")
    open(os.path.join(out, "diagnose.txt"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
