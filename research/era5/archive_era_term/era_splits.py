"""RA-29 T-B..T-E: era splits (before / after 2009-11-23) of archive results that use fixes per event. Plan: PREREGISTRATION.md.
Controls reproduce the committed overall numbers before any split is read. Bootstrap unit = season-by-era segment (season 2009 is cut
in two at the step); 2000 resamples, percentile 95%, p = two-sided bootstrap sign p. BH over the registered contrast list.
usage: era_splits.py   (writes results/splits.txt, splits_tests.csv)"""
import json, os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.abspath(os.path.join(ERA, "..", ".."))
STEP = pd.Timestamp("2009-11-23")
HF = 71.7
B = 2000
rng = np.random.default_rng(2929)
lines, TESTS = [], []
def w(s=""):
    lines.append(s); print(s, flush=True)
def bh(p):
    p = np.asarray(p, float); n = len(p); od = np.argsort(p); q = np.empty(n)
    r = p[od] * n / (np.arange(n) + 1); q[od] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1); return q


def seg_of(season, t):
    """segment id: (season, era). era = 1 when t >= STEP."""
    era = (t >= STEP).astype(int)
    return pd.Series(list(zip(season, era)), index=t.index), era


def boot_units(df, unit, stat, era_col="era", n=B):
    """Resample whole segments within each era, recompute stat(df_resampled); returns (point, array)."""
    groups = {e: {u: g for u, g in df[df[era_col] == e].groupby(unit)} for e in (0, 1)}
    keys = {e: list(groups[e]) for e in (0, 1)}
    out = []
    for _ in range(n):
        parts = []
        for e in (0, 1):
            pick = rng.integers(0, len(keys[e]), len(keys[e]))
            parts += [groups[e][keys[e][i]] for i in pick]
        out.append(stat(pd.concat(parts)))
    return np.array(out)


def rec(name, est, bt, extra=""):
    lo, hi = np.nanpercentile(bt, [2.5, 97.5])
    p = min(1.0, 2 * (1 + min((bt > 0).sum(), (bt < 0).sum())) / (1 + len(bt)))
    TESTS.append(dict(test=name, est=est, lo=lo, hi=hi, p=p))
    w(f"  {name}: {est:+.4f} [{lo:+.4f}, {hi:+.4f}]  p {p:.4f} {extra}")


# ============================================================ T-B  count reconciliation (archive events by fixes, per year of exposure)
w("== T-B  Archive events by HF fixes per event, by era (count_reconcile definitions: hf-lows.json, duplicate timestamps dropped, Category HF, seasons 2004-05..2025-26) ==")
d = json.load(open(os.path.join(ROOT, "docs/data/hf-lows.json")))
F = d["lowFields"]
R = []
for raw in d["lows"]:
    l = dict(zip(F, raw))
    f = pd.DataFrame(l["fixes"], columns=["date", "lat", "lon", "cat", "pres"])
    f["t"] = pd.to_datetime(f.date.astype("int64").astype(str), format="%Y%m%d%H")
    f = f.drop_duplicates("t").sort_values("t")
    h = f[f.cat == "HF"]
    if len(h) and 2004 <= l["season"] <= 2025:
        R.append(dict(basin=l["basin"], season=l["season"], n_hf=len(h), t_on=h.t.iloc[0]))
E = pd.DataFrame(R)
E["unit"], E["era"] = seg_of(E.season, E.t_on)
# exposure in years per segment: season s runs 1 Jun s .. 31 May s+1
def expo(s, e):
    a, b = pd.Timestamp(f"{s}-06-01"), pd.Timestamp(f"{s+1}-06-01")
    if s == 2009:
        return ((STEP - a).days if e == 0 else (b - STEP).days) / 365.25
    return 1.0 if (e == (1 if s > 2009 else 0)) else 0.0
rate_tab = {}
for bsn, lab in (("atl", "Atlantic"), ("pac", "Pacific")):
    X = E[E.basin == bsn]
    ex = {(s, e): expo(s, e) for s in range(2004, 2026) for e in (0, 1)}
    ex = {k: v for k, v in ex.items() if v > 0}
    # control: overall per season
    tot = {k: (X.n_hf >= k).sum() / 22 for k in (1, 2, 3)}
    w(f" {lab} control, events per season over 22 seasons: >=1 {tot[1]:.1f}, >=2 {tot[2]:.1f}, >=3 {tot[3]:.1f}  (count_reconcile.txt: Atl 46.0 / 33.8 / 22.4, Pac 38.9 / 29.1 / 20.5)")
    def per_era(df, k):
        out = []
        for e in (0, 1):
            ev = (df[df.era == e].n_hf >= k).sum()
            units = {u for u in df.unit.unique()} | set()
            out.append(ev)
        return out
    # rates need exposure of the *resampled* units: build unit table incl. empty units
    U = pd.DataFrame([dict(unit=k, era=k[1], expo=v) for k, v in ex.items()])
    cnt = {k: X.assign(ge=X.n_hf >= k).groupby("unit").ge.sum() for k in (1, 2, 3)}
    for k in (1, 2, 3):
        U[f"ge{k}"] = U.unit.map(cnt[k]).fillna(0)
    def stat_rate(df, k, e):
        s = df[df.era == e]
        return s[f"ge{k}"].sum() / s.expo.sum()
    point = {(k, e): stat_rate(U, k, e) for k in (1, 2, 3) for e in (0, 1)}
    groups = {e: [g for _, g in U[U.era == e].groupby("unit")] for e in (0, 1)}
    bt = {k: [] for k in ("r1", "r2", "r3", "s2", "s3")}
    for _ in range(B):
        parts = {}
        for e in (0, 1):
            pick = rng.integers(0, len(groups[e]), len(groups[e]))
            parts[e] = pd.concat([groups[e][i] for i in pick])
        rr = {(k, e): parts[e][f"ge{k}"].sum() / parts[e].expo.sum() for k in (1, 2, 3) for e in (0, 1)}
        bt["r1"].append(rr[(1, 1)] / rr[(1, 0)]); bt["r2"].append(rr[(2, 1)] / rr[(2, 0)]); bt["r3"].append(rr[(3, 1)] / rr[(3, 0)])
        bt["s2"].append(rr[(2, 1)] / rr[(1, 1)] - rr[(2, 0)] / rr[(1, 0)]); bt["s3"].append(rr[(3, 1)] / rr[(1, 1)] - rr[(3, 0)] / rr[(1, 0)])
    w(f" {lab}: events per year, before (5.48 yr) -> after (16.52 yr)")
    for k in (1, 2, 3):
        lo, hi = np.percentile(bt[f"r{k}"], [2.5, 97.5])
        w(f"   >= {k} fixes: {point[(k,0)]:.1f} -> {point[(k,1)]:.1f}  ratio {point[(k,1)]/point[(k,0)]:.3f} [{lo:.3f}, {hi:.3f}]")
    for k, key in ((2, "s2"), (3, "s3")):
        est = point[(k, 1)] / point[(1, 1)] - point[(k, 0)] / point[(1, 0)]
        w(f"   share with >= {k} fixes: {point[(k,0)]/point[(1,0)]:.3f} -> {point[(k,1)]/point[(1,1)]:.3f}")
        bt_ = np.array(bt[key])
        if k == 3:
            rec(f"TB {lab} share >=3 fixes, after minus before", est, bt_)
        else:
            lo, hi = np.percentile(bt_, [2.5, 97.5]); w(f"   (share >=2 after minus before {est:+.4f} [{lo:+.4f},{hi:+.4f}], descriptive)")
    rate_tab[bsn] = point

# ============================================================ T-C  atlas hours at HF per event
w("\n== T-C  Atlas: hours at HF per event (6 h x number of HF fixes), archive vs ERA5 pipeline A proxy, by era (event first HF fix before / after 2009-11-23) ==")
A = pd.read_csv(os.path.join(ERA, "climatology_atlas/results/events_archive.csv"), parse_dates=["t_on"])
P = pd.read_csv(os.path.join(ERA, "climatology_atlas/results/events_proxy.csv"), parse_dates=["t_on"])
A["src"], P["src"] = "arch", "prox"
for X in (A, P):
    X["unit"], X["era"] = seg_of(X.season, X.t_on)
for bsn, lab in (("atl", "Atlantic"), ("pac", "Pacific")):
    a, p = A[A.basin == bsn], P[P.basin == bsn]
    w(f" {lab} control (2004-05..2025-26, all events in the atlas tables): archive median {a.hf_h.median():.0f} h mean {a.hf_h.mean():.2f}; proxy median {p.hf_h.median():.0f} h mean {p.hf_h.mean():.2f}  (atlas README: archive medians 12 Atl, 18 Pac; proxy 12)")
    for e, nm in ((0, "before"), (1, "after")):
        w(f"   {nm}: archive n {int((a.era==e).sum())} mean {a[a.era==e].hf_h.mean():.2f} median {a[a.era==e].hf_h.median():.0f} | proxy n {int((p.era==e).sum())} mean {p[p.era==e].hf_h.mean():.2f} median {p[p.era==e].hf_h.median():.0f}")
    both = pd.concat([a[["unit", "era", "hf_h", "src"]], p[["unit", "era", "hf_h", "src"]]])
    def st(df, which):
        m = lambda s, e: df[(df.src == s) & (df.era == e)].hf_h.mean()
        if which == "arch": return m("arch", 1) - m("arch", 0)
        if which == "prox": return m("prox", 1) - m("prox", 0)
        return (m("arch", 1) - m("prox", 1)) - (m("arch", 0) - m("prox", 0))
    for which, lab2 in (("arch", "archive mean hours, after minus before"), ("prox", "proxy mean hours, after minus before"), ("gap", "archive minus proxy gap, after minus before")):
        est = st(both, which)
        bt = boot_units(both, "unit", lambda df: st(df, which))
        if which == "prox":
            lo, hi = np.percentile(bt, [2.5, 97.5]); w(f"  {lab} {lab2}: {est:+.3f} [{lo:+.3f},{hi:+.3f}] (descriptive)")
        else:
            rec(f"TC {lab} {lab2}", est, bt)
    # median archive hours before/after as stated in the atlas
    w(f"   share of events with a single HF fix (6 h): archive before {np.mean(a[a.era==0].hf_h==6):.3f} after {np.mean(a[a.era==1].hf_h==6):.3f}; proxy before {np.mean(p[p.era==0].hf_h==6):.3f} after {np.mean(p[p.era==1].hf_h==6):.3f}")

# ============================================================ T-E  probability tracks: archive-listed share of proxy events
w("\n== T-E  P(HF) tracks: share of pipeline A proxy events (gust >= 71.7 kt, 2004-05..2025-26) matched to an archive HF event (gust_depth/matches.csv), by era of the track's gust peak ==")
At = pd.read_csv(os.path.join(ERA, "hf_history/results/all_tracks.csv.gz"), dtype={"peak_time": str, "start": str})
At = At[(At.season >= 2004) & (At.season <= 2025)].copy()
mt = pd.read_csv(os.path.join(ERA, "gust_depth/results/matches.csv"))
At["listed"] = At.track.isin(set(mt.track.dropna().astype(int)))
At = At[At.gust800_kt >= HF].copy()
At["t"] = pd.to_datetime(At.peak_time, format="%Y%m%d%H", errors="coerce")
At["unit"], At["era"] = seg_of(At.season, At.t)
w(f" control: proxy events {len(At)}, listed {At.listed.mean():.3f}  (README: 61% peak-position rule over 2004-2025)")
for bsn, lab in (("atl", "Atlantic"), ("pac", "Pacific")):
    x = At[At.basin == bsn]
    w(f"  {lab}: before {x[x.era==0].listed.mean():.3f} (n {int((x.era==0).sum())}), after {x[x.era==1].listed.mean():.3f} (n {int((x.era==1).sum())})")
    st = lambda df: df[df.era == 1].listed.mean() - df[df.era == 0].listed.mean()
    rec(f"TE {lab} listed share, after minus before", st(x), boot_units(x, "unit", st))

# ============================================================ T-D  skill limits: ceiling HSS vs archive label by era
w("\n== T-D  HF skill limits: ERA5 hf24 label scored against the archive label 'any archive HF fix in the next 24 h' (matching as hf_skill_limits/analysis.py Q5), by era ==")
fx = pd.read_csv(os.path.join(ERA, "hf_probability_tracks/results/fix_probs_2004.csv.gz"), dtype={"time": str})
fx["dt"] = pd.to_datetime(fx.time, format="%Y%m%d%H")
Aa = []
for b, f in (("atl", "HF_Data_-_Atl.csv"), ("pac", "HF_Data_-_Pac.csv")):
    a = pd.read_csv(os.path.join(ROOT, "data/hf_lows", f)); a["basin"] = b; Aa.append(a)
Aa = pd.concat(Aa)
Aa["dt"] = pd.to_datetime(Aa.date.astype(str), format="%Y%m%d%H", errors="coerce")
Aa = Aa[(Aa.Category == "HF") & Aa.dt.notna()]
Aa["season"] = np.where(Aa.dt.dt.month >= 6, Aa.dt.dt.year, Aa.dt.dt.year - 1)
Aa = Aa[(Aa.season >= 2004) & (Aa.season <= 2025)].reset_index(drop=True)
Aa["lon360"] = Aa.Longitude % 360
def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2); dl = np.radians(lon2 - lon1)
    a_ = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 3440.065 * 1.852 * np.arcsin(np.sqrt(a_))
byt = {t: g for t, g in fx.groupby("time")}
mtr = []
for r in Aa.itertuples():
    best = (None, 1e9)
    for dh, lim in ((0, 600),) if r.dt.hour in (0, 12) else ((-6, 800), (6, 800)):
        t = (r.dt + pd.Timedelta(hours=dh)).strftime("%Y%m%d%H")
        if t in byt:
            g = byt[t]; dd = hav(r.Latitude, r.lon360, g.lat.values, g.lon.values); j = dd.argmin()
            if dd[j] < lim and dd[j] < best[1]:
                best = (g.track.values[j], dd[j])
    mtr.append(best[0])
Aa["track"] = mtr
Am = Aa[Aa.track.notna()].copy(); Am["track"] = Am.track.astype(int)
arch_by_track = {t: np.sort(g.dt.values) for t, g in Am.groupby("track")}
def arch24(tr, dt_):
    ts = arch_by_track.get(tr)
    if ts is None: return 0
    return int(((ts > np.datetime64(dt_)) & (ts <= np.datetime64(dt_ + pd.Timedelta(hours=24)))).any())
fx["arch24"] = [arch24(a, b) for a, b in zip(fx.track, fx.dt)]
fx["fc"] = fx.P_F >= fx.cut_F
fx["era"] = (fx.dt >= STEP).astype(int)
fx["unit"] = list(zip(fx.season, fx.era))
def hss(df, fc):
    fc = np.asarray(fc, bool); ob = (df.arch24.values == 1)
    a = (fc & ob).sum(); b_ = (fc & ~ob).sum(); c = (~fc & ob).sum(); dd = (~fc & ~ob).sum()
    den = (a + c) * (c + dd) + (a + b_) * (b_ + dd)
    return 2 * (a * dd - b_ * c) / den if den else np.nan
fx["lab_ceil"] = fx.hf24 == 1
fx["fcm"] = fx.fc
w(f" control: fixes {len(fx)}, ERA5 label positives {int(fx.hf24.sum())}, archive-label positives {int(fx.arch24.sum())}; HSS ceiling (ERA5 label vs archive label) {hss(fx, fx["lab_ceil"]):.3f}, PR 12 model vs archive label {hss(fx, fx.fcm):.3f}  (README: 0.657 and 0.543)")
for nm, col in (("ceiling", "lab_ceil"), ("PR 12 model", "fcm")):
    for e, en in ((0, "before"), (1, "after")):
        s = fx[fx.era == e]
        w(f"   {nm} {en}: fixes {len(s)}, archive-label positives {int(s.arch24.sum())} ({s.arch24.mean():.4f}), HSS {hss(s, s[col]):.3f}")
st = lambda df: hss(df[df.era == 1], df[df.era == 1]["lab_ceil"]) - hss(df[df.era == 0], df[df.era == 0]["lab_ceil"])
rec("TD ceiling HSS vs archive label, after minus before", st(fx), boot_units(fx, "unit", st, n=1000))
st2 = lambda df: hss(df[df.era == 1], df[df.era == 1]["fcm"]) - hss(df[df.era == 0], df[df.era == 0]["fcm"])
w("  PR 12 model vs archive label, after minus before (descriptive):")
bt2 = boot_units(fx, "unit", st2, n=1000); lo, hi = np.nanpercentile(bt2, [2.5, 97.5]); w(f"    {st2(fx):+.4f} [{lo:+.4f}, {hi:+.4f}]")

T = pd.DataFrame(TESTS)
T["q_bh"] = bh(T.p.values)
T.to_csv(os.path.join(HERE, "results", "splits_tests.csv"), index=False, float_format="%.5g")
w("\n== Registered contrasts (T-B share >=3, T-C archive hours and gap, T-D ceiling, T-E listed share), BH over all ==")
for r in T.itertuples():
    w(f"  {r.test}: {r.est:+.4f} [{r.lo:+.4f},{r.hi:+.4f}] p {r.p:.4f} q {r.q_bh:.4f}")
w(f"  {int((T.q_bh<0.05).sum())} of {len(T)} pass q < 0.05")
open(os.path.join(HERE, "results", "splits.txt"), "w").write("\n".join(lines) + "\n")
