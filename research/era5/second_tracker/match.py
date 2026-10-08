"""Match tracker M and V tracks to pipeline A HF events and compute the registered tests (PREREGISTRATION.md).
Pipeline A = research/era5/hf_history (a proxy). Trackers M, V = common.py. Seasons 2004-05..2021-22.
usage: match.py WORK_DIR OUT_DIR
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
import common as C

WORK, OUT = sys.argv[1], sys.argv[2]
HFH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_history", "results")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261008)
NB = 10000
SEAS = list(range(2004, 2022))
DKM = {"M": 500.0, "V": 700.0}
LOG = open(os.path.join(OUT, "match.txt"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


# ---------------------------------------------------------------- pipeline A inputs
cat = pd.read_csv(f"{HFH}/era5_hf_catalog.csv")
ev = cat[(cat.role == "event") & cat.season.between(SEAS[0], SEAS[-1])].copy()
fx = pd.read_csv(f"{HFH}/era5_hf_catalog_tracks.csv")
fx = fx[fx.track.isin(ev.track)]
life = pd.read_csv(f"{HFH}/lifecycle_events.csv")[["track", "tc"]]
allt = pd.read_csv(f"{HFH}/all_tracks.csv.gz")
allt = allt[allt.season.between(SEAS[0], SEAS[-1])].copy()

HF = fx[fx.g800 >= 71.7].copy()
first = HF.groupby("track").time.min().rename("first_hf")
nhf = HF.groupby("track").size().rename("n_hf")
ev = ev.merge(first, on="track").merge(nhf, on="track").merge(life, on="track", how="left")


def md(t):
    return (t // 10000) % 100 * 100 + (t // 100) % 100


ev["md"] = ev.first_hf.map(md)
ev["octapr"] = (ev.md >= 1005) | (ev.md <= 425)
say(f"A HF events 2004-05..2021-22: {len(ev)} (atl {int((ev.basin=='atl').sum())}, pac {int((ev.basin=='pac').sum())}); "
    f"first HF fix 5 Oct-25 Apr: {int(ev.octapr.sum())}")


def basin_of(lat, lon):
    return np.where((lat >= 30) & (lat <= 67) & ((lon >= 262) | (lon <= 10)), "atl",
                    np.where((lat >= 27) & (lat <= 67) & (lon >= 135) & (lon <= 240), "pac", ""))


# ---------------------------------------------------------------- tracker inputs
TR = {}
for t in ("M", "V"):
    d = pd.read_csv(f"{WORK}/out/{t}_fixes.csv.gz")
    TR[t] = d
    say(f"tracker {t}: {d.track.nunique()} tracks, {len(d)} fixes")
BYTIME = {t: {k: g[["track", "lat", "lon"]].values for k, g in TR[t].groupby("time")} for t in TR}


def match_event(t, fixes, D, need=0.5):
    """fixes: rows (time, lat, lon) of the A event to match. Returns (matched, best_track, n_hit)."""
    hits = {}
    for tm, la, lo in fixes:
        a = BYTIME[t].get(int(tm))
        if a is None:
            continue
        d = C.dist(la, lo, a[:, 1], a[:, 2])
        for k in a[d <= D, 0]:
            hits[int(k)] = hits.get(int(k), 0) + 1
    if not hits:
        return False, -1, 0
    k = max(hits, key=hits.get)
    return hits[k] >= max(1, need * len(fixes)), k, hits[k]


HFX = {k: g[["time", "lat", "lon"]].values for k, g in HF.groupby("track")}
ALLX = {k: g[["time", "lat", "lon"]].values for k, g in fx.groupby("track")}


def recall_table(t, D, rule="hf"):
    rows = []
    src = HFX if rule == "hf" else ALLX
    for r in ev.itertuples():
        m, k, n = match_event(t, src[r.track], D)
        rows.append((r.track, m, k, n))
    return pd.DataFrame(rows, columns=["track", "matched", "tr_track", "n_hit"])


def season_boot(df, col="matched"):
    """pooled recall by basin and its season-block bootstrap draws."""
    out = {}
    cnt = df.groupby(["season", "basin"]).agg(n=(col, "size"), k=(col, "sum")).reset_index()
    for b in ("atl", "pac"):
        c = cnt[cnt.basin == b].set_index("season").reindex(SEAS).fillna(0)
        n, k = c.n.values, c.k.values
        pt = k.sum() / n.sum()
        idx = rng.integers(0, len(SEAS), (NB, len(SEAS)))
        draws = k[idx].sum(1) / np.maximum(n[idx].sum(1), 1)
        out[b] = (pt, draws, int(n.sum()), int(k.sum()))
    return out


def ci(draws, level):
    a = (1 - level) / 2
    return np.quantile(draws, [a, 1 - a])


results = {}
# ---------------------------------------------------------------- primary and S4
say("\n== PRIMARY: recall of A HF events (first HF fix 5 Oct-25 Apr), 98.75% season-block interval (Bonferroni over 4) ==")
prim = []
rec = {}
for t in ("M", "V"):
    r = recall_table(t, DKM[t])
    rec[(t, "base")] = r
    e = ev.merge(r, on="track")
    e.to_csv(f"{OUT}/events_{t}.csv", index=False)
    sub = e[e.octapr]
    bs = season_boot(sub)
    for b in ("atl", "pac"):
        pt, dr, n, k = bs[b]
        lo, hi = ci(dr, 1 - 0.05 / 4)
        verdict = "agree" if lo >= 0.85 else ("disagree" if hi < 0.85 else "cannot tell")
        say(f"{t} {b}: recall {100*pt:5.1f}%  ({k}/{n})  98.75% [{100*lo:5.1f}, {100*hi:5.1f}]  -> {verdict}")
        prim.append(dict(tracker=t, basin=b, n=n, k=k, recall=pt, lo=lo, hi=hi, verdict=verdict))
pd.DataFrame(prim).to_csv(f"{OUT}/primary.csv", index=False)
vs = [p["verdict"] for p in prim]
overall = "yes" if all(v == "agree" for v in vs) else ("no" if "disagree" in vs else "cannot tell")
say(f"OVERALL (registered rule): {overall}")

say("\n== SECONDARY: M all-year (labelled secondary), 95% interval ==")
e = ev.merge(rec[("M", "base")], on="track")
bs = season_boot(e)
for b in ("atl", "pac"):
    pt, dr, n, k = bs[b]
    lo, hi = ci(dr, 0.95)
    say(f"M all-year {b}: {100*pt:5.1f}% ({k}/{n}) [{100*lo:5.1f}, {100*hi:5.1f}]")

say("\n== S4 recall sensitivities (Oct-Apr events, 95% interval) ==")
s4 = []
for t in ("M", "V"):
    variants = [("D/2", dict(D=DKM[t] / 2)), ("2D", dict(D=2 * DKM[t])), ("whole-track", dict(D=DKM[t], rule="all"))]
    for name, kw in variants:
        r = recall_table(t, kw["D"], kw.get("rule", "hf"))
        sub = ev.merge(r, on="track")
        sub = sub[sub.octapr]
        bs = season_boot(sub)
        for b in ("atl", "pac"):
            pt, dr, n, k = bs[b]
            lo, hi = ci(dr, 0.95)
            say(f"{t} {name:11s} {b}: {100*pt:5.1f}% ({k}/{n}) [{100*lo:5.1f}, {100*hi:5.1f}]")
            s4.append(dict(tracker=t, variant=name, basin=b, recall=pt, lo=lo, hi=hi, n=n))
    base = ev.merge(rec[(t, "base")], on="track")
    base = base[base.octapr]
    for mn in (2, 3):
        sub = base[base.n_hf >= mn]
        bs = season_boot(sub)
        for b in ("atl", "pac"):
            pt, dr, n, k = bs[b]
            lo, hi = ci(dr, 0.95)
            say(f"{t} >= {mn} HF fixes {b}: {100*pt:5.1f}% ({k}/{n}) [{100*lo:5.1f}, {100*hi:5.1f}]")
            s4.append(dict(tracker=t, variant=f">={mn} HF fixes", basin=b, recall=pt, lo=lo, hi=hi, n=n))
pd.DataFrame(s4).to_csv(f"{OUT}/s4_sensitivities.csv", index=False)


# ---------------------------------------------------------------- POST HOC control: chance-level recall
say("\n== POST HOC (not registered): chance control. Same rule with A's HF fixes moved 15 degrees east, or 5 days later ==")
import datetime as _dt


def shifted(src, kind):
    out = {}
    for k, a in src.items():
        a = a.copy()
        if kind == "lon":
            a[:, 2] = (a[:, 2] + 15) % 360
        else:
            a[:, 0] = [int((_dt.datetime.strptime(str(int(x)), "%Y%m%d%H") + _dt.timedelta(days=5)).strftime("%Y%m%d%H")) for x in a[:, 0]]
        out[k] = a
    return out


ctl = []
for kind in ("lon", "time"):
    sh = shifted(HFX, kind)
    for t in ("M", "V"):
        rows = [(r.track, match_event(t, sh[r.track], DKM[t])[0]) for r in ev.itertuples()]
        x = ev.merge(pd.DataFrame(rows, columns=["track", "matched"]), on="track")
        x = x[x.octapr]
        for b in ("atl", "pac"):
            y = x[x.basin == b]
            say(f"control shift {kind:4s} {t} {b}: {100*y.matched.mean():5.1f}% ({int(y.matched.sum())}/{len(y)})")
            ctl.append(dict(shift=kind, tracker=t, basin=b, share=y.matched.mean()))
pd.DataFrame(ctl).to_csv(f"{OUT}/posthoc_chance_control.csv", index=False)

# ---------------------------------------------------------------- S1 reverse direction, S2, S3
say("\n== S1 reverse: top-N most intense tracker tracks (N = A HF events in that season and basin) matched to an A HF event ==")


def season_of(t):
    y, m = t // 1000000, (t // 10000) % 100
    return y if m >= 6 else y - 1


AHF = HF.copy()
AHF_BYTIME = {k: g[["lat", "lon"]].values for k, g in AHF.groupby("time")}
nA = ev[ev.octapr].groupby(["season", "basin"]).size()
s1 = []
topsel = {}
for t in ("M", "V"):
    d = TR[t]
    if t == "M":
        pk = d.loc[d.groupby("track").val.idxmin()]
    else:
        pk = d.loc[d.groupby("track").val.idxmax()]
    pk = pk.copy()
    pk["basin"] = basin_of(pk.lat.values, pk.lon.values)
    pk["season"] = pk.time.map(season_of)
    pk["md"] = pk.time.map(md)
    pk = pk[(pk.basin != "") & pk.season.between(SEAS[0], SEAS[-1]) & ((pk.md >= 1005) | (pk.md <= 425))]
    allpk = pk.copy()
    sel = []
    for (s, b), n in nA.items():
        g = pk[(pk.season == s) & (pk.basin == b)]
        g = g.sort_values("val", ascending=(t == "M")).head(int(n))
        sel.append(g)
    sel = pd.concat(sel)
    hit = []
    for r in sel.itertuples():
        a = AHF_BYTIME.get(int(r.time))
        hit.append(bool(a is not None and (C.dist(r.lat, r.lon, a[:, 0], a[:, 1]) <= DKM[t]).any()))
    sel["matched"] = hit
    topsel[t] = (sel, allpk)
    sel.to_csv(f"{OUT}/top_tracks_{t}.csv", index=False)
    bs = season_boot(sel)
    for b in ("atl", "pac"):
        pt, dr, n, k = bs[b]
        lo, hi = ci(dr, 0.95)
        say(f"{t} top-N matched {b}: {100*pt:5.1f}% ({k}/{n}) [{100*lo:5.1f}, {100*hi:5.1f}]")
        s1.append(dict(tracker=t, basin=b, share=pt, lo=lo, hi=hi, n=n))
# ceiling: A ranked by minimum pressure
a = allt.copy()
a["md"] = a.peak_time.map(md)
a = a[((a.md >= 1005) | (a.md <= 425))]
hfids = set(ev.track)
for b in ("atl", "pac"):
    ks = []
    for s in SEAS:
        n = int(nA.get((s, b), 0))
        g = a[(a.season == s) & (a.basin == b)].sort_values("minp").head(n)
        ks.append((n, int(g.track.isin(hfids).sum())))
    ks = np.array(ks)
    pt = ks[:, 1].sum() / ks[:, 0].sum()
    idx = rng.integers(0, len(SEAS), (NB, len(SEAS)))
    dr = ks[:, 1][idx].sum(1) / ks[:, 0][idx].sum(1)
    lo, hi = ci(dr, 0.95)
    say(f"ceiling A-by-minp {b}: {100*pt:5.1f}% ({ks[:,1].sum()}/{ks[:,0].sum()}) [{100*lo:5.1f}, {100*hi:5.1f}]")
    s1.append(dict(tracker="A-minp", basin=b, share=pt, lo=lo, hi=hi, n=int(ks[:, 0].sum())))
pd.DataFrame(s1).to_csv(f"{OUT}/s1_reverse.csv", index=False)

say("\n== S2 seasonal track counts, 5 Oct-25 Apr peak date, 18 seasons ==")
s2 = []
a["basin_t"] = a.basin
for t in ("M", "V"):
    allpk = topsel[t][1]
    for b in ("atl", "pac"):
        x = allpk[allpk.basin == b].groupby("season").size().reindex(SEAS).fillna(0).values
        y = a[a.basin == b].groupby("season").size().reindex(SEAS).fillna(0).values
        pr, pp = stats.pearsonr(x, y)
        sr, sp = stats.spearmanr(x, y)
        say(f"{t} {b}: tracker mean {x.mean():.0f}, A mean {y.mean():.0f}, ratio {x.mean()/y.mean():.2f}, Pearson {pr:.2f} (p {pp:.3f}), Spearman {sr:.2f} (p {sp:.3f})")
        s2.append(dict(tracker=t, basin=b, mean_t=x.mean(), mean_A=y.mean(), pearson=pr, p_pearson=pp, spearman=sr, p_spearman=sp))
pd.DataFrame(s2).to_csv(f"{OUT}/s2_counts.csv", index=False)

say("\n== S3 weekly counts: A HF events vs tracker top-N tracks, Pearson r pooled over weeks ==")
s3 = []


def week_of(time, season):
    """whole weeks since 5 October of the season's start year"""
    t = pd.to_datetime(time.astype("int64").astype(str), format="%Y%m%d%H")
    start = pd.to_datetime(season.astype("int64").astype(str) + "-10-05")
    return (t - start).dt.days // 7


evw = ev[ev.octapr].copy()
evw["wk"] = week_of(evw.first_hf, evw.season)
for t in ("M", "V"):
    sel = topsel[t][0].copy()
    sel["wk"] = week_of(sel.time, sel.season)
    for b in ("atl", "pac"):
        ca = evw[evw.basin == b].groupby(["season", "wk"]).size()
        ct = sel[sel.basin == b].groupby(["season", "wk"]).size()
        grid = pd.MultiIndex.from_product([SEAS, range(0, 30)], names=["season", "wk"])
        df = pd.concat([ca.rename("A"), ct.rename("T")], axis=1).reindex(grid).fillna(0).reset_index()
        r = np.corrcoef(df.A, df["T"])[0, 1]
        by = {s: g for s, g in df.groupby("season")}
        bsr = []
        for _ in range(2000):
            g = pd.concat([by[s] for s in rng.choice(SEAS, len(SEAS))])
            bsr.append(np.corrcoef(g.A, g["T"])[0, 1])
        lo, hi = np.quantile(bsr, [0.025, 0.975])
        say(f"{t} {b}: weekly r = {r:.2f} [{lo:.2f}, {hi:.2f}]  (weeks {len(df)}, A events {int(df.A.sum())}, tracker top-N {int(df['T'].sum())})")
        s3.append(dict(tracker=t, basin=b, r=r, lo=lo, hi=hi, weeks=len(df)))
pd.DataFrame(s3).to_csv(f"{OUT}/s3_weekly.csv", index=False)

# ---------------------------------------------------------------- S5 who is missed
say("\n== S5 who is missed (events in the primary set; run for each tracker and basin whose recall is below 95%) ==")
# Greenland proximity from the 1.5 degree land-sea mask (coarser than the 0.25 degree mask of research/era5/highlat)
import urllib.request
from scipy import ndimage
raw = urllib.request.urlopen(f"{C.WB2}/land_sea_mask/0.0").read()
lsm = np.frombuffer(C.CODEC.decode(raw), "<f4").reshape(240, 121)   # [lon, lat]
land = lsm >= 0.5
lab, _ = ndimage.label(land, structure=np.ones((3, 3)))
gi = lab[int(round(320 / 1.5)) % 240, int(round((72 + 90) / 1.5))]
gl = lab == gi
LA2, LO2 = np.meshgrid(np.radians(C.LAT), np.radians(C.LON), indexing="xy")  # (121, 240)
xyz = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)
tree = cKDTree(xyz[gl])


def dgl(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    p = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    d, _ = tree.query(p)
    return 2 * np.arcsin(np.clip(d / 2, 0, 1)) * C.RE


s5 = []
for t in ("M", "V"):
    e = ev.merge(rec[(t, "base")], on="track")
    e = e[e.octapr].copy()
    e["north60"] = e.peak_lat >= 60
    e["gl100"] = dgl(e.peak_lat.values, e.peak_lon.values) <= 100
    for b in ("atl", "pac"):
        x = e[e.basin == b]
        rc = x.matched.mean()
        if rc >= 0.95:
            say(f"{t} {b}: recall {100*rc:.1f}% >= 95%, S5 not run (as registered)")
            continue
        for fac in ("north60", "gl100", "tc"):
            xx = x.dropna(subset=[fac]) if fac == "tc" else x
            f = xx[fac].astype(bool)
            tab = [[int((f & xx.matched).sum()), int((f & ~xx.matched).sum())], [int((~f & xx.matched).sum()), int((~f & ~xx.matched).sum())]]
            _, p = stats.fisher_exact(tab)
            say(f"{t} {b} {fac}: recall {100*xx[f].matched.mean() if f.any() else float('nan'):.1f}% with flag (n={int(f.sum())}) vs {100*xx[~f].matched.mean():.1f}% without; Fisher p {p:.4f}")
            s5.append(dict(tracker=t, basin=b, factor=fac, p=p))
        for fac in ("n_hf", "n_fix", "minp"):
            u = stats.mannwhitneyu(x[x.matched][fac], x[~x.matched][fac])
            say(f"{t} {b} {fac}: median matched {x[x.matched][fac].median():.1f} vs missed {x[~x.matched][fac].median():.1f}; Mann-Whitney p {u.pvalue:.4f}")
            s5.append(dict(tracker=t, basin=b, factor=fac, p=u.pvalue))
if s5:
    s5 = pd.DataFrame(s5)
    from statsmodels.stats.multitest import multipletests
    s5["q_bh"] = multipletests(s5.p, method="fdr_bh")[1]
    s5.to_csv(f"{OUT}/s5_missed.csv", index=False)
    say("S5 BH q (within family):\n" + s5.round(4).to_string(index=False))
LOG.close()
