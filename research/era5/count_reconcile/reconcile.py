"""Why are archive HF counts 1.5-2x the published OPC ones? Runs the checks of PREREG.md in its order.

usage: python3 -I research/era5/count_reconcile/reconcile.py   (reads committed files only; writes results/)

archive = OPC HF archive (docs/data/hf-lows.json, built from data/hf_lows/*.csv on main).
proxy   = ERA5 pipeline A catalog (research/era5/hf_history/results), a PROXY, used only as a physical benchmark.
Season = 1 June to 31 May, labelled by its start year.
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RES = os.path.join(HERE, "results")
HFH = os.path.join(ROOT, "research", "era5", "hf_history", "results")
os.makedirs(RES, exist_ok=True)
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


BN = {"atl": "Atlantic", "pac": "Pacific"}
PUB_VA = {"atl": {2001: 22, 2002: 23, 2003: 15}, "pac": {2001: 15, 2002: 22, 2003: 22}}   # Von Ahn et al. 2006
REF = {"atl": 20.0, "pac": 19.7}      # Von Ahn mean; Jelenak ~25 as the second reference
JEL = {"atl": 289, "pac": 269}        # Jelenak et al. 2013 slides, 2000-2010 (10 or 11 seasons)
FIRST, LAST = 2004, 2025
MONTHS = [10, 11, 12, 1, 2, 3, 4]


def km(la1, lo1, la2, lo2):
    p = np.pi / 180
    a = np.sin((la2 - la1) * p / 2) ** 2 + np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(a))


# ------------------------------------------------------------------ load
d = json.load(open(os.path.join(ROOT, "docs", "data", "hf-lows.json")))
F = d["lowFields"]
lows = []
for raw in d["lows"]:
    l = dict(zip(F, raw))
    f = pd.DataFrame(l["fixes"], columns=["date", "lat", "lon", "cat", "pres"])
    f["t"] = pd.to_datetime(f.date.astype("int64").astype(str), format="%Y%m%d%H")
    f = f.drop_duplicates("t").sort_values("t").reset_index(drop=True)
    b = l["basin"]
    lon = f.lon.astype(float)
    f["lon"] = lon % 360 if b == "pac" else ((lon + 180) % 360) - 180
    l["f"] = f
    l["key"] = f'{b}:{l["id"]}'
    lows.append(l)
TROP = {"TY", "TS", "TYPHOON", "TYPH"}
rows = []
for l in lows:
    f = l["f"]
    h = f[f.cat == "HF"]
    dh = f[f.cat.isin(["HF", "DHF"])]
    r = dict(key=l["key"], basin=l["basin"], season=l["season"], cls=l["cls"], split=bool(l["split"]), n=len(f),
             n_hf=len(h), n_dhf=int((f.cat == "DHF").sum()), n_warn=len(dh), trop=bool(f.cat.isin(TROP).any()))
    if len(h):
        r.update(t_on=h.t.iloc[0], lat_on=h.lat.iloc[0], lon_on=h.lon.iloc[0], mon=h.t.iloc[0].month)
    r.update(t0=f.t.iloc[0], t1=f.t.iloc[-1])
    rows.append(r)
E = pd.DataFrame(rows)
E["hf"] = E.n_hf >= 1


def per_season(mask, lo=FIRST, hi=LAST):
    s = E[mask & E.season.between(lo, hi)].groupby(["basin", "season"]).size().unstack(0)
    for b in ("atl", "pac"):
        if b not in s:
            s[b] = 0
    return s.reindex(range(lo, hi + 1)).fillna(0)


def boot(x, n=2000, seed=1):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    m = [rng.choice(x, len(x)).mean() for _ in range(n)]
    return np.percentile(m, [2.5, 97.5])


def summ(label, mask, lo=FIRST, hi=LAST):
    s = per_season(mask, lo, hi)
    out = {}
    for b in ("atl", "pac"):
        x = s[b].values
        ci = boot(x)
        out[b] = (x.mean(), ci[0], ci[1])
    say(f"{label:<58s} Atl {out['atl'][0]:5.1f} [{out['atl'][1]:4.1f},{out['atl'][2]:4.1f}]   Pac {out['pac'][0]:5.1f} [{out['pac'][1]:4.1f},{out['pac'][2]:4.1f}]")
    return {b: out[b][0] for b in out}


say("ARCHIVE HF COUNT RECONCILIATION (see PREREG.md). Seasons 2004-05 to 2025-26 unless stated; mean events per season, 95% season-block bootstrap.")
say("Archive tabs start: Pacific Feb 2002, Atlantic Sep 2003 (build_hf_lows.py), so 2002-03 (Pac) and 2003-04 (Atl) are the earliest complete seasons.")
say("")
say("== B0 baseline and check 1 (unit) ==")
base = summ("B0  >=1 HF-category fix", E.hf)
summ("1a  >=1 HF or DHF fix (warned incl. developing)", E.n_warn >= 1)
summ("1b  >=2 HF fixes", E.n_hf >= 2)
summ("1c  >=3 HF fixes (>=12 h at HF)", E.n_hf >= 3)
summ("1d  class 'low' only (drop tipjet, nocentre)", E.hf & (E.cls == "low"))
say("single-fix share of B0 events (Atl, Pac):",
    *[f"{(E[E.hf & (E.basin == b) & E.season.between(FIRST, LAST)].n_hf == 1).mean():.3f}" for b in ("atl", "pac")])

# ------------------------------------------------------------------ check 2: splits, merges, duplicates
say("")
say("== check 2: splits, merges and duplicates ==")
parent = {k: k for k in E.key}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


links = []
for b in ("atl", "pac"):
    L = [l for l in lows if l["basin"] == b]
    for i, A in enumerate(L):
        fa = A["f"]
        for B in L[i + 1:]:
            fb = B["f"]
            if abs((fa.t.iloc[0] - fb.t.iloc[0]).days) > 15:
                continue
            done = False
            for X, Y in ((fa, fb), (fb, fa)):
                gap = (Y.t.iloc[0] - X.t.iloc[-1]).total_seconds() / 3600
                if 0 <= gap <= 12 and km(X.lat.iloc[-1], X.lon.iloc[-1], Y.lat.iloc[0], Y.lon.iloc[0]) <= 500:
                    links.append((A["key"], B["key"], "sequential"))
                    done = True
                    break
            if not done:
                m = fa.merge(fb, on="t", suffixes=("a", "b"))
                if len(m) and (km(m.lata, m.lona, m.latb, m.lonb) <= 500).any():
                    links.append((A["key"], B["key"], "concurrent"))
for a, b_, k in links:
    parent[find(a)] = find(b_)
pd.DataFrame(links, columns=["id_a", "id_b", "kind"]).to_csv(os.path.join(RES, "linked_ids_for_review.csv"), index=False)
E["comp"] = E.key.map(find)
nseq = sum(1 for x in links if x[2] == "sequential")
ncon = sum(1 for x in links if x[2] == "concurrent")
say(f"links found (sequential <=12 h and <=500 km; concurrent within 500 km): {nseq} sequential, {ncon} concurrent, all seasons")
cr = pd.read_csv(os.path.join(ROOT, "data", "hf_lows", "collision_pairs.csv"))
say(f"build's own collision_pairs.csv: {len(cr)} pairs ({cr.kind.value_counts().to_dict()}); build 'split' flag: {int(E.split.sum())} lows")
comp_hf = E[E.hf].groupby("comp").agg(basin=("basin", "first"), season=("season", "min")).reset_index()
sp = comp_hf[comp_hf.season.between(FIRST, LAST)].groupby(["basin", "season"]).size().unstack(0).reindex(range(FIRST, LAST + 1)).fillna(0)
say(f"2   events after linking IDs                              Atl {sp['atl'].mean():5.1f}   Pac {sp['pac'].mean():5.1f}   (baseline {base['atl']:.1f}, {base['pac']:.1f})")
summ("2b  drop lows the build flags as split", E.hf & ~E.split)


# ------------------------------------------------------------------ check 3: domain
def inarea(df, b, dlon=0.0):
    if b == "atl":
        return (df.lon_on <= -35 + dlon) & df.lat_on.between(31, 67)
    return (df.lon_on >= 160 - dlon) & df.lat_on.between(30, 60)


say("")
say("== check 3: domain (OPC area assumed: Atl lon<=-35 and 31-67N; Pac lon>=160E and 30-60N; tested on the first HF fix) ==")
for lab, dl in (("assumed edges", 0), ("lon edge 5 deg outward (Atl 30W, Pac 155E)", 5), ("lon edge 5 deg inward (Atl 40W, Pac 165E)", -5)):
    out = {}
    for b in ("atl", "pac"):
        sub = E[E.hf & (E.basin == b) & E.season.between(FIRST, LAST)]
        ok = inarea(sub, b, dlon=dl)
        out[b] = (ok.groupby(sub.season).sum().reindex(range(FIRST, LAST + 1)).fillna(0).mean(), ok.mean())
    say(f"3   {lab:<46s} Atl {out['atl'][0]:5.1f} ({out['atl'][1]*100:4.1f}% of events)   Pac {out['pac'][0]:5.1f} ({out['pac'][1]*100:4.1f}%)")
sa = E[E.hf & (E.basin == "atl") & E.season.between(FIRST, LAST)]
sp_ = E[E.hf & (E.basin == "pac") & E.season.between(FIRST, LAST)]
say(f"    Atlantic first HF fix: east of 35W {(sa.lon_on > -35).mean()*100:.1f}%, east of 0 {(sa.lon_on > 0).mean()*100:.1f}%, "
    f"north of 67N {(sa.lat_on > 67).mean()*100:.1f}%, south of 31N {(sa.lat_on < 31).mean()*100:.1f}%")
say(f"    Pacific first HF fix:  west of 160E {(sp_.lon_on < 160).mean()*100:.1f}%, west of 150E {(sp_.lon_on < 150).mean()*100:.1f}%, "
    f"north of 60N {(sp_.lat_on > 60).mean()*100:.1f}%, south of 30N {(sp_.lat_on < 30).mean()*100:.1f}%")
bad = E[E.hf & (((E.basin == "atl") & (E.lon_on > 40)) | ((E.basin == "pac") & (E.lon_on < 100)))]
say(f"    first-HF positions in the other basin's longitudes (Atl lon>40E, Pac lon<100E): {len(bad)}")

# ------------------------------------------------------------------ check 4: months
say("")
say("== check 4: season window (month of first HF fix) ==")
summ("4a  Oct-Apr only", E.hf & E.mon.isin(MONTHS))
summ("4b  Nov-Mar only", E.hf & E.mon.isin([11, 12, 1, 2, 3]))

# ------------------------------------------------------------------ check 5: TC
say("")
say("== check 5: tropical and post-tropical ==")
cat = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
life = pd.read_csv(os.path.join(HFH, "lifecycle_events.csv"))
cat = cat[cat.role == "event"].merge(life[["track", "tc"]], on="track", how="left")
tcmap = {}
for _, r in cat.dropna(subset=["archive_events"]).iterrows():
    for k in str(r.archive_events).split(";"):
        tcmap[k.strip()] = bool(r.tc) if not pd.isna(r.tc) else False
E["matched"] = E.key.map(lambda k: k in tcmap)
E["tc_proxy"] = E.key.map(lambda k: tcmap.get(k, False))
sub = E[E.hf & E.season.between(FIRST, LAST)]
say(f"    archive events matched to a pipeline A event: Atl {sub[sub.basin == 'atl'].matched.mean()*100:.1f}%, Pac {sub[sub.basin == 'pac'].matched.mean()*100:.1f}%")
say(f"    TC-linked among matched (proxy 400 km IBTrACS flag): Atl {sub[(sub.basin == 'atl') & sub.matched].tc_proxy.mean()*100:.1f}%, "
    f"Pac {sub[(sub.basin == 'pac') & sub.matched].tc_proxy.mean()*100:.1f}%")
say(f"    archive events with a tropical category (TY/TS/TYPHOON/TYPH) in the track: {int(sub.trop.sum())}")
summ("5   drop TC-linked (matched flag) and tropical-category lows", E.hf & ~E.tc_proxy & ~E.trop)
summ("5b  also drop first HF in Jun-Oct (upper bound for TC season)", E.hf & ~E.tc_proxy & ~E.trop & ~E.mon.isin([6, 7, 8, 9, 10]))

# ------------------------------------------------------------------ check 6: era
say("")
say("== check 6: analysis era (archive counts, no trend fitted) ==")
S = per_season(E.hf, 2001, LAST)
say("per-season archive events (B0), seasons from 2001-02:")
say("season " + " ".join(f"{y:4d}" for y in S.index))
for b in ("atl", "pac"):
    say(f"{b:>6s} " + " ".join(f"{int(v):4d}" for v in S[b]))
for lab, lo, hi in (("2004-05..2009-10 (QuikSCAT era, to Nov 2009)", 2004, 2009), ("2010-11..2016-17", 2010, 2016),
                    ("2017-18..2025-26", 2017, 2025), ("2004-05..2012-13", 2004, 2012), ("2013-14..2025-26", 2013, 2025)):
    say(f"6   {lab:<46s} Atl {S.loc[lo:hi, 'atl'].mean():5.1f}   Pac {S.loc[lo:hi, 'pac'].mean():5.1f}")
pvals = {}
for b, yr in (("pac", 2013), ("atl", 2017)):
    x = S.loc[FIRST:LAST, b]
    pre, post = x[x.index < yr], x[x.index >= yr]
    t = stats.ttest_ind(post, pre, equal_var=False)
    u = stats.mannwhitneyu(post, pre)
    pvals[f"step_{b}_{yr}"] = t.pvalue
    say(f"    step test {BN[b]} {yr}-{str(yr+1)[2:]}: {pre.mean():.1f} before ({len(pre)} seasons) vs {post.mean():.1f} after ({len(post)}); "
        f"Welch t p={t.pvalue:.3f}, Mann-Whitney p={u.pvalue:.3f}")
say(f"    first complete archive seasons: Atl 2003-04 {int(S.loc[2003, 'atl'])} -> 2004-05 {int(S.loc[2004, 'atl'])}; "
    f"Pac 2002-03 {int(S.loc[2002, 'pac'])}, 2003-04 {int(S.loc[2003, 'pac'])} -> 2004-05 {int(S.loc[2004, 'pac'])}")

# ------------------------------------------------------------------ check 7: reproduction
say("")
say("== check 7: reproduction of published counts ==")
area = E.hf & (((E.basin == "atl") & (E.lon_on <= -35) & E.lat_on.between(31, 67)) | ((E.basin == "pac") & (E.lon_on >= 160) & E.lat_on.between(30, 60)))
rule = {"B0 >=1 HF fix": E.hf, "1a HF or DHF": E.n_warn >= 1, "1c >=3 HF fixes": E.n_hf >= 3,
        "assumed OPC area": area, "OPC area, Oct-Apr": area & E.mon.isin(MONTHS), "OPC area, >=2 HF fixes": area & (E.n_hf >= 2)}
tab7 = []
for name, m in rule.items():
    s = per_season(m, 2001, 2025)
    for b in ("atl", "pac"):
        covered = [2002, 2003] if b == "pac" else [2003]    # seasons the archive covers fully
        a = int(sum(s.loc[y, b] for y in covered))
        p = sum(PUB_VA[b][y] for y in covered)
        tab7.append(dict(rule=name, basin=b, seasons="/".join(str(y) for y in covered), archive=a, published=p, ratio=a / p,
                         p=stats.binomtest(a, a + p, 0.5).pvalue))
        j = JEL[b] / 10.5
        tab7.append(dict(rule=name, basin=b, seasons="2004-09 per season vs Jelenak", archive=s.loc[2004:2009, b].mean(), published=j,
                         ratio=s.loc[2004:2009, b].mean() / j, p=np.nan))
t7 = pd.DataFrame(tab7)
t7.to_csv(os.path.join(RES, "reproduction.csv"), index=False)
say(t7.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
for name in ("B0 >=1 HF fix", "assumed OPC area"):
    s = E[rule[name] & (E.t_on >= "2003-10-01") & (E.t_on < "2005-06-01")]
    say(f"    Chelton window, part the archive covers fully (Oct 2003-May 2005, 1.67 seasons), rule '{name}': {len(s)} events = {len(s)/1.67:.1f} per season "
        f"both basins (published 49 per season for Oct 2001-May 2005)")

# ------------------------------------------------------------------ check 8: proxy
say("")
say("== check 8: pipeline A proxy as a physical benchmark (not the archive) ==")
pc = pd.read_csv(os.path.join(HFH, "era5_hf_counts_by_season.csv"))
pc["y"] = pc.season.str[:4].astype(int)
for lo, hi, lab in ((2001, 2003, "2001-02..2003-04 (Von Ahn years)"), (2004, 2009, "2004-05..2009-10"), (2010, 2025, "2010-11..2025-26")):
    s = pc[pc.y.between(lo, hi)]
    say(f"8   proxy events per season {lab:<34s} Atl {s.atl_era5.mean():5.1f}   Pac {s.pac_era5.mean():5.1f}   (archive {s.atl_archive.mean():5.1f}, {s.pac_archive.mean():5.1f})")
vs = pc[pc.y.between(2001, 2003)]
say(f"    proxy vs published, Von Ahn seasons: Atl {int(vs.atl_era5.sum())} vs {sum(PUB_VA['atl'].values())}; Pac {int(vs.pac_era5.sum())} vs {sum(PUB_VA['pac'].values())}")
cat["hit"] = cat.archive_events.notna()
for lo, hi in ((2001, 2003), (2004, 2009), (2010, 2025)):
    q = cat[cat.season.between(lo, hi)]
    say(f"    share of proxy events with an archive match {lo}-{hi}: Atl {q[q.basin == 'atl'].hit.mean()*100:.0f}%, Pac {q[q.basin == 'pac'].hit.mean()*100:.0f}%")
say("    (the proxy's calibration is 2021-26 against the archive, so a proxy level in 2001-04 near the archive's later level means the published counts were"
    " not a lower physical count; it does not say what OPC counted)")


# ---- check 8b (added after the archive checks; the rules are the ones fixed in the sequential chain, run once, none re-picked)
say("")
say("== check 8b (post hoc addition, see PREREG.md deviations): the same matched rules applied to pipeline A events, which cover 1979 on without the archive's start-up ==")
ct = pd.read_csv(os.path.join(HFH, "era5_hf_catalog_tracks.csv"))
ev = cat[["track", "basin", "season", "tc"]].copy()
g = ct[ct.track.isin(ev.track) & (ct.g800 >= 71.7)].sort_values(["track", "time"])
fst = g.groupby("track").agg(n_hf=("g800", "size"), lat_on=("lat", "first"), lon_on=("lon", "first"), t_on=("time", "first")).reset_index()
ev = ev.merge(fst, on="track", how="left")
ev["mon"] = (ev.t_on // 10000 % 100)
ev["lon_on"] = ev.lon_on % 360
ev["area"] = np.where(ev.basin == "atl", (ev.lon_on >= 180) & (ev.lon_on <= 325) & ev.lat_on.between(31, 67), (ev.lon_on >= 160) & ev.lat_on.between(30, 60))
ev["oa"] = ev.mon.isin(MONTHS)
ev["tcf"] = ev.tc.fillna(False).astype(bool)
say(f"    proxy events with a computable first HF fix: {int(ev.n_hf.notna().sum())} of {len(ev)}")
chain = [("all pipeline A events", ev.n_hf.notna()), ("+ >=2 HF fixes", ev.n_hf >= 2), ("+ OPC area", (ev.n_hf >= 2) & ev.area),
         ("+ Oct-Apr", (ev.n_hf >= 2) & ev.area & ev.oa), ("+ drop TC-linked", (ev.n_hf >= 2) & ev.area & ev.oa & ~ev.tcf)]
rows8 = []
for lab, m in chain:
    for lo, hi, per in ((2001, 2003, "2001-02..2003-04"), (2004, 2009, "2004-05..2009-10"), (2010, 2025, "2010-11..2025-26")):
        x = ev[m & ev.season.between(lo, hi)]
        r = dict(rule=lab, period=per)
        for b in ("atl", "pac"):
            r[b] = len(x[x.basin == b]) / (hi - lo + 1)
        rows8.append(r)
t8 = pd.DataFrame(rows8)
t8.to_csv(os.path.join(RES, "proxy_matched.csv"), index=False)
say(t8.to_string(index=False, float_format=lambda v: f"{v:.1f}"))
for b in ("atl", "pac"):
    x = ev[chain[-1][1] & (ev.basin == b) & ev.season.between(2001, 2003)].groupby("season").size().reindex([2001, 2002, 2003]).fillna(0)
    a = int(x.sum()); p = sum(PUB_VA[b].values())
    say(f"    proxy, fully matched rules, Von Ahn seasons {BN[b]}: {list(x.astype(int))} vs published {list(PUB_VA[b].values())}; total {a} vs {p}, ratio {a/p:.2f}, binomial p={stats.binomtest(a, a + p, 0.5).pvalue:.2f}")

# ------------------------------------------------------------------ sequential attribution
say("")
say("== summary: sequential attribution on 2004-05..2025-26 against the Von Ahn reference (Atl 20.0, Pac 19.7) ==")
E["in_area"] = False
for b in ("atl", "pac"):
    ok = inarea(E, b)
    E.loc[E.basin == b, "in_area"] = ok[E.basin == b].fillna(False)


def cnt(mask):
    s = per_season(mask)
    return {b: s[b].mean() for b in ("atl", "pac")}


first_in_comp = E[E.hf].sort_values("t_on").groupby("comp").key.first()
cur = E.hf.copy()
seq = [("B0 baseline", cnt(cur))]
cur = cur & (E.n_hf >= 2)
seq.append(("+1 unit: >=2 HF fixes", cnt(cur)))
cur = cur & E.key.isin(first_in_comp)
seq.append(("+2 one event per linked group", cnt(cur)))
cur = cur & E.in_area
seq.append(("+3 OPC area (assumed edges)", cnt(cur)))
cur = cur & E.mon.isin(MONTHS)
seq.append(("+4 Oct-Apr", cnt(cur)))
cur = cur & ~E.tc_proxy & ~E.trop
seq.append(("+5 drop TC-linked / tropical category", cnt(cur)))
prev = None
for lab, c in seq:
    line = f"{lab:<42s} Atl {c['atl']:5.1f}   Pac {c['pac']:5.1f}"
    if prev:
        for b in ("atl", "pac"):
            gap = seq[0][1][b] - REF[b]
            line += f"   {b}: {c[b]-prev[b]:+5.1f} ({(prev[b]-c[b])/gap*100:4.1f}% of gap)"
    say(line)
    prev = c
fin = seq[-1][1]
for b in ("atl", "pac"):
    gap = seq[0][1][b] - REF[b]
    say(f"{BN[b]}: baseline {seq[0][1][b]:.1f}; after all definitional factors {fin[b]:.1f}; Von Ahn ref {REF[b]}, Jelenak ref 25; "
        f"gap closed {(seq[0][1][b]-fin[b])/gap*100:.0f}%; remaining ratio to Von Ahn {fin[b]/REF[b]:.2f}, to Jelenak {fin[b]/25:.2f}")
pd.DataFrame([dict(step=l, atl=c["atl"], pac=c["pac"]) for l, c in seq]).to_csv(os.path.join(RES, "sequential.csv"), index=False)
S.to_csv(os.path.join(RES, "archive_events_per_season_B0.csv"))
ps = list(pvals.items())
order = np.argsort([p for _, p in ps])
mm = len(ps)
say("BH-adjusted step-test q: " + ", ".join(f"{ps[i][0]} q={min(1, ps[i][1]*mm/(r+1)):.3f}" for r, i in enumerate(order)))
open(os.path.join(RES, "reconcile.txt"), "w").write("\n".join(OUT) + "\n")
